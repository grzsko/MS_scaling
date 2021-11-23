/*
 * =====================================================================================
 *
 *       Filename:  optimal_transport.cpp
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  11/04/19 10:55:59
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:
 *   Organization:
 *
 * =====================================================================================
 */
#include <cmath>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "Eigen/src/SparseCore/SparseMatrix.h"
#include "optimal_transport.h"

#include <sys/time.h>
#include <sys/resource.h>


using namespace Eigen;
typedef Map<VectorXd> MapT;
typedef Triplet<double> T;
typedef SparseMatrix<double> SparseM;


VectorXd proxdiv_KL(VectorXd s, VectorXd u, VectorXd p, double lambda,
                    double epsilon) {
    ArrayXd e_power = (-u.array() / (lambda + epsilon)).exp();
    ArrayXd p_s = (p.array() / s.array()).pow(lambda / (lambda + epsilon));
    return (p_s * e_power).matrix();
}

VectorXd proxdiv_TV(VectorXd s, VectorXd u, VectorXd p, double lambda,
                    double epsilon) {
    ArrayXd e_lam = ((lambda - u.array()) / epsilon).exp();
    ArrayXd e_minus_lam = (-(lambda + u.array()) / epsilon).exp();
    return e_lam.cwiseMin(e_minus_lam.cwiseMax(p.array() / s.array())).matrix();
}

VectorXd proxdiv_i(VectorXd s, VectorXd u, VectorXd p, double lambda,
                   double epsilon) {
    return (p.array() / s.array()).matrix();
}

SparseM calc_distance_matrix_sparse(VectorXd mzs1, VectorXd mzs2,
        double epsilon, double dist_upper_bound) {
    // TODO it can be done quicker
    std::vector<T> nonzeros;
    for (int i = 0; i < mzs1.size(); i++) {
        for (int j = 0; j < mzs2.size(); j++) {
            double dist = fabs(mzs1(i) - mzs2(j));
            if (dist <= dist_upper_bound) {
                nonzeros.push_back(T(i, j, dist));
            }
        }
    }
    SparseM dists(mzs1.size(), mzs2.size());
    dists.setFromTriplets(nonzeros.begin(), nonzeros.end());
    return dists;
}

SparseM calc_K_from_dists(SparseM dists, double epsilon) {
    SparseM K(dists.rows(), dists.cols());
    std::vector<T> triplets;
    for (int k = 0; k < dists.outerSize(); ++k) {
        for (SparseM::InnerIterator it(dists, k); it; ++it) {
            int i = it.row();
            int j = it.col();
            triplets.push_back(T(i, j, exp(-it.value() / epsilon)));
        }
    }
    K.setFromTriplets(triplets.begin(), triplets.end());
    return K;
}

void update_K(SparseM &K, SparseM dists, double epsilon, VectorXd u, VectorXd v) {
    for (int outer_index = 0; outer_index < dists.outerSize(); ++outer_index) {
        SparseM::InnerIterator K_it(K, outer_index);
        for (SparseM::InnerIterator dist_it(dists, outer_index); dist_it;
                ++dist_it) {
            int i = dist_it.row();
            int j = dist_it.col();
            K_it.valueRef() = exp((u(i) + v(j) - dist_it.value()) / epsilon);
            ++K_it;
        }
        assert (!K_it && "K and dists matrix unzipped!");
    }
}

SparseM multiply_cwise_aKb(VectorXd a, SparseM K, VectorXd b) {
    // computes (a_i * K_{ij} * b_j) for all {ij}
    // it turns out that asDiagonal is only a wrapper on vector, so it does use
    // only linear memory. Moreover, it works extremely faster.
    return a.asDiagonal() * K * b.asDiagonal();
}

SparseM  calc_transport_plan(VectorXd ints1, VectorXd ints2,VectorXd mzs1,
                             VectorXd mzs2, SparseM sparse_dists,
                             double lambda, double epsilon, double tol,
                             double threshold, double max_iter) {
    // Declarations and initialization
    int n1 = ints1.rows();
    int n2 = ints2.rows();

    VectorXd a = VectorXd::Ones(n1);
    VectorXd b = VectorXd::Ones(n2);

    VectorXd u = VectorXd::Zero(n1);
    VectorXd v = VectorXd::Zero(n2);

    int tick = 0;
    bool coverged = false;
    double coverging_val = 2 * tol; // initializeas too much

    SparseM K = calc_K_from_dists(sparse_dists, epsilon);
    SparseM tp = multiply_cwise_aKb(a, K, b); // tp like a transport plan
    SparseM old_tp;

    do {
        old_tp = tp;

        a = proxdiv_TV(K * b, u, ints1, lambda, epsilon); // * is matrix-vector
        b = proxdiv_TV(K.transpose() * a, v, ints2, lambda, epsilon);

        if ((((a.array().log().abs() - threshold) > 0).any()) or
            (((b.array().log().abs() - threshold) > 0).any())) {
            // Stabilizing
            u += epsilon * a.array().log().matrix();
            v += epsilon * b.array().log().matrix();
            // below K becomes exp((u_i + v_j - dists_{ij}) / eps) forall i,j
            update_K(K, sparse_dists, epsilon, u, v);
            b = VectorXd::Ones(n2);
        }

		if (tick % 10 == 0) {
            // for speeding up computations we check convergence not every loop
            tp = multiply_cwise_aKb(a, K, b);
            SparseM tp_diff = (tp - old_tp).cwiseAbs();
            coverging_val = tp_diff.coeffs().maxCoeff();
        }
        tick++;
        coverged = coverging_val < tol or tick >= max_iter;
    } while (not coverged);

    return tp;
}

double calc_distance_from_plan(SparseM transport_plan, SparseM distances,
                               VectorXd ints1, VectorXd ints2, double lambda) {
    double transport = (transport_plan.cwiseProduct(distances)).sum();

    double trash = 0;
    VectorXd ones_cols = VectorXd::Ones(transport_plan.cols());
    VectorXd ones_rows = VectorXd::Ones(transport_plan.rows());
    trash += lambda * ((transport_plan * ones_cols) - ints1).cwiseAbs().sum();
    // ^^^ rowwise sums
    trash += lambda * (((ones_rows.transpose() * transport_plan).transpose() - ints2).cwiseAbs().sum());
    // ^^ colwise sum is row, not column
    return transport + trash;
}

void saveData(std::string fileName, MatrixXd  matrix) {
    //https://eigen.tuxfamily.org/dox/structEigen_1_1IOFormat.html
    const static IOFormat CSVFormat(FullPrecision, DontAlignCols, ", ", "\n");

    std::ofstream file(fileName);
    if (file.is_open())
    {
        file << matrix.format(CSVFormat);
        file.close();
    }
}

double calc_distance(std::vector<double> & mzs1, std::vector<double> & ints1, std::vector<double> & mzs2,
                         std::vector<double> & ints2, double lambda,
                         double epsilon, double tol, double threshold,
                         int max_iter, double dist_upper_bound) {
    /*  TODO description
     *
     */
    MapT mzs1_map(mzs1.data(), mzs1.size());
    MapT mzs2_map(mzs2.data(), mzs2.size());
    MapT ints1_map(ints1.data(), ints1.size());
    MapT ints2_map(ints2.data(), ints2.size());
    SparseM dists = calc_distance_matrix_sparse(mzs1_map, mzs2_map, epsilon,
                                                dist_upper_bound);
    SparseM transport_plan = calc_transport_plan(ints1_map, ints2_map, mzs1_map,
            mzs2_map, dists, lambda, epsilon, tol, threshold, max_iter);

    double distance = calc_distance_from_plan(transport_plan, dists, ints1_map, ints2_map,
                                   lambda);
    return distance;
}

std::vector<std::tuple<int, int, double>> calc_moves(std::vector<double> & mzs1, std::vector<double> & ints1, std::vector<double> & mzs2,
                         std::vector<double> & ints2, double lambda,
                         double epsilon, double tol, double threshold,
                         int max_iter, double dist_upper_bound) {
    /*  TODO do this better
     *
     *  ints1, and ints2 should be normalized */
    MapT mzs1_map(mzs1.data(), mzs1.size());
    MapT mzs2_map(mzs2.data(), mzs2.size());
    MapT ints1_map(ints1.data(), ints1.size());
    MapT ints2_map(ints2.data(), ints2.size());
    SparseM dists = calc_distance_matrix_sparse(mzs1_map, mzs2_map, epsilon,
            dist_upper_bound);
    SparseM transport_plan = calc_transport_plan(ints1_map, ints2_map, mzs1_map,
            mzs2_map, dists, lambda, epsilon, tol, threshold, max_iter);
	std::vector<std::tuple<int, int, double>> moves;
	for (int k = 0; k < transport_plan.outerSize(); ++k) {
		for (SparseMatrix<double>::InnerIterator it(transport_plan, k); it; ++it) {
            std::tuple <int, int, double> step(it.row(), it.col(), it.value());
            moves.push_back(step);
		}
	}
    return moves;
}

