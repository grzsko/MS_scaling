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
 *         Author:  YOUR NAME (),
 *   Organization:
 *
 * =====================================================================================
 */
#include "optimal_transport_dense.h"
#include <stdlib.h>
#include <Eigen/Dense>
#include <fstream>
#include <iostream>


using namespace Eigen;
typedef Map<VectorXd> MapT;

#define KL "KL"
#define I "I"
#define TV "TV"

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

MatrixXd calc_distance_matrix(VectorXd mzs1, VectorXd mzs2) {
  return (mzs1 * VectorXd::Ones(mzs2.rows()).transpose() -
          VectorXd::Ones(mzs1.rows()) * mzs2.transpose()).cwiseAbs();
}

MatrixXd calc_transport_plan(VectorXd ints1, VectorXd ints2, MatrixXd dists,
                             double lambda, double epsilon, double tol,
                             double threshold, int max_iter,
                             const std::string method) {
  int n1 = ints1.rows();
  int n2 = ints2.rows();
  MatrixXd K_0 = (-dists / epsilon).array().exp().matrix();
  MatrixXd K = K_0;

  VectorXd a = VectorXd::Ones(n1);
  VectorXd b = VectorXd::Ones(n2);

  VectorXd u = VectorXd::Zero(n1);
  VectorXd v = VectorXd::Zero(n2);

  MatrixXd transport_plan = a.asDiagonal() * K * b.asDiagonal();
  MatrixXd transport_plan_old;

  int tick = 0;
  bool coverged = false;
  double coverging_val = 2 * tol;  // initialize as too much for stop criterion

  do {
    transport_plan_old = transport_plan;

    if (method == TV) {
      a = proxdiv_TV(K * b, u, ints1, lambda, epsilon);  // * is matrix-matrix
      b = proxdiv_TV(K.transpose() * a, v, ints2, lambda, epsilon);
    } else if (method == KL) {
      a = proxdiv_KL(K * b, u, ints1, lambda, epsilon);  // * is matrix-matrix
      b = proxdiv_KL(K.transpose() * a, v, ints2, lambda, epsilon);
    } else if (method == I) {
      a = proxdiv_i(K * b, u, ints1, lambda, epsilon);  // * is matrix-matrix
      b = proxdiv_i(K.transpose() * a, v, ints2, lambda, epsilon);
    } else {
      throw std::invalid_argument("Unknown proxdiv method: " + method);
    }
    if ((((a.array().log().abs() - threshold) > 0).any()) or
        (((b.array().log().abs() - threshold) > 0).any())) {
      // Stabilizing
      u += epsilon * a.array().log().matrix();
      v += epsilon * b.array().log().matrix();
      // below is exp((u_i + v_j - dists_{ij}) / eps) forall i,j
      K = (u / epsilon).array().exp().matrix().asDiagonal() * K_0 *
          (v / epsilon).array().exp().matrix().asDiagonal();
      b = VectorXd::Ones(n2);
    }

    if (tick % 10 ==
        0) {  // for better efficiency we calculate tp only one per 10 times
      // below is (a_i * K_{ij} * b_j)_{ij}
      transport_plan = a.asDiagonal() * K * b.asDiagonal();
      coverging_val =
          (transport_plan - transport_plan_old).cwiseAbs().maxCoeff();
    }
    tick++;
    coverged = coverging_val < tol or tick >= max_iter;
  } while (not coverged);
  return transport_plan;
}

double distance_from_plan(MatrixXd transport_plan, MatrixXd distances,
                          VectorXd ints1, VectorXd ints2, double lambda) {
  double transport = (transport_plan.array() * distances.array()).sum();
  // ^^^ coeff-wise multiplication ^^^

  double trash = 0;
  trash += lambda * ((transport_plan.rowwise().sum() - ints1).cwiseAbs().sum());
  trash +=
      lambda *
      ((transport_plan.colwise().sum().transpose() - ints2).cwiseAbs().sum());
  // ^^ colwise sum is row, not column
  return transport + trash;
}

std::vector<std::tuple<int, int, double>> steps_from_plan(
    MatrixXd transport_plan) {
  std::vector<std::tuple<int, int, double>> steps;
  for (int i = 0; i < transport_plan.cols(); i++) {
    for (int j = 0; j < transport_plan.rows(); j++) {
      double value = transport_plan(j, i);
      if (value > 0) {
        steps.push_back(std::tuple<int, int, double>{j, i, value});
      }
    }
  }
  return steps;
}

double calc_distance_dense(std::vector<double>& mzs1,
                           std::vector<double>& ints1,
                           std::vector<double>& mzs2,
                           std::vector<double>& ints2, double lambda,
                           double epsilon, double tol, double threshold,
                           int max_iter, const std::string method) {
  /*  TODO description
   *
   *  ints1, and ints2 should be normalized */
  MapT mzs1_map(mzs1.data(), mzs1.size());
  MapT mzs2_map(mzs2.data(), mzs2.size());
  MapT ints1_map(ints1.data(), ints1.size());
  MapT ints2_map(ints2.data(), ints2.size());

  MatrixXd dists = calc_distance_matrix(mzs1_map, mzs2_map);
  MatrixXd transport_plan =
      calc_transport_plan(ints1_map, ints2_map, dists, lambda, epsilon, tol,
                          threshold, max_iter, method);

  return distance_from_plan(transport_plan, dists, ints1_map, ints2_map,
                            lambda);
}

std::vector<std::tuple<int, int, double>> calc_steps_dense(
    std::vector<double>& mzs1, std::vector<double>& ints1,
    std::vector<double>& mzs2, std::vector<double>& ints2, double lambda,
    double epsilon, double tol, double threshold, int max_iter,
    const std::string method) {
  /*  TODO
   *
   *  ints1, and ints2 should be normalized */
  MapT mzs1_map(mzs1.data(), mzs1.size());
  MapT mzs2_map(mzs2.data(), mzs2.size());
  MapT ints1_map(ints1.data(), ints1.size());
  MapT ints2_map(ints2.data(), ints2.size());

  MatrixXd dists = calc_distance_matrix(mzs1_map, mzs2_map);
  MatrixXd transport_plan =
      calc_transport_plan(ints1_map, ints2_map, dists, lambda, epsilon, tol,
                          threshold, max_iter, method);

  return steps_from_plan(transport_plan);
}

double calc_distance_dense(std::vector<double>& ints1,
                           std::vector<double>& ints2, MatrixXd dists,
                           double lambda, double epsilon, double tol,
                           double threshold, int max_iter,
                           const std::string method) {
  /*  TODO description
   *
   *  ints1, and ints2 should be normalized */
  MapT ints1_map(ints1.data(), ints1.size());
  MapT ints2_map(ints2.data(), ints2.size());

  MatrixXd transport_plan =
      calc_transport_plan(ints1_map, ints2_map, dists, lambda, epsilon, tol,
                          threshold, max_iter, method);

  return distance_from_plan(transport_plan, dists, ints1_map, ints2_map,
                            lambda);
}

std::vector<std::tuple<int, int, double>> calc_steps_dense(
    std::vector<double>& ints1, std::vector<double>& ints2, MatrixXd dists,
    double lambda, double epsilon, double tol, double threshold, int max_iter,
    const std::string method) {
  /*  TODO
   *
   *  ints1, and ints2 should be normalized */
  MapT ints1_map(ints1.data(), ints1.size());
  MapT ints2_map(ints2.data(), ints2.size());

  MatrixXd transport_plan =
      calc_transport_plan(ints1_map, ints2_map, dists, lambda, epsilon, tol,
                          threshold, max_iter, method);

  return steps_from_plan(transport_plan);
}
