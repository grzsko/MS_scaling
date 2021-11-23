/*
 * =====================================================================================
 *
 *       Filename:  optimal_transport.h
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  11/04/19 10:54:46
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Grzegorz Skoraczynski
 *   Organization:  MIMUW
 *
 * =====================================================================================
 */

#include <Eigen/Dense>
#include <tuple>
#include <vector>

#ifndef optimal_transport_dense_h__
#define optimal_transport_dense_h__
double calc_distance_dense(std::vector<double> &mzs1, std::vector<double> &,
                           std::vector<double> &, std::vector<double> &, double,
                           double, double, double, int, const std::string);

Eigen::MatrixXd calc_moves_dense(
    std::vector<double> &, std::vector<double> &, std::vector<double> &,
    std::vector<double> &, double, double, double, double, int,
    const std::string);

double calc_distance_dense(std::vector<double> &mzs1, std::vector<double> &,
                           Eigen::MatrixXd, double, double, double, double, int,
                           const std::string);

Eigen::MatrixXd calc_moves_dense(
    std::vector<double> &, std::vector<double> &, Eigen::MatrixXd, double,
    double, double, double, int, const std::string);
#endif
