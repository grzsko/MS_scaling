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
#ifndef optimal_transport_h__
#define optimal_transport_h__

#include <vector>
#include <tuple>
double calc_distance(std::vector<double> & mzs1, std::vector<double> & ints1, std::vector<double> & mzs2,
                         std::vector<double> & ints2, double lambda,
                         double epsilon, double tol, double threshold,
                         int max_iter, double dist_upper_bound);

std::vector<std::tuple<int, int, double>> calc_steps(std::vector<double> & mzs1, std::vector<double> & ints1, std::vector<double> & mzs2,
                         std::vector<double> & ints2, double lambda,
                         double epsilon, double tol, double threshold,
                         int max_iter, double dist_upper_bound);
#endif
