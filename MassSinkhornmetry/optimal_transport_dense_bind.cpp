#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "optimal_transport_dense.h"

PYBIND11_MODULE(optimal_transport_dense, m) {
  m.doc() = "Compare mass spectra by optimal transport, dense implementation.";
  m.def("calc_distance_dense",
        static_cast<double (*)(std::vector<double> &, std::vector<double> &,
                               std::vector<double> &, std::vector<double> &,
                               double, double, double, double, int,
                               const std::string)>(&calc_distance_dense),
        "Distance dense implementation.");
  m.def("calc_distance_dense",
        static_cast<double (*)(std::vector<double> &, std::vector<double> &,
                               Eigen::MatrixXd, double, double, double, double,
                               int, const std::string)>(&calc_distance_dense),
        "Distance dense implementation.");
  m.def("calc_moves_dense",
        static_cast<Eigen::MatrixXd (*)(
            std::vector<double> &, std::vector<double> &, std::vector<double> &,
            std::vector<double> &, double, double, double, double, int,
            const std::string)>(&calc_moves_dense),
        "Transport plan dense implementation.");
  m.def("calc_moves_dense",
        static_cast<Eigen::MatrixXd (*)(
          std::vector<double> &, std::vector<double> &, Eigen::MatrixXd, double,
          double, double, double, int, const std::string)>(&calc_moves_dense),
      "Transport plan dense implementation.");
}
