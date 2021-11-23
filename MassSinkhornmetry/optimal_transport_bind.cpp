#include "optimal_transport.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

PYBIND11_MODULE(optimal_transport, m) {
    m.doc() = "Compare mass spectra by optimal transport. Sparse implementation.";

    m.def("calc_distance", &calc_distance, "Distance sparse implementation.");
    m.def("calc_moves", &calc_moves, "Transport plan sparse implementation.");

}
