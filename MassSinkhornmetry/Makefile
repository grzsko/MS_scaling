all: optimal_transport.so optimal_transport_dense.so

PYBIND_INCLUDE := $(shell python3 -m pybind11 --includes)

optimal_transport.so: optimal_transport.cpp optimal_transport.h optimal_transport_bind.cpp
	c++ -DEIGEN_MPL2_ONLY -O3 -Wall -msse2 -shared -std=c++11 -I eigen -fPIC $(PYBIND_INCLUDE) optimal_transport_bind.cpp optimal_transport.cpp -o optimal_transport.so

optimal_transport_dense.so: optimal_transport_dense.h optimal_transport_dense_bind.cpp optimal_transport_dense.cpp
	c++ -O3 -Wall -msse2 -ffast-math -fno-finite-math-only -shared -std=c++11 -I eigen -fPIC $(PYBIND_INCLUDE) optimal_transport_dense_bind.cpp optimal_transport_dense.cpp -o optimal_transport_dense.so
