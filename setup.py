from setuptools import find_packages, setup
from pybind11.setup_helpers import Pybind11Extension


module_dense = Pybind11Extension(
    'optimal_transport_dense',
    sources=["MassSinkhornmetry/optimal_transport_dense.cpp",
             "MassSinkhornmetry/optimal_transport_dense_bind.cpp"],
    include_dirs = ['include/eigen'],
    extra_compile_args=['-DEIGEN_MPL2_ONLY -O3 -fno-finite-math-only -undefined dynamic_lookup -msse2']
)
module_dense.cxx_std=11 # Of course, not documented in pybind11. Never try to
# manually add standard in extra_compile_args. It does not work for strange
# reasons.

module_sparse = Pybind11Extension(
    'optimal_transport',
    sources=["MassSinkhornmetry/optimal_transport.cpp",
             "MassSinkhornmetry/optimal_transport_bind.cpp"],
    include_dirs = ['include/eigen'],
    extra_compile_args=['-DEIGEN_MPL2_ONLY -O3 -fno-finite-math-only -msse2 -undefined dynamic_lookup']
)
module_sparse.cxx_std=11

# Should work in Linuxes with g++ and MacOSes with clang...
setup(
    name='MassSinkhornmetry',
    packages=find_packages(),
    version='0.8.0',
    description='Comparing MS spectra using Sinkhorn algorithm.',
    author='Grzegorz Skoraczynski',
    license='MIT',
    install_requires=["pybind11", "numpy"],
    python_requires='>=3.6, <4',
    ext_modules=[module_dense, module_sparse]
)
