from setuptools import find_packages, setup
from pybind11.setup_helpers import Pybind11Extension



module_dense = Pybind11Extension(
    'optimal_transport_dense',
    sources=["MassSinkhornmetry/optimal_transport_dense.cpp",
             "MassSinkhornmetry/optimal_transport_dense_bind.cpp"],
    include_dirs = ['include/eigen'],
    extra_compile_args=['-DEIGEN_MPL2_ONLY -O3 -fno-finite-math-only -msse2']
)

module_sparse = Pybind11Extension(
    'optimal_transport',
    sources=["MassSinkhornmetry/optimal_transport.cpp",
             "MassSinkhornmetry/optimal_transport_bind.cpp"],
    include_dirs = ['include/eigen'],
    extra_compile_args=['-DEIGEN_MPL2_ONLY -O3 -fno-finite-math-only -msse2']
)



setup(
    name='MassSinkhornmetry',
    packages=find_packages(),
    version='0.8.0',
    description='Comparing MS spectra using Sinkhorn algorithm.',
    author='Grzegorz Skoraczynski',
    license='MIT',
    install_requires=["pybind11", "numpy"],
    python_requires='>=3.6, <4',
    ext_modules=[module_dense, module_sparse],
)
# ext_modules = [
#     Pybind11Extension(
#         "python_example",
#         sorted(glob("src/*.cpp")),  # Sort source files for reproducibility
#     ),
# ]
#
# setup(
#     ...,
#     ext_modules=ext_modules
# )
