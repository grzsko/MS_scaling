from setuptools import find_packages, setup

from glob import glob
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension



module1 = Pybind11Extension('MassSinkhornmetryCppToPy',
                            sources =["MassSinkhornmetry/optimal_transport_dense.cpp"],
                            include_dirs = ['include/eigen'])


setup(
    name='MassSinkhornmetry',
    packages=find_packages(),
    version='0.8.0',
    description='Comparing MS spectra using Sinkhorn algorithm.',
    author='Grzegorz Skoraczynski',
    license='MIT',
    install_requires=["pybind11", "numpy"],
    python_requires='>=3.6, <4',
    ext_modules=[module1],
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
