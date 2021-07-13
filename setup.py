from setuptools import find_packages, setup
from distutils.core import setup, Extension

module1 = Extension('MS_scalingCppToPy',
                    sources = ['MS_scaling/optimal_transport.cpp'],
                    include_dirs = ['MS_scaling/eigen'])


setup(
    name='MS_scaling',
    packages=find_packages(),
    version='0.1.0',
    description='Comparing MS spectra using Sinkhorn algorithm.',
    author='Grzegorz Skoraczynski',
    license='MIT',
    install_requires=["pybind11"],
    python_requires='>=3.6, <4',
    ext_modules=[module1],
)
