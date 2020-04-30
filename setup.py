#!/usr/bin/python

from setuptools import setup, find_packages

setup(
    name="SNP Conversion",
    version="2.0",
    packages=find_packages(),
    install_requires=["numpy", "pandas", "python-dateutil", "scikit-allel", "toolz"]



)