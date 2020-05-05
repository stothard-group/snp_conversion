#!/usr/bin/env python3

from setuptools import setup, find_packages

setup(
    name="SNP Conversion",
    version="2.0",
    packages=find_packages(include=['lib/*']),
    install_requires=["numpy", "pandas", "python-dateutil", "scikit-allel", "toolz"],
    python_requires='>=3',



)