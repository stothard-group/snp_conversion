#!/usr/bin/env python3

from setuptools import setup, find_packages

setup(
    name="snp_conversion_ekherman",
    author="Emily Herman",
    description="Tools for manipulating and converting SNP information",
    url="https://github.com/stothard-group/snp_conversion",
    version="2.0",
    packages=find_packages(include=['lib/*']),
    install_requires=["numpy", "pandas", "python-dateutil", "dask==2.14.0", "scikit-allel", "toolz"],
    python_requires='>=3.6',
    classifiers=["Programming Language :: Python :: 3",
                 "Operating System :: Unix",
                 "License :: OSI Approved :: GNU General Public License v3 (GPLv3)"],

)