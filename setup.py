#!/usr/bin/env python

from setuptools import setup, Extension

setup(
    author="Ajay Vibhute",
    author_email="ajayvibhute@gmail.com",
    name="find_objects",
    version="1.0.0",
    url="https://github.com/ajayvibhute/pyfind_objects",
    description="Python package to find sources in the radio image",
    packages=[],
    setup_requires=["astropy","matplotlib","math","time","numpy"],
    install_requires=["astropy","matplotlib","math","time","numpy"],
    python_requires=">=3.6"
)


