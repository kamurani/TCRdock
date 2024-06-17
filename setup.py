"""Setuptools configuration for the project."""

import os
from setuptools import find_packages, setup

PROJECT_ROOT = os.path.dirname(os.path.realpath(__file__)) 

setup(
    name="tcrdock",
    version="2.0.0",
    description="Python tools for TCR:peptide-MHC modeling and analysis",
    long_description="Python tools for TCR:peptide-MHC modeling and analysis",
    author="Philip Bradley",
    author_email="pbradley@fredhutch.org", 
    
    packages=find_packages(),

    include_package_data=True,
    install_requires=[
        "biopython==1.79",
        "click",
        "numpy",
        "pandas",
        "scipy",
        "matplotlib",
    ],
    entry_points={
        "console_scripts": [
            "tcrdock = tcrdock.cli:main",
            "tcrd = tcrdock.cli:main",
        ],
    },
)