#!/usr/bin/env python3
"""
Setup script for PicoMol - Molecular Visualization and Bioinformatics Suite
"""

from setuptools import setup, find_packages
import os

# Read the README file for long description
def read_readme():
    with open("README.md", "r", encoding="utf-8") as fh:
        return fh.read()

# Read requirements from requirements.txt
def read_requirements():
    with open("requirements.txt", "r", encoding="utf-8") as fh:
        return [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="picomol",
    version="0.0.3",
    author="Jack Magson",
    author_email="your.email@example.com",  # Update with actual email
    description="Molecular Visualization and Bioinformatics Suite",
    long_description=read_readme(),
    long_description_content_type="text/markdown",
    url="https://github.com/ultimafounding/PicoMol",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Education",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Visualization",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Operating System :: OS Independent",
        "Environment :: X11 Applications :: Qt",
    ],
    python_requires=">=3.7",
    install_requires=read_requirements(),
    extras_require={
        "dev": [
            "pytest>=6.0",
            "black>=21.0",
            "flake8>=3.8",
        ],
    },
    entry_points={
        "console_scripts": [
            "picomol=picomol:main",
        ],
    },
    include_package_data=True,
    package_data={
        "": ["*.md", "*.txt", "*.py"],
        "blast_utils": ["*.py"],
    },
    keywords="molecular visualization bioinformatics protein structure PDB BLAST",
    project_urls={
        "Bug Reports": "https://github.com/ultimafounding/PicoMol/issues",
        "Source": "https://github.com/ultimafounding/PicoMol",
        "Documentation": "https://github.com/ultimafounding/PicoMol/blob/main/README.md",
    },
)