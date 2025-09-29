#!/usr/bin/env python3
"""
Setup script for NGS Pipeline
"""

from setuptools import setup, find_packages
import os

# Read the README file for long description
def read_readme():
    readme_path = os.path.join(os.path.dirname(__file__), 'README.md')
    if os.path.exists(readme_path):
        with open(readme_path, 'r', encoding='utf-8') as f:
            return f.read()
    return "NGS Analysis Pipeline"

# Read requirements from file if it exists
def read_requirements():
    req_path = os.path.join(os.path.dirname(__file__), 'requirements.txt')
    if os.path.exists(req_path):
        with open(req_path, 'r', encoding='utf-8') as f:
            return [line.strip() for line in f.readlines() if line.strip() and not line.startswith('#')]
    return [
        'numpy>=1.24.0',
        'pandas>=2.0.0',
        'scipy>=1.10.0',
        'matplotlib>=3.7.0',
        'seaborn>=0.12.0',
        'scikit-learn>=1.3.0',
        'typer>=0.9.0',
        'rich>=13.0.0',
        'click>=8.1.0',
        'pyyaml>=6.0.0',
        'jinja2>=3.1.0',
        'requests>=2.31.0',
        'plotly>=5.17.0',
        'kaleido>=0.2.1',
    ]

setup(
    name="ngs_pipeline",
    version="1.0.0",
    author="NGS Pipeline Team",
    author_email="pipeline@example.com",
    description="A comprehensive NGS analysis pipeline for multi-condition experiments",
    long_description=read_readme(),
    long_description_content_type="text/markdown",
    url="https://github.com/example/ngs_pipeline",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: R",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Artificial Intelligence",
    ],
    python_requires=">=3.11",
    install_requires=read_requirements(),
    extras_require={
        'dev': [
            'pytest>=7.4.0',
            'pytest-cov>=4.1.0',
            'black>=23.7.0',
            'flake8>=6.0.0',
            'mypy>=1.5.0',
        ],
        'notebook': [
            'jupyter>=1.0.0',
            'jupyterlab>=4.0.0',
            'ipykernel>=6.25.0',
        ],
        'bio': [
            # Bioinformatics packages that work on Windows
            'pysam>=0.21.0',
            'pyBigWig>=0.3.22',
            'pybedtools>=0.9.1',
        ]
    },
    entry_points={
        'console_scripts': [
            'ngs_pipeline=ngs_pipeline.cli:app',
        ],
    },
    include_package_data=True,
    package_data={
        'ngs_pipeline': [
            'templates/*.html',
            'templates/*.css',
            'templates/*.js',
            'data/*.yml',
            'data/*.json',
        ],
    },
    project_urls={
        "Bug Reports": "https://github.com/example/ngs_pipeline/issues",
        "Source": "https://github.com/example/ngs_pipeline",
        "Documentation": "https://github.com/example/ngs_pipeline/blob/main/README.md",
    },
)