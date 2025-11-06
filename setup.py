from setuptools import setup, find_packages

setup(
    name="cmbm",
    version="0.1.0",
    description="Cloud-by-cloud Multiphase Bayesian ionization Modeling (CMBM) for quasar absorption systems",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="Sameer",
    author_email="sameeresque@gmail.com",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.24.4",
        "scipy>=1.10.1", 
        "pandas>=2.0.3",
        "astropy>=5.2.2",
        "ultranest>=4.4.0",
        "VoigtFit>=3.21.3", 
        "mpi4py>=3.1.5",
        "pyyaml>=5.0",
        "roman>=3.0",
    ],
    extras_require={
        "dev": ["pytest", "pytest-cov", "black", "flake8"],
        "docs": ["sphinx", "sphinx-rtd-theme"],
    },
    python_requires=">=3.8",
    include_package_data=True,
    package_data={
        "cmbm": ["config/*.yaml"],
    },
    entry_points={
        "console_scripts": [
            "cmbm-fit=cmbm.cli:main",
        ],
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
    keywords="astronomy spectroscopy ionization-modeling bayesian nested-sampling absorption-lines multiphase cmbm",
    project_urls={
        "Bug Reports": "https://github.com/sameeresque/cmbm/issues",
        "Source": "https://github.com/sameeresque/cmbm",
    },
)
