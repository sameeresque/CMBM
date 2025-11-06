"""
CMBM: Cloud-by-cloud Multiphase Bayesian ionization Modeling

A Python package for fitting quasar absorption line systems using
Cloudy photoionization models and nested sampling.

CMBM enables users to model complex multi-phase absorption systems by
fitting individual cloud components with physical parameters derived
from photoionization modeling.
"""

__version__ = "0.1.0"
__author__ = "sameer"
__email__ = "sameeresque@gmail.com"

# Import main classes for easy access
from .fitter import CMBMFitter

# Define what gets imported with "from cmbm import *"
__all__ = [
    "CMBMFitter"
]


