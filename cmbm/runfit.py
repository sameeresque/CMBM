import cmbm
import subprocess
import os
import argparse

parser = argparse.ArgumentParser(description='Run CMBM fit')

parser.add_argument('--config', required=True, help='Path to config YAML file')
args = parser.parse_args()

fitter = cmbm.CMBMFitter(args.config)
fitter.load_data()

results = fitter.run_fit()
