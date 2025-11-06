#!/usr/bin/env python
"""
Basic usage example for CMBM (Cloud-by-cloud Multiphase Bayesian ionization Modeling)

This example demonstrates how to use CMBM to fit a simple two-phase absorption system
using both direct parameters and configuration files.

Requirements:
- VoigtFit dataset (.hdf5)
- VoigtFit initial fit results (.pickle) 
- Velocity masks (.pkl)
- Cloudy photoionization grids
- Atomic transition data file

Author: Your Name
Date: 2024
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Import CMBM
try:
    import cmbm
except ImportError:
    print("CMBM not installed. Install with: pip install cmbm")
    sys.exit(1)

# ============================================================================
# EXAMPLE 1: Basic usage with direct parameters
# ============================================================================

def example_1_direct_parameters():
    """
    Example 1: Fit a simple two-phase system using direct parameters.
    
    This example fits a system with:
    - Low ionization phase: SiII and MgII 
    - High ionization phase: CIV and SiIV
    """
    
    print("="*60)
    print("EXAMPLE 1: Direct parameter specification")
    print("="*60)
    
    # Define your data paths (UPDATE THESE FOR YOUR DATA)
    data_paths = {
        'dataset_path': '/path/to/your/dataset.hdf5',
        'voigtfit_path': '/path/to/your/bestfit_result.pickle',
        'mask_path': '/path/to/your/masks.pkl',
        'transition_library_path': '/path/to/your/atomdata.dat',
        'grid_paths': {
            'gridpath': '/path/to/your/cloudy/grids/',
            'gridpath_thick': '/path/to/your/cloudy/grids/'
        }
    }
    
    # Check if files exist (for demonstration)
    for key, path in data_paths.items():
        if key != 'grid_paths' and not os.path.exists(path):
            print(f"WARNING: {key} file not found: {path}")
            print("Please update the paths in this example to point to your data")
            return
    
    # Define the absorption system structure
    phases = {
        'Low_ionization': ['SiII_0', 'MgII_1'],   # Low-ion clouds
        'High_ionization': ['CIV_0', 'SiIV_1']   # High-ion clouds  
    }
    
    # Select spectral lines to fit
    use_lines = [
        'SiII_1190',    # Low ionization
        'MgII_2796',
        'CIV_1548',     # High ionization  
        'SiIV_1393',
        'HI_1215'       # Neutral hydrogen
    ]
    
    # Create CMBM fitter
    print("Creating CMBM fitter...")
    fitter = cmbm.CMBMFitter(
        **data_paths,
        phases=phases,
        use_lines=use_lines,
        grid_z='0.2',
        output_dir='example_1_results'
    )
    
    # Load data
    print("Loading data...")
    fitter.load_data()
    
    # Print setup summary
    fitter.summary()
    
    # Run the fit
    print("Running nested sampling fit...")
    print("This may take 1-4 hours depending on system complexity...")
    
    results = fitter.run_fit(
        min_num_live_points=400,  # Use more for higher precision
        dlogz=0.5
    )
    
    # Get organized results
    summary = fitter.get_results()
    
    print(f"\nFitting completed!")
    print(f"Log evidence: {summary['log_evidence']:.2f} ± {summary['log_evidence_error']:.2f}")
    print(f"Number of samples: {summary['n_samples']}")
    
    # Print parameter results
    print("\nBest-fit parameters:")
    for param, stats in summary['parameter_summaries'].items():
        median = stats['median']
        lower = median - stats['percentiles']['16th']  
        upper = stats['percentiles']['84th'] - median
        print(f"  {param:20s}: {median:8.4f} +{upper:6.4f} -{lower:6.4f}")
    
    # Get best-fit model
    model_profiles, data_profiles = fitter.get_best_fit_model()
    
    print(f"\nResults saved to: example_1_results/")
    print("Check the output directory for detailed results and plots")
    
    return fitter, results


# ============================================================================
# EXAMPLE 2: Using a configuration file
# ============================================================================

def example_2_config_file():
    """
    Example 2: Fit using a YAML configuration file.
    
    This demonstrates the recommended way to use CMBM for complex setups.
    """
    
    print("\n" + "="*60)
    print("EXAMPLE 2: Configuration file usage")
    print("="*60)
    
    # First, create an example configuration file
    config_filename = 'example_config.yaml'
    create_example_config(config_filename)
    
    print(f"Created example configuration: {config_filename}")
    print("You should edit this file with your actual data paths\n")
    
    # Create fitter from config file
    try:
        fitter = cmbm.CMBMFitter(config_filename)
        print("Configuration loaded successfully!")
        
        # You can still override specific parameters
        fitter.config['output_dir'] = 'example_2_results'
        fitter.config['min_num_live_points'] = 600
        
        print("Configuration summary:")
        fitter.summary()
        
        # Load data (will fail if paths aren't updated)
        # fitter.load_data()
        # results = fitter.run_fit()
        
        print("To run this example:")
        print("1. Edit example_config.yaml with your data paths")
        print("2. Run: fitter.load_data()")
        print("3. Run: results = fitter.run_fit()")
        
    except FileNotFoundError as e:
        print(f"Error: {e}")
        print("Please check your data paths in the configuration file")


# ============================================================================
# EXAMPLE 3: Advanced usage with custom settings
# ============================================================================

def example_3_advanced_usage():
    """
    Example 3: Advanced usage with custom phase structures and high precision.
    
    This example shows:
    - Complex multi-component system
    - High-precision sampling settings
    - Custom line selection
    """
    
    print("\n" + "="*60)
    print("EXAMPLE 3: Advanced multi-component system")
    print("="*60)
    
    # Complex phase structure with multiple velocity components
    complex_phases = {
        'Component_A': ['SiIII_0', 'CII_1', 'HI_2'],
        'Component_B': ['CIV_3', 'SiIV_4', 'NV_5'], 
        'Component_C': ['SiII_6', 'FeII_7', 'MgI_8']
    }
    
    # Extended line list for comprehensive analysis
    extended_lines = [
        # Low ionization
        'SiII_1190', 'SiII_1193', 'SiII_1260',
        'MgII_2796', 'MgII_2803',
        'FeII_1144', 'FeII_1608', 'FeII_2374', 'FeII_2600',
        'MgI_2852',
        
        # Intermediate ionization
        'SiIII_1206', 'SIII_1012', 'SIII_1190',
        'CII_1334', 'CII_1036',
        
        # High ionization  
        'CIV_1548', 'CIV_1550',
        'SiIV_1393', 'SiIV_1402',
        'NV_1238', 'NV_1242',
        'OVI_1031', 'OVI_1037',
        
        # Hydrogen Lyman series
        'HI_1215', 'HI_1025', 'HI_972', 'HI_949'
    ]
    
    print(f"Phase structure: {len(complex_phases)} phases")
    for phase, components in complex_phases.items():
        print(f"  {phase}: {components}")
    
    print(f"\nLines to fit: {len(extended_lines)} transitions")
    
    # This would create a system with 9 components × 5 parameters = 45 parameters!
    n_components = sum(len(comps) for comps in complex_phases.values())
    n_parameters = n_components * 5
    print(f"Total parameters: {n_parameters}")
    
    # High-precision sampling settings for complex system
    high_precision_settings = {
        'min_num_live_points': 2000,  # More live points for complex system
        'dlogz': 0.1,                 # Tighter convergence
        'min_ess': 1000,              # Larger effective sample size
        'nsteps_factor': 3,           # More steps per iteration
    }
    
    print(f"\nHigh-precision sampling settings:")
    for key, value in high_precision_settings.items():
        print(f"  {key}: {value}")
    
    print(f"\nThis would be a very complex fit!")
    print(f"Estimated runtime: 12-24 hours with MPI parallelization")
    print(f"Recommended: Start with simpler 2-3 phase systems first")


# ============================================================================
# Helper function to create example configuration
# ============================================================================

def create_example_config(filename):
    """Create an example configuration file."""
    
    config_content = """
# Example CMBM Configuration File
# Update the paths below with your actual data locations

# ==============================================================================
# DATA PATHS - UPDATE THESE FOR YOUR DATA
# ==============================================================================
dataset_path: "/path/to/your/dataset.hdf5"
voigtfit_path: "/path/to/your/bestfit_result.pickle"  
mask_path: "/path/to/your/masks.pkl"
transition_library_path: "/path/to/your/atomdata.dat"

grid_paths:
  gridpath: "/path/to/your/cloudy/grids/"
  gridpath_thick: "/path/to/your/cloudy/grids/"

grid_z: "0.2"

# ==============================================================================
# SYSTEM DEFINITION
# ==============================================================================

# Two-phase system: low and high ionization
phases:
  Low_ionization:
    - "SiII_0"
    - "MgII_1" 
    - "FeII_2"
  High_ionization:
    - "CIV_0"
    - "SiIV_1"
    - "OVI_2"

# Lines to include in fit
use_lines:
  - "SiII_1190"
  - "MgII_2796" 
  - "FeII_1144"
  - "CIV_1548"
  - "SiIV_1393"
  - "OVI_1031"
  - "HI_1215"

# ==============================================================================
# SAMPLING SETTINGS
# ==============================================================================
min_num_live_points: 800
dlogz: 0.3
min_ess: 400

output_dir: "my_absorption_system"
resume: true
"""
    
    with open(filename, 'w') as f:
        f.write(config_content.strip())


# ============================================================================
# Plotting function for results
# ============================================================================

def plot_results(fitter, output_dir='plots'):
    """
    Create basic plots of the fitting results.
    
    Parameters:
    -----------
    fitter : cmbm.CMBMFitter
        Fitted CMBM instance
    output_dir : str
        Directory to save plots
    """
    
    if fitter.results is None:
        print("No results to plot - run fit first")
        return
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Get best-fit model
    model_profiles, data_profiles = fitter.get_best_fit_model()
    
    # Plot each line
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()
    
    plot_idx = 0
    for ion in model_profiles:
        for trans in model_profiles[ion]:
            if plot_idx >= len(axes):
                break
                
            ax = axes[plot_idx]
            
            # Get data
            vel, model = model_profiles[ion][trans]
            flux, error, mask = data_profiles[ion][trans]
            
            # Calculate wavelength for plotting
            l0 = fitter.species[ion][trans][0]
            wave = l0 * (1 + fitter.dataset.redshift) * (1 + vel / 299792.458)
            
            # Plot data points (only unmasked)
            masked_indices = np.where(mask)
            ax.errorbar(vel[masked_indices], flux[masked_indices], 
                       yerr=error[masked_indices], fmt='o', alpha=0.7, 
                       markersize=3, label='Data')
            
            # Plot model
            ax.plot(vel, model, 'r-', linewidth=2, label='Best fit')
            
            ax.set_xlabel('Velocity (km/s)')
            ax.set_ylabel('Normalized Flux')
            ax.set_title(f'{ion} {trans}')
            ax.legend()
            ax.grid(True, alpha=0.3)
            
            plot_idx += 1
    
    # Remove empty subplots
    for i in range(plot_idx, len(axes)):
        fig.delaxes(axes[i])
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/line_fits.png', dpi=150, bbox_inches='tight')
    plt.show()
    
    print(f"Plots saved to {output_dir}/")


# ============================================================================
# Main execution
# ============================================================================

if __name__ == "__main__":
    print("CMBM Basic Usage Examples")
    print("=" * 60)
    
    # Run examples
    try:
        # Example 1: Direct parameters (will need real data paths)
        # fitter1, results1 = example_1_direct_parameters()
        
        # Example 2: Configuration file  
        example_2_config_file()
        
        # Example 3: Advanced usage
        example_3_advanced_usage()
        
        print("\n" + "="*60)
        print("Examples completed!")
        print("="*60)
        print("\nTo run with real data:")
        print("1. Update the file paths in the examples")
        print("2. Prepare your VoigtFit dataset and results")
        print("3. Set up your Cloudy model grids")
        print("4. Run the fitting examples")
        print("\nFor questions, see: https://github.com/yourusername/cmbm")
        
    except Exception as e:
        print(f"Error running examples: {e}")
        print("This is expected if data paths are not set up correctly")
        print("Please update the file paths and try again")
