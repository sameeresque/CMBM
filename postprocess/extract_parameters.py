#!/usr/bin/env python
"""
Comprehensive post-processing and analysis of CMBM results.

This script:
1. Loads posterior samples and computes derived quantities
2. Computes quantile statistics (16th, 50th, 84th percentiles)
3. Computes ion-specific properties with total column densities
4. Generates corner plots for all clouds
5. Generates model comparison plots overlaying data
6. Saves all results in organized format

Each step can be skipped if already completed.
All plotting configuration (zgal, velocity range, line exclusions) is read from config file.
Figure size and layout are automatically determined based on number of lines.
Lines are sorted: HI first (by f*λ), then metals (alphabetically, by ionization, by f*λ).
Velocities are computed relative to zgal from config (or z_sys if not specified).
"""

import argparse
import json
import numpy as np
import pandas as pd
import pickle
from pathlib import Path
from collections import OrderedDict
import multiprocessing as mp
import sys

# Add parent directory to path to import cmbm
sys.path.insert(0, str(Path(__file__).parent.parent))
import cmbm

def parse_args():
    parser = argparse.ArgumentParser(description='Analyze CMBM results')
    parser.add_argument('--config', required=True, help='Config file used for fitting')
    parser.add_argument('--output_dir', required=True, help='CMBM output directory')
    parser.add_argument('--save_dir', required=True, help='Directory to save analysis results')
    parser.add_argument('--n_processes', type=int, default=12, help='Number of parallel processes')
    
    # Skip options
    parser.add_argument('--skip_derived', action='store_true', help='Skip computing derived quantities (load from file)')
    parser.add_argument('--skip_stats', action='store_true', help='Skip computing statistics (load from file)')
    parser.add_argument('--skip_corner', action='store_true', help='Skip corner plot generation')
    
    # Plot options
    parser.add_argument('--make_plots', action='store_true', help='Generate model comparison plots')
    
    return parser.parse_args()

def load_posterior(output_dir):
    """Load posterior samples from UltraNest output."""
    output_path = Path(output_dir)
    
    # Load equal weighted posterior
    post_file = output_path / 'chains' / 'equal_weighted_post.txt'
    if not post_file.exists():
        raise FileNotFoundError(f"Posterior file not found: {post_file}")
    
    post = pd.read_csv(post_file, sep='\s+').dropna(axis=1)
    print(f"  Loaded {len(post)} posterior samples")
    
    # Load weighted posterior for statistics
    weighted_file = output_path / 'chains' / 'weighted_post.txt'
    weighted = pd.read_csv(weighted_file, sep='\s+')
    weights = weighted['weight']
    logl = weighted['logl']
    
    # Load results for MLE
    results_file = output_path / 'info' / 'results.json'
    if not results_file.exists():
        raise FileNotFoundError(f"Results file not found: {results_file}")
    
    with open(results_file, 'r') as f:
        post_summary = json.load(f)
    
    return post, weights, logl, post_summary

def btherm(te, mass):
    """Calculate thermal broadening parameter."""
    return (0.0166289252 * 10**te / mass)**0.5

def binfo(b_turb, T, mass):
    """Calculate net b-value from turbulent and thermal components."""
    return np.sqrt(b_turb**2 + btherm(T, mass)**2)

def findV(zAbs, z):
    """Calculate velocity offset."""
    from astropy import constants as const
    v = const.c.cgs.value * 1e-5 * (z - zAbs) / (1.0 + zAbs)
    return v

# Processing functions for parallel computation
def process_pie_L(Z, nH, NHI, grids):
    """Compute cloud size for PIE model."""
    if NHI <= 16.0:
        try:
            N = grids['Nipiethin']['logNHtot']([Z, nH])[0] + (NHI - 14.0)
            return np.log10((10**N / 10**nH) / (3.086e21))
        except:
            return np.nan
    else:
        try:
            N = grids['Nipiethick']['logNHtot']([Z, nH, NHI])[0]
            return np.log10((10**N / 10**nH) / (3.086e21))
        except:
            return np.nan

def process_pie_T(Z, nH, NHI, grids):
    """Compute temperature for PIE model."""
    if NHI <= 16.0:
        try:
            return grids['Nipiethin']['logT']([Z, nH])[0]
        except:
            return np.nan
    else:
        try:
            return grids['Nipiethick']['logT']([Z, nH, NHI])[0]
        except:
            return np.nan

def process_pie_NH(Z, nH, NHI, grids):
    """Compute total hydrogen column density for PIE model."""
    if NHI <= 16.0:
        try:
            return grids['Nipiethin']['logNHtot']([Z, nH])[0] + (NHI - 14.0)
        except:
            return np.nan
    else:
        try:
            return grids['Nipiethick']['logNHtot']([Z, nH, NHI])[0]
        except:
            return np.nan

def process_pie_ion(Z, nH, NHI, ion, grids):
    """Compute ion column density for PIE model."""
    if NHI <= 16.0:
        try:
            return grids['Nipiethin'][ion]([Z, nH])[0] + (NHI - 14.0)
        except:
            return np.nan
    else:
        try:
            return grids['Nipiethick'][ion]([Z, nH, NHI])[0]
        except:
            return np.nan

def compute_derived_quantities(post, fitter, n_processes=12):
    """Compute derived quantities using multiprocessing."""
    phases = fitter.config['phases']
    interest_ions = list(set([line.split('_')[0] for line in fitter.config['use_lines']]))
    interest_ions = [ion for ion in interest_ions if ion in fitter.config['cloudy_ions']]
    
    # Flatten phases to get list of clouds
    clouds = []
    for phase_name, components in phases.items():
        clouds.extend(components)
    
    print(f"  Computing derived quantities for {len(clouds)} clouds...")
    
    # Get zgal from config for velocity computation
    zgal = fitter.config.get('zgal', None)
    if zgal is None:
        zgal = fitter.dataset.redshift
        print(f"    Using z_sys = {zgal:.6f} for velocity reference (zgal not specified in config)")
    else:
        print(f"    Using zgal = {zgal:.6f} for velocity reference")
    
    with mp.Pool(processes=n_processes) as pool:
        for cloud in clouds:
            print(f"    Processing cloud: {cloud}")
            
            # Prepare input grids
            pie_grid = [list(post[f'{cloud}_Z']), list(post[f'{cloud}_nH']), list(post[f'{cloud}_NHI'])]
            
            # Compute size
            L_args = [(Z, nH, NHI, fitter.grids) for Z, nH, NHI in zip(*pie_grid)]
            post[f'{cloud}_L'] = pool.starmap(process_pie_L, L_args)
            
            # Compute temperature
            T_args = [(Z, nH, NHI, fitter.grids) for Z, nH, NHI in zip(*pie_grid)]
            post[f'{cloud}_T'] = pool.starmap(process_pie_T, T_args)
            
            # Compute total NH
            NH_args = [(Z, nH, NHI, fitter.grids) for Z, nH, NHI in zip(*pie_grid)]
            post[f'{cloud}_NH'] = pool.starmap(process_pie_NH, NH_args)
            
            # Compute ion column densities
            for ion in interest_ions:
                ion_args = [(Z, nH, NHI, ion, fitter.grids) for Z, nH, NHI in zip(*pie_grid)]
                post[f'{cloud}_{ion}'] = pool.starmap(process_pie_ion, ion_args)
            
            # Compute velocity relative to zgal from config (or z_sys if not specified)
            V_args = [(zgal, z) for z in post[f'{cloud}_z']]
            post[f'{cloud}_V'] = pool.starmap(findV, V_args)
            
            # Compute thermal b-value for HI
            mass_HI = fitter.species['HI'][list(fitter.species['HI'].keys())[0]][2]
            btherm_args = [(T, mass_HI) for T in post[f'{cloud}_T']]
            post[f'{cloud}_btherm_HI'] = pool.starmap(btherm, btherm_args)
            
            # Compute net b-value for HI
            post[f'{cloud}_bnet_HI'] = np.sqrt(
                np.array(post[f'{cloud}_btherm_HI'])**2 + 
                np.array(post[f'{cloud}_bturb'])**2
            )
    
    # Compute total ion column densities across all clouds
    print("  Computing total ion column densities...")
    for ion in interest_ions:
        post[f'{ion}_total'] = post.apply(
            lambda x: np.log10(sum([10**x[f'{cloud}_{ion}'] for cloud in clouds])), 
            axis=1
        )
    
    return post

def _quantile(x, q, weights=None):
    """
    Compute (weighted) quantiles from an input set of samples.
    
    Parameters
    ----------
    x : array_like
        Input samples
    q : array_like
        Quantiles to compute from [0., 1.]
    weights : array_like, optional
        Sample weights
        
    Returns
    -------
    quantiles : array_like
        The weighted sample quantiles computed at q
    """
    x = np.atleast_1d(x)
    q = np.atleast_1d(q)
    
    if np.any(q < 0.0) or np.any(q > 1.0):
        raise ValueError("Quantiles must be between 0. and 1.")
    
    if weights is None:
        return np.percentile(x, list(100.0 * q))
    else:
        weights = np.atleast_1d(weights)
        if len(x) != len(weights):
            raise ValueError("Dimension mismatch: len(weights) != len(x).")
        idx = np.argsort(x)
        sw = weights[idx]
        cdf = np.cumsum(sw)[:-1]
        cdf /= cdf[-1]
        cdf = np.append(0, cdf)
        quantiles = np.interp(q, cdf, x[idx]).tolist()
        return quantiles

def sigma_arr(sigma):
    """Get quantile array for given sigma level."""
    if sigma == 1:
        return [(50.0-68.27/2.0)/100.0, 0.5, (50.0+68.27/2.0)/100.0]
    elif sigma == 2:
        return [(50.0-95.45/2.0)/100.0, 0.5, (50.0+95.45/2.0)/100.0]
    else:  # sigma == 3
        return [(50.0-99.73/2.0)/100.0, 0.5, (50.0+99.73/2.0)/100.0]

def compute_statistics(post, weights, fitter):
    """Compute quantile statistics for all parameters."""
    phases = fitter.config['phases']
    clouds = []
    for phase_name, components in phases.items():
        clouds.extend(components)
    
    interest_ions = list(set([line.split('_')[0] for line in fitter.config['use_lines']]))
    interest_ions = [ion for ion in interest_ions if ion in fitter.config['cloudy_ions']]
    
    print("  Computing quantile statistics...")
    
    # Sort by log-likelihood for weighted quantiles
    ilocs = np.argsort(weights.index)
    
    # Structure: pars[cloud][param][sigma] = (median, -err, +err)
    pars_basic = OrderedDict()
    for cloud in clouds:
        pars_basic[cloud] = OrderedDict()
        for param in ['z', 'V', 'Z', 'nH', 'T', 'NHI', 'NH', 'L', 'bturb', 'btherm_HI', 'bnet_HI']:
            pars_basic[cloud][param] = OrderedDict()
            for sigma in [1, 2, 3]:
                ql, qm, qh = _quantile(
                    post[f'{cloud}_{param}'].iloc[ilocs].values, 
                    sigma_arr(sigma), 
                    weights=weights.iloc[ilocs].values
                )
                q_minus, q_plus = qm - ql, qh - qm
                pars_basic[cloud][param][sigma] = (qm, q_minus, q_plus)
    
    # Ion-specific properties
    pars_ions = OrderedDict()
    for cloud in clouds:
        pars_ions[cloud] = OrderedDict()
        for ion in interest_ions:
            pars_ions[cloud][ion] = OrderedDict()
            
            # Velocity (same for all ions in a cloud)
            ql, qm, qh = _quantile(
                post[f'{cloud}_V'].iloc[ilocs].values,
                sigma_arr(1),
                weights=weights.iloc[ilocs].values
            )
            pars_ions[cloud][ion]['V'] = (qm, qm - ql, qh - qm)
            
            # Redshift
            ql, qm, qh = _quantile(
                post[f'{cloud}_z'].iloc[ilocs].values,
                sigma_arr(1),
                weights=weights.iloc[ilocs].values
            )
            pars_ions[cloud][ion]['z'] = (qm, qm - ql, qh - qm)
            
            # b-value (ion-specific)
            if ion == 'HI':
                b_col = f'{cloud}_bnet_HI'
            else:
                # Compute b for this specific ion
                mass = fitter.species[ion][list(fitter.species[ion].keys())[0]][2]
                post[f'{cloud}_{ion}_b'] = post.apply(
                    lambda x: binfo(x[f'{cloud}_bturb'], x[f'{cloud}_T'], mass),
                    axis=1
                )
                b_col = f'{cloud}_{ion}_b'
            
            ql, qm, qh = _quantile(
                post[b_col].iloc[ilocs].values,
                sigma_arr(1),
                weights=weights.iloc[ilocs].values
            )
            pars_ions[cloud][ion]['b'] = (qm, qm - ql, qh - qm)
            
            # Column density
            if ion == 'HI':
                N_col = f'{cloud}_NHI'
            else:
                N_col = f'{cloud}_{ion}'
            
            ql, qm, qh = _quantile(
                post[N_col].iloc[ilocs].values,
                sigma_arr(1),
                weights=weights.iloc[ilocs].values
            )
            pars_ions[cloud][ion]['logN'] = (qm, qm - ql, qh - qm)
            
            # Total column density (across all clouds)
            ql, qm, qh = _quantile(
                post[f'{ion}_total'].iloc[ilocs].values,
                sigma_arr(1),
                weights=weights.iloc[ilocs].values
            )
            pars_ions[cloud][ion]['Total'] = (qm, qm - ql, qh - qm)
    
    return pars_basic, pars_ions

def make_corner_plots(post, weights, fitter, save_dir):
    """Generate corner plots for all clouds."""
    try:
        import corner
        import matplotlib.pyplot as plt
    except ImportError:
        print("  Warning: corner or matplotlib not installed. Skipping corner plots.")
        print("  Install with: pip install corner matplotlib")
        return
    
    phases = fitter.config['phases']
    clouds = []
    for phase_name, components in phases.items():
        clouds.extend(components)
    
    print("  Generating corner plots...")
    
    plots_dir = Path(save_dir) / 'corner_plots'
    plots_dir.mkdir(exist_ok=True)
    
    for cloud in clouds:
        print(f"    Creating corner plot for {cloud}...")
        
        # Select columns for this cloud
        cols = [f'{cloud}_Z', f'{cloud}_nH', f'{cloud}_NHI', f'{cloud}_bturb', 
                f'{cloud}_T', f'{cloud}_L', f'{cloud}_NH', f'{cloud}_z']
        
        # Filter to existing columns
        cols = [c for c in cols if c in post.columns]
        
        use_df = post[cols]
        
        # Create corner plot
        fig = corner.corner(
            use_df,
            weights=weights,
            labels=cols,
            title_fmt='.2f',
            show_titles=True,
            use_math_text=True
        )
        
        plot_file = plots_dir / f'{cloud}_corner.pdf'
        plt.savefig(plot_file, bbox_inches='tight')
        plt.close(fig)
        print(f"      Saved: {plot_file}")

def generate_model_profiles(pars_basic, fitter):
    """
    Generate median model profiles for each cloud and combined profile.
    
    Parameters
    ----------
    pars_basic : dict
        Parameter statistics from compute_statistics
    fitter : CMBMFitter
        Fitter object with loaded data
        
    Returns
    -------
    model_profile : dict
        Individual cloud profiles
    model_profile_comb : dict
        Combined profile across all clouds
    infos : dict
        Cloud information dictionaries
    clouds : list
        List of cloud names
    """
    from VoigtFit import voigt
    from scipy.signal import fftconvolve
    from scipy.signal.windows import gaussian
    from astropy import constants as const
    
    def eval_profile(x, info, z_sys, ion, trans, species, kernel, sampling=3, kernel_nsub=1):
        """Evaluate Voigt profile for single cloud."""
        if isinstance(kernel, float):
            dx = np.mean(np.diff(x))
            xmin = np.log10(x.min() - 50*dx)
            xmax = np.log10(x.max() + 50*dx)
            N = int(sampling * len(x))
            profile_wl = np.logspace(xmin, xmax, N)
            pxs = np.diff(profile_wl)[0] / profile_wl[0] * const.c.cgs.value*1e-5
            kernel = kernel / pxs / 2.35482
        elif isinstance(kernel, np.ndarray):
            N = int(kernel_nsub * len(x))
            if kernel_nsub > 1:
                profile_wl = np.linspace(x.min(), x.max(), N)
            else:
                profile_wl = x.copy()
        else:
            raise TypeError(f"Invalid type of kernel: {type(kernel)}")
        
        tau = np.zeros_like(profile_wl)
        l0, f, gam = species[ion][trans][0], species[ion][trans][1], species[ion][trans][3]
        l_center = l0 * (z_sys + 1.)
        
        velcen = []
        for component in info:
            z = info[component][ion][0]
            if x.min() < l0*(z+1.0) < x.max() or ion == 'HI':
                b = info[component][ion][1]
                N_col = info[component][ion][2]
                tau += voigt.Voigt(profile_wl, l0, f, float(N_col), 1.e5*float(b), gam, z=z)
                if b != 0.00001:
                    velcen.append(-findV(z, z_sys))
        
        profile = np.exp(-tau)
        
        if isinstance(kernel, float):
            LSF = gaussian(10*int(kernel) + 1, kernel)
            LSF = LSF/LSF.sum()
            profile_broad = fftconvolve(profile, LSF, 'same')
            profile_obs = np.interp(x, profile_wl, profile_broad)
        else:
            profile_broad = voigt.convolve(profile, kernel)
            if kernel_nsub > 1:
                profile_obs = np.interp(x, profile_wl, profile_broad)
            else:
                profile_obs = profile_broad
        
        vel = const.c.cgs.value*1e-5*(x-l_center)/(l_center)
        return vel, profile_obs, velcen, tau
    
    def eval_comb_profile(x, infos, z_sys, ion, trans, species, kernel, sampling=3, kernel_nsub=1):
        """Evaluate combined Voigt profile for all clouds."""
        if isinstance(kernel, float):
            dx = np.mean(np.diff(x))
            xmin = np.log10(x.min() - 50*dx)
            xmax = np.log10(x.max() + 50*dx)
            N = int(sampling * len(x))
            profile_wl = np.logspace(xmin, xmax, N)
            pxs = np.diff(profile_wl)[0] / profile_wl[0] * const.c.cgs.value*1e-5
            kernel = kernel / pxs / 2.35482
        elif isinstance(kernel, np.ndarray):
            N = int(kernel_nsub * len(x))
            if kernel_nsub > 1:
                profile_wl = np.linspace(x.min(), x.max(), N)
            else:
                profile_wl = x.copy()
        else:
            raise TypeError(f"Invalid type of kernel: {type(kernel)}")
        
        l0, f, gam = species[ion][trans][0], species[ion][trans][1], species[ion][trans][3]
        l_center = l0 * (z_sys + 1.)
        
        tau_list = []
        for info in infos:
            tau = np.zeros_like(profile_wl)
            for component in info:
                z = info[component][ion][0]
                if x.min() < l0*(z+1.0) < x.max() or ion == 'HI':
                    b = info[component][ion][1]
                    N_col = info[component][ion][2]
                    tau += voigt.Voigt(profile_wl, l0, f, float(N_col), 1.e5*float(b), gam, z=z)
            tau_list.append(tau)
        
        net_tau = [sum(lista) for lista in zip(*tau_list)]
        profile = np.exp(-np.asarray(net_tau))
        
        if isinstance(kernel, float):
            LSF = gaussian(10*int(kernel) + 1, kernel)
            LSF = LSF/LSF.sum()
            profile_broad = fftconvolve(profile, LSF, 'same')
            profile_obs = np.interp(x, profile_wl, profile_broad)
        else:
            profile_broad = voigt.convolve(profile, kernel)
            if kernel_nsub > 1:
                profile_obs = np.interp(x, profile_wl, profile_broad)
            else:
                profile_obs = profile_broad
        
        vel = const.c.cgs.value*1e-5*(x-l_center)/(l_center)
        return vel, profile_obs
    
    # Get clouds list
    phases = fitter.config['phases']
    clouds = []
    for phase_name, components in phases.items():
        clouds.extend(components)
    
    interest_ions = list(set([line.split('_')[0] for line in fitter.config['use_lines']]))
    interest_ions = [ion for ion in interest_ions if ion in fitter.config['cloudy_ions']]
    
    # Build info dict for each cloud using median values
    infos = {}
    for num, cloud in enumerate(clouds):
        infos[num] = OrderedDict()
        infos[num]['0'] = OrderedDict()
        
        # Get median values
        z_med = pars_basic[cloud]['z'][1][0]  # median from 1-sigma
        Z_med = pars_basic[cloud]['Z'][1][0]
        nH_med = pars_basic[cloud]['nH'][1][0]
        NHI_med = pars_basic[cloud]['NHI'][1][0]
        bturb_med = pars_basic[cloud]['bturb'][1][0]
        T_med = pars_basic[cloud]['T'][1][0]
        
        # Compute properties for each ion
        for ion in interest_ions:
            mass = fitter.species[ion][list(fitter.species[ion].keys())[0]][2]
            
            # Get column density
            if NHI_med <= 16.0:
                N = fitter.grids['Nipiethin'][ion]([Z_med, nH_med])[0] + (NHI_med - 14.0)
            else:
                N = fitter.grids['Nipiethick'][ion]([Z_med, nH_med, NHI_med])[0]
            
            # Compute b-value
            b = np.sqrt(bturb_med**2 + btherm(T_med, mass)**2)
            
            # Store as (z, b, N_linear)
            infos[num]['0'][ion] = (z_med, b, 10**N)
    
    # Generate model profiles for individual clouds
    z_sys = fitter.dataset.redshift
    model_profile = {}
    
    print("  Generating model profiles for individual clouds...")
    for i in infos:
        model_profile[i] = OrderedDict()
        for ion in interest_ions:
            model_profile[i][ion] = OrderedDict()
            for trans in fitter.species[ion]:
                line_id = f'{ion}_{trans}'
                if line_id in fitter.dataset.all_lines:
                    line_data = fitter.dataset.find_line(line_id)[0]
                    x = line_data.wl
                    kernel = line_data.kernel
                    kernel_nsub = line_data.kernel_nsub
                    
                    model_profile[i][ion][trans] = eval_profile(
                        x, infos[i], z_sys, ion, trans, fitter.species, 
                        kernel, sampling=3, kernel_nsub=kernel_nsub
                    )
    
    # Generate combined profile
    print("  Generating combined model profile...")
    model_profile_comb = OrderedDict()
    for ion in interest_ions:
        model_profile_comb[ion] = OrderedDict()
        for trans in fitter.species[ion]:
            line_id = f'{ion}_{trans}'
            if line_id in fitter.dataset.all_lines:
                line_data = fitter.dataset.find_line(line_id)[0]
                x = line_data.wl
                kernel = line_data.kernel
                kernel_nsub = line_data.kernel_nsub
                
                model_profile_comb[ion][trans] = eval_comb_profile(
                    x, tuple([infos[i] for i in infos]), z_sys, ion, trans, 
                    fitter.species, kernel, sampling=3, kernel_nsub=kernel_nsub
                )
    
    return model_profile, model_profile_comb, infos, clouds

def plot_model_comparison(model_profile, model_profile_comb, infos, clouds, post, pars_basic, fitter, save_dir):
    """
    Create comprehensive plot comparing model to data with parameter distributions.
    
    All configuration (zgal, velocity range, line exclusions) is read from config file.
    Figure size and layout are automatically determined based on number of lines.
    Lines are sorted: HI first (by f*λ), then metals (alphabetically, by ionization, by f*λ).
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec
        from matplotlib.ticker import MultipleLocator
        import random
    except ImportError:
        print("  Warning: matplotlib not installed. Skipping model comparison plot.")
        return
    
    def reject_outliers(data, m=3):
        """Remove outliers beyond m standard deviations."""
        return data[abs(data - np.mean(data)) < m * np.std(data)]
    
    def sort_ions_and_transitions_for_plotting(ions, species):
        """
        Sort ions and their transitions for optimal plotting order.
        - HI lines come first, sorted by f*lambda descending (strongest first)
        - Metal lines sorted alphabetically by element, then by ionization state
        - Within each ion, transitions sorted by f*lambda descending
        
        Parameters
        ----------
        ions : list
            List of ion names (e.g., ['CII', 'HI', 'SiIII'])
        species : dict
            Species dictionary with transition information
            species[ion][trans] = (wavelength, f_value, mass, gamma)
            
        Returns
        -------
        sorted_ion_trans_list : list of tuples
            List of (ion, transition) tuples with HI first, then metals alphabetically
        """
        import roman
        
        def parse_ion(ion):
            """Parse ion name into element and ionization state."""
            # Find where roman numeral starts
            for i in range(1, len(ion)+1):
                try:
                    ionization = roman.fromRoman(ion[i:])
                    element = ion[:i]
                    return element, ionization
                except:
                    continue
            return ion, 0
        
        hi_pairs = []
        metal_ions = {}
        
        # Separate HI from metal ions
        for ion in ions:
            if ion in species:
                if ion == 'HI':
                    # Collect all HI transitions with f*lambda
                    for trans in species[ion]:
                        wavelength = species[ion][trans][0]
                        f_value = species[ion][trans][1]
                        f_lambda = f_value * wavelength
                        hi_pairs.append((ion, trans, f_lambda))
                else:
                    # Group metal ions for sorting
                    element, ionization = parse_ion(ion)
                    if element not in metal_ions:
                        metal_ions[element] = []
                    
                    # Collect all transitions for this ion with f*lambda
                    ion_trans = []
                    for trans in species[ion]:
                        wavelength = species[ion][trans][0]
                        f_value = species[ion][trans][1]
                        f_lambda = f_value * wavelength
                        ion_trans.append((ion, trans, f_lambda))
                    
                    # Sort transitions by f*lambda descending (strongest first)
                    ion_trans.sort(key=lambda x: x[2], reverse=True)
                    
                    metal_ions[element].append((ionization, ion_trans))
        
        # Sort HI by f*lambda descending (strongest first: Lyα, Lyβ, Lyγ, etc.)
        hi_pairs.sort(key=lambda x: x[2], reverse=True)
        
        # Build result: HI first
        result = [(ion, trans) for ion, trans, _ in hi_pairs]
        
        # Then metals: alphabetically by element, then by ionization state
        for element in sorted(metal_ions.keys()):
            # Sort by ionization state (I, II, III, IV, etc.)
            metal_ions[element].sort(key=lambda x: x[0])
            
            # Add all transitions for this element in order
            for ionization, ion_trans_list in metal_ions[element]:
                for ion, trans, _ in ion_trans_list:
                    result.append((ion, trans))
        
        return result
    
    interest_ions = list(set([line.split('_')[0] for line in fitter.config['use_lines']]))
    interest_ions = [ion for ion in interest_ions if ion in fitter.config['cloudy_ions']]
    
    z_sys = fitter.dataset.redshift
    
    # Get zgal from config, default to z_sys if not specified
    zgal = fitter.config.get('zgal', None)
    if zgal is None:
        zgal = z_sys
        print(f"  Using z_sys = {z_sys:.4f} for velocity reference")
    else:
        print(f"  Using zgal = {zgal:.4f} for velocity reference (z_sys = {z_sys:.4f})")
    
    # Get plotting configuration from config file
    plotting_config = fitter.config.get('plotting', {})
    deactivate_lines = plotting_config.get('deactivate_lines', [])
    masked = plotting_config.get('masked_lines', [])
    velocity_range = plotting_config.get('velocity_range', [-30, 390])
    
    vmin, vmax = velocity_range
    print(f"  Velocity range: [{vmin}, {vmax}] km/s")
    
    if deactivate_lines:
        print(f"  Excluding {len(deactivate_lines)} lines from plot")
    if masked:
        print(f"  Showing {len(masked)} lines as fully masked")
    
    # Get all ion-transition pairs sorted: HI first, then by alphabetical/ionization
    sorted_ion_trans = sort_ions_and_transitions_for_plotting(interest_ions, fitter.species)
    
    # Filter to only lines in dataset and not deactivated
    lines_to_plot = []
    for ion, trans in sorted_ion_trans:
        line_id = f'{ion}_{trans}'
        if line_id not in deactivate_lines and line_id in fitter.dataset.all_lines:
            lines_to_plot.append((ion, trans, line_id))
    
    n_lines_to_plot = len(lines_to_plot)
    
    if n_lines_to_plot == 0:
        print("  Warning: No lines to plot after applying filters!")
        return
    
    print(f"  Plotting {n_lines_to_plot} absorption lines")
    print(f"  Line ordering: HI first (by f*λ), then metals (alphabetically, by ionization, by f*λ)")
    
    # Print HI lines if present
    hi_lines = [(ion, trans, line_id) for ion, trans, line_id in lines_to_plot if ion == 'HI']
    if hi_lines:
        print(f"  HI lines ({len(hi_lines)}) - sorted by f*λ:")
        for i, (ion, trans, line_id) in enumerate(hi_lines, 1):
            wavelength = fitter.species[ion][trans][0]
            f_lambda = fitter.species[ion][trans][0] * fitter.species[ion][trans][1]
            print(f"    {i}. {line_id}: λ = {wavelength:.1f} Å, f*λ = {f_lambda:.2f}")
    
    # Print metal lines
    non_hi_lines = [(ion, trans, line_id) for ion, trans, line_id in lines_to_plot if ion != 'HI']
    if non_hi_lines:
        print(f"  Metal lines ({len(non_hi_lines)}) - alphabetical, by ionization, by f*λ:")
        for i, (ion, trans, line_id) in enumerate(non_hi_lines[:15], 1):
            f_lambda = fitter.species[ion][trans][0] * fitter.species[ion][trans][1]
            print(f"    {i}. {line_id}: f*λ = {f_lambda:.2f}")
        if len(non_hi_lines) > 15:
            print(f"    ... and {len(non_hi_lines) - 15} more lines")
    
    # Calculate grid dimensions
    rows = 6  # Fixed: number of parameters to display
    n_cols_lines = int(np.ceil(n_lines_to_plot / rows))  # Columns needed for lines
    n_cols_total = n_cols_lines + 1  # Add 1 column for parameter distributions
    
    # Calculate figure size dynamically
    # Width: 3 inches per line column + 4 inches for parameter column
    # Height: 2.3 inches per row (for 6 rows)
    fig_width = 3.0 * n_cols_lines + 4.0
    fig_height = 2.3 * rows
    
    print(f"  Figure layout: {rows} rows × {n_cols_total} columns ({n_cols_lines} line columns + 1 param column)")
    print(f"  Figure size: {fig_width:.1f} × {fig_height:.1f} inches")
    
    # Generate colors for clouds
    number_of_colors = len(infos)
    random.seed(0)
    colora = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
              for i in range(number_of_colors)]
    
    plt.rc('font', family='serif')
    plt.rcParams["axes.edgecolor"] = "black"
    
    fig = plt.figure(figsize=(fig_width, fig_height))
    gs = gridspec.GridSpec(rows, n_cols_total)
    gs.update(wspace=0.0, hspace=0.0)
    k = 0
    l = 0
    
    # Plot absorption lines in sorted order
    ax0 = None
    for ion, trans, line_id in lines_to_plot:
        # Create subplot
        if k == 0 and l == 0:
            ax = fig.add_subplot(gs[k, l])
            ax0 = ax
        else:
            ax = fig.add_subplot(gs[k, l], sharex=ax0, sharey=ax0)
        
        # Hide labels appropriately
        if l != 0:
            plt.setp(ax.get_yticklabels(), visible=False)
        if k != rows-1:  # Not in last row
            if l != n_cols_lines-1 or k < rows//2:  # Not in last column or upper half
                plt.setp(ax.get_xticklabels(), visible=False)
        
        # Plot individual cloud models
        for numb, i in enumerate(infos):
            vel, flux, velcen, tau = model_profile[i][ion][trans]
            ax.plot(vel + findV(zgal, z_sys), flux, 
                   color=colora[i], label=f'{numb}', 
                   linewidth=1.5, linestyle='-', zorder=1)
            for m in velcen:
                ax.axvline(x=m + findV(zgal, z_sys), ymin=0.8, ymax=0.95,
                          color=colora[i], linewidth=1.5, linestyle='-')
        
        # Plot combined model
        try:
            vel_comb, flux_comb = model_profile_comb[ion][trans]
            ax.plot(vel_comb + findV(zgal, z_sys), flux_comb, 
                   'black', linewidth=1.75, linestyle='-', zorder=3)
        except:
            print(f"  {line_id}: no combined model")
        
        # Annotate
        ax.annotate(f'{ion} {trans}', xy=(0.03, 0.10), 
                   xycoords='axes fraction', color='navy', fontsize=14)
        
        # Plot data
        try:
            line_data = fitter.dataset.find_line(line_id)[0]
            vel_data = vel_comb + findV(zgal, z_sys)
            
            ax.errorbar(vel_data, line_data.flux, yerr=line_data.err,
                       color='gray', fmt='.', ls='none')
            
            mask_data = line_data.mask
            if line_id not in masked:
                ax.fill_between(vel_data, 0, 1, where=(mask_data == 0),
                               facecolor='grey', alpha=0.5)
            else:
                ax.fill_between(vel_data, 0, 1, 
                               facecolor='grey', alpha=0.5)
        except:
            print(f"  {line_id}: no data")
        
        # Formatting - use velocity range from config
        ax.set_xlim([vmin, vmax])
        ax.set_ylim([-0.1, 1.3])
        ax.tick_params(axis='both', direction='in', which='major', 
                      length=10, width=2, labelsize=16)
        ax.tick_params(axis='both', direction='in', which='minor',
                      length=5, width=1, labelsize=16)
        ax.yaxis.set_major_locator(MultipleLocator(0.5))
        ax.yaxis.set_minor_locator(MultipleLocator(0.1))
        
        # Dynamic tick spacing based on velocity range
        vel_span = vmax - vmin
        if vel_span <= 100:
            major_tick = 20
            minor_tick = 5
        elif vel_span <= 300:
            major_tick = 50
            minor_tick = 10
        elif vel_span <= 600:
            major_tick = 100
            minor_tick = 25
        else:
            major_tick = 150
            minor_tick = 50
        
        ax.xaxis.set_major_locator(MultipleLocator(major_tick))
        ax.xaxis.set_minor_locator(MultipleLocator(minor_tick))
        
        k = k + 1
        if k == rows:
            k = 0
            l = l + 1
    
    # Add parameter distribution panels
    print("  Adding parameter distributions...")
    ax_params = []
    for i in range(rows):
        # Parameter column is the last column
        ax_param = fig.add_subplot(gs[i, n_cols_lines:n_cols_total])
        ax_params.append(ax_param)
        if i > 0:
            ax_param.sharex(ax_params[0])
    
    # Format parameter axes - use velocity range from config
    for axis in ax_params:
        axis.yaxis.set_label_position("right")
        axis.yaxis.tick_right()
        axis.set_xlim([vmin, vmax])
        axis.tick_params(axis='both', direction='in', which='major',
                        length=10, width=2, labelsize=16)
        axis.tick_params(axis='both', direction='in', which='minor',
                        length=5, width=1, labelsize=16)
        
        # Use same dynamic tick spacing
        axis.xaxis.set_major_locator(MultipleLocator(major_tick))
        axis.xaxis.set_minor_locator(MultipleLocator(minor_tick))
        
        if axis != ax_params[-1]:
            plt.setp(axis.get_xticklabels(), visible=False)
    
    # Extract posterior distributions
    postZ = post[post.columns[post.columns.str.endswith('_Z')]]
    postnH = post[post.columns[post.columns.str.endswith('_nH')]]
    postNHI = post[post.columns[post.columns.str.endswith('_NHI')]]
    postT = post[post.columns[post.columns.str.endswith('_T')]]
    postL = post[post.columns[post.columns.str.endswith('_L')]]
    postbturb = post[post.columns[post.columns.str.endswith('_bturb')]]
    
    # Plot violin plots for each parameter
    param_data = [
        (ax_params[0], postZ, 'Z', r'log $\frac{Z}{Z_{\odot}}$'),
        (ax_params[1], postnH, 'nH', r'log $\frac{n_{H}}{cm^{-3}}$'),
        (ax_params[2], postNHI, 'NHI', r'log $\frac{N(\rm HI)}{cm^{-2}}$'),
        (ax_params[3], postT, 'T', r'log $\frac{T}{K}$'),
        (ax_params[4], postL, 'L', r'log $\frac{L}{kpc}$'),
        (ax_params[5], postbturb, 'bturb', r'b$_{nt}$ [km/s]'),
    ]
    
    for ax, post_param, param_name, ylabel in param_data:
        # Get velocities for x-axis - use velocities from pars_basic which are now relative to zgal
        velocities = [pars_basic[cloud]['V'][1][0] for cloud in clouds]
        
        # Dynamic violin width based on velocity range
        violin_width = (vmax - vmin) / 10
        
        # 1-sigma violin plots
        parts1 = ax.violinplot(
            [reject_outliers(post_param[f'{cloud}_{param_name}'].dropna().values, 1) 
             for cloud in clouds],
            velocities,
            showmeans=False, showmedians=False, showextrema=False, widths=violin_width
        )
        
        j = 0
        for pc in parts1['bodies']:
            pc.set_facecolor(colora[j])
            pc.set_edgecolor('black')
            pc.set_alpha(0.9)
            j = j + 1
        
        # 3-sigma violin plots  
        parts2 = ax.violinplot(
            [reject_outliers(post_param[f'{cloud}_{param_name}'].dropna().values, 3)
             for cloud in clouds],
            velocities,
            showmeans=False, showmedians=False, showextrema=False, widths=violin_width
        )
        
        j = 0
        for pc in parts2['bodies']:
            pc.set_facecolor('none')
            pc.set_edgecolor('black')
            pc.set_linestyle('dashed')
            pc.set_alpha(0.7)
            j = j + 1
        
        # Plot median points
        medians = [pars_basic[cloud][param_name][1][0] for cloud in clouds]
        ax.scatter(velocities, medians, s=50, marker='*', color='black')
        
        ax.set_ylabel(ylabel, size=18)
    
    # Add axis labels
    fig.text(0.5, 0.08, r'Relative Velocity [km s$^{-1}$]', 
            ha='center', va='center', fontsize=18)
    fig.text(0.08, 0.5, 'Normalized Flux', 
            ha='center', va='center', rotation='vertical', fontsize=18)
    
    # Save figure
    plot_file = Path(save_dir) / 'model_comparison.pdf'
    plt.savefig(plot_file, bbox_inches='tight')
    plt.close(fig)
    
    print(f"  ✓ Model comparison plot saved: {plot_file}")

def main():
    args = parse_args()
    
    print("="*60)
    print("CMBM Results Analysis")
    print("="*60)
    
    print("\n[1/7] Loading configuration and data...")
    fitter = cmbm.CMBMFitter(args.config)
    fitter.load_data()
    
    save_path = Path(args.save_dir)
    save_path.mkdir(parents=True, exist_ok=True)
    
    print("\n[2/7] Loading posterior samples...")
    post, weights, logl, post_summary = load_posterior(args.output_dir)
    
    # Step 3: Compute derived quantities
    if args.skip_derived:
        print("\n[3/7] Loading derived quantities from file...")
        posterior_file = save_path / 'posterior_full.pkl'
        if not posterior_file.exists():
            raise FileNotFoundError(f"Cannot skip - file not found: {posterior_file}")
        with open(posterior_file, 'rb') as f:
            post = pickle.load(f)
        print(f"  ✓ Loaded from: {posterior_file}")
    else:
        print("\n[3/7] Computing derived quantities...")
        post = compute_derived_quantities(post, fitter, args.n_processes)
    
    # Step 4: Compute statistics
    if args.skip_stats:
        print("\n[4/7] Loading statistics from file...")
        stats_file = save_path / 'parameter_statistics.pkl'
        ions_file = save_path / 'ion_properties.pkl'
        
        if not stats_file.exists() or not ions_file.exists():
            raise FileNotFoundError(f"Cannot skip - files not found in {save_path}")
        
        with open(stats_file, 'rb') as f:
            pars_basic = pickle.load(f)
        with open(ions_file, 'rb') as f:
            pars_ions = pickle.load(f)
        print(f"  ✓ Loaded statistics from: {save_path}")
    else:
        print("\n[4/7] Computing statistics...")
        pars_basic, pars_ions = compute_statistics(post, weights, fitter)
    
    # Step 5: Save results (only if not skipped)
    if not args.skip_derived or not args.skip_stats:
        print("\n[5/7] Saving results...")
        
        if not args.skip_derived:
            # Save full posterior with derived quantities
            posterior_file = save_path / 'posterior_full.pkl'
            with open(posterior_file, 'wb') as f:
                pickle.dump(post, f, protocol=2)
            print(f"  ✓ Full posterior: {posterior_file}")
            
            # Save merged dataframe
            df_merge = pd.concat([weights.reset_index(drop=True), logl.reset_index(drop=True), post.reset_index(drop=True)], axis=1)
            merged_file = save_path / 'posterior_merged.txt'
            np.savetxt(merged_file, df_merge.values, fmt='%e')
            print(f"  ✓ Merged posterior: {merged_file}")
        
        if not args.skip_stats:
            # Save basic parameter statistics
            stats_file = save_path / 'parameter_statistics.pkl'
            with open(stats_file, 'wb') as f:
                pickle.dump(pars_basic, f, protocol=2)
            print(f"  ✓ Basic statistics: {stats_file}")
            
            # Save ion-specific properties
            ions_file = save_path / 'ion_properties.pkl'
            with open(ions_file, 'wb') as f:
                pickle.dump(pars_ions, f, protocol=2)
            print(f"  ✓ Ion properties: {ions_file}")
    else:
        print("\n[5/7] Skipping save (all steps were loaded from file)")
    
    # Step 6: Corner plots
    if not args.skip_corner:
        print("\n[6/7] Creating corner plots...")
        make_corner_plots(post, weights, fitter, save_path)
    else:
        print("\n[6/7] Skipping corner plots (--skip_corner flag set)")
    
    # Step 7: Model comparison plots
    if args.make_plots:
        print("\n[7/7] Generating model comparison plots...")
        model_profile, model_profile_comb, infos, clouds = generate_model_profiles(pars_basic, fitter)
        plot_model_comparison(model_profile, model_profile_comb, infos, clouds, post, pars_basic, fitter, save_path)
    else:
        print("\n[7/7] Skipping model comparison plots (use --make_plots to enable)")
    
    print("\n" + "="*60)
    print("✓ Analysis complete!")
    print(f"  Results saved to: {save_path}")
    print("="*60)

if __name__ == '__main__':
    main()