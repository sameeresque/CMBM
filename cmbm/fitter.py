"""
CMBM: Cloud-by-cloud Multiphase Bayesian ionization Modeling

This module contains the main CMBMFitter class that wraps your original 
spectral fitting functionality into a clean interface.
"""

import time
import sys
import os
import numpy as np
import pickle
import roman
from collections import OrderedDict
import scipy
import scipy.signal
import multiprocessing as mp
from scipy.signal import fftconvolve
from scipy.signal.windows import gaussian
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import LinearNDInterpolator
from astropy import constants as const
from numpy import convolve
import VoigtFit
from VoigtFit import voigt
from VoigtFit import show_transitions
import pandas as pd
import re
import math
import json
import itertools
import ultranest
import ultranest.stepsampler
import mpi4py
import yaml
from pathlib import Path

ALPHA_ELEMENT_MAP = {
    'Mgalpha': 'Mg',
    'Nalpha':  'N',
    'Calpha':  'C',
}


# Skip warnings when we add the -np.inf log likelihood value
np.seterr(invalid='ignore')


class CMBMFitter:
    """
    Cloud-by-cloud Multiphase Bayesian ionization Modeling (CMBM) fitter.
    
    This class provides a clean interface to fit quasar absorption line systems
    using Cloudy photoionization models and nested sampling.
    
    Parameters:
    -----------
    config_file : str, optional
        Path to YAML configuration file
    **kwargs : dict
        Direct parameter specification (overrides config file values)
        
    Key parameters:
    - dataset_path : str, path to VoigtFit dataset (.hdf5)
    - voigtfit_path : str, path to VoigtFit results (.pickle)
    - mask_path : str, path to velocity masks (.pkl)  
    - grid_paths : dict, paths to Cloudy model grids
    - transition_library_path : str, path to atomic data file
    - phases : dict, phase structure (e.g., {'Ph1': ['SiIII_0', 'HI_1']})
    - use_lines : list, spectral lines to fit
    - output_dir : str, output directory for results
    
    Examples:
    ---------
    >>> # Method 1: Using config file
    >>> fitter = CMBMFitter('my_config.yaml')
    >>> results = fitter.run_fit()
    
    >>> # Method 2: Direct parameters  
    >>> fitter = CMBMFitter(
    ...     dataset_path='data.hdf5',
    ...     voigtfit_path='bestfit.pickle',
    ...     phases={'Ph1': ['SiIII_0', 'HI_1']},
    ...     use_lines=['CII_1334', 'HI_1215', 'MgII_2796']
    ... )
    >>> results = fitter.run_fit()
    """
    
    def __init__(self, config_file=None, **kwargs):
        """Initialize the CMBM fitter."""
        # Load configuration
        self.config = self._load_config(config_file, **kwargs)
        
        # Initialize data containers
        self.dataset = None
        self.voigt_data = None
        self.species = None
        self.grids = None
        self.masks = None
        
        # Results
        self.sampler = None
        self.results = None
        
        print("CMBMFitter initialized")
    
    def _load_config(self, config_file, **kwargs):
        """Load configuration from user-supplied file."""
        if not config_file:
            raise ValueError("A config file must be supplied.")
        
        with open(config_file, 'r') as f:
            config = yaml.safe_load(f)
        
        # Override with direct parameters if supplied
        for key, value in kwargs.items():
            config[key] = value
        
        return config
    
    
    def load_data(self):
        """Load all required data files."""
        self._load_grids()
        print("Loading data files...")
        
        # Load dataset
        self.dataset = VoigtFit.load_dataset(self.config['dataset_path'])
        self._fix_dataset_errors()
        print(f"Loaded dataset: {len(self.dataset.all_lines)} lines")
        
        # Load VoigtFit results
        with open(self.config['voigtfit_path'], 'rb') as f:
            self.voigt_data = pickle.load(f, encoding='bytes')
        print("Loaded VoigtFit results")
        
        # Load transition library and create species dictionary
        transition_library = pd.read_table(
            self.config['transition_library_path'], 
            sep=r'\s+', header=None, comment="#"
        )
        transition_library = np.asarray(transition_library)
        
        all_ions = [i.split('_')[0] for i in self.dataset.all_lines]
        interest_ions = list(set(all_ions))
        interest_ions = [ion for ion in interest_ions if ion in self.grids['Nipiethin']]

        self.species = self._create_species_dict(interest_ions, transition_library)
        print(f"Loaded species data for {len(interest_ions)} ions")
        
        # Load grids
        self._load_grids()
        print("Loaded Cloudy grids")
        
        # Load and apply masks
        self.masks = pd.read_pickle(self.config['mask_path'])

        # Add custom masks from config if specified
        custom_masks = self.config.get('masks', {}).get('custom_masks', {})
        for line_id, regions in custom_masks.items():
            if line_id not in self.masks:
                self.masks[line_id] = []
            else:
                # Convert numpy array to list if needed
                if isinstance(self.masks[line_id], np.ndarray):
                    self.masks[line_id] = self.masks[line_id].tolist()
            
            # Now extend with custom masks
            self.masks[line_id].extend(regions)

        self.masks = self._update_velocity_masks(self.masks)
        self._apply_masks()
        print("Applied velocity masks")
        
        print("All data loaded successfully")
    
    def _fix_dataset_errors(self):
        """Fix negative errors in dataset."""
        for line_name in self.dataset.all_lines:
            line_obj = self.dataset.find_line(line_name)[0]
            df = pd.DataFrame({'A': line_obj.err})
            positive_errors = line_obj.err[line_obj.err > 0]
            if len(positive_errors) > 0:
                median_error = np.nanmedian(positive_errors)
                df[df < 0] = median_error
                line_obj.err = df['A'].values
    
    def _create_species_dict(self, plot_ions, transition_library, choose=None):
        """Create species dictionary"""
        species = {}
        strongest_trans = OrderedDict()
        second_trans = OrderedDict()
        
        for search_trans in plot_ions:
            species[search_trans] = OrderedDict()
            max_oscillator = 0.
            m1 = m2 = float('-inf')
            oldwav1 = oldwav2 = " "
            
            for transition in transition_library:
                atom, tranwav = transition[0].split("_")
                current_oscillator = transition[3]
                gamma = transition[4]
                atomic_mass = transition[5]
                
                if atom == search_trans:
                    species[atom][tranwav] = (
                        transition[2], current_oscillator, atomic_mass, gamma
                    )
                    
                    # Find strongest and 2nd strongest transitions
                    if transition[3] > m2:
                        if transition[3] >= m1:
                            m1, m2 = transition[3], m1
                            oldwav1, oldwav2 = tranwav, oldwav1
                            strongest_trans[search_trans] = tranwav
                            second_trans[search_trans] = oldwav1
                        else:
                            m2 = transition[3]
                            oldwav2 = tranwav
                            second_trans[search_trans] = tranwav
            
            if search_trans in second_trans:
                second_trans[search_trans] = oldwav2
        
        if choose is not None:
            for var in choose:
                if var in species:
                    for transition in species[var].copy().keys():
                        if transition not in choose[var]:
                            species[var].pop(transition)
        
        return species
    
    def _load_grids(self):
        with open(self.config['grids']['thin'], 'rb') as f:
            Nipiethin = pickle.load(f)
        with open(self.config['grids']['thick'], 'rb') as f:
            Nipiethick = pickle.load(f)
        
        self.grids = {
            'Nipiethin':  Nipiethin,
            'Nipiethick': Nipiethick,
        }
    
    def _update_velocity_masks(self, mask_dict, lowlim=-160, uplim=190, maxvel=500):
        """Update velocity masks"""
        # Get velocity limits from config
        vel_limits = self.config.get('masks', {}).get('velocity_limits', {})
        lowlim = vel_limits.get('lowlim', -160)
        uplim = vel_limits.get('uplim', 190)
        maxvel = vel_limits.get('maxvel', 500)

        updated_dict = {}
        
        for key, masks in mask_dict.items():
            if len(masks) == 0:
                updated_dict[key] = [[-maxvel, lowlim], [uplim, maxvel]]
                continue
            
            mask_array = np.array(masks)
            if len(mask_array.shape) == 1:
                mask_array = mask_array.reshape(-1, 2)
            
            new_masks = []
            for mask in masks:
                if isinstance(mask, np.ndarray):
                    mask = mask.tolist()
                
                if mask[0] < lowlim:
                    new_masks.append([-maxvel, lowlim])
                    if mask[1] > lowlim:
                        if mask[1] <= uplim:
                            new_masks.append([lowlim, mask[1]])
                        else:
                            new_masks.append([lowlim, uplim])
                            new_masks.append([uplim, maxvel])
                elif mask[0] < uplim:
                    if mask[1] <= uplim:
                        new_masks.append([mask[0], mask[1]])
                    else:
                        new_masks.append([mask[0], uplim])
                        new_masks.append([uplim, maxvel])
                else:
                    new_masks.append([uplim, maxvel])
            
            if not any(mask[0] <= -maxvel for mask in new_masks):
                new_masks.append([-maxvel, lowlim])
            if not any(mask[1] >= maxvel for mask in new_masks):
                new_masks.append([uplim, maxvel])
            
            new_masks.sort(key=lambda x: x[0])
            merged = []
            for mask in new_masks:
                if not merged or merged[-1][1] < mask[0]:
                    merged.append(mask)
                else:
                    merged[-1][1] = max(merged[-1][1], mask[1])
            
            updated_dict[key] = merged
        
        return updated_dict
    
    def _apply_masks(self):
        """Apply masks to dataset."""
        interest_ions = list(set([line.split('_')[0] for line in self.config['use_lines']]))
        interest_ions = [ion for ion in interest_ions if ion in self.grids['Nipiethin']]
        # First mask all lines
        for ion in interest_ions:
            if ion not in self.species:
                continue
            for trans in self.species[ion]:
                line_id = f'{ion}_{trans}'
                if line_id in self.dataset.all_lines:
                    self.dataset.find_line(line_id)[-1].mask[:] = True
        
        # Then unmask specified regions
        for line, regions in self.masks.items():
            if (line in self.dataset.all_lines and 
                line.split('_')[0] in interest_ions):
                
                line_obj = self.dataset.find_line(line)[0]
                ion, trans = line.split('_')
                
                if ion in self.species and trans in self.species[ion]:
                    l0 = self.species[ion][trans][0]
                    velr = self._getVel(line_obj.wl, l0, self.dataset.redshift)
                    
                    for region in regions:
                        loc = np.where((velr >= region[0]) & (velr <= region[1]))
                        line_obj.mask[loc] = False
    
    def _getVel(self, wave, l0, zabs):
        """Calculate velocity"""
        wave_center = l0 * (zabs + 1.)
        vel = const.c.cgs.value * 1e-5 * (wave - wave_center) / wave_center
        return vel

    def _speciesname(self,ion):
        """Determine the species"""
        romannum=['I','II','III','IV','V','VI','VII','VIII','IX','X']
        if ion[1:] not in romannum:
            ionname=ion[0:2]
            trans=ion[2:]
        else:
            ionname=ion[0:1]
            trans=ion[1:]
        return ionname,roman.fromRoman(trans)

    
    def _btherm(self, te, mass):
        """Calculate thermal broadening"""
        return (0.0166289252 * 10**te / mass)**0.5
    
    def _Nxpiethin(self, Z_val, Hden_val, NHI_val, ion):
        """Interpolate thin grid"""
        try:
            N = self.grids['Nipiethin'][ion]([Z_val, Hden_val])[0] + (NHI_val - 14.0)
            T = self.grids['Nipiethin']['logT']([Z_val, Hden_val])[0]
            return N, T
        except:
            return np.nan, np.nan
    
    def _Nxpiethick(self, Z_val, Hden_val, NHI_val, ion):
        """Interpolate thick grid"""
        try:
            N = self.grids['Nipiethick'][ion]([Z_val, Hden_val, NHI_val])[0]
            T = self.grids['Nipiethick']['logT']([Z_val, Hden_val, NHI_val])[0]
            return N, T
        except:
            return np.nan, np.nan
    
    def _eval_comb_profile(self, x, infos, z_sys, ion, trans, species, kernel, 
                          sampling=3, kernel_nsub=1):
        """Evaluate combined profile"""
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
            raise TypeError(f"Invalid type of `kernel`: {type(kernel)}")

        l0, f, gam = species[ion][trans][0], species[ion][trans][1], species[ion][trans][3]
        l_center = l0*(z_sys + 1.)

        tau_list = []
        for info in infos:
            tau = np.zeros_like(profile_wl)
            for component in info:
                z = info[component][ion][0]
                if x.min() < l0*(z+1.0) < x.max():
                    b = info[component][ion][1]
                    N = info[component][ion][2]
                    tau += voigt.Voigt(profile_wl, l0, f, float(N), 1.e5*float(b), gam, z=z)
                elif ion == 'HI':
                    b = info[component][ion][1]
                    N = info[component][ion][2]
                    tau += voigt.Voigt(profile_wl, l0, f, float(N), 1.e5*float(b), gam, z=z)
                else:
                    continue
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


    def _apply_alpha(self, N, item, alpha_vals):
        """Apply abundance offset if item's element matches an alpha parameter."""
        if not alpha_vals:
            return N
        elem = self._speciesname(item)[0]   # e.g. 'Mg', 'N', 'C', 'Si', ...
        for alpha_name, alpha_val in alpha_vals.items():
            if ALPHA_ELEMENT_MAP.get(alpha_name) == elem:
                return N + alpha_val
        return N


    
    def _cmodel(self, pars):
        """Cloud model evaluation"""
        phases = self.config['phases']
        use_lines = self.config['use_lines']
        interest_ions = list(set([line.split('_')[0] for line in use_lines]))
        interest_ions = [ion for ion in interest_ions if ion in self.grids['Nipiethin']]
        alpha_params = self.config.get('alpha_params', {})
        allphases = []
        param_idx = 0
        
        for phase_name, components in phases.items():
            n_clouds = len(components)
            optim_ions = []
            for j in range(n_clouds):
                optim_ion = components[j].split('_')[0]
                optim_ion_loc = components[j].split('_')[1]
                optim_ions.append((optim_ion, optim_ion_loc))
            
            allinfo = OrderedDict()
            for cloud in range(n_clouds):
                s_cloud = str(cloud)
                vcloud = str(optim_ions[cloud][1])
                component = components[cloud]
                allinfo[s_cloud] = OrderedDict()
                
                # Extract parameters for this cloud
                Z = pars[param_idx]
                nH = pars[param_idx + 1]
                NHI = pars[param_idx + 2]
                bturb = pars[param_idx + 3]
                z = pars[param_idx + 4]
                param_idx += 5


                # --- Alpha parameters for this cloud (0 or more) ---
                cloud_alphas = alpha_params.get(component, [])
                alpha_vals = {}
                for alpha_name in cloud_alphas:
                    alpha_vals[alpha_name] = pars[param_idx]
                    param_idx += 1


                
                for item in interest_ions:
                    weight = self.species[item][list(self.species[item].keys())[0]][2]
                    
                    # Use appropriate grid based on NHI
                    if NHI <= 16.0:
                        N, te = self._Nxpiethin(Z, nH, NHI, item)
                    else:
                        N, te = self._Nxpiethick(Z, nH, NHI, item)
                    
                    b = ((bturb)**2 + (self._btherm(te, weight))**2)**0.5
                    # Apply alpha offset if element matches

                    N_final = self._apply_alpha(N, item, alpha_vals)

                    allinfo[s_cloud][item] = z, b, 10**N_final


            
            allphases.append(allinfo)

        model_profile_comb = OrderedDict()  
        data_profile = OrderedDict()

        for ion in interest_ions:
            if ion not in self.species:
                continue
            model_profile_comb[ion] = OrderedDict()  
            data_profile[ion] = OrderedDict()
            
            for trans in self.species[ion]:
                line_id = f'{ion}_{trans}'
                if line_id in use_lines and line_id in self.dataset.all_lines:
                    line_obj = self.dataset.find_line(line_id)[-1]
                    x = line_obj.wl
                    kernel = line_obj.kernel
                    
                    vel, profile = self._eval_comb_profile(
                        x, tuple(allphases), self.dataset.redshift, ion, trans, 
                        self.species, kernel, sampling=3, 
                        kernel_nsub=getattr(line_obj, 'kernel_nsub', 1)
                    )
                    
                    model_profile_comb[ion][trans] = (vel, profile)
                    data_profile[ion][trans] = (line_obj.flux, line_obj.err, line_obj.mask)
        
        return model_profile_comb, data_profile
    
    def _generate_parameters(self):
        """Generate parameter names and bounds."""
        parameters = []
        minval = []
        maxval = []

        alpha_params = self.config.get('alpha_params', {})
        alpha_bounds = self.config.get('alpha_bounds', {})  # e.g. {'Mgalpha': [-1.0, 1.0]}

        prior_settings = self.config.get('prior_settings', {})
        z_sigma     = prior_settings.get('z_sigma', 2)
        bturb_sigma = prior_settings.get('bturb_sigma', 2)
        NHI_min = prior_settings['NHI_min']
        NHI_max = prior_settings['NHI_max']
        
        for phase_name, components in self.config['phases'].items():
            for component in components:
                ion, comp_id = component.split('_')
                base_name = f"{ion}_{comp_id}"
                
                # Add parameter names
                param_names = [f"{base_name}_Z", f"{base_name}_nH", f"{base_name}_NHI", 
                             f"{base_name}_bturb", f"{base_name}_z"]
                parameters.extend(param_names)
                
                # Add bounds based on your original logic
                for param_name in param_names:
                    param_type = param_name.split('_')[-1]
                    
                    if param_type == 'z':
                        z_val = self.voigt_data[ion][comp_id]['z'][0]
                        z_err = self.voigt_data[ion][comp_id]['z'][1]
                        minval.append(z_val - z_sigma * z_err)
                        maxval.append(z_val + z_sigma * z_err)
                    
                    elif param_type == 'Z':
                        Z_axis = self.grids['Nipiethin']['logT'].grid[0]
                        minval.append(float(Z_axis.min()))
                        maxval.append(float(Z_axis.max()))

                    elif param_type == 'nH':
                        nH_axis = self.grids['Nipiethin']['logT'].grid[1]
                        minval.append(float(nH_axis.min()))
                        maxval.append(float(nH_axis.max()))
                    
                    elif param_type == 'NHI':
                        minval.append(NHI_min)
                        maxval.append(NHI_max)
                    
                    elif param_type == 'bturb':
                        b_val = self.voigt_data[ion][comp_id]['b'][0]
                        b_err = self.voigt_data[ion][comp_id]['b'][1]
                        minval.append(0)
                        maxval.append(b_val + bturb_sigma * b_err)


                for alpha in alpha_params.get(component, []):
                    parameters.append(f"{base_name}_{alpha}")
                    bounds = alpha_bounds.get(alpha, [-1.0, 1.0])
                    minval.append(bounds[0])
                    maxval.append(bounds[1])                

        
        return parameters, minval, maxval
    
    def _prior(self, pars, minval, maxval):
        """Prior transform"""
        params = pars.copy()
        for i in range(len(pars)):
            params[i] = minval[i] + (maxval[i] - minval[i]) * pars[i]
        return params
    
    def _ln_likelihood(self, pars):
        """Log likelihood"""
        try:
            model_data = self._cmodel(pars)
            resid = []
            
            for specie in model_data[0]:
                for trans in model_data[0][specie]:
                    model = model_data[0][specie][trans][1]
                    ydata = model_data[1][specie][trans][0]
                    ysigma = model_data[1][specie][trans][1]
                    loc_mask = np.where(model_data[1][specie][trans][2])
                    resi = (ydata[loc_mask] - model[loc_mask]) / ysigma[loc_mask]
                    resid.append(resi)
            
            residue = list(itertools.chain.from_iterable(resid))
            lnprobval = -0.5 * np.dot(residue, residue)

            if math.isinf(lnprobval) or math.isnan(lnprobval):
                return -1e90
            else:
                return lnprobval
        except Exception as e:
            print(f"Error in likelihood calculation: {e}")
            return -1e90
    
    def run_fit(self, output_dir=None, **kwargs):
        """
        Run the nested sampling fit.
        
        Parameters:
        -----------
        output_dir : str, optional
            Output directory for results (overrides config)
        **kwargs : dict
            Sampling parameters to override:
            - min_num_live_points : int
            - dlogz : float
            - min_ess : int
            - nsteps_factor : int
        
        Returns:
        --------
        dict
            UltraNest results dictionary
            
        Examples:
        ---------
        >>> # Basic usage
        >>> results = fitter.run_fit()
        
        >>> # Custom output directory and sampling
        >>> results = fitter.run_fit(
        ...     output_dir='my_results',
        ...     min_num_live_points=800,
        ...     dlogz=0.3
        ... )
        """
        if self.dataset is None:
            raise ValueError("Data must be loaded first. Call load_data().")
        
        # Generate parameters and bounds
        parameters, minval, maxval = self._generate_parameters()
        ndim = len(parameters)
        
        print(f"Starting nested sampling with {ndim} parameters:")
        for i, param in enumerate(parameters):
            print(f"  {param}: [{minval[i]:.3f}, {maxval[i]:.3f}]")
        
        # Setup output directory
        if output_dir is None:
            output_dir = self.config['output_dir']
        
        # Sampling parameters
        min_num_live_points = kwargs.get('min_num_live_points', self.config['min_num_live_points'])
        dlogz = kwargs.get('dlogz', self.config['dlogz']) + 0.1 * ndim
        min_ess = kwargs.get('min_ess', self.config['min_ess'])
        nsteps_factor = kwargs.get('nsteps_factor', self.config['nsteps_factor'])
        
        print(f"Sampling parameters:")
        print(f"  Live points: {min_num_live_points}")
        print(f"  Target dlogz: {dlogz:.2f}")
        print(f"  Min ESS: {min_ess}")
        print(f"  Output: {output_dir}")
        
        # Create prior and likelihood functions
        def prior_transform(pars):
            return self._prior(pars, minval, maxval)
        
        def log_likelihood(pars):
            return self._ln_likelihood(pars)
        
        # Create sampler
        self.sampler = ultranest.ReactiveNestedSampler(
            parameters, 
            log_likelihood, 
            prior_transform,
            log_dir=output_dir, 
            resume=self.config['resume']
        )
        
        # Setup step sampler
        self.sampler.stepsampler = ultranest.stepsampler.SliceSampler(
            nsteps=nsteps_factor * ndim,
            generate_direction=ultranest.stepsampler.generate_mixture_random_direction,
        )
        
        # Run sampling
        print("Starting nested sampling...")
        self.results = self.sampler.run(
            min_num_live_points=min_num_live_points,
            dlogz=dlogz,
            min_ess=min_ess,
            update_interval_volume_fraction=self.config['update_interval_volume_fraction'],
            max_num_improvement_loops=self.config['max_num_improvement_loops']
        )
        
        print("Fitting completed!")
        print(f"Log evidence: {self.results['logz']:.2f} ± {self.results['logzerr']:.2f}")
        print(f"Samples: {len(self.results['samples'])}")


        self.sampler.plot_run()
        self.sampler.plot_trace()
        self.sampler.plot_corner()
        
        return self.results
    
    def get_results(self):
        """
        Get organized results with parameter summaries.
        
        Returns:
        --------
        dict
            Dictionary with results and parameter statistics
        """
        if self.results is None:
            raise ValueError("No results available. Run fit first.")
        
        parameters, _, _ = self._generate_parameters()
        samples = self.results['samples']
        
        # Calculate parameter summaries
        param_summaries = {}
        for i, param in enumerate(parameters):
            param_samples = samples[:, i]
            param_summaries[param] = {
                'median': np.median(param_samples),
                'mean': np.mean(param_samples),
                'std': np.std(param_samples),
                'percentiles': {
                    '16th': np.percentile(param_samples, 16),
                    '84th': np.percentile(param_samples, 84)
                }
            }
        
        return {
            'sampler_results': self.results,
            'parameters': parameters,
            'log_evidence': self.results['logz'],
            'log_evidence_error': self.results['logzerr'],
            'parameter_summaries': param_summaries,
            'n_samples': len(samples)
        }
    
    def get_best_fit_model(self):
        """
        Get the best-fit model profiles.
        
        Returns:
        --------
        tuple
            (model_profiles, data_profiles) for the best-fit parameters
        """
        if self.results is None:
            raise ValueError("No results available. Run fit first.")
        
        # Find best-fit parameters (highest likelihood)
        max_like_idx = np.argmax(self.results['logl'])
        best_params = self.results['samples'][max_like_idx]
        
        # Evaluate model at best-fit parameters
        model_profiles, data_profiles = self._cmodel(best_params)
        
        return model_profiles, data_profiles
    
    def summary(self):
        """Print a summary of the fitter setup and results."""
        print("\n" + "="*60)
        print("CMBM Fitter Summary")
        print("="*60)
        
        # Configuration
        print("\nConfiguration:")
        print(f"  Dataset: {self.config.get('dataset_path', 'Not set')}")
        print(f"  VoigtFit results: {self.config.get('voigtfit_path', 'Not set')}")
        print(f"  Mask file: {self.config.get('mask_path', 'Not set')}")
        print(f"  Atomic data: {self.config.get('transition_library_path', 'Not set')}")
        
        # Grid paths
        grid_paths = self.config.get('grid_paths', {})
        print(f"  Grid (thin): {grid_paths.get('gridpath', 'Not set')}")
        print(f"  Grid (thick): {grid_paths.get('gridpath_thick', 'Not set')}")
        print(f"  Grid redshift: {self.config.get('grid_z', 'Not set')}")
        
        print(f"  System redshift: {self.dataset.redshift:.6f}" if self.dataset else "  System redshift: Not loaded")
        
        phases = self.config['phases']
        total_clouds = sum(len(components) for components in phases.values())



        print(f"  Phases: {total_clouds}")
        for phase, components in phases.items():
            print(f"    {phase}: {components}")
        
        use_lines = self.config['use_lines']
        print(f"  Lines to fit: {len(use_lines)}")
        print(f"  Line list: {', '.join(use_lines)}")

        # Compute parameters
        print(f"\nModel Parameters:")
        parameters, _, _ = self._generate_parameters() if self.dataset else ([], [], [])
        ndim = len(parameters)
        print(f"  Number of parameters (ndim): {ndim}")
        if ndim > 0:
            print(f"  Parameters: {', '.join(parameters)}")
        
        # Nested sampling settings
        print(f"\nNested Sampling Settings:")
        print(f"  Live points: {self.config.get('min_num_live_points', 'Not set')}")
        print(f"  Target dlogz: {self.config.get('dlogz', 'Not set')}")
        print(f"  Minimum ESS: {self.config.get('min_ess', 'Not set')}")
        nsteps_factor = self.config.get('nsteps_factor', 'Not set')
        if ndim > 0 and nsteps_factor != 'Not set':
            print(f"  Step sampler steps: {nsteps_factor} × {ndim} = {nsteps_factor * ndim}")
        else:
            print(f"  Step sampler factor: {nsteps_factor}")
        print(f"  Output directory: {self.config.get('output_dir', 'Not set')}")
        print(f"  Resume: {self.config.get('resume', 'Not set')}")
        
        # System information
        print(f"\nSystem Information:")
        n_cpu = mp.cpu_count()
        print(f"  Available CPUs: {n_cpu}")
        try:
            from mpi4py import MPI
            comm = MPI.COMM_WORLD
            n_mpi = comm.Get_size()
            rank = comm.Get_rank()
            print(f"  MPI processes: {n_mpi}")
            if n_mpi > 1:
                print(f"  Current MPI rank: {rank}")
        except:
            print(f"  MPI processes: 1 (MPI not initialized)")


        # SLURM Job Settings
        slurm_config = self.config.get('slurm', {})
        if slurm_config:
            print(f"\nSLURM Job Settings:")
            job_name = slurm_config.get('job_name', 'Not set')
            partition = slurm_config.get('partition', 'Not set')
            container = slurm_config.get('container', 'Not set')
            nodes = slurm_config.get('nodes', 'Not set')
            ntasks = slurm_config.get('ntasks_per_node', 'Not set')
            mem = slurm_config.get('mem_per_cpu', 'Not set')
            walltime = slurm_config.get('time', 'Not set')
            
            print(f"  Job name: {job_name}")
            print(f"  Partition: {partition}")
            print(f"  Container: {container}")
            if nodes != 'Not set' and ntasks != 'Not set':
                print(f"  Total MPI tasks: {nodes} × {ntasks} = {nodes * ntasks}")
            print(f"  Nodes: {nodes}")
            print(f"  Tasks per node: {ntasks}")
            print(f"  Memory per CPU: {mem}")
            print(f"  Wall time: {walltime}")

        
        # Results
        if self.results is not None:
            print(f"\nResults:")
            print(f"  Log evidence: {self.results['logz']:.2f} ± {self.results['logzerr']:.2f}")
            print(f"  Samples: {len(self.results['samples'])}")
            
            # Parameter summary
            results = self.get_results()
            print(f"\nBest-fit parameters (median ± 1σ):")
            for param, stats in results['parameter_summaries'].items():
                lower = stats['median'] - stats['percentiles']['16th']
                upper = stats['percentiles']['84th'] - stats['median']
                print(f"  {param:15s}: {stats['median']:8.4f} +{upper:6.4f} -{lower:6.4f}")
        else:
            print(f"\nNo results yet - run fit first")
        
        print("="*60 + "\n")
    
    def save_config(self, filename):
        """
        Save current configuration to YAML file.
        
        Parameters:
        -----------
        filename : str
            Output filename for configuration
        """
        import yaml
        with open(filename, 'w') as f:
            yaml.dump(self.config, f, default_flow_style=False, indent=2)
        print(f"Configuration saved to {filename}")


