# Antarctic-Sea-ice-thickness
_Authors: Magata Mangatane , Marcello Vichi_

 This is a repository containing scripts to estimate Antarctic Sea-ice thickness with three algorithms. The algorithms have been used for the data analysis in Mangatane and Vichi (Submitted to The Cryosphere)

## Brief description
This collection of scripts reconstructs the Antarctic circumpolar seasonal Sea ice thickness with three innovative algorithms, namely. the improved One-Layer Method (OLMi), improved Buoyancy equation (BOC), and the freeboard differencing method (Diff method). The ice freeboard data used are from the ICESat-2 and CryoSat-2 satellites (see reference in the README.md)

# Workflow
The workflow is divided into the preprocessing of the data and the analyses. The raw data and preprocessed data are not provided due to their size but would be located in the data/ directory of this repository. The raw data are publicly available for download. 

 The scripts are provided in the scripts/ directory. This directory is divided into two directories, preprocessing/ and analysis/ to preprocess the raw data and to estimate the Sea-ice thickness, respectively. 

# Analysis

Drifter diagnostics (Drifter_diagnostics.py)
Drift velocity
Meander coefficient (MC)
Single-buoy absolute dispersion
Trajectory shape
Drifter response to wind and oceanic forcing (Drifter_diagnostics_metE5.py)
Wind factor, turning angle, correlation (Pearson and vector), and residual current
Power Spectral Analysis of the drifter and ERA5 wind velocity
Wavelet analysis and Butterworth High-pass filter (Drifter_wavelet_analysis.py)
NOTE: Read in Functions_wavelet.py file before
Absolute(single-particle) dispersion of the buoy cluster (Cluster_abs_dispersion.py)
Relative (two-particle) dispersion of the buoy cluster (Cluster_rel_dispersion.py)
NOTE: Read in Functions_rel_dispersion.py file before
Includes the computation of the deformation proxy and its spectral analysis

# Data source
