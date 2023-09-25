# Antarctic-Sea-ice-thickness
_Authors: Magata Mangatane , Marcello Vichi_
This is a repository containing scripts to estimate Antarctic Sea-ice thickness with three algorithms. The algorithms have been used for the data analysis in Mangatane and Vichi (Submitted to The Cryosphere)

#Brief description
This collection of scripts processes drifter data from the Antarctic MIZ, in the Atlantic and Indian Ocean sectors. It has been tested on drifter data collected during the SA Agulhas II expeditions, and public data from the AWI catalogue (see reference in the README.md). These scripts compute the buoys’ drift diagnostics, the drift response to ERA5 atmospheric forcing, the spectral analysis and wavelet analysis of the drift velocity, the cluster absolute and relative dispersion statistics, and the spectral analysis of the deformation proxy.  

Atmospheric data will need to be downloaded separately from ERA5 reanalysis (hourly on single levels from 1940 to present). Variables should include mean sea level pressure, wind velocity components (u, v), and 2-m air temperature for the period and region of the buoys’ drift.

#Workflow
The workflow is divided into the preprocessing of the data and the analyses. The raw data and preprocessed data are both provided in the data/ directory of this repository.

The scripts and functions are provided in the scripts/ directory.

#Analysis

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

#Data source
