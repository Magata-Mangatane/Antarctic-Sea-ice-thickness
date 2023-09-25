# Antarctic-Sea-ice-thickness: a repository of algorithms to reconstruct Antarctic Sea-ice thickness with satellite observations
_Authors: Magata J Mangatane , Marcello Vichi_

 This is a repository containing scripts to estimate Antarctic Sea-ice thickness with three algorithms. The algorithms have been used for the data analysis in Mangatane and Vichi (Submitted to The Cryosphere)

## Brief description
This collection of scripts reconstructs the Antarctic circumpolar seasonal Sea-ice thickness (SIT) with three innovative algorithms, namely, the improved One-Layer Method (OLMi), the improved Buoyancy equation (BOC), and the freeboard differencing method (Diff method). The ice freeboard data used are from the ICESat-2 and CryoSat-2 satellites covering the period 2019-2022 and 2019 only, respectively (see reference in the README.md)

## Workflow
The workflow is divided into the preprocessing of the data and the analyses. The raw data and preprocessed data are not provided due to their size but would be located in the data/ directory of this repository. The raw data are publicly available for download. 

 The scripts are provided in the scripts/ directory. This directory is divided into two sub-directories, preprocessing/ and analysis/ to preprocess the raw data and to estimate the Sea-ice thickness, respectively. 

## Preprocessing of data
The preprocessing scripts are found in scripts/preprocessing/

 This is the sequence of operations:
 1. Extract ICESat-2 freeboard data from HDF files to NetCDF file format [extract_IS2_freeboards.py](scripts/analysis/preprocessing/extract_IS2_freeboards.py)
 2. Extract AMSR Sea-ice concentration data from He5 to NetCDF file format [extract_amsr_sic_to_nc.py](scripts/analysis/preprocessing/extract_amsr_sic_to_nc.py)
 3. Grid CryoSat-2 data to a daily 25 km polar stereographic grid [grid_CS2_freeboards.py](scripts/analysis/preprocessing/grid_CS2_freeboards.py)
 4. Weigh CryoSat-2 data with AMSR data [weigh_cs2_freeboards.py](scripts/analysis/preprocessing/weigh_cs2_freeboards.py)

## Analysis
The scripts to estimate Sea-ice thickness are found in scripts/analysis/
 
 These are the different algorithms used:
 1. Estimation of SIT with the improved One-Layer Method is done here [OLMi_sit_estimation.py](scripts/analysis/OLMi_sit_estimation.py)
 2. Estimation of SIT with the improved Buoyancy equation is done here [BOC_sit_estimation.py](scripts/analysis/BOC_sit_estimation.py)
 3. Estimation of SIT with the freeboard differencing method is done here [diff_method_sit_estimation.py](scripts/analysis/diff_method_sit_estimation.py)

## Data source
The raw data required to complete the reconstruction of Antarctic Sea-ice thickness demonstrated here are publicly available as follows:
* Kwok, R., A. A. Petty, G. Cunningham, T. Markus, D. Hancock, A. Ivanoff, J. Wimert, M. Bagnardi, N. Kurtz, and  the ICESat-2 Science Team. (2021). ATLAS/ICESat-2 L3A Sea Ice Height, Version 5 [Data Set]. Boulder, Colorado USA. NASA National Snow and Ice Data Center Distributed Active Archive Center. https://doi.org/10.5067/ATLAS/ATL07.005. Date Accessed 09-25-2023.
* European Space Agency, 2019, L2 SAR Precise Orbit. Baseline D. https://doi.org/10.5270/CR2-53hztdl
* Markus, T., J. C. Comiso, and W. N. Meier. (2018). AMSR-E/AMSR2 Unified L3 Daily 25 km Brightness Temperatures & Sea Ice Concentration Polar Grids, Version 1 [Data Set]. Boulder, Colorado USA. NASA National Snow and Ice Data Center Distributed Active Archive Center. https://doi.org/10.5067/TRUIAL3WPAUP. Date Accessed 09-24-2023
 
