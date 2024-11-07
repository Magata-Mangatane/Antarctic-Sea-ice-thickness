# Antarctic-Sea-ice-thickness: a repository of algorithms to reconstruct Antarctic Sea-ice thickness with satellite observations
_Authors: Magata J Mangatane , Marcello Vichi_

 This repository contains scripts to estimate Antarctic Sea-ice thickness with six satellite methods. The methods are used in Mangatane and Vichi manuscript, submitted to the Journal of Geophysical Research, and is currently under review.

## Brief description
This collection of scripts reconstructs the Antarctic circumpolar seasonal Sea-ice thickness (SIT) with six methods from 2019 to 2023. Four of the methods only require ICESat-2 sea-ice freeboard data, namely, the Zero sea-ice freeboard method (**ZIF**), the Empirical relationship method (**ERM**),the Buoyancy equation and empirical relationship (**BERM**), and the One-layer Method (**OLM**). The Sea Ice Climate Change Initiative (**SICC**) additionally requires independent snow thickness data, while the Freeboard differencing method (**FDM**) requires ICESat-2 and CryoSat-2 freeboard data. (see Data sources below)

## Workflow
The workflow is divided into the preprocessing of the data and the analysis. The raw data and preprocessed data are not provided due to their size but would be located in the [data/](data/) directory of this repository. These data are publicly available for download. 

 The scripts are provided in the [scripts/](scripts/) directory. This directory is divided into two sub-directories, [preprocessing/](preprocessing/) and [analysis/](analysis/) to preprocess the raw data and to estimate the Sea-ice thickness, respectively. We provide the processed data on zenodo.org


## Preprocessing of data
The preprocessing scripts are found in [scripts/preprocessing/](scripts/preprocessing/)

 Below are the preprocessing steps:
 1. Extract ICESat-2 freeboard data from HDF files to NetCDF file format ([extract_IS2_freeboards.py](scripts/preprocessing/extract_IS2_freeboards.py))
 2. Extract AMSR Sea-ice concentration data from He5 to NetCDF file format ([extract_amsr_sic_to_nc.py](scripts/preprocessing/extract_amsr_sic_to_nc.py))
 3. Grid CryoSat-2 data to a daily 25 km polar stereographic grid ([grid_CS2_freeboards.py](scripts/preprocessing/grid_CS2_freeboards.py))
 4. Weigh CryoSat-2 data with AMSR data ([weigh_cs2_freeboards.py](scripts/preprocessing/weigh_cs2_freeboards.py))
 5. Regrid AMSR snow depths onto the 25 km polar stereographic grid ([regrid_snow_depth_to_is2_grid.py](scripts/preprocessing/regrid_snow_depth_to_is2_grid.py))

## Analysis
The scripts to estimate Sea-ice thickness are found in [scripts/analysis/](scripts/analysis/)
 
 These are the different algorithms used:
 1. Estimation of SIT with the Zero sea-ice freeboard method ([ZIF.py](scripts/analysis/ZIF.py))
 2. Estimation of SIT with the Empirical relationship method ([ERM_and_OLM.py](scripts/analysis/ERM_and_OLM.py))
 3. Estimation of SIT with the Buoyancy equation and empirical relationship ([BERM.py](scripts/analysis/BERM.py))
 4. Estimation of SIT with the One-Layer Method ([ERM_and_OLM.py](scripts/analysis/ERM_and_OLM.py))
 5. Estimation of SIT with the Sea Ice Climate Change Initiative method ([SICC.py](scripts/analysis/SICC.py))
 6. Estimation of SIT with the freeboard differencing method ([FDM.py](scripts/analysis/FDM.py))

# Key figures
The scripts to plot key figures from the manuscript are found in the ([figures](scripts/figures)) folder

## Data sources
The raw data required to complete the reconstruction of Antarctic Sea-ice thickness demonstrated here are publicly available as follows:
* Kwok, R., A. A. Petty, G. Cunningham, T. Markus, D. Hancock, A. Ivanoff, J. Wimert, M. Bagnardi, N. Kurtz, and  the ICESat-2 Science Team. (2021). ATLAS/ICESat-2 L3A Sea Ice Height, Version 5 [Data Set]. Boulder, Colorado USA. NASA National Snow and Ice Data Center Distributed Active Archive Center. https://doi.org/10.5067/ATLAS/ATL07.005. Date Accessed 09-25-2023.
* European Space Agency, 2019, L2 SAR Precise Orbit. Baseline D. https://doi.org/10.5270/CR2-53hztdl
* Meier, W. N., Markus, T. & Comiso, J. C. (2018). AMSR-E/AMSR2 Unified L3 Daily 12.5 km Brightness Temperatures, Sea Ice Concentration, Motion & Snow Depth Polar Grids. (AU_SI12, Version 1). [Data Set]. Boulder, Colorado USA. NASA National Snow and Ice Data Center Distributed Active Archive Center. https://doi.org/10.5067/RA1MIJOYPK3P. [describe subset used if applicable]. Date Accessed 11-07-2024.

## References
Below are references to public repositories that were adopted to the methods presented in this repository. 
* Anthony Arendt, Ben Smith, David Shean, Amy Steiker, Alek Petty, Fernando Perez, Scott Henderson, Fernando Paolo, Johan Nilsson, Maya Becker, Susheel Adusumilli, Daniel Shapero, Bruce Wallin, Axel Schweiger, Suzanne Dickinson, Nicholas Hoschuh, Matthew Siegfried, Thomas Neumann. (2019). ICESAT-2HackWeek/ICESat2_hackweek_tutorials (Version 1.0). Zenodo. http://doi.org/10.5281/zenodo.3360994
* Jack Landy (2022) “jclandy/CryoSat2_Summer_SIT: CryoSat2_Summer_SIT_2022”. Zenodo. doi: 10.5281/zenodo.6558483.

 
