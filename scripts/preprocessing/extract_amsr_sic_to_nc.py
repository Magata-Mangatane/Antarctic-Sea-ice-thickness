# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 17:17:21 2022

@author: mngma
"""

#script to extract amsr unified sic data to nc for weighing cs2 freeboards
#set directory to sea ice concentration data 
import os
path = 'D:/manuscript_1/data/raw_data/amsre_sic/'
#make necessary imports
import h5py
import numpy as np
import xarray as xr
import glob
import pandas as pd

#set month names
months = ['may','june','july','august','september','october','november','december']

#extract monthly data
monthly_data = []
for j in range(len(months)):
    os.chdir(path+months[j])
    filenames = glob.glob('*.he5')
    daily_data = []
    for i in range(len(filenames)):
        #extract date from filename
        date = os.path.splitext(os.path.basename(filenames[i]))[0].split('_')[5] 
        time = pd.date_range(date, periods=1)
        #extract data from file into a dataset
        f = h5py.File(filenames[i], mode='r')
        data_fields = f['/HDFEOS/GRIDS/SpPolarGrid25km']
        sic = data_fields['Data Fields/SI_25km_SH_ICECON_DAY']
        lats = data_fields['lat'][:]
        lons = data_fields['lon'][:]
        units = sic.attrs['units']
        long_name = sic.attrs['long_name']
        comment = sic.attrs['comment']
        # Handle FillValue
        sic = np.ma.masked_equal(sic[:], 120) #mask land
        sic = np.ma.masked_equal(sic[:], 0)   #mask open water
        sic = np.ma.masked_equal(sic[:], 110)  # mask missing data
        sic = np.expand_dims(sic,2)
        sic = xr.Dataset(data_vars=dict(concentration=(["x", "y","time"], sic),),coords=dict(lons=(["x", "y"], lons),lats=(["x", "y"], lats),time = time,),attrs=dict(description=str(long_name),comment=str(comment)),)
        daily_data.append(sic)
    monthly_data.append(xr.merge(daily_data))
        

#save monthly datasets in netcdf format
processed_data_path = 'D:/manuscript_1/data/processed_data/amsr_sic/'
os.chdir(processed_data_path)

for j in range(len(months)):
    monthly_data[j].to_netcdf(months[j] + '_monthly_25km_sea-ice_concentration.nc')
    
    
    
    
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    