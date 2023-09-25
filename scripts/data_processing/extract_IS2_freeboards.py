# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 13:25:43 2023

@author: mngma
"""

#script to extract and save ICESat-2 ATL20 freeboard data to netcdf file for processing 
import os
path = 'E:/manuscript_1/data/raw_data/icesat_2_freeboards/'
os.chdir(path+'/')
#make necessary imports
import h5py
import numpy as np
import xarray as xr
import pandas as pd
import glob
import warnings

#run the script for each month by changing the month name and directory below
#first define empty lists to hold dates and filenames
month_dates = [[],[],[],[]]
month_names = [[],[],[],[]]
file_dates = [[],[],[],[]]

#define the year and month names
years = ["2019","2020","2021","2022"]
nums = ["01","02","03","04","05","06","07","08","09","10","11","12"]
names = ["january","february","march","april","may","june","july","august","september","october","november","december"]

#sort the month names according to the years
for i in range(0,4):
    month_date = month_dates[i]
    month_name = month_names[i]
    file_date = file_dates[i]
    year = years[i]
    for j in range(0,12):
        month_date.append(str(year)+"-"+str(nums[j])+"-")
        month_name.append(str(year)+"_"+str(names[j]))
        file_date.append(str(year)+str(nums[j]))

#load the filenames
filenames = []
for i in range(0,4):
    year = years[i]
    date = file_dates[i]
    os.chdir(path+str(year)+'/')
    filenames.append(glob.glob('ATL20-02_'+str(date)+'*.h5'))


#extract freeboard data into datasets for later
yearly_ds_list = []   #empty list to hold the datasets

for k in range(0,4):
    os.chdir(path+'/')
    monthly_ds_list = []
    month_date = month_dates[k]
    month_name = month_names[k]
    year = years[k]
    os.chdir(path+'/'+str(year)+'/')
    for i in range(0,12):
        data = []
        with h5py.File(filenames[k][i], mode='r') as file:
            for x in range(1,32):
                grid_lon = file['/grid_lon'][:]
                grid_lat = file['/grid_lat'][:]
                if x < 10:
                    try:
                        dset_name = '/daily/day0'+str(x)+'/mean_fb'
                        datavar = file[dset_name]
                        dset = datavar[:]
                        sd_name = '/daily/day0'+str(x)+'/sigma'
                        sdvar = file[sd_name]
                        sd = sdvar[:]
                        units = datavar.attrs['units']
                        long_name = datavar.attrs['long_name']
                        _FillValue = datavar.attrs['_FillValue']
                        time = pd.date_range(month_date[i] + str(x), periods=1)
                        dset = np.expand_dims(dset,2)
                        sd = np.expand_dims(sd,2)
                        dset[dset == _FillValue] = np.nan
                        dset = np.ma.masked_where(np.isnan(dset), dset)    
                        
                        ds = xr.Dataset(
                            data_vars=dict(
                                freeboard=(["x", "y","time"], dset),
                                standard_deviation=(["x", "y","time"], sd),
                                ),
                            coords=dict(
                                lons=(["x", "y"], grid_lon),
                                lats=(["x", "y"], grid_lat),
                                time = time,
                                ),
                            attrs=dict(description='IS-2 '+str(month_name[i])+' freeboard in meters'),
                            )   
                        data.append(ds)
                    except:
                        warnings.warn("Failed at i = "+str(i)+' and at x = '+ str(x))
                else:
                    try:
                        dset_name = '/daily/day'+str(x)+'/mean_fb'
                        datavar = file[dset_name]
                        dset = datavar[:]
                        sd_name = '/daily/day'+str(x)+'/sigma' 
                        sdvar = file[sd_name]
                        sd = sdvar[:]
                        units = datavar.attrs['units']
                        long_name = datavar.attrs['long_name']
                        _FillValue = datavar.attrs['_FillValue']
                        time = pd.date_range(month_date[i] + str(x), periods=1)
                        dset = np.expand_dims(dset,2)
                        sd = np.expand_dims(sd,2)
                        dset[dset == _FillValue] = np.nan
                        dset = np.ma.masked_where(np.isnan(dset), dset)    
                        
                        ds = xr.Dataset(
                            data_vars=dict(
                                freeboard=(["x", "y","time"], dset),
                                standard_deviation=(["x", "y","time"], sd),
                                ),
                            coords=dict(
                                lons=(["x", "y"], grid_lon),
                                lats=(["x", "y"], grid_lat),
                                time = time,
                                ),
                            attrs=dict(description='IS-2 '+str(month_name[i])+' freeboard in meters'),
                            )
                        data.append(ds)
                    except:
                        warnings.warn("Failed at i ="+str(i) +' and at x = '+ str(x))
        mds = xr.merge(data)
        monthly_ds_list.append(mds)
    yearly_ds_list.append(monthly_ds_list)
    print(year)
    

#save the monthly data to netcdf files
processed_data_path = 'E:/manuscript_1/data/processed_data/icesat_2_freeboards/'
for k in range(0,4):
    year = years[k]
    month = month_names[k]
    os.chdir(processed_data_path+str(year)+'/')
    ds_list = yearly_ds_list[k]
    for i in range(0,12):
        ds_list[i].to_netcdf(str(month[i])+'_IS-2_daily_gridded_freeboard.nc')
        


































