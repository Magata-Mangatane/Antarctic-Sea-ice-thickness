# -*- coding: utf-8 -*-
"""
Created on Thu May 23 16:57:54 2024

@author: mngma
"""

#estimate sea-ice thickness with assumption of zero-ice freeboard
import os
path = 'processed_data/IS-2_freeboards/'
destination_dir = 'processed_data/ZIF/'

os.chdir(path+'/')

#make necessary imports
import h5py
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import glob
import warnings
import cartopy.feature as cfeature

#first define empty lists to hold dates and filenames for the 5 years
month_names = [[],[],[],[],[]]
years = ["2019","2020","2021","2022","2023"]
nums = ["01","02","03","04","05","06","07","08","09","10","11","12"]
names = ["january","february","march","april","may","june","july","august","september","october","november","december"]

#sort the month names according to the years and only computing 11 months in 2023 due to lack of december data 
for i in range(0,5):
    month_name = month_names[i]
    year = years[i]
    if i ==4:
        for j in range(0,11):
            month_name.append(str(year)+"_"+str(names[j]))
    else:
        for j in range(0,12):
            month_name.append(str(year)+"_"+str(names[j]))

#load freeboard data  
yearly_freeboard_list = []   #empty list to hold the freeboard datasets

for k in range(0,5):
    year = years[k]
    os.chdir(path+str(year)+'/')
    monthly_list = []
    if k==4:
        for i in range(0,11):
            ds = xr.open_dataset(str(year)+'_'+str(names[i])+'_IS-2_gridded_daily_freeboard.nc')
            monthly_list.append(ds)
        yearly_freeboard_list.append(monthly_list)
    else:
        for i in range(0,12):
            ds = xr.open_dataset(str(year)+'_'+str(names[i])+'_IS-2_gridded_daily_freeboard.nc')
            monthly_list.append(ds)
        yearly_freeboard_list.append(monthly_list)


#convert freeboard to thickness
#define densities 
rhoi = 915.1 #kg/m^3
rhos = 300 #kg/m^3
rhow = 1023.9 #kg/m^3

#define uncertainities
ierror = 15 #kg/m^3
serror = 100 #kg/m^3
werror = 3 #kg/m^3



for k in range(0,5):
    month_name = month_names[k]
    year = years[k]
    os.chdir(path+'/'+str(year)+'/')
    if k == 4:
        for j in range(0,11):
            #os.chdir(path+'/'+str(year)+'/')
            ds = yearly_freeboard_data[k][j]
            stdfr = ds.standard_deviation 
            T = ds.freeboard*(rhos/(rhow-rhoi))
            T = T.where(T < 12) ; T = T.rename('SIT_ZRIF')
            freeboard = ds.freeboard
            ferror = (stdfr**2)*((rhos/(rhow-rhoi))**2)
            srho = (serror**2)*((freeboard/(rhow-rhoi))**2)
            wrho = (werror**2)*(-(freeboard*rhos/(rhow-rhoi)**2)**2)
            irho = (ierror**2)*((freeboard*rhos/(rhow-rhoi)**2)**2)
            variance = ferror+srho+wrho+irho
            uncertainity = variance**(1/2) ; uncertainity = uncertainity.rename('uncertainity')
            mds = xr.merge([T,uncertainity,freeboard,stdfr])
            os.chdir(destination_dir+'/'+str(year)+'/')
            mds.to_netcdf(str(month_names[k][j])+'_IS-2_gridded_daily_sit_from_ZIF.nc')
    else:
        for j in range(0,12):
            #os.chdir(path+'/'+str(year)+'/')
            ds = yearly_freeboard_data[k][j]
            stdfr = ds.standard_deviation 
            T = ds.freeboard*(rhos/(rhow-rhoi))
            T = T.where(T < 12) ; T = T.rename('SIT_ZRIF')
            freeboard = ds.freeboard
            ferror = (stdfr**2)*((rhos/(rhow-rhoi))**2)
            srho = (serror**2)*((freeboard/(rhow-rhoi))**2)
            wrho = (werror**2)*(-(freeboard*rhos/(rhow-rhoi)**2)**2)
            irho = (ierror**2)*((freeboard*rhos/(rhow-rhoi)**2)**2)
            variance = ferror+srho+wrho+irho
            uncertainity = variance**(1/2) ; uncertainity = uncertainity.rename('uncertainity')
            mds = xr.merge([T,uncertainity,freeboard,stdfr])
            os.chdir(destination_dir+'/'+str(year)+'/')
            mds.to_netcdf(str(month_names[k][j])+'_IS-2_gridded_daily_sit_from_ZIF.nc')









