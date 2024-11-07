# -*- coding: utf-8 -*-
"""
Created on Fri May 24 09:08:06 2024

@author: mngma
"""

#script to regrid AMSR-E/2 snow depth data onto IS-2 grid following the link below
#https://www.icesat-2-sea-ice-state.info/content/7_gridded_data_wrangling.html

import os
base_dir = 'raw_data/AMSRE_snow_depth/'
destination_dir = 'processed_data/AMSRE_snow_depth/'


#make necessary imports
import h5py
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import glob
import cartopy.feature as cfeature
import pyproj
import scipy

#set dates, running a single month for each year. the example below shows december 2023
file_root = '2023_december'
year_and_month = "2023-12-" 

#regrid data onto IS-2 grid 
os.chdir(base_dir+'/25km_grid/')

grid_var = xr.open_dataarray('extract_IS2_grid.nc')
out_lons = grid_var.lons.values
out_lats = grid_var.lats.values

# Initialize map projection and project data to it
out_proj = 'EPSG:3412'

mapProj = pyproj.Proj("+init=" + out_proj)
xptsIS2, yptsIS2 = mapProj(out_lons, out_lats)


os.chdir(base_dir+'12.5km_grid/2023')
filenames = glob.glob('AMSR_U2_L3_SeaIce12km_B04_202312*.he5')

data_list = []
for x in range(0,len(filenames)):
    if x < 9: 
        time = pd.date_range(str(year_and_month)+"0" + str(x+1), periods=1)
        file = filenames[x]
        f = h5py.File(file, mode='r')
        data_fields = f['/HDFEOS/GRIDS/SpPolarGrid12km']
        lats = data_fields['lat'][:]
        lons = data_fields['lon'][:]
        xdim = data_fields['XDim'][:]
        ydim = data_fields['YDim'][:]
        snow_depth = data_fields['Data Fields/SI_12km_SH_SNOWDEPTH_5DAY']
        units = snow_depth.attrs['units']
        long_name = snow_depth.attrs['long_name']
        comment = snow_depth.attrs['comment']
        snow_depth = snow_depth[:]
        # Handle FillValue
        snow_depth = np.ma.masked_equal(snow_depth, 120) ;snow_depth = np.ma.masked_equal(snow_depth, 140)
        snow_depth = np.ma.masked_equal(snow_depth, 130) ;snow_depth = np.ma.masked_equal(snow_depth, 150)
        snow_depth = np.ma.masked_equal(snow_depth, 110) ;snow_depth = np.ma.masked_equal(snow_depth, 160)
        xptsNEW, yptsNEW = mapProj(lons,lats)
        snow_depth = np.expand_dims(snow_depth,2)
        snow_depth = xr.Dataset(data_vars=dict(snow_thickness=(["x", "y","time"], snow_depth),),coords=dict(lons=(["x", "y"], lons),lats=(["x", "y"], lats),time = time,),attrs=dict(description=str(long_name),comment=str(comment)),)
        snow_depth = scipy.interpolate.griddata((xptsNEW.flatten(),yptsNEW.flatten()), snow_depth.snow_thickness.mean('time').values.flatten(), (xptsIS2, yptsIS2), method = 'nearest')
        snow_depth = np.expand_dims(snow_depth,2)
        snow_depth = xr.Dataset(data_vars=dict(snow_thickness=(["x", "y","time"], snow_depth),),coords=dict(lons=(["x", "y"], out_lons),lats=(["x", "y"], out_lats),time = time,),attrs=dict(description=str(long_name),comment=str(comment)),)
        data_list.append(snow_depth)
    else:
        time = pd.date_range(str(year_and_month) + str(x+1), periods=1)
        file = filenames[x]
        f = h5py.File(file, mode='r')
        data_fields = f['/HDFEOS/GRIDS/SpPolarGrid12km']
        lats = data_fields['lat'][:]
        lons = data_fields['lon'][:]
        xdim = data_fields['XDim'][:]
        ydim = data_fields['YDim'][:]
        snow_depth = data_fields['Data Fields/SI_12km_SH_SNOWDEPTH_5DAY']
        units = snow_depth.attrs['units']
        long_name = snow_depth.attrs['long_name']
        snow_depth = snow_depth[:]
        # Handle FillValue
        snow_depth = np.ma.masked_equal(snow_depth, 120) ;snow_depth = np.ma.masked_equal(snow_depth, 140)
        snow_depth = np.ma.masked_equal(snow_depth, 130) ;snow_depth = np.ma.masked_equal(snow_depth, 150)
        snow_depth = np.ma.masked_equal(snow_depth, 110) ;snow_depth = np.ma.masked_equal(snow_depth, 160)
        xptsNEW, yptsNEW = mapProj(lons,lats)
        snow_depth = np.expand_dims(snow_depth,2)
        snow_depth = xr.Dataset(data_vars=dict(snow_thickness=(["x", "y","time"], snow_depth),),coords=dict(lons=(["x", "y"], lons),lats=(["x", "y"], lats),time = time,),attrs=dict(description=str(long_name),comment=str(comment)),)
        snow_depth = scipy.interpolate.griddata((xptsNEW.flatten(),yptsNEW.flatten()), snow_depth.snow_thickness.mean('time').values.flatten(), (xptsIS2, yptsIS2), method = 'nearest')
        snow_depth = np.expand_dims(snow_depth,2)
        snow_depth = xr.Dataset(data_vars=dict(snow_thickness=(["x", "y","time"], snow_depth),),coords=dict(lons=(["x", "y"], out_lons),lats=(["x", "y"], out_lats),time = time,),attrs=dict(description=str(long_name),comment=str(comment)),)
        data_list.append(snow_depth)

#merge the monthly data
mds = xr.merge(data_list)

#quick look at the snow depths
snow_thickness = mds.snow_thickness.mean('time')


fig = plt.figure(dpi=300)
ax = fig.add_subplot(1,1,1,projection=ccrs.SouthPolarStereo())
ax.set_extent([-180, 180, -90,  -55], ccrs.PlateCarree())
plt.pcolormesh(snow_thickness.lons[:], snow_thickness.lats[:], snow_thickness.values,shading='nearest', transform=ccrs.PlateCarree())
ax.coastlines()
gl = ax.gridlines(draw_labels=True)
ax.add_feature(cfeature.LAND, color='tab:gray')
ax.set_facecolor("white")
plt.colorbar()
plt.show()

#save the monthly data to a nc file 
os.chdir(destination+'/25km_grid/2023/')
mds.to_netcdf(str(file_root)+'_amsr_snow_depth.nc')





    




