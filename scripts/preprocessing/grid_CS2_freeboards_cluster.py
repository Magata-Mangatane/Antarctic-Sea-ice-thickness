# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 09:36:04 2023

@author: mngma
"""

#This script grids cryosat-2 along-tract freeboards data onto a 25 km grids and weighs them with AMSRE Sea-ice concentration
#set path to data and change into the folder 
path = '/scratch/mngmag002/cs2_freeboards/'

#Make other necessary imports  
import os
import xarray as xr 
import glob
import pyproj
import numpy as np
import pandas as pd
from scipy import stats


#run the script for each month by changing the month name and directory below
#selected start changes the current month and associated date 
month_names = ['january','february','march','april','may','june','july','august','september','october','november','december']
month_dates = ["201901","201902","201903","201904","201905","201906","201907","201908","201909","201910","201911","201912"]

#define function to generate 25 km grid coordinates
#define grid for data
def create_grid_nsidc(epsg_string='3976', nx=316, ny=332, leftx=-3950000, dxRes=25000, uppery=4350000, dyRes=25000):
    """ Use pyproj to generate the NSIDC South Polar Stereographic grid covering the given domain (defined by the projection and the corner lat/lons)"""

    crs = pyproj.CRS.from_string("epsg:"+epsg_string)
    p=pyproj.Proj(crs)
    
    print(dxRes, dyRes)

    x=leftx+dxRes*np.indices((ny,nx),np.float32)[1]
    y=uppery-dxRes*np.indices((ny,nx),np.float32)[0]
    lons, lats = p(x, y, inverse=True)
    
    return x, y, lats, lons, p

xptsG, yptsG, latG, lonG, proj = create_grid_nsidc()

#define function to grid the data onto the generated coordinates 
def bin_data(xpts, ypts, var, xptsG, yptsG, dx):
    """ Bin data using numpy histogram2d
    
    Adapted for the NSIDC grid which has its orgin in the top left corner.

    """
    # Need to flip the arrays because the origin is in the top left but the histogram2d function needs x/y increasing.
    xptsG2=np.flipud(xptsG)
    yptsG2=np.flipud(yptsG)
    # get bin edges by subtracting half a grid-width and adding on another bin in both directions
    xbins=xptsG2[0]-(dx/2)
    ybins=yptsG2[:, 0]-(dx/2)
    xbins2=np.append(xbins, xbins[-1]+dx)
    ybins2=np.append(ybins, ybins[-1]+dx)
    print('binning..')
    counts, xedges, yedges = np.histogram2d(xpts, ypts,bins=(xbins2, ybins2))
    z, _, _ = np.histogram2d(xpts, ypts,bins=(xbins2, ybins2), weights=var)
    varG = z / counts
    #do same for std
    std =  stats.binned_statistic_2d(xpts, ypts, var, statistic='std', bins=[xbins2, ybins2])
    std = np.flipud(std.statistic.T)
    std = np.ma.masked_invalid(std)
    #count
    stdcounts =  stats.binned_statistic_2d(xpts, ypts, var, statistic='count', bins=[xbins2, ybins2])
    stdcounts = np.flipud(stdcounts.statistic.T)
    stdcounts = np.ma.masked_invalid(stdcounts)
    # Need to re-flip the arrays then transpose because of how histogram2d works across columns then rows.
    varG=np.flipud(varG.T)
    counts=np.flipud(counts.T)
    return varG,std,stdcounts

#Define a projection 
mapProj = pyproj.Proj("+init=EPSG:3976")

#define empty array for holding monthly data then merge into one file for the month
monthly_data = []

for i in range(0,12):
    data = []
    month_date = month_dates[i] ; month_name = month_names[i]
    os.chdir(path+month_name)
    for x in range(1,32):
        if x < 10:
            try:
                files = glob.glob('CS_LTA__SIR_SAR_2__'+month_date+'0'+ str(x) + '*.nc')
                ds_list = [xr.open_dataset(file) for file in files]
                ds_list = [ds_list[i].freeboard_20_ku for i in range(len(ds_list))]
                ds = xr.merge(ds_list)
                ds = ds.freeboard_20_ku
                ds = ds.where(ds > 0, drop=True)
                ds = ds.rename("freeboard")
                x10, y10=mapProj(ds.lon_poca_20_ku, ds.lat_poca_20_ku.values)
                varG=bin_data(x10, y10, ds.values, xptsG, yptsG, 25000)
                fr = varG[0]
                std = varG[1]
                countsN = varG[2]
                #expand dimension to add day
                std = np.expand_dims(std,2)
                countsN = np.expand_dims(countsN,2)
                fr1 = np.expand_dims(fr,2)
                time = pd.date_range(month_date +'0'+ str(x), periods=1)
                #create dataset
                ds = xr.Dataset(
                    data_vars=dict(
                        freeboard=(["x", "y","time"], fr1),
                        standard_dev=(["x", "y","time"], std),
                        N_counts=(["x", "y","time"], countsN),
                        ),
                    coords=dict(
                        lons=(["x", "y"], lonG),
                        lats=(["x", "y"], latG),
                        time = time,
                        ),
                    attrs=dict(description='CS-2 2019 gridded freeboard',comment='data gridded on 25 km NSIDC South stereographic grid'),
                    )
                data.append(ds)
            except:
                  pass
        else:
            try:
                files = glob.glob('CS_LTA__SIR_SAR_2__'+month_date+ str(x) + '*.nc')
                ds_list = [xr.open_dataset(file) for file in files]
                ds_list = [ds_list[i].freeboard_20_ku for i in range(len(ds_list))]
                ds = xr.merge(ds_list)
                ds = ds.freeboard_20_ku
                ds = ds.where(ds > 0, drop=True)
                ds = ds.rename("freeboard")
                x10, y10=mapProj(ds.lon_poca_20_ku, ds.lat_poca_20_ku.values)
                varG=bin_data(x10, y10, ds.values, xptsG, yptsG, 25000)
                fr = varG[0]
                std = varG[1]
                countsN = varG[2]
                #expand dimension to add day
                std = np.expand_dims(std,2)
                countsN = np.expand_dims(countsN,2)
                fr1 = np.expand_dims(fr,2)
                time = pd.date_range(month_date + str(x), periods=1)
                #create dataset
                ds = xr.Dataset(
                    data_vars=dict(
                        freeboard=(["x", "y","time"], fr1),
                        standard_dev=(["x", "y","time"], std),
                        N_counts=(["x", "y","time"], countsN),
                        ),
                    coords=dict(
                        lons=(["x", "y"], lonG),
                        lats=(["x", "y"], latG),
                        time = time,
                        ),
                    attrs=dict(description='CS-2 2019 gridded freeboard',comment='data gridded on 25 km NSIDC South stereographic grid'),
                    )
                data.append(ds)
            except:
                  pass
    monthly_data.append(xr.merge(data)) ; os.chdir(path+'/gridded/') ; monthly_data[i].to_netcdf('cs2_'+month_names[i]+'_circum_gridded_freeboard.nc')


#save monthly data to a single nc file
#os.chdir(path+'/gridded/')

#for i in range(0,2):
#    monthly_data[i].to_netcdf('cs2_'+month_names[i]+'_circum_gridded_freeboard.nc')


