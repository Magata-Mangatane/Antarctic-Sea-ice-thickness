# -*- coding: utf-8 -*-
"""
Created on Sun Sep 24 21:20:33 2023

@author: MNGMA
"""

#script to weigh cryosat-2 freeboard data with amsr sic data
import os
processed_data_path = 'D:/manuscript_1/data/processed_data/'

#Make other necessary imports  
import xarray as xr 
import pandas as pd
import numpy as np
import glob 

#selected start changes the current month and associated date 
months = ['may','june','july','august','september','october','november','december']

#extract data from he5 file to an xarray dataset
weighted_monthly = []
for i in range(len(months)):
    sic = xr.open_dataset(processed_data_path+'/amsr_sic/'+months[i]+'_monthly_25km_sea-ice_concentration.nc')
    fr = xr.open_dataset(processed_data_path+'cryosat_2_freeboards/'+months[i]+'_circum_gridded_freeboard.nc')
    weighted_daily_fr = []
    valid_dates = fr.time
    filenames = glob.glob('D:/manuscript_1/data/raw_data/amsre_sic/'+months[i]+'/*.he5')
    for x in range(len(valid_dates)):
        try:
            date = valid_dates[i]
            dailysic = sic.sel(time=date)
            dailyfr = fr.sel(time=date)
            weighted_fr = dailyfr.freeboard*(dailysic.concentration/100)
            weighted_fr = np.expand_dims(weighted_fr,2)
            time = pd.date_range(pd.to_datetime(valid_dates[x].values), periods=1)
            flag = np.zeros([332,316, 1])   #weighted freeboard
            ds = xr.Dataset(
                data_vars=dict(
                    freeboard=(["x", "y","time"], weighted_fr),
                    flags = (["x", "y","time"],flag),
                    ),
                coords=dict(
                    lons=(["x", "y"], fr.lons.values),
                    lats=(["x", "y"], fr.lats.values),
                    time = time,
                    ),
                attrs=dict(description='CS-2 2019 gridded data',comment='data gridded on 25 km NSIDC South stereographic grid'),
                )
            weighted_daily_fr.append(ds)
        except:            
            date = valid_dates[i]
            meansic = sic.concentration.mean('time') 
            meanweighted_fr = fr.sel(time=date)
            meanweighted_fr = meanweighted_fr.freeboard*(meansic/100)
            meanweighted_fr = np.expand_dims(meanweighted_fr,2)
            flag = np.ones([332,316, 1])   #freeboard weighted with the mean
            time = pd.date_range(pd.to_datetime(valid_dates[x].values), periods=1)
            ds = xr.Dataset(
                data_vars=dict(
                    freeboard=(["x", "y","time"], meanweighted_fr),
                    flags = (["x", "y","time"],flag),
                    ),
                coords=dict(
                    lons=(["x", "y"], fr.lons.values),
                    lats=(["x", "y"], fr.lats.values),
                    time = time,
                    ),
                attrs=dict(description='CS-2 2019 gridded data',comment='data gridded on 25 km NSIDC South stereographic grid'),
                )
            weighted_daily_fr.append(ds)
    weighted_monthly.append(xr.merge(weighted_daily_fr))


#save datasets in netcdf format
os.chdir(processed_data_path+'/amsr_weighted_cs2_freeboards')
for i in range(len(months)):
    weighted_monthly[i].to_netcdf(months[i]+'_cs2_amsr_weighted_freeboard.nc')

