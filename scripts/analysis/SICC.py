# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 14:35:29 2024

@author: mngma
"""

#script to estimate sea-ice thickness with the SICC Method using ICESat-2 freeboards
#import packages
import os
import xarray as xr

path = 'processed_data/IS-2_freeboards/'
destination_dir = 'processed_data/SICC/'

os.chdir(path)

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
yearly_ds_list = []   #empty list to hold the freeboard datasets

for k in range(0,5):
    year = years[k]
    os.chdir(path+str(year)+'/')
    monthly_list = []
    if k==4:
        for i in range(0,11):
            ds = xr.open_dataset(str(year)+'_'+str(names[i])+'_IS-2_gridded_daily_freeboard.nc')
            monthly_list.append(ds)
        yearly_ds_list.append(monthly_list)
    else:
        for i in range(0,12):
            ds = xr.open_dataset(str(year)+'_'+str(names[i])+'_IS-2_gridded_daily_freeboard.nc')
            monthly_list.append(ds)
        yearly_ds_list.append(monthly_list)


#load snow_depth data  
yearly_snow_list = []   #empty list to hold the datasets

for k in range(0,5):
    year = years[k]
    os.chdir(destination_dir+'AMSRE_snow_depth/25km_grid/'+str(year)+'/')
    monthly_list = []
    if k==4:
        for i in range(0,11):
            ds = xr.open_dataset(str(year)+'_'+str(names[i])+'_amsr_snow_depth.nc')
            monthly_list.append(ds)
        yearly_snow_list.append(monthly_list)
    else:
        for i in range(0,12):
            ds = xr.open_dataset(str(year)+'_'+str(names[i])+'_amsr_snow_depth.nc')
            monthly_list.append(ds)
        yearly_snow_list.append(monthly_list)

#convert freeboard to thickness
#define densities 
rhoi = 915.1 #kg/m^3
rhos = 300 #kg/m^3
rhow = 1023.9 #kg/m^3

#define uncertainities
ierror = 15 #kg/m^3
serror = 100 #kg/m^3
werror = 3 #kg/m^3


yearly_sit_list = []

for k in range(0,5):
    month_list = yearly_ds_list[k]
    snow_list = yearly_snow_list[k]
    sit_list = []
    if k ==4:
        for j in range(0,11):
            sn = snow_list[j].snow_thickness/100
            ds = month_list[j]
            mds = xr.merge([ds,sn])
            mds = mds.where(mds.freeboard > 0) ; mds = mds.where(mds.snow_thickness > 0)
            stdfr = mds.standard_deviation 
            #where F > S
            T_1 = mds.where(mds.freeboard > mds.snow_thickness)
            T1 = T_1.freeboard*(rhow/(rhow-rhoi)) + (T_1.snow_thickness)*((rhos-rhow)/(rhow-rhoi))
            #compute uncertainity
            error1 = (T_1.standard_deviation*((rhow)/(rhow-rhoi)))**2+((T_1.snow_thickness*0.3)*((rhos-rhow)/(rhow-rhoi)))**2+(ierror*((rhow*T_1.freeboard+rhos*T_1.snow_thickness-rhow*T_1.snow_thickness)/(rhow-rhoi)**2))**2+(serror*(T_1.snow_thickness)/(rhow-rhoi))**2
            error1 = error1**(0.5); error1 = error1.rename('uncertainity')
            #where F<= S
            T_2 = mds.where(mds.freeboard <= mds.snow_thickness)
            T2 = T_2.freeboard*(rhos/(rhow-rhoi))
            #compute uncertainity
            error2 = (T_2.standard_deviation*(rhow/(rhow-rhoi)))**2+(serror*(T_2.freeboard/(rhow-rhoi)))**2+(ierror*(rhow*T_2.freeboard/(rhow-rhoi)**2))**2
            error2 = error2**(0.5) ; error2 = error2.rename('uncertainity')
            T1 = T1.where(T1 < 12) ; T2 = T1.where(T2 < 12)
            mds = xr.merge([T1.rename('SIT_SICC'),T2.rename('SIT_SICC'),T1.rename('SIT_F_greater'),error1,error2])
            sit_list.append(mds)
        yearly_sit_list.append(sit_list)
    else:
        for j in range(0,12):
            sn = snow_list[j].snow_thickness/100
            ds = month_list[j]
            mds = xr.merge([ds,sn])
            mds = mds.where(mds.freeboard > 0) ; mds = mds.where(mds.snow_thickness > 0)
            stdfr = mds.standard_deviation 
            #where F > S
            T_1 = mds.where(mds.freeboard > mds.snow_thickness)
            T1 = T_1.freeboard*(rhow/(rhow-rhoi)) + (T_1.snow_thickness)*((rhos-rhow)/(rhow-rhoi))
            #compute uncertainity
            error1 = (T_1.standard_deviation*((rhow)/(rhow-rhoi)))**2+((T_1.snow_thickness*0.3)*((rhos-rhow)/(rhow-rhoi)))**2+(ierror*((rhow*T_1.freeboard+rhos*T_1.snow_thickness-rhow*T_1.snow_thickness)/(rhow-rhoi)**2))**2+(serror*(T_1.snow_thickness)/(rhow-rhoi))**2
            error1 = error1**(0.5); error1 = error1.rename('uncertainity')
            #where F<= S
            T_2 = mds.where(mds.freeboard <= mds.snow_thickness)
            T2 = T_2.freeboard*(rhos/(rhow-rhoi))
            #compute uncertainity
            error2 = (T_2.standard_deviation*(rhow/(rhow-rhoi)))**2+(serror*(T_2.freeboard/(rhow-rhoi)))**2+(ierror*(rhow*T_2.freeboard/(rhow-rhoi)**2))**2
            error2 = error2**(0.5) ; error2 = error2.rename('uncertainity')
            T1 = T1.where(T1 < 12) ; T2 = T1.where(T2 < 12)
            mds = xr.merge([T1.rename('SIT_SICC'),T2.rename('SIT_SICC'),T1.rename('SIT_F_greater'),error1,error2])
            sit_list.append(mds)
        yearly_sit_list.append(sit_list)

#save monthly data to a single nc file
for i in range(0,5):
    month_name = month_names[i]
    year = years[i]
    ds_list = yearly_sit_list[i]
    os.chdir(destination_dir+'SICC/'+str(year)+'/')
    if i == 4:
        for j in range(0,11):
            ds_list[j].to_netcdf(str(month_name[j])+'_SICC_daily_gridded_sit.nc')
    else:
        for j in range(0,12):
            ds_list[j].to_netcdf(str(month_name[j])+'_SICC_daily_gridded_sit.nc')
        

