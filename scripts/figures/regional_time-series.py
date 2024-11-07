# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 14:22:34 2024

@author: mngma
"""

#regional time-series of the six lidar methods 
#set directory
import os
data_dir = 'processed_data/'

#make necessary imports
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt


#define empty lists for holding data and define months
month_names_19 = [] ; month_names_20 = [] ; month_names_21 = [] ; month_names_22 = [] ; month_names_23 = []
month_names = [month_names_19,month_names_20,month_names_21,month_names_22,month_names_23]
nums = ["01","02","03","04","05","06","07","08","09","10","11","12"]
years = ["2019","2020","2021","2022","2023"]
names = ["january","february","march","april","may","june","july","august","september","october","november","december"]


for i in range(0,5):
    month_name = month_names[i]
    year = years[i]
    if i ==4:
        for j in range(0,11):
            month_name.append(str(year)+"_"+str(names[j]))
    else:
        for j in range(0,12):
            month_name.append(str(year)+"_"+str(names[j]))


filenames = [[],[],[],[],[]]
for i in range(0,5):
    os.chdir(data_dir+'OLM/')
    year = years[i]
    os.chdir(data_dir+'OLM/'+str(year))
    if i ==4:
        for j in range(0,11):
            filenames[i].append(str(year)+'_'+str(names[j])+'_OLMi_daily_gridded_sit.nc')
    else:
        for j in range(0,12):
            filenames[i].append(str(year)+'_'+str(names[j])+'_OLMi_daily_gridded_sit.nc')

yearly_sit_list = [] ; yearly_OLM_uncert_list = [] ; yearly_cons_list = []
ERM_yearly_sit_list = [] ; yearly_erm_uncert_list = []

for i in range(0,5):
    os.chdir(data_dir+'OLM/')
    year = years[i]
    os.chdir(data_dir+'OLM/'+str(year))
    sit_list = [] ; uncert_list = [] ; uncert_con_list = [] ; ERM_list = [] ; erm_uncert_list = []
    filename = filenames[i]
    if i ==4:
        for j in range(0,11):
            ds = xr.open_dataset(filename[j])
            uncert = ds.sit_error_OLMi_std ; erm_uncert = ds.ERM_uncertainity
            ERM_SIT = ds.I
            ds = ds.SIT_OLMi
            sit_list.append(ds); ERM_list.append(ERM_SIT) ; uncert_list.append(uncert) ; erm_uncert_list.append(erm_uncert)
        yearly_sit_list.append(sit_list) ;ERM_yearly_sit_list.append(ERM_list); yearly_OLM_uncert_list.append(uncert_list) ; yearly_erm_uncert_list.append(erm_uncert_list)
    else:
        for j in range(0,12):
            ds = xr.open_dataset(filename[j])
            uncert = ds.sit_error_OLMi_std ; erm_uncert = ds.ERM_uncertainity
            ERM_SIT = ds.I
            ds = ds.SIT_OLMi
            sit_list.append(ds); ERM_list.append(ERM_SIT) ; uncert_list.append(uncert); erm_uncert_list.append(erm_uncert)
        yearly_sit_list.append(sit_list) ; ERM_yearly_sit_list.append(ERM_list); yearly_OLM_uncert_list.append(uncert_list) ; yearly_erm_uncert_list.append(erm_uncert_list)
            

#load BOC data
filenames = [[],[],[],[],[]]
for i in range(0,5):
    os.chdir(data_dir+'BOC/')
    year = years[i]
    os.chdir(data_dir+'BOC/'+str(year))
    if i == 4:
        for j in range(0,11):
            filenames[i].append(str(year)+'_'+str(names[j])+'_BOC_daily_gridded_sit.nc')
    else:
        for j in range(0,12):
            filenames[i].append(str(year)+'_'+str(names[j])+'_BOC_daily_gridded_sit.nc')
        
BOC_yearly_sit_list = [] ; boc_yearly_uncert_list = []

for i in range(0,5):
    os.chdir(data_dir+'BOC/')
    year = years[i]
    os.chdir(data_dir+'BOC/'+str(year))
    sit_list = [] ; uncert_list = []
    filename = filenames[i]
    if i ==4:
        for j in range(0,11):
            mds = xr.open_dataset(filename[j])
            ds = mds.BOC_SIT ; uncert = mds.sit_error_BOC_std
            sit_list.append(ds) ; uncert_list.append(uncert)
        BOC_yearly_sit_list.append(sit_list)  ; boc_yearly_uncert_list.append(uncert_list)
    else:
        for j in range(0,12):
            mds = xr.open_dataset(filename[j])
            ds = mds.BOC_SIT ; uncert = mds.sit_error_BOC_std
            sit_list.append(ds) ; uncert_list.append(uncert)
        BOC_yearly_sit_list.append(sit_list) ; boc_yearly_uncert_list.append(uncert_list)
            
#load ZIF
filenames = [[],[],[],[],[]]
for i in range(0,5):
    os.chdir(data_dir+'zero_ice_freeboard/')
    year = years[i]
    os.chdir(data_dir+'zero_ice_freeboard/'+str(year))
    if i ==4:
        for j in range(0,11):
            filenames[i].append(str(year)+'_'+str(names[j])+'_IS-2_gridded_daily_sit_from_ZIF.nc')
    else:
        for j in range(0,12):
            filenames[i].append(str(year)+'_'+str(names[j])+'_IS-2_gridded_daily_sit_from_ZIF.nc')
        
ZIF_yearly_sit_list = [] ; ZIF_yearly_uncert_list = []

for i in range(0,5):
    os.chdir(data_dir+'zero_ice_freeboard/')
    year = years[i]
    os.chdir(data_dir+'zero_ice_freeboard/'+str(year))
    sit_list = [] ; uncert_list = []
    filename = filenames[i]
    if i ==4:
        for j in range(0,11):
            mds = xr.open_dataset(filename[j])
            ds = mds.SIT_ZRIF ; uncert = mds.uncertainity
            sit_list.append(ds) ; uncert_list.append(uncert) 
        ZIF_yearly_sit_list.append(sit_list) ; ZIF_yearly_uncert_list.append(uncert_list)
    else:
        for j in range(0,12):
            mds = xr.open_dataset(filename[j])
            ds = mds.SIT_ZRIF ; uncert = mds.uncertainity
            sit_list.append(ds) ; uncert_list.append(uncert)
        ZIF_yearly_sit_list.append(sit_list)  ; ZIF_yearly_uncert_list.append(uncert_list)

#SICC
filenames = [[],[],[],[],[]]
for i in range(0,5):
    os.chdir(data_dir+'SICC/')
    year = years[i]
    os.chdir(data_dir+'SICC/'+str(year))
    if i ==4:
        for j in range(0,11):
            filenames[i].append(str(year)+'_'+str(names[j])+'_SICC_daily_gridded_sit.nc')
    else:
        for j in range(0,12):
            filenames[i].append(str(year)+'_'+str(names[j])+'_SICC_daily_gridded_sit.nc')

yearly_sicc_list = [] ; yearly_sicc_uncert_list = []


for i in range(0,5):
    os.chdir(data_dir+'SICC/')
    year = years[i]
    os.chdir(data_dir+'SICC/'+str(year))
    sit_list = [] ; uncert_list = [] ; uncert_con_list = [] 
    filename = filenames[i]
    if i ==4:
        for j in range(0,11):
            ds = xr.open_dataset(filename[j])
            uncert = ds.uncertainity
            ds = ds.SIT_SICC
            sit_list.append(ds);  uncert_list.append(uncert) 
        yearly_sicc_list.append(sit_list) ; yearly_sicc_uncert_list.append(uncert_list) 
    else:
        for j in range(0,12):
            ds = xr.open_dataset(filename[j])
            uncert = ds.uncertainity
            ds = ds.SIT_SICC
            sit_list.append(ds);  uncert_list.append(uncert)
        yearly_sicc_list.append(sit_list) ; yearly_sicc_uncert_list.append(uncert_list) 

#FDM
filenames = [[],[],[],[],[]]
for i in range(0,5):
    os.chdir(data_dir+'FDM/')
    year = years[i]
    os.chdir(data_dir+'FDM/'+str(year))
    if i ==4:
        for j in range(0,11):
            filenames[i].append(str(year)+'_'+str(names[j])+'_circum_sit_from_FDM_approach_25km.nc')
    else:
        for j in range(0,12):
            filenames[i].append(str(year)+'_'+str(names[j])+'_circum_sit_from_FDM_approach_25km.nc')

yearly_fdm_list = [] ; yearly_fdm_uncert_list = []


for i in range(0,5):
    os.chdir(data_dir+'FDM/')
    year = years[i]
    os.chdir(data_dir+'FDM/'+str(year))
    sit_list = [] ; uncert_list = [] ; uncert_con_list = [] 
    filename = filenames[i]
    if i ==4:
        for j in range(0,11):
            ds = xr.open_dataset(filename[j])
            uncert = ds.uncertainity
            ds = ds.sea_ice_thickness
            sit_list.append(ds);  uncert_list.append(uncert) 
        yearly_fdm_list.append(sit_list) ; yearly_fdm_uncert_list.append(uncert_list) 
    else:
        for j in range(0,12):
            ds = xr.open_dataset(filename[j])
            uncert = ds.uncertainity
            ds = ds.sea_ice_thickness
            sit_list.append(ds);  uncert_list.append(uncert)
        yearly_fdm_list.append(sit_list) ; yearly_fdm_uncert_list.append(uncert_list) 

#CS2WFA
altimeter_path = '/CS2WFA/'
os.chdir(altimeter_path)
output_dirs = ['2019', '2020', '2021']

months = ['January','February','March','April','May','June','July','August','September','October','November','December']
month_numbers_2019 = [201901,201902,201903,201904,201905,201906,201907,201908,201909,201910,201911,201912]
month_numbers_2020 = [202001,202002,202003,202004,202005,202006,202007,202008,202009,202010,202011,202012]
month_numbers_2021 = [202101,202102,202103,202104,202105,202106,202107,202108]

month_numbers = [month_numbers_2019,month_numbers_2020,month_numbers_2021]

monthly_means = []
for i in range(0,3):
    os.chdir(altimeter_path+'/'+output_dirs[i])
    mean_values = []
    if i==2:
        for j in range(0,7):
            ds = xr.open_dataset(months[j]+'_CS2WFA_25km_'+str(month_numbers[i][j])+'.nc')
            ds = ds.rename({'lat': 'lats', 'lon': 'lons'})
            ds = ds.assign_coords(lats=(("y", "x"), ds['lats'].values), lons=(("y", "x"), ds['lons'].values))
            ds = ds.sea_ice_thickness.mean('time')
            mean_values.append(ds)
        monthly_means.append(mean_values)
    else:
        for j in range(0,12):
            ds = xr.open_dataset(months[j]+'_CS2WFA_25km_'+str(month_numbers[i][j])+'.nc')
            ds = ds.rename({'lat': 'lats', 'lon': 'lons'})
            ds = ds.assign_coords(lats=(("y", "x"), ds['lats'].values), lons=(("y", "x"), ds['lons'].values))
            ds = ds.sea_ice_thickness.mean('time')
            mean_values.append(ds)
        monthly_means.append(mean_values)


###
#plot with 2023
#plot time series of years following each other
y = ["J","F","M","A","M","J","J","A","S","O","N","D"
     ,"J","F","M","A","M","J","J","A","S","O","N","D"
     ,"J","F","M","A","M","J","J","A","S","O","N","D"
     ,"J","F","M","A","M","J","J","A","S","O","N","D"
     ,"J","F","M","A","M","J","J","A","S","O","N"]


os.chdir('scripts/')

import time_series_means
#the original does not include CS2WFA, we use this after 2021 because those data are not computing past the year
import time_series_means_original
import time_series_uncertainity_means


#compute yearly means using time series means function 

#load uncertainity values 
year_value = 0
year_2019_uncertainity = time_series_uncertainity_means.year_monthly_means(yearly_OLM_uncert_list,boc_yearly_uncert_list,yearly_erm_uncert_list,ZIF_yearly_uncert_list,yearly_sicc_uncert_list,yearly_fdm_uncert_list,year_value)
year_value = 1
year_2020_uncertainity = time_series_uncertainity_means.year_monthly_means(yearly_OLM_uncert_list,boc_yearly_uncert_list,yearly_erm_uncert_list,ZIF_yearly_uncert_list,yearly_sicc_uncert_list,yearly_fdm_uncert_list,year_value)
year_value = 2
year_2021_uncertainity = time_series_uncertainity_means.year_monthly_means(yearly_OLM_uncert_list,boc_yearly_uncert_list,yearly_erm_uncert_list,ZIF_yearly_uncert_list,yearly_sicc_uncert_list,yearly_fdm_uncert_list,year_value)
year_value = 3
year_2022_uncertainity = time_series_uncertainity_means.year_monthly_means(yearly_OLM_uncert_list,boc_yearly_uncert_list,yearly_erm_uncert_list,ZIF_yearly_uncert_list,yearly_sicc_uncert_list,yearly_fdm_uncert_list,year_value)
year_value = 4
year_2023_uncertainity = time_series_uncertainity_means.year_monthly_means(yearly_OLM_uncert_list,boc_yearly_uncert_list,yearly_erm_uncert_list,ZIF_yearly_uncert_list,yearly_sicc_uncert_list,yearly_fdm_uncert_list,year_value)


year_value = 0
year_2019 = time_series_means.year_monthly_means(yearly_sit_list,BOC_yearly_sit_list,ERM_yearly_sit_list,ZIF_yearly_sit_list,yearly_sicc_list,yearly_fdm_list,monthly_means,year_value)
year_value = 1
year_2020 = time_series_means.year_monthly_means(yearly_sit_list,BOC_yearly_sit_list,ERM_yearly_sit_list,ZIF_yearly_sit_list,yearly_sicc_list,yearly_fdm_list,monthly_means,year_value)
year_value = 2
year_2021 = time_series_means.year_monthly_means(yearly_sit_list,BOC_yearly_sit_list,ERM_yearly_sit_list,ZIF_yearly_sit_list,yearly_sicc_list,yearly_fdm_list,monthly_means,year_value)
year_value = 3
year_2022 = time_series_means_original.year_monthly_means(yearly_sit_list,BOC_yearly_sit_list,ERM_yearly_sit_list,ZIF_yearly_sit_list,yearly_sicc_list,yearly_fdm_list,year_value)
year_value = 4
year_2023 = time_series_means_original.year_monthly_means(yearly_sit_list,BOC_yearly_sit_list,ERM_yearly_sit_list,ZIF_yearly_sit_list,yearly_sicc_list,yearly_fdm_list,year_value)


###################


#plot monthly sea ice thickness distributions
#region_sn_old = ['KH','WS','EA','AB','RS','SH']
region_sn = ['AB','WS','KH','EA','RS','SH']
seasons = ['March-May','June-Aug','Sept-Nov']
#region_names_old = ['King Haakon VII','Weddell Sea','East Antarctic','Bellingshausen/Amundsen','Ross Sea', 'Circumpolar']
region_names = ['Amundsen/Bellingshausen','Weddell Sea','King Haakon VII','East Antarctic','Ross Sea','Circumpolar']


#plot regions on the same figure
region_names = ['Amundsen/Bellingshausen','Weddell Sea','King Haakon VII','East Antarctic','Ross Sea','Circumpolar']

fig = plt.figure(figsize=(18, 30),dpi=300)
for i in range(0,5):
    ax = fig.add_subplot(5,1,i+1)
    plt.plot(np.arange(0,12,1),year_2019[0][i],'purple',marker='.',label='ZIF')
    plt.plot(np.arange(12,24,1),year_2020[0][i],'purple',marker='.')
    plt.plot(np.arange(24,36,1),year_2021[0][i],'purple',marker='.')
    plt.plot(np.arange(36,48,1),year_2022[0][i],'purple',marker='.')
    plt.plot(np.arange(0,12,1),year_2019[1][i],'b',marker='.',label='ERM')
    plt.plot(np.arange(12,24,1),year_2020[1][i],'b',marker='.')
    plt.plot(np.arange(24,36,1),year_2021[1][i],'b',marker='.')
    plt.plot(np.arange(36,48,1),year_2022[1][i],'b',marker='.')
    plt.plot(np.arange(0,12,1),year_2019[2][i],'r',marker='.',label='BERM')
    plt.plot(np.arange(12,24,1),year_2020[2][i],'r',marker='.')
    plt.plot(np.arange(24,36,1),year_2021[2][i],'r',marker='.')
    plt.plot(np.arange(36,48,1),year_2022[2][i],'r',marker='.')
    plt.plot(np.arange(0,12,1),year_2019[3][i],'g',marker='.',label='OLM')
    plt.plot(np.arange(12,24,1),year_2020[3][i],'g',marker='.')
    plt.plot(np.arange(24,36,1),year_2021[3][i],'g',marker='.')
    plt.plot(np.arange(36,48,1),year_2022[3][i],'g',marker='.')
    plt.plot(np.arange(0,12,1),year_2019[4][i],'k',marker='.',label='SICC')
    plt.plot(np.arange(12,24,1),year_2020[4][i],'k',marker='.')
    plt.plot(np.arange(24,36,1),year_2021[4][i],'k',marker='.')
    plt.plot(np.arange(36,48,1),year_2022[4][i],'k',marker='.')
    plt.text(1.01, 0.5, region_names[i],transform=ax.transAxes,fontsize=16, va='center', ha='left', rotation=90)
    plt.legend(loc="center right",bbox_to_anchor=(1.1,0.5))
    #plt.xlabel('Month',fontsize='xx-large')
    if i==0:
        plt.text(5,2.8,'2019',fontsize='xx-large')
        plt.text(18,2.8,'2020',fontsize='xx-large')
        plt.text(30,2.8,'2021',fontsize='xx-large')
        plt.text(42,2.8,'2022',fontsize='xx-large')
    if i == 4:
        plt.xlabel('Month',fontsize='xx-large')
      #  plt.yticks(ticks=np.arange(0,4,0.2),fontsize='xx-large')
    plt.xticks(ticks=np.arange(0,59,1),labels=y,fontsize='xx-large')
    plt.ylabel('Sea-ice thickness (m)',fontsize='xx-large')
    plt.xticks(ticks=np.arange(0,48,1),labels=y[:-11],fontsize='xx-large')
    plt.yticks(fontsize='xx-large')
    plt.margins(0.025,0.025)

#plot with uncertainity information 
fig = plt.figure(figsize=(18, 30),dpi=300)
for i in range(0,6):
    ax = fig.add_subplot(6,1,i+1)
    #SICC
    plt.plot(np.arange(0,12,1),year_2019[4][i],'k',marker='.',label='SICC')
    plt.fill_between(x=np.arange(0,12,1),y1=np.add(np.array(year_2019[4][i]), np.array(year_2019_uncertainity[4][i])),y2=np.subtract(np.array(year_2019[4][i]), np.array(year_2019_uncertainity[4][i])),color='gray', alpha=0.2)
    plt.plot(np.arange(12,24,1),year_2020[4][i],'k',marker='.')
    plt.fill_between(x=np.arange(12,24,1),y1=np.add(np.array(year_2020[4][i]), np.array(year_2020_uncertainity[4][i])),y2=np.subtract(np.array(year_2020[4][i]), np.array(year_2020_uncertainity[4][i])),color='gray', alpha=0.2)
    plt.plot(np.arange(24,36,1),year_2021[4][i],'k',marker='.')
    plt.fill_between(x=np.arange(24,36,1),y1=np.add(np.array(year_2021[4][i]), np.array(year_2021_uncertainity[4][i])),y2=np.subtract(np.array(year_2021[4][i]), np.array(year_2021_uncertainity[4][i])),color='gray', alpha=0.2)
    plt.plot(np.arange(36,48,1),year_2022[4][i],'k',marker='.')
    plt.fill_between(x=np.arange(36,48,1),y1=np.add(np.array(year_2022[4][i]), np.array(year_2022_uncertainity[4][i])),y2=np.subtract(np.array(year_2022[4][i]), np.array(year_2022_uncertainity[4][i])),color='gray', alpha=0.2)
    #FDM
    plt.plot(np.arange(0,12,1),year_2019[5][i],'tan',marker='.',label='FDM')
    plt.fill_between(x=np.arange(0,12,1),y1=np.add(np.array(year_2019[5][i]), np.array(year_2019_uncertainity[5][i])),y2=np.subtract(np.array(year_2019[5][i]), np.array(year_2019_uncertainity[5][i])),color='tan', alpha=0.2)
    plt.plot(np.arange(12,24,1),year_2020[5][i],'tan',marker='.')
    plt.fill_between(x=np.arange(12,24,1),y1=np.add(np.array(year_2020[5][i]), np.array(year_2020_uncertainity[5][i])),y2=np.subtract(np.array(year_2020[5][i]), np.array(year_2020_uncertainity[5][i])),color='tan', alpha=0.2)
    plt.plot(np.arange(24,36,1),year_2021[5][i],'tan',marker='.')
    plt.fill_between(x=np.arange(24,36,1),y1=np.add(np.array(year_2021[5][i]), np.array(year_2021_uncertainity[5][i])),y2=np.subtract(np.array(year_2021[5][i]), np.array(year_2021_uncertainity[5][i])),color='tan', alpha=0.2)
    plt.plot(np.arange(36,48,1),year_2022[5][i],'tan',marker='.')
    plt.fill_between(x=np.arange(36,48,1),y1=np.add(np.array(year_2022[5][i]), np.array(year_2022_uncertainity[5][i])),y2=np.subtract(np.array(year_2022[5][i]), np.array(year_2022_uncertainity[5][i])),color='tan', alpha=0.2)
    #BERM
    plt.plot(np.arange(0,12,1),year_2019[2][i],'r',marker='.',label='BERM')
    plt.fill_between(x=np.arange(0,12,1),y1=np.add(np.array(year_2019[2][i]), np.array(year_2019_uncertainity[2][i])),y2=np.subtract(np.array(year_2019[2][i]), np.array(year_2019_uncertainity[2][i])),color='red', alpha=0.2)
    plt.plot(np.arange(12,24,1),year_2020[2][i],'r',marker='.')
    plt.fill_between(x=np.arange(12,24,1),y1=np.add(np.array(year_2020[2][i]), np.array(year_2020_uncertainity[2][i])),y2=np.subtract(np.array(year_2020[2][i]), np.array(year_2020_uncertainity[2][i])),color='red', alpha=0.2)
    plt.plot(np.arange(24,36,1),year_2021[2][i],'r',marker='.')
    plt.fill_between(x=np.arange(24,36,1),y1=np.add(np.array(year_2021[2][i]), np.array(year_2021_uncertainity[2][i])),y2=np.subtract(np.array(year_2021[2][i]), np.array(year_2021_uncertainity[2][i])),color='red', alpha=0.2)
    #ZIF
    plt.plot(np.arange(36,48,1),year_2022[2][i],'r',marker='.')
    plt.fill_between(x=np.arange(36,48,1),y1=np.add(np.array(year_2022[2][i]), np.array(year_2022_uncertainity[2][i])),y2=np.subtract(np.array(year_2022[2][i]), np.array(year_2022_uncertainity[2][i])),color='red', alpha=0.2)
    plt.plot(np.arange(0,12,1),year_2019[0][i],'purple',marker='.',label='ZIF')
    plt.fill_between(x=np.arange(0,12,1),y1=np.add(np.array(year_2019[0][i]), np.array(year_2019_uncertainity[0][i])),y2=np.subtract(np.array(year_2019[0][i]), np.array(year_2019_uncertainity[0][i])),color='purple', alpha=0.2)
    plt.plot(np.arange(12,24,1),year_2020[0][i],'purple',marker='.')
    plt.fill_between(x=np.arange(12,24,1),y1=np.add(np.array(year_2020[0][i]), np.array(year_2020_uncertainity[0][i])),y2=np.subtract(np.array(year_2020[0][i]), np.array(year_2020_uncertainity[0][i])),color='purple', alpha=0.2)
    plt.plot(np.arange(24,36,1),year_2021[0][i],'purple',marker='.')
    plt.fill_between(x=np.arange(24,36,1),y1=np.add(np.array(year_2021[0][i]), np.array(year_2021_uncertainity[0][i])),y2=np.subtract(np.array(year_2021[0][i]), np.array(year_2021_uncertainity[0][i])),color='purple', alpha=0.2)
    plt.plot(np.arange(36,48,1),year_2022[0][i],'purple',marker='.')
    plt.fill_between(x=np.arange(36,48,1),y1=np.add(np.array(year_2022[0][i]), np.array(year_2022_uncertainity[0][i])),y2=np.subtract(np.array(year_2022[0][i]), np.array(year_2022_uncertainity[0][i])),color='purple', alpha=0.2)
    #ERM
    plt.plot(np.arange(0,12,1),year_2019[1][i],'b',marker='.',label='ERM')
    plt.fill_between(x=np.arange(0,12,1),y1=np.add(np.array(year_2019[1][i]), np.array(year_2019_uncertainity[1][i])),y2=np.subtract(np.array(year_2019[1][i]), np.array(year_2019_uncertainity[1][i])),color='blue', alpha=0.2)
    plt.plot(np.arange(12,24,1),year_2020[1][i],'b',marker='.')
    plt.fill_between(x=np.arange(12,24,1),y1=np.add(np.array(year_2020[1][i]), np.array(year_2020_uncertainity[1][i])),y2=np.subtract(np.array(year_2020[1][i]), np.array(year_2020_uncertainity[1][i])),color='blue', alpha=0.2)
    plt.plot(np.arange(24,36,1),year_2021[1][i],'b',marker='.')
    plt.fill_between(x=np.arange(24,36,1),y1=np.add(np.array(year_2021[1][i]), np.array(year_2021_uncertainity[1][i])),y2=np.subtract(np.array(year_2021[1][i]), np.array(year_2021_uncertainity[1][i])),color='blue', alpha=0.2)
    plt.plot(np.arange(36,48,1),year_2022[1][i],'b',marker='.')
    plt.fill_between(x=np.arange(36,48,1),y1=np.add(np.array(year_2022[1][i]), np.array(year_2022_uncertainity[1][i])),y2=np.subtract(np.array(year_2022[1][i]), np.array(year_2022_uncertainity[1][i])),color='blue', alpha=0.2)
    #OLM
    plt.plot(np.arange(0,12,1),year_2019[3][i],'g',marker='.',label='OLM')
    plt.fill_between(x=np.arange(0,12,1),y1=np.add(np.array(year_2019[3][i]), np.array(year_2019_uncertainity[3][i])),y2=np.subtract(np.array(year_2019[3][i]), np.array(year_2019_uncertainity[3][i])),color='green', alpha=0.2)
    plt.plot(np.arange(12,24,1),year_2020[3][i],'g',marker='.')
    plt.fill_between(x=np.arange(12,24,1),y1=np.add(np.array(year_2020[3][i]), np.array(year_2020_uncertainity[3][i])),y2=np.subtract(np.array(year_2020[3][i]), np.array(year_2020_uncertainity[3][i])),color='green', alpha=0.2)
    plt.plot(np.arange(24,36,1),year_2021[3][i],'g',marker='.')
    plt.fill_between(x=np.arange(24,36,1),y1=np.add(np.array(year_2021[3][i]), np.array(year_2021_uncertainity[3][i])),y2=np.subtract(np.array(year_2021[3][i]), np.array(year_2021_uncertainity[3][i])),color='green', alpha=0.2)
    plt.plot(np.arange(36,48,1),year_2022[3][i],'g',marker='.')
    plt.fill_between(x=np.arange(36,48,1),y1=np.add(np.array(year_2022[3][i]), np.array(year_2022_uncertainity[3][i])),y2=np.subtract(np.array(year_2022[3][i]), np.array(year_2022_uncertainity[3][i])),color='green', alpha=0.2)
    #CS2WFA
    #ax = fig.add_subplot(6,1,i+1)
    plt.plot(np.arange(0,12,1),year_2019[6][i],'orange',marker='.',label='CS2WFA')
    plt.plot(np.arange(12,24,1),year_2020[6][i],'orange',marker='.')
    plt.plot(np.arange(24,31,1),year_2021[6][i],'orange',marker='.')
    plt.text(1.01, 0.5, region_names[i],transform=ax.transAxes,fontsize=16, va='center', ha='left', rotation=90)
    plt.legend(loc="center right",bbox_to_anchor=(1.15,0.5))
    #plt.xlabel('Month',fontsize='xx-large')
    if i==0:
        plt.text(5,6,'2019',fontsize='xx-large')
        plt.text(18,6,'2020',fontsize='xx-large')
        plt.text(30,6,'2021',fontsize='xx-large')
        plt.text(42,6,'2022',fontsize='xx-large')
    if i == 4:
        plt.xlabel('Month',fontsize='xx-large')
      #  plt.yticks(ticks=np.arange(0,4,0.2),fontsize='xx-large')
    plt.xticks(ticks=np.arange(0,59,1),labels=y,fontsize='xx-large')
    plt.ylabel('Sea-ice thickness (m)',fontsize='xx-large')
    plt.xticks(ticks=np.arange(0,48,1),labels=y[:-11],fontsize='xx-large')
    plt.yticks(fontsize='xx-large')
    plt.margins(0.025,0.025)
    #plt.show()
    #plt.savefig('Regional_time_series.png')






