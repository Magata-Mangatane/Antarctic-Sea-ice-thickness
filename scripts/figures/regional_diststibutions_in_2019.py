# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 10:49:13 2024

@author: mngma
"""
#script to intercompare sit estimated with ICESat-2 freeboards 
#set directory
import os
base_dir = 'processed_data/'
os.chdir(base_dir)

#import packages
import xarray as xr 
import matplotlib.pyplot as plt
import numpy as np
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import seaborn as sns
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib as mpl
import scipy.stats as stats
#load data from the different algorithms 

#OLM and ERM
os.chdir(base_dir+'OLM/2019')


month_date = ["2019-03-01","2019-04-01","2019-05-01","2019-06-01","2019-07-01","2019-08-01","2019-09-01","2019-10-01","2019-11-01"]
month_nms = ["2019_march","2019_april","2019_may","2019_june","2019_july","2019_august","2019_september","2019_october","2019_november"]

olm_w_list = [] ; xoc_w_list = []
olm_s_list = [] ; xoc_s_list = []
olm_a_list = [] ; xoc_a_list = []

for i in range(0,9):
    if i < 3:
        ds = xr.open_dataset(str(month_nms[i])+'_OLMi_daily_gridded_sit.nc')
        olm_ds = ds.SIT_OLMi ; xoc_ds = ds.I
        olm_ds = olm_ds.where((olm_ds > 0) & (olm_ds < 12))
        xoc_ds = xoc_ds.where((xoc_ds > 0) & (xoc_ds < 12))
        olm_a_list.append(olm_ds) ; xoc_a_list.append(xoc_ds)
    elif 2 < i < 6:
        ds = xr.open_dataset(str(month_nms[i])+'_OLMi_daily_gridded_sit.nc')
        olm_ds = ds.SIT_OLMi ; xoc_ds = ds.I
        olm_ds = olm_ds.where((olm_ds > 0) & (olm_ds < 12))  
        xoc_ds = xoc_ds.where((xoc_ds > 0) & (xoc_ds < 12))
        olm_w_list.append(olm_ds)  ; xoc_w_list.append(xoc_ds)
    else:
        ds = xr.open_dataset(str(month_nms[i])+'_OLMi_daily_gridded_sit.nc')
        olm_ds = ds.SIT_OLMi ; xoc_ds = ds.I
        olm_ds = olm_ds.where((olm_ds > 0) & (olm_ds < 12))
        xoc_ds = xoc_ds.where((xoc_ds > 0) & (xoc_ds < 12))
        olm_s_list.append(olm_ds) ; xoc_s_list.append(xoc_ds)


#merge into seasons
#OLMi
olm_win = xr.merge(olm_w_list[:])
olm_win = olm_win.SIT_OLMi.mean('time')
olm_spr = xr.merge(olm_s_list[:])
olm_spr = olm_spr.SIT_OLMi.mean('time')
olm_aut = xr.merge(olm_a_list[:])
olm_aut = olm_aut.SIT_OLMi.mean('time')


#xoc/ERM
xoc_win = xr.merge(xoc_w_list[:])
xoc_win = xoc_win.I.mean('time')
xoc_spr = xr.merge(xoc_s_list[:])
xoc_spr = xoc_spr.I.mean('time')
xoc_aut = xr.merge(xoc_a_list[:])
xoc_aut = xoc_aut.I.mean('time')


#boc sit
os.chdir(base_dir+'BOC/2019')

boc_w_list = []
boc_s_list = []
boc_a_list = []

for i in range(0,9):
    if i < 3:
        ds = xr.open_dataset(str(month_nms[i])+'_BERM_daily_gridded_sit.nc')
        ds = ds.BOC_SIT
        ds = ds.where((ds > 0) & (ds < 12))
        boc_a_list.append(ds)
    elif 2 < i < 6:
        ds = xr.open_dataset(str(month_nms[i])+'_BERM_daily_gridded_sit.nc')
        ds = ds.BOC_SIT
        ds = ds.where((ds > 0) & (ds < 12))
        boc_w_list.append(ds)
    else:
        ds = xr.open_dataset(str(month_nms[i])+'_BERM_daily_gridded_sit.nc')
        ds = ds.BOC_SIT
        ds = ds.where((ds > 0) & (ds < 12))
        boc_s_list.append(ds)


#merge into seasons
boc_win = xr.merge(boc_w_list[:])
boc_win = boc_win.BOC_SIT.mean('time')
boc_spr = xr.merge(boc_s_list[:])
boc_spr = boc_spr.BOC_SIT.mean('time')
boc_aut = xr.merge(boc_a_list[:])
boc_aut = boc_aut.BOC_SIT.mean('time')

#zero ice freeboard sit
os.chdir(base_dir+'zero_ice_freeboard/2019')

zrif_w_list = []
zrif_s_list = []
zrif_a_list = []

for i in range(0,9):
    if i < 3:
        ds = xr.open_dataset(str(month_nms[i])+'_IS-2_gridded_daily_sit_from_ZIF.nc')
        ds = ds.SIT_ZRIF
        ds = ds.where((ds > 0) & (ds < 12))
        zrif_a_list.append(ds)
    elif 2 < i < 6:
        ds = xr.open_dataset(str(month_nms[i])+'_IS-2_gridded_daily_sit_from_ZIF.nc')
        ds = ds.SIT_ZRIF
        ds = ds.where((ds > 0) & (ds < 12))
        zrif_w_list.append(ds)
    else:
        ds = xr.open_dataset(str(month_nms[i])+'_IS-2_gridded_daily_sit_from_ZIF.nc')
        ds = ds.SIT_ZRIF
        ds = ds.where((ds > 0) & (ds < 12))
        zrif_s_list.append(ds)


#merge into seasons
zrif_win = xr.merge(zrif_w_list[:])
zrif_win = zrif_win.SIT_ZRIF.mean('time')
zrif_spr = xr.merge(zrif_s_list[:])
zrif_spr = zrif_spr.SIT_ZRIF.mean('time')
zrif_aut = xr.merge(zrif_a_list[:])
zrif_aut = zrif_aut.SIT_ZRIF.mean('time')

#FDM
os.chdir(base_dir+'FDM/2019')


month_date = ["2019-03-01","2019-04-01","2019-05-01","2019-06-01","2019-07-01","2019-08-01","2019-09-01","2019-10-01","2019-11-01"]
month_nms = ["2019_march","2019_april","2019_may","2019_june","2019_july","2019_august","2019_september","2019_october","2019_november"]


fdm_w_list = [] 
fdm_s_list = [] 
fdm_a_list = [] 

for i in range(0,9):
    if i < 3:
        ds = xr.open_dataset(str(month_nms[i])+'_circum_sit_from_FDM_approach_25km.nc')
        fdm_ds = ds.sea_ice_thickness -0.32
        fdm_ds = fdm_ds.where((fdm_ds > 0) & (fdm_ds < 12))
        fdm_a_list.append(fdm_ds) 
    elif 2 < i < 6:
        ds = xr.open_dataset(str(month_nms[i])+'_circum_sit_from_FDM_approach_25km.nc')
        fdm_ds = ds.sea_ice_thickness -0.32
        fdm_ds = fdm_ds.where((fdm_ds > 0) & (fdm_ds < 12))
        fdm_w_list.append(fdm_ds)  
    else:
        ds = xr.open_dataset(str(month_nms[i])+'_circum_sit_from_FDM_approach_25km.nc')
        fdm_ds = ds.sea_ice_thickness -0.32
        fdm_ds = fdm_ds.where((fdm_ds > 0) & (fdm_ds < 12))
        fdm_s_list.append(fdm_ds) 


#merge into seasons
#fdm
fdm_win = xr.merge(fdm_w_list[:])
fdm_win = fdm_win.sea_ice_thickness.mean('time')
fdm_spr = xr.merge(fdm_s_list[:])
fdm_spr = fdm_spr.sea_ice_thickness.mean('time')
fdm_aut = xr.merge(fdm_a_list[:])
fdm_aut = fdm_aut.sea_ice_thickness.mean('time')

#sicc sit
os.chdir(base_dir+'SICC/2019')

sicc_w_list = []
sicc_s_list = []
sicc_a_list = []

for i in range(0,9):
    if i < 3:
        ds = xr.open_dataset(str(month_nms[i])+'_SICC_daily_gridded_sit.nc')
        ds = ds.SIT_SICC
        ds = ds.where((ds > 0) & (ds < 12))
        sicc_a_list.append(ds)
    elif 2 < i < 6:
        ds = xr.open_dataset(str(month_nms[i])+'_SICC_daily_gridded_sit.nc')
        ds = ds.SIT_SICC
        ds = ds.where((ds > 0) & (ds < 12))
        sicc_w_list.append(ds)
    else:
        ds = xr.open_dataset(str(month_nms[i])+'_SICC_daily_gridded_sit.nc')
        ds = ds.SIT_SICC
        ds = ds.where((ds > 0) & (ds < 12))
        sicc_s_list.append(ds)


#merge into seasons
sicc_win = xr.merge(sicc_w_list[:])
sicc_win = sicc_win.SIT_SICC.mean('time')
sicc_spr = xr.merge(sicc_s_list[:])
sicc_spr = sicc_spr.SIT_SICC.mean('time')
sicc_aut = xr.merge(sicc_a_list[:])
sicc_aut = sicc_aut.SIT_SICC.mean('time')


dates = ['March-May','June-Aug','Sept-Nov','March-May','June-Aug','Sept-Nov','March-May','June-Aug','Sept-Nov']


#define empty arrays to hold climatological statistics
KH = [-15,70]   # King Haakon VII
WS = [-60,-15]  # Weddell Sea
EA = [70,170]   # East Antarctica
AB = [-130,-60] # Amundsen/Bellingshausen Seas
#separate the Ross Sea into two to avoid an issue related to the longitude format (i.e. -180 to 180)
RS_1 = [-180,-130] # Ross Sea
RS_2 = [160,180]
SH = [-180,180] # Southern Hemisphere


region_old = [KH,WS,EA,AB,RS_1,RS_2,SH]

region = [AB,WS,KH,EA,RS_1,RS_2,SH]


#slice thickness data into regions
olm_list = [olm_aut,olm_win,olm_spr]
sit_list = []

for sit_ds in olm_list:
    ds_list = []
    for i in range(0,7):
        ds =  (sit_ds.where((sit_ds.lons >region[i][0])&(sit_ds.lons <region[i][1])))
        ds_list.append(ds)
    sit_list.append(ds_list) 

RS = []
for i in range(0,3):
    rs_ds = xr.merge([sit_list[i][4],sit_list[i][5]])
    RS.append(rs_ds.SIT_OLMi)
    
for i in range(0,3):
    sit_list[i].pop(4)
    sit_list[i].pop(4)
    sit_list[i].insert(4,RS[i])

#BOC regions
#slice thickness data into regions
boc_list = [boc_aut,boc_win,boc_spr]
sit_boc_list = []

for sit_ds in boc_list:
    ds_list = []
    for i in range(0,7):
        ds =  (sit_ds.where((sit_ds.lons >region[i][0])&(sit_ds.lons <region[i][1])))
        ds_list.append(ds)
    sit_boc_list.append(ds_list) 

RS = []
for i in range(0,3):
    rs_ds = xr.merge([sit_boc_list[i][4],sit_boc_list[i][5]])
    RS.append(rs_ds.BOC_SIT)
    
for i in range(0,3):
    sit_boc_list[i].pop(4)
    sit_boc_list[i].pop(4)
    sit_boc_list[i].insert(4,RS[i])


#xoc regions
xoc_list = [xoc_aut,xoc_win,xoc_spr]
sit_xoc_list = []

for sit_ds in xoc_list:
    ds_list = []
    for i in range(0,7):
        ds =  (sit_ds.where((sit_ds.lons >region[i][0])&(sit_ds.lons <region[i][1])))
        ds_list.append(ds)
    sit_xoc_list.append(ds_list) 

RS = []
for i in range(0,3):
    rs_ds = xr.merge([sit_xoc_list[i][4],sit_xoc_list[i][5]])
    RS.append(rs_ds.I)
    
for i in range(0,3):
    sit_xoc_list[i].pop(4)
    sit_xoc_list[i].pop(4)
    sit_xoc_list[i].insert(4,RS[i])

#xoc regions
zrif_list = [zrif_aut,zrif_win,zrif_spr]
sea_ice_thickness_list = []

for sit_ds in zrif_list:
    ds_list = []
    for i in range(0,7):
        ds =  (sit_ds.where((sit_ds.lons >region[i][0])&(sit_ds.lons <region[i][1])))
        ds_list.append(ds)
    sea_ice_thickness_list.append(ds_list) 

RS = []
for i in range(0,3):
    rs_ds = xr.merge([sea_ice_thickness_list[i][4],sea_ice_thickness_list[i][5]])
    RS.append(rs_ds.SIT_ZRIF)
    
for i in range(0,3):
    sea_ice_thickness_list[i].pop(4)
    sea_ice_thickness_list[i].pop(4)
    sea_ice_thickness_list[i].insert(4,RS[i])

#compute climatological statistics for each month 
sit_df_list_OLM = []
means_list_OLM = []
mins_list_OLM = []
medians_list_OLM = []
pers_list_OLM = []
stds_list_OLM = []

for i in range(0,3):
    ds = sit_list[i]
    df_list = [] ;  means = [] ; mins = [] ; medians = [] ; pers = [] ; stds = []
    for j in range(0,6):
        df = ds[j].rename('Sea-ice thickness'); df = df.to_dataframe(); df = df['Sea-ice thickness']
        df_list.append(df.dropna())
        means.append(round(np.nanmean(ds[j]),1))
        mins.append(round(np.nanmin(ds[j]),1))
        medians.append(round(np.nanmedian(ds[j]),1))
        stds.append(round(np.nanstd(ds[j]),1))
        pers.append(round(np.nanpercentile(ds[j], 95),1))
    sit_df_list_OLM.append(df_list) ; means_list_OLM.append(means); mins_list_OLM.append(mins)
    medians_list_OLM.append(medians) ; pers_list_OLM.append(pers) ; stds_list_OLM.append(stds)


#ERM statistics
#compute climatological statistics for each month 
sit_df_list_boc = []
means_list_boc = []
mins_list_boc = []
medians_list_boc = []
pers_list_boc = []
stds_list_boc = []

for i in range(0,3):
    ds = sit_boc_list[i]
    df_list = [] ;  means = [] ; mins = [] ; medians = [] ; pers = [] ; stds = []
    for j in range(0,6):
        df = ds[j].rename('Sea-ice thickness'); df = df.to_dataframe(); df = df['Sea-ice thickness']
        df_list.append(df.dropna())
        means.append(round(np.nanmean(ds[j]),1))
        mins.append(round(np.nanmin(ds[j]),1))
        medians.append(round(np.nanmedian(ds[j]),1))
        stds.append(round(np.nanstd(ds[j]),1))
        pers.append(round(np.nanpercentile(ds[j], 95),1))
    sit_df_list_boc.append(df_list) ; means_list_boc.append(means); mins_list_boc.append(mins)
    medians_list_boc.append(medians) ; pers_list_boc.append(pers) ; stds_list_boc.append(stds)


#ERM/XOC statistics 
sit_df_list_xoc = []
means_list_xoc = []
mins_list_xoc = []
medians_list_xoc = []
pers_list_xoc = []
stds_list_xoc = []

for i in range(0,3):
    ds = sit_xoc_list[i]
    df_list = [] ;  means = [] ; mins = [] ; medians = [] ; pers = [] ; stds = []
    for j in range(0,6):
        df = ds[j].rename('Sea-ice thickness'); df = df.to_dataframe(); df = df['Sea-ice thickness']
        df_list.append(df.dropna())
        means.append(round(np.nanmean(ds[j]),1))
        mins.append(round(np.nanmin(ds[j]),1))
        medians.append(round(np.nanmedian(ds[j]),1))
        stds.append(round(np.nanstd(ds[j]),1))
        pers.append(round(np.nanpercentile(ds[j], 95),1))
    sit_df_list_xoc.append(df_list) ; means_list_xoc.append(means); mins_list_xoc.append(mins)
    medians_list_xoc.append(medians) ; pers_list_xoc.append(pers) ; stds_list_xoc.append(stds)


#Zero-ice freeboard statistics
#compute climatological statistics for each month 
sit_df_list_zrif = []
means_list_zrif = []
mins_list_zrif = []
medians_list_zrif = []
pers_list_zrif = []
stds_list_zrif = []

for i in range(0,3):
    ds = sea_ice_thickness_list[i]
    df_list = [] ;  means = [] ; mins = [] ; medians = [] ; pers = [] ; stds = []
    for j in range(0,6):
        df = ds[j].rename('sea_ice_thickness'); df = df.to_dataframe(); df = df['sea_ice_thickness']
        df_list.append(df.dropna())
        means.append(round(np.nanmean(ds[j]),1))
        mins.append(round(np.nanmin(ds[j]),1))
        medians.append(round(np.nanmedian(ds[j]),1))
        stds.append(round(np.nanstd(ds[j]),1))
        pers.append(round(np.nanpercentile(ds[j], 95),1))
    sit_df_list_zrif.append(df_list) ; means_list_zrif.append(means); mins_list_zrif.append(mins)
    medians_list_zrif.append(medians) ; pers_list_zrif.append(pers) ; stds_list_zrif.append(stds)


#slice thickness data into regions
fdm_list = [fdm_aut,fdm_win,fdm_spr]
sit_fdm_list = []

for sit_ds in fdm_list:
    ds_list = []
    for i in range(0,7):
        ds =  (sit_ds.where((sit_ds.lons >region[i][0])&(sit_ds.lons <region[i][1])))
        ds_list.append(ds)
    sit_fdm_list.append(ds_list) 

RS = []
for i in range(0,3):
    rs_ds = xr.merge([sit_fdm_list[i][4],sit_fdm_list[i][5]])
    RS.append(rs_ds.sea_ice_thickness)
    
for i in range(0,3):
    sit_fdm_list[i].pop(4)
    sit_fdm_list[i].pop(4)
    sit_fdm_list[i].insert(4,RS[i])

#sicc regions
#slice thickness data into regions
sicc_list = [sicc_aut,sicc_win,sicc_spr]
sit_sicc_list = []

for sit_ds in sicc_list:
    ds_list = []
    for i in range(0,7):
        ds =  (sit_ds.where((sit_ds.lons >region[i][0])&(sit_ds.lons <region[i][1])))
        ds_list.append(ds)
    sit_sicc_list.append(ds_list) 

RS = []
for i in range(0,3):
    rs_ds = xr.merge([sit_sicc_list[i][4],sit_sicc_list[i][5]])
    RS.append(rs_ds.SIT_SICC)
    
for i in range(0,3):
    sit_sicc_list[i].pop(4)
    sit_sicc_list[i].pop(4)
    sit_sicc_list[i].insert(4,RS[i])

#compute climatological statistics for each month 
sit_df_list_sicc = []
means_list_sicc = []
mins_list_sicc = []
medians_list_sicc = []
pers_list_sicc = []
stds_list_sicc = []

for i in range(0,3):
    ds = sit_sicc_list[i]
    df_list = [] ;  means = [] ; mins = [] ; medians = [] ; pers = [] ; stds = []
    for j in range(0,6):
        df = ds[j].rename('Sea-ice thickness'); df = df.to_dataframe(); df = df['Sea-ice thickness']
        df_list.append(df.dropna())
        means.append(round(np.nanmean(ds[j]),1))
        mins.append(round(np.nanmin(ds[j]),1))
        medians.append(round(np.nanmedian(ds[j]),1))
        stds.append(round(np.nanstd(ds[j]),1))
        pers.append(round(np.nanpercentile(ds[j], 95),1))
    sit_df_list_sicc.append(df_list) ; means_list_sicc.append(means); mins_list_sicc.append(mins)
    medians_list_sicc.append(medians) ; pers_list_sicc.append(pers) ; stds_list_sicc.append(stds)


#FDM statistics
#compute climatological statistics for each month 
sit_df_list_fdm = []
means_list_fdm = []
mins_list_fdm = []
medians_list_fdm = []
pers_list_fdm = []
stds_list_fdm = []

for i in range(0,3):
    ds = sit_fdm_list[i]
    df_list = [] ;  means = [] ; mins = [] ; medians = [] ; pers = [] ; stds = []
    for j in range(0,6):
        df = ds[j].rename('Sea-ice thickness'); df = df.to_dataframe(); df = df['Sea-ice thickness']
        df_list.append(df.dropna())
        means.append(round(np.nanmean(ds[j]),2))
        mins.append(round(np.nanmin(ds[j]),2))
        medians.append(round(np.nanmedian(ds[j]),2))
        stds.append(round(np.nanstd(ds[j]),2))
        pers.append(round(np.nanpercentile(ds[j], 95),2))
    sit_df_list_fdm.append(df_list) ; means_list_fdm.append(means); mins_list_fdm.append(mins)
    medians_list_fdm.append(medians) ; pers_list_fdm.append(pers) ; stds_list_fdm.append(stds)
    
    

#plot monthly sea ice thickness distributions
region_sn_old = ['KH','WS','EA','AB','RS','SH']
region_sn = ['AB','WS','KH','EA','RS','SH']
seasons = ['March-May','June-Aug','Sept-Nov']


# Create a figure with subplots and a shared y-axis at the hemispheric scale
fig, axs= plt.subplots(1, 3, sharey=True,figsize=(13, 5),dpi=300)
for i in range(0,3):
    axs[i].set_xlim([0, 6])
    #plt.xticks(ticks=np.arange(0,5,1))
    axs[i].set_ylim([0, 1.7])
    axs[i].set_title(seasons[i],fontsize=12)
    axs[i].set_xlabel('Sea-ice thickness (m)')
    sns.kdeplot(sit_df_list_zrif[i][5],ax=axs[i],linestyle="-",fill=False,label='ZIF'+'      '+ str(means_list_zrif[i][5])+'     '+ str(medians_list_zrif[i][5])+'     '+ str(stds_list_zrif[i][5]))
    sns.kdeplot(sit_df_list_xoc[i][5],ax=axs[i],linestyle="-",fill=False,label='ERM'+'      '+ str(means_list_xoc[i][5])+'     '+ str(medians_list_xoc[i][5])+'     '+ str(stds_list_xoc[i][5]))
    sns.kdeplot(sit_df_list_boc[i][5],ax=axs[i],linestyle="-",fill=False,label='BERM'+'     '+ str(means_list_boc[i][5])+'      '+ str(medians_list_boc[i][5])+'     '+ str(stds_list_boc[i][5]))
    sns.kdeplot(sit_df_list_OLM[i][5],ax=axs[i],linestyle="-",fill=False,label='OLM'+'    '+ str(means_list_OLM[i][5])+'      '+ str(medians_list_OLM[i][5])+'     '+ str(stds_list_OLM[i][5]))
    sns.kdeplot(sit_df_list_sicc[i][5],ax=axs[i],linestyle="-",fill=False,label='SICC'+'      '+ str(means_list_sicc[i][5])+'     '+ str(medians_list_sicc[i][5])+'     '+ str(stds_list_sicc[i][5]))
    sns.kdeplot(sit_df_list_fdm[i][5],ax=axs[i],linestyle="-",fill=False,label='FDM'+'    '+ str(means_list_fdm[i][5])+'      '+ str(medians_list_fdm[i][5])+'     '+ str(stds_list_fdm[i][5]))
    #sns.kdeplot(sit_df_list_zrif[i][5],ax=axs[i],linestyle="-",fill=False,label='SICC'+'      '+ str(means_list_zrif[i][5])+'     '+ str(medians_list_zrif[i][5])+'     '+ str(stds_list_zrif[i][5]))
    axs[i].legend(title= '  Method  ' + ' Mean '+ ' Median ' + '  S.D',handlelength=0.8,loc='best',fontsize=10,frameon=True)
    plt.tight_layout()


#plot regions on the same figure
# Assuming you have 5 regions, 3 seasons, and the data is organized accordingly
seasons = ['Autumn', 'Winter', 'Spring']
region_names = ['Amundsen/Bellingshausen','Weddell Sea','King Haakon VII','East Antarctic','Ross Sea','Circumpolar']

# Create the subplots with 6 rows and 3 columns
fig, axs = plt.subplots(6, 3, sharey=True, sharex=False, figsize=(15, 21), dpi=300)

# Loop through each region and season
for j in range(0, 6):
    for i in range(3):       
        # Plot the KDE plots for each method
        sns.kdeplot(sit_df_list_zrif[i][j], ax=axs[j, i], linestyle="-", fill=False,
                    label='ZIF'+'        '+ str(medians_list_zrif[i][j])+'         '+ str(means_list_zrif[i][j])+'         '+ str(stds_list_zrif[i][j]))
        sns.kdeplot(sit_df_list_OLM[i][j], ax=axs[j, i], linestyle="-", fill=False,
                    label='OLM'+'       '+ str(medians_list_OLM[i][j])+'         '+ str(means_list_OLM[i][j])+'         '+ str(stds_list_OLM[i][j]))
        sns.kdeplot(sit_df_list_boc[i][j], ax=axs[j, i], linestyle="-", fill=False,
                    label='BERM'+'        '+ str(medians_list_boc[i][j])+'       '+ str(means_list_boc[i][j])+'       '+ str(stds_list_boc[i][j]))
        sns.kdeplot(sit_df_list_xoc[i][j], ax=axs[j, i], linestyle="-", fill=False,
                    label='ERM'+'       '+ str(medians_list_xoc[i][j])+'         '+ str(means_list_xoc[i][j])+'        '+ str(stds_list_xoc[i][j]))
        sns.kdeplot(sit_df_list_sicc[i][j], ax=axs[j, i], linestyle="-", fill=False,
                    label='SICC'+'       '+ str(medians_list_sicc[i][j])+'         '+ str(means_list_sicc[i][j])+'         '+ str(stds_list_sicc[i][j]))
        sns.kdeplot(sit_df_list_fdm[i][j], ax=axs[j, i], linestyle="-", fill=False,
                    label='FDM'+'       '+ str(medians_list_fdm[i][j])+'       '+ str(means_list_fdm[i][j])+'       '+ str(stds_list_fdm[i][j]))
        
        # Set legend font size to 15
        axs[j, i].legend(title='  Method  ' + ' Median ' + ' Mean ' + '  S.D',
                         handlelength=0.8,title_fontsize=14, loc='best', fontsize=12, frameon=True)

        # Set y-axis label font size to 15
        axs[j, i].set_ylabel('Density', fontsize=15)
        
        # Increase y-tick label size to 15
        axs[j, i].tick_params(axis='y', labelsize=15)

        axs[j, i].set_xlim([0, 4])
        axs[j, i].tick_params(axis='x', labelsize=15)
        axs[j, i].set_xlabel(None)
        
        # Add title for the first plot of each row
        if j == 0:
            axs[j, i].set_title(seasons[i], fontsize=15)

        # Add region name labels to the rightmost plots
        if i == 2:
            axs[j, i].text(1.05, 0.5, region_names[j], transform=axs[j, i].transAxes, fontsize=15,
                           va='center', ha='left', rotation=90)

        # Set x-axis labels for the bottom row
        if j == 5:
            axs[j, i].set_xlabel('Sea-ice thickness (m)', fontsize=15)

# Adjust layout for better spacing
plt.tight_layout()
plt.show()
