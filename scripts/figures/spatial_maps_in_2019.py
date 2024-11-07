# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 10:03:53 2024

@author: mngma
"""

#script to plot spatial maps of the SIT from the six satellite methods in 2019 
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
os.chdir(base_dir+'BERM/2019')

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
os.chdir(base_dir+'FDM/2019/')


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

#define a list with all seasonal data of each method for plotting
data_list = [zrif_aut,zrif_win,zrif_spr,xoc_aut,xoc_win,xoc_spr,boc_aut,boc_win,boc_spr,olm_aut,olm_win,olm_spr,sicc_aut,sicc_win,sicc_spr,fdm_aut,fdm_win,fdm_spr]

dates = ['March-May','June-Aug','Sept-Nov','March-May','June-Aug','Sept-Nov','March-May','June-Aug','Sept-Nov']
#sea ice thickness circumpolar maps 

#conditions for plotting
import cmocean
bounds_olm = np.arange(0.1,2.8,0.2)
bounds = bounds_olm
#cmap = mpl.cm.gnuplot2
cmap = cmocean.cm.thermal
cmap_diff = cmap
norm = mpl.colors.BoundaryNorm(bounds, cmap_diff.N+1, extend='both')
norm_olm = mpl.colors.BoundaryNorm(bounds_olm, cmap.N+1, extend='both')

import matplotlib.path as mpath
theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)

#plot each month, specifying colorbar locations using if,elif, and else statements.
fig = plt.figure(figsize=(15, 30),dpi=300)
for i in range(0,18):
    if (i == 0):
        ax = fig.add_subplot(6,3,i+1,projection=ccrs.SouthPolarStereo())
        ax.set_extent([-180, 180, -90,  -52], ccrs.PlateCarree())
        ax.set_boundary(circle, transform=ax.transAxes)
        plt.pcolormesh(data_list[i].lons[:], data_list[i].lats[:], data_list[i].values, cmap=cmap,norm=norm_olm,shading='nearest', transform=ccrs.PlateCarree())
       # plt.contour(sic_list[i].lons,sic_list[i].lats,sic_list[i],[15],transform=ccrs.PlateCarree()) 
        ax.coastlines()
        plt.title(dates[i],fontsize=16)
        plt.text(-0.08, 0.5, 'ZIF', ha='left', va='center', rotation='vertical',transform=ax.transAxes,fontsize=16)
        #ax.add_feature(cfeature.LAND, color='tab:gray')
       # ax.add_feature(cfeature.OCEAN, facecolor=("black"))
        #ax.coastlines(resolution='110m')
    elif i==3:
        ax = fig.add_subplot(6,3,i+1,projection=ccrs.SouthPolarStereo())
        ax.set_extent([-180, 180, -90,  -52], ccrs.PlateCarree())
        ax.set_boundary(circle, transform=ax.transAxes)
        plt.pcolormesh(data_list[i].lons[:], data_list[i].lats[:], data_list[i].values, cmap=cmap,norm=norm_olm,shading='nearest', transform=ccrs.PlateCarree())
       # plt.contour(sic_list[i].lons,sic_list[i].lats,sic_list[i],[15],transform=ccrs.PlateCarree()) 
        ax.coastlines()
        plt.text(-0.08, 0.5, 'ERM', ha='left', va='center', rotation='vertical',transform=ax.transAxes,fontsize=16)
        #ax.add_feature(cfeature.LAND, color='tab:gray')
        #ax.add_feature(cfeature.OCEAN, facecolor=("black"))
    elif (i == 1) or (i == 2):
        ax = fig.add_subplot(6,3,i+1,projection=ccrs.SouthPolarStereo())
        ax.set_extent([-180, 180, -90,  -52], ccrs.PlateCarree())
        ax.set_boundary(circle, transform=ax.transAxes)
        plt.pcolormesh(data_list[i].lons[:], data_list[i].lats[:], data_list[i].values, cmap=cmap,norm=norm_olm,shading='nearest', transform=ccrs.PlateCarree())
        #plt.contour(sic_list[i].lons,sic_list[i].lats,sic_list[i],[15],transform=ccrs.PlateCarree()) 
        ax.coastlines()
        plt.title(dates[i],fontsize=16)
        #ax.add_feature(cfeature.LAND, color='tab:gray')
       # ax.add_feature(cfeature.OCEAN, facecolor=("black"))
    elif i==4 or i==5:
        ax = fig.add_subplot(6,3,i+1,projection=ccrs.SouthPolarStereo())
        ax.set_extent([-180, 180, -90,  -52], ccrs.PlateCarree())
        ax.set_boundary(circle, transform=ax.transAxes)
        plt.pcolormesh(data_list[i].lons[:], data_list[i].lats[:], data_list[i].values, cmap=cmap,norm=norm_olm,shading='nearest', transform=ccrs.PlateCarree())
        #plt.contour(sic_list[i].lons,sic_list[i].lats,sic_list[i],[15],transform=ccrs.PlateCarree()) 
        ax.coastlines()
        #ax.add_feature(cfeature.LAND, color='tab:gray')
        #ax.add_feature(cfeature.OCEAN, facecolor=("black"))
    elif i==6:
        ax = fig.add_subplot(6,3,i+1,projection=ccrs.SouthPolarStereo())
        ax.set_extent([-180, 180, -90,  -52], ccrs.PlateCarree())
        ax.set_boundary(circle, transform=ax.transAxes)
        plt.pcolormesh(data_list[i].lons[:], data_list[i].lats[:], data_list[i].values, cmap=cmap_diff,norm=norm,shading='nearest', transform=ccrs.PlateCarree())
       # plt.contour(sic_list[i].lons,sic_list[i].lats,sic_list[i],[15],transform=ccrs.PlateCarree()) 
        ax.coastlines()
        plt.text(-0.08, 0.5, 'BERM', ha='left', va='center', rotation='vertical',transform=ax.transAxes,fontsize=16)
        #ax.add_feature(cfeature.LAND, color='tab:gray')
        #ax.add_feature(cfeature.OCEAN, facecolor=("black"))
    elif i==9:
        ax = fig.add_subplot(6,3,i+1,projection=ccrs.SouthPolarStereo())
        ax.set_extent([-180, 180, -90,  -52], ccrs.PlateCarree())
        ax.set_boundary(circle, transform=ax.transAxes)
        plt.pcolormesh(data_list[i].lons[:], data_list[i].lats[:], data_list[i].values, cmap=cmap_diff,norm=norm,shading='nearest', transform=ccrs.PlateCarree())
       # plt.contour(sic_list[i].lons,sic_list[i].lats,sic_list[i],[15],transform=ccrs.PlateCarree()) 
        ax.coastlines()
        plt.text(-0.08, 0.5, 'OLM', ha='left', va='center', rotation='vertical',transform=ax.transAxes,fontsize=16)
        #ax.add_feature(cfeature.LAND, color='tab:gray')
        #ax.add_feature(cfeature.OCEAN, facecolor=("black"))
    elif i==12:
        ax = fig.add_subplot(6,3,i+1,projection=ccrs.SouthPolarStereo())
        ax.set_extent([-180, 180, -90,  -52], ccrs.PlateCarree())
        ax.set_boundary(circle, transform=ax.transAxes)
        plt.pcolormesh(data_list[i].lons[:], data_list[i].lats[:], data_list[i].values, cmap=cmap_diff,norm=norm,shading='nearest', transform=ccrs.PlateCarree())
       # plt.contour(sic_list[i].lons,sic_list[i].lats,sic_list[i],[15],transform=ccrs.PlateCarree()) 
        ax.coastlines()
        plt.text(-0.08, 0.5, 'SICC', ha='left', va='center', rotation='vertical',transform=ax.transAxes,fontsize=16)
        #ax.add_feature(cfeature.LAND, color='tab:gray')
        #ax.add_feature(cfeature.OCEAN, facecolor=("black"))
    elif i==15:
        ax = fig.add_subplot(6,3,i+1,projection=ccrs.SouthPolarStereo())
        ax.set_extent([-180, 180, -90,  -52], ccrs.PlateCarree())
        ax.set_boundary(circle, transform=ax.transAxes)
        plt.pcolormesh(data_list[i].lons[:], data_list[i].lats[:], data_list[i].values, cmap=cmap_diff,norm=norm,shading='nearest', transform=ccrs.PlateCarree())
       # plt.contour(sic_list[i].lons,sic_list[i].lats,sic_list[i],[15],transform=ccrs.PlateCarree()) 
        ax.coastlines()
        plt.text(-0.08, 0.5, 'FDM', ha='left', va='center', rotation='vertical',transform=ax.transAxes,fontsize=16)
        #ax.add_feature(cfeature.LAND, color='tab:gray')
        #ax.add_feature(cfeature.OCEAN, facecolor=("black"))
    elif i ==16:
        ax = fig.add_subplot(6,3,i+1,projection=ccrs.SouthPolarStereo())
        ax.set_extent([-180, 180, -90,  -52], ccrs.PlateCarree())
        ax.set_boundary(circle, transform=ax.transAxes)
        plt.pcolormesh(data_list[i].lons[:], data_list[i].lats[:], data_list[i].values, cmap=cmap_diff,norm=norm,shading='nearest', transform=ccrs.PlateCarree())
        #plt.contour(sic_list[i].lons,sic_list[i].lats,sic_list[i],[15],transform=ccrs.PlateCarree()) 
        ax.coastlines()
        #ax.add_feature(cfeature.LAND, color='tab:gray')
        #ax.add_feature(cfeature.OCEAN, facecolor=("black"))
        axins1 = inset_axes(ax,width = "250%", height="6%", loc='center', bbox_to_anchor=(0, -0.6, 1, 1), bbox_transform=ax.transAxes, borderpad=0)  
        cb=plt.colorbar(mpl.cm.ScalarMappable(norm=norm_olm, cmap=cmap),orientation="horizontal",extend='both',label='Sea-ice thickness (m)',cax=axins1) 
        #axins1.yaxis.tick_left()
    else:
       ax = fig.add_subplot(6,3,i+1,projection=ccrs.SouthPolarStereo())
       ax.set_extent([-180, 180, -90,  -52], ccrs.PlateCarree())
       ax.set_boundary(circle, transform=ax.transAxes)
       plt.pcolormesh(data_list[i].lons[:], data_list[i].lats[:], data_list[i].values, cmap=cmap_diff,norm=norm,shading='nearest', transform=ccrs.PlateCarree())
       #plt.contour(sic_list[i].lons,sic_list[i].lats,sic_list[i],[15],transform=ccrs.PlateCarree()) 
       ax.coastlines()
       ##ax.add_feature(cfeature.LAND, color='tab:gray')
       #ax.add_feature(cfeature.OCEAN, facecolor=("black"))
fig.subplots_adjust(wspace=0,hspace=0)
plt.show()

















