# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 17:44:41 2024

@author: mngma
"""

#script to reconstruct sea-ice thickness with the BERM

import os
path = 'processed_data/IS-2_freeboards/'
destination_dir = 'processed_data/BERM/'

os.chdir(path)

#make necessary imports
import xarray as xr

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


#define densities 
rhoi = 915.1 #kg/m^3
rhos = 300 #kg/m^3
rhow = 1023.9 #kg/m^3

#define uncertainities
ierror = 15 #kg/m^3
serror = 100 #kg/m^3
werror = 3 #kg/m^3
cerror = 1/100 #m

WW = [-60,-45] # Southern Hemisphere
WE = [-45,20]   # East Antarctica
IO = [20,90]   # King Haakon VII
PA = [90,160]  # Weddell Sea
AB = [-130,-60] # Amundsen/Bellingshausen Seas
#separate the Ross Sea into two to avoid an issue related to the longitude format (i.e. -180 to 180)
RS_1 = [-180,-130] # Ross Sea
RS_2 = [160,180]


region = [WW,WE,IO,PA,RS_1,RS_2,AB]

WW_pars = [22/100,2.34,0.9/100,0.88] ; WE_pars = [17.8/100,2.76,-1/100,0.87]
IO_pars = [26/100,3.5,-0.2/100,0.83] ; PA_pars = [26/100,3.5,-0.2/100,0.83] 
RS_pars = [20.9/100,2.45,-0.5/100,1.05] ; AB_pars = [16.9/100,2.79,-1.4/100,0.95]

pars_list = [WW_pars,WE_pars,IO_pars,PA_pars,RS_pars,RS_pars,AB_pars]

#define a list to hold the SIT data and loop over each month in each year computing the SIT and uncertainity

yearly_sit_list = []

for i in range(0,5):
    month_list = yearly_ds_list[i]
    sit_list = []
    if i ==4:
        for j in range(0,11):
            ds = month_list[j]
            ds_list = []
            for i in range(0,7):
                c = pars_list[i][2] ; d = pars_list[i][3] ; derror = 0.3*d
                fr_ds =  ds.where((ds.lons >region[i][0])&(ds.lons <region[i][1]))
                S = c+d*fr_ds.freeboard ; S = S.rename('snow_depth')
                fr_ds = fr_ds.where(S > 0)
                fr = fr_ds.freeboard
                stdfr = fr_ds.standard_deviation 
                deno = rhow-rhoi
                T = rhow*fr-(rhow-rhos)*(c+d*fr) ; T = T/deno
                T = T.where((T > 0) & (T < 12))
                thickness = T.rename("BERM_SIT") 
                rhowd = c*(rhoi-rhos)+(d-1)*fr*rhoi-d*fr*rhos ; rhowd = rhowd/(deno**2)
                rhoid = c*(rhos-rhow)+fr*(-d*rhow+d*rhos+rhow) ; rhoid = rhoid/(deno**2)
                rhosd = (c+d*fr)/deno
                frd = -d*rhow+d*rhos+rhow ; frd = frd/deno
                cd =  (rhos-rhow)/deno
                dd = fr*(rhos-rhow) ; dd = dd/deno
                varrhow = (rhowd**2)*(werror**2) ; varrhoi = (rhoid**2)*(ierror**2) ; varrhos = (rhosd**2)*(serror**2) ; varfr = (frd**2)*(stdfr**2)
                varc = (cd**2)*(cerror**2) ; vard = (dd**2)*(derror**2)
                std_Ti = varrhow + varrhoi + varrhos + varfr + varc + vard
                std_Ti = (std_Ti)**(1/2)
                std_Ti = std_Ti.rename('sit_error_BERM_std')
                #calculate ucnertainity with a constant freeboard uncertainity = 0.05 m 
                varfr_cons = (frd**2)*((0.05)**2)
                std_Ti_constant = varrhow + varrhoi + varrhos + varfr_cons + varc + vard
                varrhow=varrhow.rename('rhow')  
                varrhoi=varrhoi.rename('rhoi') 
                varrhos=varrhos.rename('rhos') 
                varfr_cons = fr+varfr_cons ; varfr_cons= varfr_cons-fr ; varfr_cons = varfr_cons.rename('varfr_cons') 
                varc = fr+varc ; varc= varc-fr ; varc = varc.rename('varc') 
                vard=vard.rename('vard')
                varfr=varfr.rename('varfr')
                std_Ti_constant = (std_Ti_constant)**(1/2)
                std_Ti_constant = std_Ti_constant.rename('berm_uncert_constant_fr_error')
                m_ds = xr.merge([thickness,S,fr,std_Ti])
                ds_list.append(m_ds)
            sit_list.append(xr.merge(ds_list))
        yearly_sit_list.append(sit_list)
    else:
        for j in range(0,12):
            ds = month_list[j]
            ds_list = []
            for i in range(0,7):
                c = pars_list[i][2] ; d = pars_list[i][3] ; derror = 0.3*d
                fr_ds =  ds.where((ds.lons >region[i][0])&(ds.lons <region[i][1]))
                S = c+d*fr_ds.freeboard ; S = S.rename('snow_depth')
                fr_ds = fr_ds.where(S > 0)
                fr = fr_ds.freeboard
                stdfr = fr_ds.standard_deviation #; stdfr = stdfr**(0.5)
                deno = rhow-rhoi
                T = rhow*fr-(rhow-rhos)*(c+d*fr) ; T = T/deno
                T = T.where((T > 0) & (T < 12))
                thickness = T.rename("BERM_SIT") 
                rhowd = c*(rhoi-rhos)+(d-1)*fr*rhoi-d*fr*rhos ; rhowd = rhowd/(deno**2)
                rhoid = c*(rhos-rhow)+fr*(-d*rhow+d*rhos+rhow) ; rhoid = rhoid/(deno**2)
                rhosd = (c+d*fr)/deno
                frd = -d*rhow+d*rhos+rhow ; frd = frd/deno
                cd =  (rhos-rhow)/deno
                dd = fr*(rhos-rhow) ; dd = dd/deno
                varrhow = (rhowd**2)*(werror**2) ; varrhoi = (rhoid**2)*(ierror**2) ; varrhos = (rhosd**2)*(serror**2) ; varfr = (frd**2)*(stdfr**2)
                varc = (cd**2)*(cerror**2) ; vard = (dd**2)*(derror**2)
                std_Ti = varrhow + varrhoi + varrhos + varfr + varc + vard
                std_Ti = (std_Ti)**(1/2)
                std_Ti = std_Ti.rename('sit_error_BERM_std')
                #calculate ucnertainity with a constant freeboard uncertainity = 0.05 m 
                varfr_cons = (frd**2)*((0.05)**2)
                std_Ti_constant = varrhow + varrhoi + varrhos + varfr_cons + varc + vard
                varrhow=varrhow.rename('rhow')  
                varrhoi=varrhoi.rename('rhoi') 
                varrhos=varrhos.rename('rhos') 
                varfr_cons = fr+varfr_cons ; varfr_cons= varfr_cons-fr ; varfr_cons = varfr_cons.rename('varfr_cons') 
                varc = fr+varc ; varc= varc-fr ; varc = varc.rename('varc') 
                vard=vard.rename('vard')
                varfr=varfr.rename('varfr')
                std_Ti_constant = (std_Ti_constant)**(1/2)
                std_Ti_constant = std_Ti_constant.rename('berm_uncert_constant_fr_error')
                m_ds = xr.merge([thickness,S,fr,std_Ti])
                ds_list.append(m_ds)
            sit_list.append(xr.merge(ds_list))
        yearly_sit_list.append(sit_list)



#save monthly data to a single nc file in folders of each year
for i in range(0,5):
    month_name = month_names[i]
    year = years[i]
    ds_list = yearly_sit_list[i]
    os.chdir(destination_dir+'/'+str(year)+'/')
    if i == 4:
        for j in range(0,11):
            ds_list[j].to_netcdf(str(month_name[j])+'_BERM_daily_gridded_sit.nc')
    else:
        for j in range(0,12):
            ds_list[j].to_netcdf(str(month_name[j])+'_BERM_daily_gridded_sit.nc')
        


