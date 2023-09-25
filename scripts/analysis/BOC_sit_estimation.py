# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 15:58:46 2023

@author: mngma
"""

#estimation of SIT and uncertainity with the improved buoyancy equation (BOC)
import os
processed_data_path = 'data/processed_data/icesat_2_freeboards/'
os.chdir(processed_data_path)
#make necessary imports
import xarray as xr

#run the script for each month by changing the month name and directory below
#first define empty lists to hold dates and filenames
month_names = [[],[],[],[]]

#define the year and month names
years = ["2019","2020","2021","2022"]
nums = ["01","02","03","04","05","06","07","08","09","10","11","12"]
names = ["january","february","march","april","may","june","july","august","september","october","november","december"]

#sort the month names according to the years
for i in range(0,4):
    month_name = month_names[i]
    year = years[i]
    for j in range(0,12):
        month_name.append(str(year)+"_"+str(names[j]))

#load freeboard data  
yearly_ds_list = []   #empty list to hold the datasets

for k in range(0,4):
    year = years[k]
    os.chdir(processed_data_path+str(year)+'/')
    monthly_list = []
    for i in range(0,12):
        ds = xr.open_dataset(str(year)+'_'+str(names[j])+'_IS-2_daily_gridded_freeboard.nc')
        monthly_list.append(ds)
    yearly_ds_list.append(monthly_list)
    
##convert freeboard to thickness
#define densities 
rhoi = 915.1 #kg/m^3
rhos = 300 #kg/m^3
rhow = 1023.9 #kg/m^3

#define uncertainities
ierror = 15 #kg/m^3
serror = 100 #kg/m^3
werror = 3 #kg/m^3
aerror = 10/100 #m
cerror = 1/100 #m

#define sectors
WW = [-60,-45] # Weddell West 
WE = [-45,20]   # Weddell East
IO = [20,90]   # Indian Ocean
PA = [90,160]  # Pacific Ocean
RS_1 = [-180,-130] # Ross Sea
RS_2 = [160,180] # Ross Sea
AB = [-130,-60] # Amundsen/Bellingshausen Seas

region = [WW,WE,IO,PA,RS_1,RS_2,AB]

WW_pars = [22/100,2.34,0.9/100,0.88] ; WE_pars = [17.8/100,2.76,-1/100,0.87]
IO_pars = [26/100,3.5,-0.2/100,0.83] ; PA_pars = [26/100,3.5,-0.2/100,0.83] 
RS_pars = [20.9/100,2.45,-0.5/100,1.05] ; AB_pars = [16.9/100,2.79,-1.4/100,0.95]

pars_list = [WW_pars,WE_pars,IO_pars,PA_pars,RS_pars,RS_pars,AB_pars]

#estimate sea ice thickness for each year and month 
year_list = yearly_ds_list[0]            #demostrating with 2019 only
sit_list = []

for j in range(0,12):
    ds = year_list[j]
    ds_list = []
    for i in range(0,7):
        #select the coefficients according to sector and estimate their error
        c = pars_list[i][2] ; d = pars_list[i][3] ; derror = 0.3*d
        fr_ds =  ds.where((ds.lons >region[i][0])&(ds.lons <region[i][1]))
        #estimate the snow depth
        S = c+d*fr_ds.freeboard ; S = S.rename('snow_depth')
        fr_ds = fr_ds.where(S > 0)      # taking locations with postive snow depth estimations only
        fr = fr_ds.freeboard
        stdfr = fr_ds.standard_deviation ; stdfr = stdfr**(0.5)
        #estimate sea ice thickness
        deno = rhow-rhoi
        T = rhow*fr-(rhow-rhos)*(c+d*fr) ; T = T/deno
        thickness = T.where((T > 0) & (T < 12))
        thickness = T.rename("BOC_SIT") 
        #estimate the uncertainity with freeboard uncertainity equal std and then 0.05 m
        rhowd = c*(rhoi-rhos)+(d-1)*fr*rhoi-d*fr*rhos ; rhowd = rhowd/(deno**2)
        rhoid = c*(rhos-rhow)+fr*(-d*rhow+d*rhos+rhow) ; rhoid = rhoid/(deno**2)
        rhosd = (c+d*fr)/deno
        frd = -d*rhow+d*rhos+rhow ; frd = frd/deno
        cd =  (rhos-rhow)/deno
        dd = fr*(rhos-rhow) ; dd = dd/deno
        varrhow = (rhowd**2)*(werror**2) ; varrhoi = (rhoid**2)*(ierror**2) 
        varrhos = (rhosd**2)*(serror**2) ; varfr = (frd**2)*(stdfr**2) 
        varc = (cd**2)*(cerror**2) ; vard = (dd**2)*(derror**2)
        std_Ti = varrhow + varrhoi + varrhos + varfr + varc + vard
        std_Ti = (std_Ti)**(1/2) ; std_Ti = std_Ti.rename('sit_error_BOC_std')
        varfr_cons = (frd**2)*((0.05)**2)
        std_Ti_constant = varrhow + varrhoi + varrhos + varfr_cons + varc + vard
        varrhow=varrhow.rename('rhow')  ; varrhoi=varrhoi.rename('rhoi'); varrhos=varrhos.rename('rhos') 
        varfr_cons = fr+varfr_cons ; varfr_cons= varfr_cons-fr 
        varfr_cons = varfr_cons.rename('varfr_cons') 
        varc = fr+varc ; varc= varc-fr ; varc = varc.rename('varc') 
        vard=vard.rename('vard') ; varfr=varfr.rename('varfr')
        std_Ti_constant = (std_Ti_constant)**(1/2)
        std_Ti_constant = std_Ti_constant.rename('boc_uncert_constant_fr_error')
        m_ds = xr.merge([thickness,S,fr,std_Ti,std_Ti_constant,varrhow, varrhoi,varfr,
                         varrhos,varc,vard,varfr_cons])
        ds_list.append(m_ds)
    sit_list.append(xr.merge(ds_list))


#save monthly data to a single nc file
for i in range(0,1):
    month_name = month_names[i]
    year = years[i]
    os.chdir('data/processed_data/BOC_sit/'+str(year))
    for j in range(0,12):
        sit_list[j].to_netcdf(str(month_name[j])+'_BOC_daily_gridded_sit.nc')




        
