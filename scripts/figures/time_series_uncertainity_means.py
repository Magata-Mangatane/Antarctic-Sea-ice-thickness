# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 14:22:09 2024

@author: mngma
"""

#script to load mean uncertainities from six altimeter methods
import xarray as xr
import numpy as np

def year_monthly_means(yearly_OLM_uncert_list,boc_yearly_uncert_list,yearly_erm_uncert_list,ZIF_yearly_uncert_list,yearly_sicc_uncert_list,yearly_fdm_uncert_list,year_value):
    #define empty arrays to hold climatological statistics
    KH = [-15,70]   # King Haakon VII
    WS = [-60,-15]  # Weddell Sea
    EA = [70,170]   # East Antarctica
    AB = [-110,-60] # Amundsen/Bellingshausen Seas
    RS_1 = [-180,-110] # Ross Sea
    RS_2 = [170,180]
    SH = [-180,180] # Southern Hemisphere


    #region_old = [KH,WS,EA,AB,RS_1,RS_2,SH]

    region = [AB,WS,KH,EA,RS_1,RS_2,SH]


    #slice thickness data into regions
    sit_list = []

    for sit_ds in yearly_OLM_uncert_list[year_value]:
        ds_list = []
        for i in range(0,7):
            ds =  (sit_ds.where((sit_ds.lons >region[i][0])&(sit_ds.lons <region[i][1])))
            ds_list.append(ds)
        sit_list.append(ds_list) 

    RS = []
    for i in range(0,12):
        rs_ds = xr.merge([sit_list[i][4],sit_list[i][5]])
        RS.append(rs_ds.sit_error_OLMi_std)
        
    for i in range(0,12):
        sit_list[i].pop(4)
        sit_list[i].pop(4)
        sit_list[i].insert(4,RS[i])

    #BOC regions
    #slice thickness data into regions
    sit_boc_list = []

    for sit_ds in boc_yearly_uncert_list[year_value]:
        ds_list = []
        for i in range(0,7):
            ds =  (sit_ds.where((sit_ds.lons >region[i][0])&(sit_ds.lons <region[i][1])))
            ds_list.append(ds)
        sit_boc_list.append(ds_list) 

    RS = []
    for i in range(0,12):
        rs_ds = xr.merge([sit_boc_list[i][4],sit_boc_list[i][5]])
        RS.append(rs_ds.sit_error_BOC_std)
        
    for i in range(0,12):
        sit_boc_list[i].pop(4)
        sit_boc_list[i].pop(4)
        sit_boc_list[i].insert(4,RS[i])


    #xoc regions
    sit_xoc_list = []

    for sit_ds in yearly_erm_uncert_list[year_value]:
        ds_list = []
        for i in range(0,7):
            ds =  (sit_ds.where((sit_ds.lons >region[i][0])&(sit_ds.lons <region[i][1])))
            ds_list.append(ds)
        sit_xoc_list.append(ds_list) 

    RS = []
    for i in range(0,12):
        rs_ds = xr.merge([sit_xoc_list[i][4],sit_xoc_list[i][5]])
        RS.append(rs_ds.ERM_uncertainity)
        
    for i in range(0,12):
        sit_xoc_list[i].pop(4)
        sit_xoc_list[i].pop(4)
        sit_xoc_list[i].insert(4,RS[i])

    #xoc regions
    sea_ice_thickness_list = []

    for sit_ds in ZIF_yearly_uncert_list[year_value]:
        ds_list = []
        for i in range(0,7):
            ds =  (sit_ds.where((sit_ds.lons >region[i][0])&(sit_ds.lons <region[i][1])))
            ds_list.append(ds)
        sea_ice_thickness_list.append(ds_list) 

    RS = []
    for i in range(0,12):
        rs_ds = xr.merge([sea_ice_thickness_list[i][4],sea_ice_thickness_list[i][5]])
        RS.append(rs_ds.uncertainity)
        
    for i in range(0,12):
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

    for i in range(0,12):
        ds = sit_list[i]
        df_list = [] ;  means = [] ; mins = [] ; medians = [] ; pers = [] ; stds = []
        for j in range(0,6):
            df = ds[j].rename('Sea-ice thickness'); df = df.to_dataframe(); df = df['Sea-ice thickness']
            df_list.append(df.dropna())
            means.append(round(np.nanmean(ds[j]),2))
            mins.append(round(np.nanmin(ds[j]),2))
            medians.append(round(np.nanmedian(ds[j]),2))
            stds.append(round(np.nanstd(ds[j]),2))
            pers.append(round(np.nanpercentile(ds[j], 95),2))
        sit_df_list_OLM.append(df_list) ; means_list_OLM.append(means); mins_list_OLM.append(mins)
        medians_list_OLM.append(medians) ; pers_list_OLM.append(pers) ; stds_list_OLM.append(stds)

    KH = []
    WS = []
    EA = []
    AB = []
    RS = []
    SH = []

    # Rearrange the data
    for i in range(0,12):
        AB.append(means_list_OLM[i][0])
        WS.append(means_list_OLM[i][1])
        KH.append(means_list_OLM[i][2])
        EA.append(means_list_OLM[i][3])
        RS.append(means_list_OLM[i][4])
        SH.append(means_list_OLM[i][5])

    # Final rearranged list
    means_list_OLM = [AB,WS,KH,EA,RS,SH]
    #ERM statistics
    #compute climatological statistics for each month 
    sit_df_list_boc = []
    means_list_boc = []
    mins_list_boc = []
    medians_list_boc = []
    pers_list_boc = []
    stds_list_boc = []

    for i in range(0,12):
        ds = sit_boc_list[i]
        df_list = [] ;  means = [] ; mins = [] ; medians = [] ; pers = [] ; stds = []
        for j in range(0,6):
            df = ds[j].rename('Sea-ice thickness'); df = df.to_dataframe(); df = df['Sea-ice thickness']
            df_list.append(df.dropna())
            means.append(round(np.nanmean(ds[j]),2))
            mins.append(round(np.nanmin(ds[j]),2))
            medians.append(round(np.nanmedian(ds[j]),2))
            stds.append(round(np.nanstd(ds[j]),2))
            pers.append(round(np.nanpercentile(ds[j], 95),2))
        sit_df_list_boc.append(df_list) ; means_list_boc.append(means); mins_list_boc.append(mins)
        medians_list_boc.append(medians) ; pers_list_boc.append(pers) ; stds_list_boc.append(stds)

    KH = []
    WS = []
    EA = []
    AB = []
    RS = []
    SH = []

    # Rearrange the data
    for i in range(0,12):
        AB.append(means_list_boc[i][0])
        WS.append(means_list_boc[i][1])
        KH.append(means_list_boc[i][2])
        EA.append(means_list_boc[i][3])
        RS.append(means_list_boc[i][4])
        SH.append(means_list_boc[i][5])

    # Final rearranged list
    means_list_boc = [AB,WS,KH,EA,RS,SH]
    #ERM/XOC statistics 
    sit_df_list_xoc = []
    means_list_xoc = []
    mins_list_xoc = []
    medians_list_xoc = []
    pers_list_xoc = []
    stds_list_xoc = []

    for i in range(0,12):
        ds = sit_xoc_list[i]
        df_list = [] ;  means = [] ; mins = [] ; medians = [] ; pers = [] ; stds = []
        for j in range(0,6):
            df = ds[j].rename('Sea-ice thickness'); df = df.to_dataframe(); df = df['Sea-ice thickness']
            df_list.append(df.dropna())
            means.append(round(np.nanmean(ds[j]),2))
            mins.append(round(np.nanmin(ds[j]),2))
            medians.append(round(np.nanmedian(ds[j]),2))
            stds.append(round(np.nanstd(ds[j]),2))
            pers.append(round(np.nanpercentile(ds[j], 95),2))
        sit_df_list_xoc.append(df_list) ; means_list_xoc.append(means); mins_list_xoc.append(mins)
        medians_list_xoc.append(medians) ; pers_list_xoc.append(pers) ; stds_list_xoc.append(stds)

    KH = []
    WS = []
    EA = []
    AB = []
    RS = []
    SH = []

    # Rearrange the data
    for i in range(0,12):
        AB.append(means_list_xoc[i][0])
        WS.append(means_list_xoc[i][1])
        KH.append(means_list_xoc[i][2])
        EA.append(means_list_xoc[i][3])
        RS.append(means_list_xoc[i][4])
        SH.append(means_list_xoc[i][5])
        
    # Final rearranged list
    means_list_xoc = [AB,WS,KH,EA,RS,SH]

    #Zero-ice freeboard statistics
    #compute climatological statistics for each month 
    sit_df_list_zrif = []
    means_list_zrif = []
    mins_list_zrif = []
    medians_list_zrif = []
    pers_list_zrif = []
    stds_list_zrif = []

    for i in range(0,12):
        ds = sea_ice_thickness_list[i]
        df_list = [] ;  means = [] ; mins = [] ; medians = [] ; pers = [] ; stds = []
        for j in range(0,6):
            df = ds[j].rename('sea_ice_thickness'); df = df.to_dataframe(); df = df['sea_ice_thickness']
            df_list.append(df.dropna())
            means.append(round(np.nanmean(ds[j]),2))
            mins.append(round(np.nanmin(ds[j]),2))
            medians.append(round(np.nanmedian(ds[j]),2))
            stds.append(round(np.nanstd(ds[j]),2))
            pers.append(round(np.nanpercentile(ds[j], 95),2))
        sit_df_list_zrif.append(df_list) ; means_list_zrif.append(means); mins_list_zrif.append(mins)
        medians_list_zrif.append(medians) ; pers_list_zrif.append(pers) ; stds_list_zrif.append(stds)

    #rearrange data 
    KH = []
    WS = []
    EA = []
    AB = []
    RS = []
    SH = []

    # Rearrange the data
    for i in range(0,12):
        AB.append(means_list_zrif[i][0])
        WS.append(means_list_zrif[i][1])
        KH.append(means_list_zrif[i][2])
        EA.append(means_list_zrif[i][3])
        RS.append(means_list_zrif[i][4])
        SH.append(means_list_zrif[i][5])

    # Final rearranged list
    means_list_zrif = [AB,WS,KH,EA,RS,SH]

    #sicc regions
    #slice thickness data into regions
    sit_sicc_list = []

    for sit_ds in yearly_sicc_uncert_list[year_value]:
        ds_list = []
        for i in range(0,7):
            ds =  (sit_ds.where((sit_ds.lons >region[i][0])&(sit_ds.lons <region[i][1])))
            ds_list.append(ds)
        sit_sicc_list.append(ds_list) 

    RS = []
    for i in range(0,12):
        rs_ds = xr.merge([sit_sicc_list[i][4],sit_sicc_list[i][5]])
        RS.append(rs_ds.uncertainity)
        
    for i in range(0,12):
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

    for i in range(0,12):
        ds = sit_sicc_list[i]
        df_list = [] ;  means = [] ; mins = [] ; medians = [] ; pers = [] ; stds = []
        for j in range(0,6):
            df = ds[j].rename('Sea-ice thickness'); df = df.to_dataframe(); df = df['Sea-ice thickness']
            df_list.append(df.dropna())
            means.append(round(np.nanmean(ds[j]),2))
            mins.append(round(np.nanmin(ds[j]),2))
            medians.append(round(np.nanmedian(ds[j]),2))
            stds.append(round(np.nanstd(ds[j]),2))
            pers.append(round(np.nanpercentile(ds[j], 95),2))
        sit_df_list_sicc.append(df_list) ; means_list_sicc.append(means); mins_list_sicc.append(mins)
        medians_list_sicc.append(medians) ; pers_list_sicc.append(pers) ; stds_list_sicc.append(stds)

     
    #rearrange data 
    KH = []
    WS = []
    EA = []
    AB = []
    RS = []
    SH = []

    # Rearrange the data
    for i in range(0,12):
        AB.append(means_list_sicc[i][0])
        WS.append(means_list_sicc[i][1])
        KH.append(means_list_sicc[i][2])
        EA.append(means_list_sicc[i][3])
        RS.append(means_list_sicc[i][4])
        SH.append(means_list_sicc[i][5])

    # Final rearranged list
    means_list_sicc = [AB,WS,KH,EA,RS,SH]
    #fdm regions
    #slice thickness data into regions
    sit_fdm_list = []

    for sit_ds in yearly_fdm_uncert_list[year_value]:
        ds_list = []
        for i in range(0,7):
            ds =  (sit_ds.where((sit_ds.lons >region[i][0])&(sit_ds.lons <region[i][1])))
            ds_list.append(ds)
        sit_fdm_list.append(ds_list) 

    RS = []
    for i in range(0,12):
        rs_ds = xr.merge([sit_fdm_list[i][4],sit_fdm_list[i][5]])
        RS.append(rs_ds.uncertainity)
        
    for i in range(0,12):
        sit_fdm_list[i].pop(4)
        sit_fdm_list[i].pop(4)
        sit_fdm_list[i].insert(4,RS[i])

    #compute climatological statistics for each month 
    sit_df_list_fdm = []
    means_list_fdm = []
    mins_list_fdm = []
    medians_list_fdm = []
    pers_list_fdm = []
    stds_list_fdm = []

    for i in range(0,12):
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

     
    #rearrange data 
    KH = []
    WS = []
    EA = []
    AB = []
    RS = []
    SH = []

    # Rearrange the data
    for i in range(0,12):
        AB.append(means_list_fdm[i][0])
        WS.append(means_list_fdm[i][1])
        KH.append(means_list_fdm[i][2])
        EA.append(means_list_fdm[i][3])
        RS.append(means_list_fdm[i][4])
        SH.append(means_list_fdm[i][5])

    # Final rearranged list
    means_list_fdm = [AB,WS,KH,EA,RS,SH]
    return (means_list_zrif,means_list_xoc,means_list_boc,means_list_OLM,means_list_sicc,means_list_fdm)