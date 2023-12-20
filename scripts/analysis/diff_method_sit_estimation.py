# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 00:02:27 2023

@author: MNGMA
"""

#script to estimation sea ice thickness with the freeboard differencing method
import os
processed_data_path = 'data/processed_data/'

#Make other necessary imports 
import xarray as xr 
import pandas as pd
import numpy as np
from sklearn.neighbors import KDTree
import pyproj

#define functions to be used for finding the nearest neighbours and gridding the data
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
    #print(xbins2.shape)
    #print(ybins2.shape)
    counts, xedges, yedges = np.histogram2d(xpts, ypts,bins=(xbins2, ybins2))
    z, _, _ = np.histogram2d(xpts, ypts,bins=(xbins2, ybins2), weights=var)
    varG = z / counts
    
    # Need to re-flip the arrays then transpose because of how histogram2d works across columns then rows.
    varG=np.flipud(varG.T)
    counts=np.flipud(counts.T)
    
    return varG

# Define a projection 
mapProj = pyproj.Proj("+init=EPSG:3976")
dx = 25000

def get_nearest(src_points, candidates, k_neighbors=1):
    """Find nearest neighbors for all source points from a set of candidate points"""

    # Create tree from the candidate points
    tree = KDTree(candidates, leaf_size=15, metric='euclidean')

    # Find closest points and distances
    distances, indices = tree.query(src_points, k=k_neighbors)

    # Transpose to get distances and indices into arrays
    distances = distances.transpose()
    indices = indices.transpose()

    # Get closest indices and distances (i.e. array at index 0)
    # note: for the second closest points, you would take index 1, etc.
    closest = indices[0]
    closest_dist = distances[0]

    # Return indices and distances
    return (closest, closest_dist)


def nearest_neighbor(left_gdf, right_gdf, return_dist=False):
    """
    For each point in left_gdf, find closest point in right GeoDataFrame and return them.
    """

    # Ensure that index in right gdf is formed of sequential numbers
    right = right_gdf.copy().reset_index(drop=True)
    left =  left_gdf.copy().reset_index(drop=True)
    # Parse coordinates from points and insert them into a numpy array as RADIANS
    left_xy = mapProj(left['lon'],left['lat'])
    right_xy = mapProj(right['lons'],right['lats'])
    
    left_radians = np.column_stack([list(left_xy[0].flatten()),list(left_xy[1].flatten())])
    right_radians = np.column_stack([list(right_xy[0].flatten()),list(right_xy[1].flatten())])

    # Find the nearest points
    # -----------------------
    # closest ==> index in right_gdf that corresponds to the closest point
    # dist ==> distance between the nearest neighbors (in meters)

    closest, dist = get_nearest(src_points=left_radians, candidates=right_radians)

    # Return points from right GeoDataFrame that are closest to points in left GeoDataFrame
    closest_points = right.loc[closest]

    # Ensure that the index corresponds the one in left_gdf
    closest_points = closest_points.reset_index(drop=True)

    # Add distance if requested
    if return_dist:
        
        closest_points['distance'] = dist 

    return closest_points

#define month names and empty lists to hold the monthly data
month_names = ['january','february','march','april','may','june','july','august','september','october','november','december']

df1_list = []   # this list holds the cs2 freeboards
df2_list = []   # this list holds the is2 freeboards
joined_list = [] # this list holds estimated sea ice thickness data 
dates_list = []

os.chdir('data/processed_data/amsr_weighted_cs2_freeboards/')
for i in range(0,12):
    if  i == 0:
        ds1 = xr.open_dataset(month_names[i]+'_cs2_amsr_weighted_SAR_SIN_freeboard.nc')
        ds1 = ds1.to_dataframe()
        ds1.rename(columns={'freeboard':'CS2_freeboard'},inplace= True) ;
        ds1 = ds1.dropna()
        extend = xr.open_dataset(str('february')+ '_cs2_amsr_weighted_SAR_SIN_freeboard.nc')
        extend = extend.to_dataframe()
        extend.rename(columns={'freeboard':'CS2_freeboard'},inplace= True) ; extend = extend.dropna()
        extend_1 = xr.open_dataset(str('december')+ '_cs2_amsr_weighted_SAR_SIN_freeboard.nc')
        extend_1 = extend_1.to_dataframe()
        extend_1.rename(columns={'freeboard':'CS2_freeboard'},inplace= True) ; extend_1 = extend_1.dropna()
        ds1 = pd.concat([ds1,extend,extend_1]).reset_index().drop(columns=['x','y'])
        df1_list.append(ds1)
    if  i == 11:
        ds1 = xr.open_dataset(month_names[i]+'_cs2_amsr_weighted_SAR_SIN_freeboard.nc')
        ds1 = ds1.to_dataframe()
        ds1.rename(columns={'freeboard':'CS2_freeboard'},inplace= True) ; ds1 = ds1.dropna()
        extend = xr.open_dataset(str('january')+ '_cs2_amsr_weighted_SAR_SIN_freeboard.nc')
        extend = extend.to_dataframe()
        extend.rename(columns={'freeboard':'CS2_freeboard'},inplace= True) ; extend = extend.dropna()
        extend_1 = xr.open_dataset(str('november')+ '_cs2_amsr_weighted_SAR_SIN_freeboard.nc')
        extend_1 = extend_1.to_dataframe()
        extend_1.rename(columns={'freeboard':'CS2_freeboard'},inplace= True) ; extend_1 = extend_1.dropna()
        ds1 = pd.concat([ds1,extend,extend_1]).reset_index().drop(columns=['x','y'])
        df1_list.append(ds1)
    if 0 < i < 11:
        ds1 = xr.open_dataset(month_names[i]+'_cs2_amsr_weighted_SAR_SIN_freeboard.nc')
        ds1 = ds1.to_dataframe()
        ds1.rename(columns={'freeboard':'CS2_freeboard'},inplace= True) ; ds1 = ds1.dropna()
        extend = xr.open_dataset(str(month_names[i-1])+ '_cs2_amsr_weighted_SAR_SIN_freeboard.nc')
        extend = extend.to_dataframe()
        extend.rename(columns={'freeboard':'CS2_freeboard'},inplace= True) ; extend = extend.dropna()
        extend_1 = xr.open_dataset(str(month_names[i+1])+ '_cs2_amsr_weighted_SAR_SIN_freeboard.nc')
        extend_1 = extend_1.to_dataframe()
        extend_1.rename(columns={'freeboard':'CS2_freeboard'},inplace= True) ; extend_1 = extend_1.dropna()
        ds1 = pd.concat([ds1,extend,extend_1]).reset_index().drop(columns=['x','y'])
        df1_list.append(ds1)

os.chdir('data/processed_data/OLMi_sea_ice_thickness/')
for i in range(0,12):
    ds2 = xr.open_dataset('2019_'+str(month_names[i])+'_IS-2_gridded_daily_sit_from_OLM_and_OLMi.nc')
    ds2 = ds2.to_dataframe().dropna().reset_index().drop(columns=['x','y'])
    ds2.rename(columns={'freeboard':'IS2_freeboard'},inplace= True)
    ds2 = ds2.rename(columns={"lons" : "lon","lats" : "lat", 'time' : 'IS2_time'})
    df2_list.append(ds2)
    

rhow = 1024 # sea water density in kg/m^3
rhoi = 917 # sea ice density in kg/m^3
rhos = 320 # snow density in kg/m^3

eta = (1 + (0.51*rhos)/1000)**(1.5)  # refractive index

ierror = 15 # sea ice density error in kg/m^3
serror = 100 # snow density error in kg/m^3
werror = 3 # sea water density error in kg/m^3


for i in range(0,12):
    closest_grids = nearest_neighbor(df2_list[i], df1_list[i], return_dist=True)
    df_joined = closest_grids.join(df2_list[i])
    df_joined = df_joined.where(df_joined.distance < 36000.0)   # limit separation in space to the eight neighbouring grid boxes
    df_joined['snow_depth'] = (df_joined['IS2_freeboard'] - df_joined['CS2_freeboard'])/eta   # estimate the snow depth 
    df_joined['time_separation'] = (df_joined['IS2_time'] - df_joined['time']).dt.days
    df_joined = df_joined.where(df_joined.time_separation <= 10)  # limit time separation to 10 days
    df_joined = df_joined.where(df_joined.time_separation >= -10) # limit time separation to 10 days
    df_joined = df_joined.dropna()
    df_joined['thickness'] = (rhow/(rhow-rhoi))*df_joined['IS2_freeboard'] + ((rhos-rhow)/(rhow-rhoi))*df_joined['snow_depth']  # estimate the ice thickness
    # below is the uncertainity estimation (see S1 for details)
    df_joined['varfi'] = (rhow-rhos)/(eta*(rhow-rhoi))
    df_joined['varf'] = (rhos+(eta-1)*rhow)/(eta*(rhow-rhoi))
    df_joined['varrhow'] = ((df_joined['CS2_freeboard']*(rhos-rhoi))-(df_joined['IS2_freeboard']*((eta-1)*rhoi+rhos)))/(eta*(rhow-rhoi)**(2))
    df_joined['varrhos'] = (df_joined['IS2_freeboard'] - df_joined['CS2_freeboard'])/(eta*(rhow-rhoi))
    df_joined['varrhoi'] = ((df_joined['CS2_freeboard']*(rhow-rhos))+(df_joined['IS2_freeboard']*(rhos+(eta-1)*rhow)))/(eta*(rhow-rhoi)**(2))
    df_joined['errorfi'] = (df_joined['varfi']**2)*(df_joined['standard_dev']**2)
    df_joined['errorf'] = (df_joined['varf']**2)*(df_joined['fr_standard_deviation']**2)
    df_joined['errorrhow'] = (df_joined['varrhow']**2)*(werror**2)
    df_joined['errorrhos'] = (df_joined['varrhos']**2)*(serror**2)
    df_joined['errorrhoi'] = (df_joined['varrhoi']**2)*(ierror**2)
    df_joined['variance'] = df_joined['errorfi']+df_joined['errorf']+df_joined['errorrhow']+df_joined['errorrhos']+df_joined['errorrhoi']
    df_joined['uncertainity'] = df_joined['variance']**(0.5)
    dates = df_joined['IS2_time'].groupby(df_joined['IS2_time'].dt.floor('d')).size().reset_index(name='count')
    dates_list.append(dates)
    joined_list.append(df_joined)



ds_list = []
for i in range(0,12):
    dsd_list = []
    df = joined_list[i]
    dates = dates_list[i]
    for x in range(0,(dates['IS2_time'].size)):
        daydata = df[(df['IS2_time'] == str(dates['IS2_time'][x]))]
        x10, y10=mapProj(daydata['lon'], daydata['lat'])
        varU=bin_data(x10, y10, daydata['uncertainity'], xptsG, yptsG, dx)  #uncertainity variable
        varS = bin_data(x10, y10, daydata['snow_depth'], xptsG, yptsG, dx)  #snow depth variable
        varT = bin_data(x10, y10, daydata['thickness'], xptsG, yptsG, dx)   #thickness variable
        #expand dimension to add day
        sn_depth = np.expand_dims(varS,2)
        sit = np.expand_dims(varT,2)
        date = pd.date_range(dates['IS2_time'][x],periods=1)
        uncert = np.expand_dims(varU,2)
        #create dataset
        ds = xr.Dataset(
            data_vars=dict(
                sea_ice_thickness=(["x", "y","time"], sit),
                snow_depth=(["x", "y","time"], sn_depth),
                uncertainity=(["x", "y","time"], uncert),
                ),
            coords=dict(
                lons=(["x", "y"], lonG),
                lats=(["x", "y"], latG),
                time = date,
                ),
            attrs=dict(description='sea_ice_thickness estimated with the freeboard differencing method',comment='data gridded on 25 km NSIDC South stereographic grid'),
            )
        dsd_list.append(ds)    
    mds = xr.merge(dsd_list)
    ds_list.append(mds)    


#quick look at a monthly map
ds1 = ds_list[3].mean('time')
sit = ds1.sea_ice_thickness
plt.pcolormesh(sit.lons,sit.lats,sit); plt.colorbar()

os.chdir('data/processed_data/diff_method_sit/')

for i in range(0,12):
    ds = ds_list[i]
    ds.to_netcdf(str(month_names[i])+'_circum_sit_from_fr_diff_approach_25km.nc')


















