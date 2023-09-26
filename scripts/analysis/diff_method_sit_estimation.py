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

#define functions to be used later
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


month_names = ['january','february','march','april','may','june','july','august','september','october','november','december']

df1_list = []
df2_list = []
joined_list = []
dates_list = []

os.chdir(processed_data_path+'/amsr_weighted_cs2_freeboards')
filenames = []
for i in range(0,12):
    filenames.append(month_names[i]+'_cs2_amsr_weighted_freeboard.nc')
    
for i in range(0,12):
    if  i == 0:
        ds1 = xr.open_mfdataset([filenames[11],filenames[0],filenames[1]])
        ds1 = ds1.freeboard.rename('CS2_freeboard').to_dataframe().dropna().reset_index().drop(columns=['x','y'])
        df1_list.append(ds1)
    if  i == 11:
        ds1 = xr.open_mfdataset([filenames[10],filenames[11],filenames[0]])
        ds1 = ds1.freeboard.rename('CS2_freeboard').to_dataframe().dropna().reset_index().drop(columns=['x','y'])
        df1_list.append(ds1)
    if 0 < i < 11:
        ds1 = xr.open_mfdataset([filenames[i-1],filenames[i],filenames[i+1]])
        ds1 = ds1.freeboard.rename('CS2_freeboard').to_dataframe().dropna().reset_index().drop(columns=['x','y'])
        df1_list.append(ds1)

os.chdir(processed_data_path+'/icesat_2_freeboards/2019/')
for i in range(0,12):
    ds2 = xr.open_dataset('2019_'+str(month_names[i])+'_IS-2_daily_gridded_freeboard.nc')
    ds2 = ds2.freeboard.rename('IS2_freeboard').to_dataframe().dropna().reset_index().drop(columns=['x','y'])
    ds2 = ds2.rename(columns={"lons" : "lon","lats" : "lat", 'time' : 'IS2_time'})
    df2_list.append(ds2)
    

rhow = 1024 #kg.m^-3
rhoi = 917 #kg.m^-3
rhos = 320 #kg.m^-3
eta = (1 + (0.51*rhos)/1000)**(1.5)

for i in range(0,12):
    closest_grids = nearest_neighbor(df2_list[i], df1_list[i], return_dist=True)
    df_joined = closest_grids.join(df2_list[i])
    df_joined = df_joined.where(df_joined.distance < 36000.0)
    df_joined['snow_depth'] = (df_joined['IS2_freeboard'] - df_joined['CS2_freeboard'])/eta
    df_joined['time_separation'] = (df_joined['IS2_time'] - df_joined['time']).dt.days
    df_joined = df_joined.where(df_joined.time_separation <= 10)
    df_joined = df_joined.where(df_joined.time_separation >= -10)
    df_joined = df_joined.dropna()
    df_joined['thickness'] = (rhow/(rhow-rhoi))*df_joined['IS2_freeboard'] + ((rhos-rhow)/(rhow-rhoi))*df_joined['snow_depth']
    dates = df_joined['IS2_time'].groupby(df_joined['IS2_time'].dt.floor('d')).size().reset_index(name='count')
    dates_list.append(dates)
    joined_list.append(df_joined)


flag_list = []
for x in range(0,12):
    x_r = []
    dsd = np.array(joined_list[x]['distance'])
    dse = np.array(joined_list[x]['time_separation'])
    for i in range(0,(joined_list[x]['distance'].size)):
        distance = dsd[i]
        time_sep = dse[i]
        if (distance ==0) & (-1 <= time_sep <= 1):
                x_r.append(0)
        if (distance ==0) & (1 < time_sep <= 10):
                x_r.append(1)
        if (distance ==0) & (-10 <= time_sep < -1):
                x_r.append(1)
        if (distance ==0) & ( 10 < time_sep <= 15):
                x_r.append(2)
        if (distance ==0) & ( -15 <= time_sep < -10):
                x_r.append(2)
        if (distance >=0) & (-1 <= time_sep <= 1):
                x_r.append(3)
        if (distance >=0) & (1 < time_sep <= 10):
                x_r.append(4)
        if (distance >=0) & (-10 <= time_sep < -1):
                x_r.append(4)
        if (distance >=0) & (10 < time_sep <= 15):
                x_r.append(5)
        if (distance >=0) & (-15 <= time_sep < -10):
                x_r.append(5)      
    flag_list.append(x_r)    


ds_list = []
for i in range(0,12):
    dsd_list = []
    df = joined_list[i]
    df['flags'] = flag_list[i]
    dates = dates_list[i]
    for x in range(0,(dates['IS2_time'].size)):
        daydata = df[(df['IS2_time'] == str(dates['IS2_time'][x]))]
        x10, y10=mapProj(daydata['lon'], daydata['lat'])
        varS = bin_data(x10, y10, daydata['snow_depth'], xptsG, yptsG, dx)  #snow depth variable
        varT = bin_data(x10, y10, daydata['thickness'], xptsG, yptsG, dx)   #thickness variable
        varFg = bin_data(x10, y10, daydata['flags'], xptsG, yptsG, dx)      #flags
        #expand dimension to add day
        sn_depth = np.expand_dims(varS,2)
        fr1 = np.expand_dims(varT,2)
        date = pd.date_range(dates['IS2_time'][x],periods=1)
        flgs = np.expand_dims(varFg,2)
        #create dataset
        ds = xr.Dataset(
            data_vars=dict(
                sea_ice_thickness=(["x", "y","time"], fr1),
                snow_depth=(["x", "y","time"], sn_depth),
                flags=(["x", "y","time"], flgs),
                ),
            coords=dict(
                lons=(["x", "y"], lonG),
                lats=(["x", "y"], latG),
                time = date,
                ),
            attrs=dict(description='CS-2 2019 gridded data',comment='data gridded on 25 km NSIDC South stereographic grid'),
            )
        dsd_list.append(ds)    
    mds = xr.merge(dsd_list)
    ds_list.append(mds)    


os.chdir(processed_data_path+'diff_method_sit')

for i in range(0,12):
    ds = ds_list[i]
    ds.to_netcdf(str(month_names[i])+'_circum_sit_from_fr_diff_approach_25km.nc')


















