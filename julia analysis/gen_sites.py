# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 06:27:28 2024

@author: Aaron
"""


import matplotlib.pyplot as plt
from shapely.geometry import Point
import SF_shapefiles as SF
import geopandas as gpd
import utils
import numpy as np
from datetime import datetime, date
import os
import pickle
import pandas as pd


plt.rcParams['figure.dpi'] = 500
SF_boundary_file_path = 'C:/Users/Aaron/Documents/GitHub/Hitchhiking/Hitchhiking/SF Boundary/s7d02x.shp'

SF_boundary = SF.load_SF_boundaries(SF_boundary_file_path)

fig, ax = utils.init_fig(figsize=(10,10))
fig, ax = utils.plot_boundary(fig, ax, SF_boundary)


lat_start=37.707040
lat_end=37.80913
lon_start=-122.518
lon_end=-122.350



total_sites = 10000
site_ls = []
lat_ls = []
lon_ls = []
geometry_ls = []
i=1

while i < total_sites+1:

    #gen a random site
    rnd_site_lat=np.random.uniform(lat_start, lat_end)
    rnd_site_lon=np.random.uniform(lon_start, lon_end)
    
    if (rnd_site_lat>37.80) & (rnd_site_lon>-122.375): #remove any points on treasure island
        continue

    site = Point(rnd_site_lon, rnd_site_lat)

    #check to see if the site is within the 7by7 boundary
    print(SF_boundary.contains(site)[0])
    if SF_boundary.contains(site)[0] == True:
        
        site_ls.append(i)
        lat_ls.append(rnd_site_lat)
        lon_ls.append(rnd_site_lon)
        geometry_ls.append(site)
        print(i)
        i+=1
    


d = {'Site': site_ls, 'Lat': lat_ls, 'Lon': lon_ls, 'geometry': geometry_ls}

site_gdf = gpd.GeoDataFrame(d, crs='EPSG:4326')

site_gdf.plot(ax=ax)


#convert geodataframe to dataframe
site_df = pd.DataFrame(site_gdf[['Site', 'Lat', 'Lon']])

#write the geodataframe to a .dat file    
current_date = date.today().strftime('%Y-%m-%d') 
base_path = 'C:/Users/Aaron/Documents/GitHub/Hitchhiking/Hitchhiking/'
current_date = current_date + '_sites'

# Create a folder with the formatted date if it doesn't exist
folder_path = os.path.join(base_path, current_date)
if not os.path.exists(folder_path):
    os.makedirs(folder_path)


# Define the file path for saving
saveFile = os.path.join(folder_path, 'sites.txt')

# with open(saveFile, 'wb') as file:
#     pickle.dump(site_df, file)

site_df.to_csv(saveFile, header=True, index=False, sep=' ')
