# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 09:24:02 2024

@author: Aaron
"""


import pandas as pd
import glob
import pickle
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Point
import SF_shapefiles as SF
import utils




plt.rcParams['figure.dpi'] = 500

#point to the map files for plotting
muni_file_path = 'C:/Users/Aaron/Documents/GitHub/Hitchhiking/Hitchhiking/muni_routes/geo_export_6dee9e27-b549-4312-94b4-5ad68950d5fe.shp'
SF_boundary_file_path = 'C:/Users/Aaron/Documents/GitHub/Hitchhiking/Hitchhiking/SF Boundary/s7d02x.shp'
SF_Zoning_file_path = 'C:/Users/Aaron/Documents/GitHub/Hitchhiking/Hitchhiking/SF Zoning/geo_export_486f35c6-390f-4739-8f97-6171581968e5.shp'

#load the map files
#muni = SF.load_muni(muni_file_path)
SF_boundary = SF.load_SF_boundaries(SF_boundary_file_path)
#SF_boundary_main = gpd.GeoDataFrame(SF_boundary.iloc[0])
SF_zoning = SF.load_SF_zoning(SF_Zoning_file_path)



#set up Depot geodataframes
depot = {'depot': [1]}
df = pd.DataFrame(depot)
geometry = [Point(-122.4012, 37.744045)]
depot_config_parcel_gdf = gpd.GeoDataFrame(df, crs="EPSG:4326", geometry=geometry)

#depot_config_pharm
depot = {'Store': ['Walgreens', 'Walgreens', 'CVS', 'Walgreens']}
df = pd.DataFrame(depot)
geometry = [Point(-122.42371, 37.78675), Point(-122.42225, 37.74267), Point(-122.47634, 37.72684), Point(-122.47597, 37.78084)]
depot_config_pharm_gdf = gpd.GeoDataFrame(df, crs="EPSG:4326", geometry=geometry)

#depot_config_food
depot = {'depot': [1]}
df = pd.DataFrame(depot)
geometry = [SF_boundary['geometry'][0].centroid]
depot_config_food_gdf = gpd.GeoDataFrame(df, crs="EPSG:4326", geometry=geometry)

#plot the hitchhiking drone travel time graphic
fig, ax = utils.init_fig(figsize=(10,10))
fig, ax = utils.plot_boundary(fig, ax, SF_boundary)
fig, ax, handles, previous_labels = utils.plot_zoning(fig, ax, SF_zoning)
fig, ax = utils.plot_depots(fig, ax, depot_config_parcel_gdf, handles, previous_labels, marker='p', color='darkorange', legend_kwds={'label':'Parcel Facility'})
# fig, ax = utils.plot_depots(fig, ax, depot_config_pharm_gdf, color='red', legend_kwds={'label':'Pharma'})
# fig, ax = utils.plot_depots(fig, ax, depot_config_food_gdf, color='green', legend_kwds={'label': 'Restaurant'})