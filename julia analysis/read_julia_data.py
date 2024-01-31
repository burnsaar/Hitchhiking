# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 08:56:00 2023

@author: Aaron
"""

import numpy as np
import pandas as pd
import glob
import pickle
import re
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Point
import SF_shapefiles as SF

def data_paths(fp):
    files = glob.glob(fp)
    return files

def load_d_and_s(files):
    df_compiled = pd.DataFrame()
    iteration = 1
    for file in files:
        df = pd.read_csv(file, delimiter='\t', header=None, names=['Lat', 'Lon'])
        df['Iteration'] = str(iteration)
        iteration += 1
        # with open(file, 'rb') as temp:
        #     res.append(pickle.load(temp))
        #     temp.close() 
        df_compiled = pd.concat([df_compiled, df])
    return df_compiled

def load_trips(files):
    df_ls = []
    iteration = 1
    for file in files:
        df = pd.read_csv(file, delimiter='\t', header=None, names=['vertex', 'full time', '1-way time', 'Lat', 'Lon'])
        df['Iteration'] = str(iteration)
        iteration += 1
        # with open(file, 'rb') as temp:
        #     res.append(pickle.load(temp))
        #     temp.close() 
        df_ls.append(df)
    return df_ls

def compile_depots(depots_df):
    df = depots_df #.drop_duplicates(subset=['Lat']) #subset=['Lat', 'Lon']
    geometry = [Point(x,y) for x,y in zip(df['Lon'], df['Lat'])]
    geo_df = gpd.GeoDataFrame(df, crs="EPSG:4326", geometry=geometry)
    
    return df, geo_df

def compile_trip_time(trips):
    #setup a list to store each of the series in
    ser_ls = []
    for trip in trips:
        #find the index of the delivery vertex
        site_idx = trip.index[trip['vertex'].str.contains(pat = 's')]
        #pull the row based on the delivery index and add to the compiled dataframe
        ser = trip.loc[site_idx[0]]

        ser_ls.append(ser)
        
    trip_time_df = pd.DataFrame(ser_ls)
    trip_time_df.reset_index(inplace = True)
    trip_time_df.rename(columns = {'Iteration': 'trip', 'index': 'Vertex Index'}, inplace = True)
    trip_time_df['full time (min)'] = trip_time_df['full time'] / 60
    
    #create GeoDataFrame as well
    geometry = [Point(x,y) for x,y in zip(trip_time_df['Lon'], trip_time_df['Lat'])]
    geo_df = gpd.GeoDataFrame(trip_time_df, crs="EPSG:4326", geometry=geometry)
        
    return trip_time_df, geo_df





if __name__ == '__main__':
    
    #set the results parent folder
    file_path = 'C:/Users/Aaron/AppData/Local/Programs/Julia-1.6.7/MultiAgentAllocationTransit.jl/results/2024-01-30'
    
    #load the depot locations
    d_file_paths = data_paths(file_path + '/depots/*.dat')
    depots = load_d_and_s(d_file_paths)
    
    #compile the depot locations
    depots_df, depots_gdf = compile_depots(depots)
    
    #load the delivery site locations
    s_file_paths = data_paths(file_path + '/sites/*.dat')
    sites = load_d_and_s(s_file_paths)
    
    #load the states for each trip, t
    t_file_path = data_paths(file_path + '/states/*.dat')
    trips = load_trips(t_file_path)
    
    #compile each of the travel time from the depot to the site
    trip_time_df, trips_gdf = compile_trip_time(trips)
    
    
    
    
    # #plotting
    # #point to the map files for plotting
    # muni_file_path = 'C:/Users/Aaron/Documents/GitHub/Hitchhiking/Hitchhiking/muni_routes/geo_export_6dee9e27-b549-4312-94b4-5ad68950d5fe.shp'
    # SF_boundary_file_path = 'C:/Users/Aaron/Documents/GitHub/Hitchhiking/Hitchhiking/SF Boundary/s7d02x.shp'
    # SF_Zoning_file_path = 'C:/Users/Aaron/Documents/GitHub/Hitchhiking/Hitchhiking/SF Zoning/geo_export_486f35c6-390f-4739-8f97-6171581968e5.shp'
    
    # #load the map files
    # muni = SF.load_muni(muni_file_path)
    # SF_boundary = SF.load_SF_boundaries(SF_boundary_file_path)
    # SF_zoning = SF.load_SF_zoning(SF_Zoning_file_path)
    
    # #load the bus data
    # bus_only = SF.muni_bus_only(muni)
    # bus_freq, bus_freq_full = SF.add_bus_freq(bus_only)
    
    
    # SF.plot_muni(muni, SF_boundary, SF_zoning)
    # SF.plot_muni_freq(bus_freq, SF_boundary, SF_zoning)
    
    # #plot the travel time graphic
    # fig, ax = init_fig(figsize=(10,10))
    # fig, ax = plot_boundary(SF_boundary)
    # fig, ax = plot_zoning(SF_zoning)
    # fix, ax = plot_bus_routes(fig, ax, bus_freq)
    # fig, ax = plot_travel_time(fig, ax, trips_gdf)
    # fig, ax = plot_depots(fig, ax, depots_gdf)
    # #fig, ax = plot_drone_rad(fig, ax, depots_gdf)
    
    
    
    # fig, ax = plt.subplots(figsize=(10,10))
    # #plot SF boundary
    
    # #plot depots
    
    # #plot SF bus routes?
    
    # geometry = [Point(x,y) for x,y in zip(trip_time_df['Lon'], trip_time_df['Lat'])]
    # geo_df = gpd.GeoDataFrame(trip_time_df, crs="EPSG:4326", geometry=geometry)
    # #geo_df.plot(column = 'full time')
    # #plt.plot(geometry)
    # geo_df.plot(column='full time (min)', 
    #             ax=ax, 
    #             legend=True,
    #             legend_kwds={'label':'Trip Time (minutes)'})
    # plt.legend()



# df_actions = pd.read_csv("C:/Users/Aaron/Documents/GitHub/Hitchhiking/Hitchhiking/julia output/actions.txt", 
#                  delimiter = "\t", header=None, names = ['raw_action', 'distance (km)']) 

# df_states = pd.read_csv("C:/Users/Aaron/Documents/GitHub/Hitchhiking/Hitchhiking/julia output/states.txt", 
#                  delimiter = "\t", header=None, names = ['vertex', 'time', 'latlon']) 






# f = open("C:/Users/Aaron/Documents/GitHub/Hitchhiking/Hitchhiking/julia output/test.txt", "r")
# print(f.read())

# f = open("C:/Users/Aaron/Documents/GitHub/Hitchhiking/Hitchhiking/julia output/test_2.txt", "r")
# print(f.read())

# df = pd.read_fwf("C:/Users/Aaron/Documents/GitHub/Hitchhiking/Hitchhiking/julia output/test.txt")



#df2 = pd.read_fwf("C:/Users/Aaron/Documents/GitHub/Hitchhiking/Hitchhiking/julia output/states.txt", delimter = "\t", header=None) 
