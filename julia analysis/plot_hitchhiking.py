# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 15:27:33 2024

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
import read_julia_data as read_jl
#import distance_matrix_API as API


def init_fig(figsize):
    fig, ax = plt.subplots(figsize = figsize)
    plt.xlim([-122.52480, -122.33907])
    plt.ylim([37.70, 37.84])
    plt.xlabel('longitude')
    plt.ylabel('latitude')
    plt.legend()
    return fig, ax

def plot_boundary(SF_boundary):
    SF_boundary.boundary.plot(ax=ax, edgecolor = 'k')    
    return fig, ax

def plot_zoning(SF_zoning):
    SF_zoning.plot(ax=ax, 
                   legend = True,
                   edgecolor='k')
    return fig, ax

def plot_travel_time(fig, ax, gdf):
    gdf.plot(column='full time (min)', 
                ax=ax,
                cmap='plasma',
                vmax=60,
                legend=True,
                legend_kwds={'label':'Trip Time (minutes)'}) #, 'orientation': 'horizontal'
    plt.title('Hitchhiking Isochrone')
    return fig, ax

def plot_depots(fig, ax, gdf):
    gdf.plot(ax=ax,
             legend=True,
             markersize=1000,
             color='aqua',
             edgecolor='k',
             marker='*',
             legend_kwds={'label':'Depot'})
    return fig, ax

def plot_bus_routes(fig, ax, gdf):
    gdf.plot(ax=ax,
             color='slategrey',
             alpha=0.5)
    return fig, ax

def plot_drone_rad(fig, ax, gdf): #TODO:need to implement, not working yet
    for i in range(len(gdf)):
        radius = gdf.iloc[i]['geometry'].buffer(250).boundary
        radius = gpd.GeoSeries(radius)
        radius.plot(ax=ax,
                    color='green')
    
    return fig, ax

def plot_driving_time(fig, ax, gdf):
    gdf_driving = gdf[gdf['mode'] == 'driving']
    gdf_driving.plot(column='duration (min)', 
                ax=ax,
                cmap='plasma',
                vmax=60,
                legend=True,
                legend_kwds={'label':'Driving Time (minutes)'})
    plt.title('Driving Isochrone')
    
    return fig, ax

def plot_bicycling_time(fig, ax, gdf):
    gdf_bicycling = gdf[gdf['mode'] == 'bicycling']
    gdf_bicycling.plot(column='duration (min)', 
                ax=ax,
                cmap='plasma',
                vmax=60,
                legend=True,
                legend_kwds={'label':'Bicycling Time (minutes)'})
    plt.title('Bicycling Isochrone')
    
    return fig, ax

def plot_delta_driving(fig, ax, gdf):
    gdf_driving = gdf[gdf['mode']=='driving']
    gdf_driving.plot(column='delta_to_hitch (perc)', 
                ax=ax,
                cmap='plasma',
                vmin=-200,
                vmax=100,
                legend=True,
                legend_kwds={'label':'Percent Difference between Hitchhiking and Driving Time (minutes)'})
    plt.title('Percent Difference Hitchhiking and Driving Isochrone')
    
    return fig, ax

def plot_delta_bicycling(fig, ax, gdf):
    gdf_bicycling = gdf[gdf['mode']=='bicycling']
    gdf_bicycling.plot(column='delta_to_hitch (perc)', 
                ax=ax,
                cmap='plasma',
                vmin=-200,
                vmax=100,
                legend=True,
                legend_kwds={'label':'Percent Difference between Hitchhiking and Bicycling Time (minutes)'})
    plt.title('Percent Difference Hitchhiking and Bicycling Isochrone')
    
    return fig, ax




if __name__ == '__main__':

    #----------------------------Read in the Julia Data------------------------
    #set the results parent folder
    file_path = 'C:/Users/Aaron/AppData/Local/Programs/Julia-1.6.7/MultiAgentAllocationTransit.jl/results/2024-01-30 (d_100_s_100_iter_20)'
    
    #load the depot locations
    d_file_paths = read_jl.data_paths(file_path + '/depots/*.dat')
    depots = read_jl.load_d_and_s(d_file_paths)
    
    #compile the depot locations
    depots_df, depots_gdf = read_jl.compile_depots(depots)
    
    #load the delivery site locations
    s_file_paths = read_jl.data_paths(file_path + '/sites/*.dat')
    sites = read_jl.load_d_and_s(s_file_paths)
    
    #load the states for each trip, t
    t_file_path = read_jl.data_paths(file_path + '/states/*.dat')
    trips = read_jl.load_trips(t_file_path)
    
    #compile each of the travel time from the depot to the site
    trip_time_df, trips_gdf = read_jl.compile_trip_time(trips)
    
    
    #-------------------------Read in the API data-----------------------------
    API_file = 'C:/Users/Aaron/Documents/GitHub/Hitchhiking/Hitchhiking/google maps API/results/2024-01-30/API_results.dat'
    with open(API_file, 'rb') as temp:
        API_res = pickle.load(temp)
        temp.close() 
        
    API_df = pd.DataFrame(API_res)
    API_df['duration (min)'] = API_df['duration']/60
    
    #create GeoDataFrame as well
    geometry = [Point(x,y) for x,y in zip(API_df['Lon'], API_df['Lat'])]
    API_gdf = gpd.GeoDataFrame(API_df, crs="EPSG:4326", geometry=geometry)
    
    
    #-------------------------Combine the datasets-----------------------------
    combined_gdf = pd.merge(API_gdf, trip_time_df, how='left')
    combined_gdf['delta_to_hitch'] = combined_gdf['duration (min)'] - combined_gdf['full time (min)']
    combined_gdf['delta_to_hitch (perc)'] = (combined_gdf['delta_to_hitch'] / combined_gdf['duration (min)'])*100
    
    #removing Julia direct drone flights out over the water, temporary fix for now
    
    
    #---------------------------Plotting---------------------------------------
    plt.rcParams['figure.dpi'] = 500
    
    #point to the map files for plotting
    muni_file_path = 'C:/Users/Aaron/Documents/GitHub/Hitchhiking/Hitchhiking/muni_routes/geo_export_6dee9e27-b549-4312-94b4-5ad68950d5fe.shp'
    SF_boundary_file_path = 'C:/Users/Aaron/Documents/GitHub/Hitchhiking/Hitchhiking/SF Boundary/s7d02x.shp'
    SF_Zoning_file_path = 'C:/Users/Aaron/Documents/GitHub/Hitchhiking/Hitchhiking/SF Zoning/geo_export_486f35c6-390f-4739-8f97-6171581968e5.shp'
    
    #load the map files
    muni = SF.load_muni(muni_file_path)
    SF_boundary = SF.load_SF_boundaries(SF_boundary_file_path)
    SF_zoning = SF.load_SF_zoning(SF_Zoning_file_path)
    
    #load the bus data
    bus_only = SF.muni_bus_only(muni)
    bus_freq, bus_freq_full = SF.add_bus_freq(bus_only)
    
    
    SF.plot_muni(muni, SF_boundary, SF_zoning)
    SF.plot_muni_freq(bus_freq, SF_boundary, SF_zoning)
    
    #plot the hitchhiking drone travel time graphic
    fig, ax = init_fig(figsize=(10,10))
    fig, ax = plot_boundary(SF_boundary)
    fig, ax = plot_zoning(SF_zoning)
    fix, ax = plot_bus_routes(fig, ax, bus_freq)
    fig, ax = plot_travel_time(fig, ax, trips_gdf)
    fig, ax = plot_depots(fig, ax, depots_gdf)
    #fig, ax = plot_drone_rad(fig, ax, depots_gdf)
    
    
    #plot the API data for driving
    fig, ax = init_fig(figsize=(10,10))
    fig, ax = plot_boundary(SF_boundary)
    fig, ax = plot_zoning(SF_zoning)
    fix, ax = plot_bus_routes(fig, ax, bus_freq)
    fig, ax = plot_driving_time(fig, ax, API_gdf)
    fig, ax = plot_depots(fig, ax, depots_gdf)
    
    
    #plot the API data for bicycling
    fig, ax = init_fig(figsize=(10,10))
    fig, ax = plot_boundary(SF_boundary)
    fig, ax = plot_zoning(SF_zoning)
    fix, ax = plot_bus_routes(fig, ax, bus_freq)
    fig, ax = plot_bicycling_time(fig, ax, API_gdf)
    fig, ax = plot_depots(fig, ax, depots_gdf)
    
    
    #plot the delta between hitchhiking and cars
    fig, ax = init_fig(figsize=(10,10))
    fig, ax = plot_boundary(SF_boundary)
    fig, ax = plot_zoning(SF_zoning)
    fix, ax = plot_bus_routes(fig, ax, bus_freq)
    fig, ax = plot_delta_driving(fig, ax, combined_gdf)
    fig, ax = plot_depots(fig, ax, depots_gdf)
    
    
    #plot the delta between hitchhiking and bicycling
    fig, ax = init_fig(figsize=(10,10))
    fig, ax = plot_boundary(SF_boundary)
    fig, ax = plot_zoning(SF_zoning)
    fix, ax = plot_bus_routes(fig, ax, bus_freq)
    fig, ax = plot_delta_bicycling(fig, ax, combined_gdf)
    fig, ax = plot_depots(fig, ax, depots_gdf)