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
from matplotlib.colors import TwoSlopeNorm
import trip_recharging as recharging
import math
#import contextily as cx
#import geoplot as gplt
#import distance_matrix_API as API

def roundup(x, base=5): #https://stackoverflow.com/questions/26454649/python-round-up-to-the-nearest-ten
    return base * math.ceil(x/base) #https://stackoverflow.com/questions/2272149/round-to-5-or-other-number-in-python

def init_fig(figsize):
    #fig, ax = plt.subplots(figsize = figsize)
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

def plot_travel_time(fig, ax, gdf, vmin, vmax):
    gdf.plot(column='full time (min)', 
                ax=ax,
                cmap='viridis',
                legend=True,
                legend_kwds={'label':'Trip Time (minutes)'},
                vmin=vmin,
                vmax=vmax
                ) #, 'orientation': 'horizontal'
    plt.title('Hitchhiking Isochrone')
    return fig, ax

def plot_CO2_emissions(fig, ax, gdf, vmin, vmax):
    gdf.plot(column='CO2 Emissions', 
                ax=ax,
                cmap='viridis',
                legend=True,
                legend_kwds={'label':'CO2 Emissions (kg)'},
                vmin=vmin,
                vmax=vmax
                ) #, 'orientation': 'horizontal'
    plt.title('Hitchhiking Isochrone')
    return fig, ax

def plot_hitchhiking_time_boundaries(fig, ax, gdf, vmin, vmax):
    
    polys = []
    mag = []
    for item in gdf['dur boundary'].unique():
        subset = gdf[gdf['dur boundary'] == item] #get a subset of the datapoints within a driving range
        polygon = subset.unary_union.convex_hull #find the boundary of this subset of datapoints
        polys.append(polygon) #save the polygon
        mag.append(item) #save the magnitude of the driving distance
        
    boundary_dict = {'hitchhiking time': mag, 'geometry': polys}
    gdf_boundaries = gpd.GeoDataFrame(boundary_dict)
    
    gdf_boundaries.plot(column='hitchhiking time',
                        ax=ax,
                        cmap='viridis',
                        legend=True,
                        legend_kwds={'label':'Hitchhiking Time (minutes)'},
                        vmin=vmin,
                        vmax=vmax,
                        alpha = 0.75)
        
    return fig, ax, gdf_boundaries

def plot_depots(fig, ax, gdf):
    gdf.to_crs(7131)
    gdf.plot(ax=ax,
             legend=True,
             markersize=1500,
             color='darkorange',
             edgecolor='k',
             marker='p',
             legend_kwds={'label':'Depot'})
    
    #added to help prject the figure to 7131
    radius = gdf.iloc[0]['geometry'].buffer(.1/111).boundary
    radius = gpd.GeoSeries(radius, crs='7131' )
    
    #radius.to_crs(4326)
    radius.plot(ax=ax,
                color='darkorange',
                linewidth=4,
                linestyle='--',
                label="Direct Drone-Only",
                alpha=0.8)
    
    return fig, ax

def plot_test_points(fig, ax):
    #geometry = [Point(x,y) for x,y in zip(API_df['Lon'], API_df['Lat'])]
    #API_gdf = gpd.GeoDataFrame(API_df, crs="EPSG:7131", geometry=geometry)
    
    data = {'name': ['test point 1', 'test point 2', 'test point 3', 'test point 4']}
    df = pd.DataFrame(data)
    #geometry = [Point(-122.40099, 37.71257), Point(-122.36172, 37.73979), Point(-122.40162, 37.77550), Point(-122.44097, 37.74424)]
    geometry = [Point(-122.40099, 37.71257), Point(-122.36967, 37.73979), Point(-122.40162, 37.77550), Point(-122.43273, 37.74424)]
    
    gdf = gpd.GeoDataFrame(df, crs='EPSG:7131', geometry=geometry)
    gdf.plot(ax=ax,
              markersize=500)
    
    return fig, ax

# def plot_test_points(fig, ax):
#     #geometry = [Point(x,y) for x,y in zip(API_df['Lon'], API_df['Lat'])]
#     #API_gdf = gpd.GeoDataFrame(API_df, crs="EPSG:7131", geometry=geometry)
    
#     data = {'name': ['test point 1', 'test point 2', 'test point 3', 'test point 4']}
#     df = pd.DataFrame(data)
#     geometry = [Point(-13625615.87939257, 4538900.789086343), Point(-122.36172, 37.73979), Point(-122.40162, 37.77550), Point(-122.44097, 37.74424)]
    
#     gdf = gpd.GeoDataFrame(df, crs='EPSG:3857', geometry=geometry)
#     gdf.plot(ax=ax,
#              markersize=500)
    
#     return fig, ax


def plot_bus_routes(fig, ax, gdf):
    gdf.plot(ax=ax,
             color='slategrey',
             alpha=0.5)
    return fig, ax

def plot_drone_rad(fig, ax, gdf, delivery_radius, figsize):
    #gdf.to_crs(7131, inplace=True)
    gdf.drop_duplicates(subset=['geometry'], inplace=True)

    for i in range(len(gdf)):
        depot_loc = gdf.iloc[i]['geometry']
        
        radius = gdf.iloc[i]['geometry'].buffer(delivery_radius/111).boundary
        radius = gpd.GeoSeries(radius, crs='7131' )
        
        #radius.to_crs(4326)
        radius.plot(ax=ax,
                    color='k',
                    linewidth=4,
                    linestyle='--',
                    figsize=figsize,
                    label="Direct Drone-Only",
                    alpha=0.8)
    
    return fig, ax, radius

# def plot_drone_rad(fig, ax, gdf, delivery_radius, figsize):
#     gdf.to_crs(3857, inplace=True)
#     gdf.drop_duplicates(subset=['geometry'], inplace=True)

#     for i in range(len(gdf)):
#         depot_loc = gdf.iloc[i]['geometry']
        
#         r = gdf.iloc[i]['geometry'].buffer(delivery_radius*1000).boundary
#         radius = gpd.GeoSeries(r, crs="EPSG:3857" )
#         radius = radius.to_crs(4326)
#         #radius.to_crs(4326)
#         radius.plot(ax=ax,
#                     color='k',
#                     linewidth=4,
#                     linestyle='--',
#                     figsize=figsize,
#                     label="Direct Drone-Only",
#                     alpha=0.8)
    
#     return fig, ax, radius

# def plot_drone_rad(fig, ax, gdf, delivery_radius, figsize):
#     #gdf.to_crs(3171, inplace=True)
#     gdf.drop_duplicates(subset=['geometry'], inplace=True)

#     for i in range(len(gdf)):
#         depot_loc = gdf.iloc[i]['geometry']
        
#         r = gdf.iloc[i]['geometry'].buffer(delivery_radius/88).boundary
#         radius = gpd.GeoSeries(r, crs="EPSG:4326" )
#         radius = radius.to_crs(4326)
#         #radius.to_crs(4326)
#         radius.plot(ax=ax,
#                     color='k',
#                     linewidth=4,
#                     linestyle='--',
#                     figsize=figsize,
#                     label="Direct Drone-Only",
#                     alpha=0.8)
    
#     return fig, ax, radius


# def plot_drone_rad(fig, ax, gdf, delivery_radius, figsize):
#     #gdf.to_crs(3171, inplace=True)
#     gdf.drop_duplicates(subset=['geometry'], inplace=True)

#     for i in range(len(gdf)):
#         depot_loc = gdf.iloc[i]['geometry']
        
#         r = gdf.iloc[i]['geometry']
#         ellipse = ((r.x, r.y), (delivery_radius/88, delivery_radius/111), 0)
        
#         #r = gdf.iloc[i]['geometry'].buffer(delivery_radius/88).boundary
#         radius = gpd.GeoSeries(r, crs="EPSG:4326" )
#         radius = radius.to_crs(4326)
#         #radius.to_crs(4326)
#         radius.plot(ax=ax,
#                     color='k',
#                     linewidth=4,
#                     linestyle='--',
#                     figsize=figsize,
#                     label="Direct Drone-Only",
#                     alpha=0.8)
    
#     return fig, ax, radius


def plot_recharging_stations(fig, ax, gdf, delivery_radius, figsize):
    
    gdf.plot(ax=ax,
             legend=True,
             markersize=1500,
             color='lightgreen',
             edgecolor='k',
             marker='P',
             legend_kwds={'label':'Recharge Locations'})
    
    for i in range(len(gdf)):
        recharge_loc = gdf.iloc[i]['geometry']
        radius = gdf.iloc[i]['geometry'].buffer(delivery_radius/111).boundary
        radius = gpd.GeoSeries(radius, crs='7131' )
        radius.plot(ax=ax,
                    color='k',
                    linewidth=4,
                    linestyle='--',
                    figsize=figsize,
                    label="Direct Drone-Only",
                    alpha=0.8)
        
    #plt.title('Depot and Recharging Locations in San Francisco', fontsize=24)
    
    return fig, ax

def add_legend_recharging(fig, ax):
    ax.legend(['Depot Location', 'Recharging Site'], loc='center left', bbox_to_anchor=(1, 0.5), fontsize=20)
    plt.title('Depot and Recharging Locations in San Francisco', fontsize=24)
    
    return fig, ax

def plot_recharing_trips(fig, ax, gdf, vmin, vmax):
    gdf.plot(column='one-way trip time (min)', 
                ax=ax,
                cmap='viridis',
                legend=True,
                legend_kwds={'label':'Trip Time (minutes)'},
                vmin=vmin,
                vmax=vmax) #, 'orientation': 'horizontal'
    plt.title('Recharging Isochrone')
    
    return fig, ax

def plot_recharing_CO2(fig, ax, gdf, vmin, vmax):
    gdf.plot(column='CO2 Emissions Recharging', 
                ax=ax,
                cmap='viridis',
                legend=True,
                legend_kwds={'label':'CO2 Emissions(kg))'},
                vmin=vmin,
                vmax=vmax) #, 'orientation': 'horizontal'
    plt.title('Recharging Isochrone')
    
    return fig, ax

def plot_driving_time(fig, ax, gdf, vmin, vmax):
    gdf_driving = gdf[gdf['mode'] == 'driving']
    gdf_driving.plot(column='duration in traffic (min)', 
                ax=ax,
                cmap='viridis',
                legend=True,
                legend_kwds={'label':'Driving Time (minutes)'},
                vmin=vmin,
                vmax=vmax)
    plt.title('Driving Isochrone (with traffic)')
    
    return fig, ax

def plot_driving_CO2_emissions(fig, ax, gdf, vmin, vmax):
    gdf_driving = gdf[gdf['mode'] == 'driving']
    gdf_driving.plot(column='CO2 Emissions Driving', 
                ax=ax,
                cmap='viridis',
                vmin=vmin,
                vmax=vmax,
                legend=True,
                legend_kwds={'label':'CO2 Emissions (kg)'}
                )
    plt.title('CO2 Emissions Isochrone - Driving')
    
    return fig, ax

def plot_driving_time_boundaries(fig, ax, gdf, vmin, vmax):
    gdf_driving = gdf[gdf['mode'] == 'driving']
    
    polys = []
    mag = []
    for item in gdf_driving['dur boundary'].unique():
        subset = gdf_driving[gdf_driving['dur boundary'] == item] #get a subset of the datapoints within a driving range
        polygon = subset.unary_union.convex_hull #find the boundary of this subset of datapoints
        polys.append(polygon) #save the polygon
        mag.append(item) #save the magnitude of the driving distance
        
    boundary_dict = {'driving time': mag, 'geometry': polys}
    gdf_boundaries = gpd.GeoDataFrame(boundary_dict)
    
    gdf_boundaries.plot(column='driving time',
                        ax=ax,
                        cmap='viridis',
                        legend=True,
                        legend_kwds={'label':'Driving Time (minutes)'},
                        vmin=vmin,
                        vmax=vmax,
                        alpha = 1)
        
    return fig, ax, gdf_boundaries

def plot_bicycling_time(fig, ax, gdf, vmin, vmax):
    gdf_bicycling = gdf[gdf['mode'] == 'bicycling']
    gdf_bicycling.plot(column='duration no traffic (min)', 
                ax=ax,
                cmap='viridis',
                legend=True,
                legend_kwds={'label':'Bicycling Time (minutes)'},
                vmin=vmin,
                vmax=vmax)
    plt.title('Bicycling Isochrone (no traffic)')
    
    return fig, ax

def plot_bicycling_time_contour(fig, ax, gdf, vmin, vmax):
    gdf_bicycling = gdf[gdf['mode'] == 'bicycling']
    gplt.kdeplot(gdf, )
    gdf_bicycling.plot(column='duration (min)', 
                ax=ax,
                cmap='viridis',
                legend=True,
                legend_kwds={'label':'Bicycling Time (minutes)'},
                vmin=vmin,
                vmax=vmax)
    plt.title('Bicycling Isochrone')
    
    return fig, ax

def plot_bicycling_CO2_emissions(fig, ax, gdf, vmin, vmax):
    gdf_bicycling = gdf[gdf['mode'] == 'bicycling']
    gdf_bicycling.plot(column='CO2 Emissions Bicycling', 
                ax=ax,
                cmap='viridis',
                vmin=vmin,
                vmax=vmax,
                legend=True,
                legend_kwds={'label':'CO2 Emissions (kg)'}
                )
    plt.title('CO2 Emissions Isochrone - Bicycling')
    
    return fig, ax

def plot_bicycling_time_boundaries(fig, ax, gdf, vmin, vmax):
    gdf_bicycling = gdf[gdf['mode'] == 'bicycling']
    
    polys = []
    mag = []
    for item in gdf_bicycling['dur boundary'].unique():
        subset = gdf_bicycling[gdf_bicycling['dur boundary'] == item] #get a subset of the datapoints within a driving range
        polygon = subset.unary_union.convex_hull #find the boundary of this subset of datapoints
        polys.append(polygon) #save the polygon
        mag.append(item) #save the magnitude of the driving distance
        
    boundary_dict = {'bicycling time': mag, 'geometry': polys}
    gdf_boundaries = gpd.GeoDataFrame(boundary_dict)
    
    gdf_boundaries.plot(column='bicycling time',
                        ax=ax,
                        cmap='viridis',
                        legend=True,
                        legend_kwds={'label':'bicycling Time (minutes)'},
                        vmin=vmin,
                        vmax=vmax,
                        alpha = 0.5)
        
    return fig, ax, gdf_boundaries

def plot_delta_recharging(fig, ax, gdf, vmin, vmax):
    vmin, vmax, vcenter = -160, vmax, 0  #https://gis.stackexchange.com/questions/330008/center-normalize-choropleth-colors-in-geopandas
    norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
    cmap='RdBu'
    cbar=plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    
    gdf.plot(column='delta_to_hitch (perc)', ax=ax, cmap=cmap, norm=norm, legend=True)
    
    plt.title('Percent Difference in Trip Time Isochrone (Recharging - Hitchhiking)')
    return fig, ax

def plot_delta_CO2_recharging(fig, ax, gdf, vmin, vmax):
    vmin, vmax, vcenter = -175, vmax, 0
    norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
    cmap='RdBu'
    cbar=plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    
    gdf.plot(column='delta_to_hitch CO2 (perc)', ax=ax, cmap=cmap, norm=norm, legend=True) #norm=norm, 

    plt.title('Percent Difference in CO2 Isochrone (Recharging - Hitchhiking)')

    return fig, ax

def plot_delta_driving(fig, ax, gdf, vmin, vmax):
    gdf_driving = gdf[gdf['mode']=='driving']
    vmin, vmax, vcenter = -160, vmax, 0  #https://gis.stackexchange.com/questions/330008/center-normalize-choropleth-colors-in-geopandas
    norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
    cmap='RdBu'
    cbar=plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    
    gdf_driving.plot(column='delta_to_hitch time driving (perc)', ax=ax, cmap=cmap, norm=norm, legend=True)

    plt.title('Percent Difference in Trip Time Isochrone (Driving - Hitchhiking)')
    
    return fig, ax

def plot_delta_CO2_driving(fig, ax, gdf, vmin, vmax):
    gdf_driving = gdf[gdf['mode']=='driving']
    vmin, vmax, vcenter = -175, vmax, 0  #https://gis.stackexchange.com/questions/330008/center-normalize-choropleth-colors-in-geopandas
    norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
    cmap='RdBu'
    cbar=plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    
    gdf_driving.plot(column='delta_to_hitch CO2 driving (perc)', ax=ax, cmap=cmap, norm=norm, legend=True) #norm=norm, 

    plt.title('Percent Difference in CO2 Isochrone (Driving - Hitchhiking)')
    
    return fig, ax

def plot_delta_bicycling(fig, ax, gdf, vmin, vmax):
    gdf_bicycling = gdf[gdf['mode']=='bicycling']
    
    vmin, vmax, vcenter = -160, vmax, 0
    norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
    cmap='RdBu'
    cbar=plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    
    gdf_bicycling.plot(column='delta_to_hitch time bicycling (perc)', ax=ax, cmap=cmap, norm=norm, legend=True)
    
    plt.title('Percent Difference in Trip Time Isochrone (Bicycling - Hitchhiking)')

    return fig, ax

def plot_delta_CO2_bicycling(fig, ax, gdf, vmin, vmax):
    gdf_bicycling = gdf[gdf['mode']=='bicycling']
    vmin, vmax, vcenter = -175, vmax, 0
    norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
    cmap='RdBu'
    cbar=plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    
    gdf_bicycling.plot(column='delta_to_hitch CO2 bicycling (perc)', ax=ax, cmap=cmap, norm=norm, legend=True) #norm=norm, 

    plt.title('Percent Difference in CO2 Isochrone (Bicycling - Hitchhiking)')

    return fig, ax



if __name__ == '__main__':

    #----------------------------Read in the Julia Data------------------------
    #set the results parent folder for the hitchhiking data
    file_path = 'C:/Users/Aaron/AppData/Local/Programs/Julia-1.6.7/MultiAgentAllocationTransit.jl/results/2024-02-16 (d_100_s_100_iter_100_2281_sites)' #01-30 (d_100_s_100_iter_20)
    
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
    
    #load the actions for each trip, t
    a_file_path = read_jl.data_paths(file_path + '/actions/*.dat')
    actions = read_jl.load_actions(a_file_path)
    
    #compile each of the travel time from the depot to the site
    trip_time_df, trips_gdf = read_jl.compile_trip_time(trips)
    
    #compile the flight distance for the trip
    trip_dist_df = read_jl.compile_trip_dist(trip_time_df, actions)
    
    #combine trip time and dist dataframes
    trip_main_gdf = read_jl.combine_trip_time_dist(trips_gdf, trip_dist_df)
    
    #get the vmin and vmax values for trip time across modes
    trip_main_feas_gdf = trip_main_gdf[trip_main_gdf['Total Dist'] <= 7]
    #trip_main_feas_gdf.to_crs(7131, inplace=True)
    
    #add CO2 emissions
    trip_main_feas_gdf['CO2 Emissions'] = trip_main_feas_gdf['Total Dist'] * 0.012 #kg/km  #(old source, see paper for new data) 70 grams/km #https://www.sciencedirect.com/science/article/pii/S2666389922001805#:~:text=Our%20model%20shows%20that%20an,package%20in%20the%20United%20States.
    
    #------------------------Parameters----------------------------------------
    delivery_radius = 3.5 #km
    
    
    #-------------------------Read in the API data-----------------------------
    API_file = 'C:/Users/Aaron/Documents/GitHub/Hitchhiking/Hitchhiking/google maps API/results/2024-02-16 (2281 sites)_datetime1000/API_results.dat' #10am driving traffic
    #API_file = 'C:/Users/Aaron/Documents/GitHub/Hitchhiking/Hitchhiking/google maps API/results/2024-03-27 (2281 sites)_datetime1700/API_results.dat' #5pm driving traffic
    #API_file = 'C:/Users/Aaron/Documents/GitHub/Hitchhiking/Hitchhiking/google maps API/results/2024-03-27 (2281 sites)_UTCseconds1700/API_results.dat' #5pm driving traffic
    API_file = 'C:/Users/Aaron/Documents/GitHub/Hitchhiking/Hitchhiking/google maps API/results/2024-03-28 (2281 sites)_datetime1700_pessimistic/API_results.dat' #5pm driving traffic
    
    with open(API_file, 'rb') as temp:
        API_res = pickle.load(temp)
        temp.close() 
        
    API_df = pd.DataFrame(API_res)
    #need to correct the old API data, was pulling the "duration", but not the "duration-in-traffic", need this for peak hour analysis
    API_df.rename(columns={"duration":"duration no traffic"}, inplace=True)
    
    #split out the driving API runs so that I can adjust duration to duration in traffic (not a parameter for biking)
    API_driving = API_df[API_df['mode'] == 'driving']
    API_driving.reset_index(inplace=True)
    API_bicycling = API_df[API_df['mode'] == 'bicycling']
    API_bicycling.reset_index(inplace=True)
    
    #pull out the duration in traffic data from the API response
    dur_traffic_ls = [API_driving['response'][i]['rows'][0]['elements'][0]['duration_in_traffic']['value'] for i in range(len(API_driving))]
    #add the duration in traffic data to the dataframe
    API_driving['duration in traffic'] = dur_traffic_ls
        
    #scale the duration in traffic appropriately
    API_driving['duration in traffic (min)'] = API_driving['duration in traffic']/60
    
    #still scale the bicycling and driving duration as well
    API_bicycling['duration no traffic (min)'] = API_bicycling['duration no traffic']/60
    #API_driving['duration no traffic (min)'] = API_driving['duration no traffic']/60
    
    #add in the CO2 emissions estimates
    #API_driving['Fuel Consumption (gal)'] = API_driving['distance']/1000*0.621371 / 24.4 #https://afdc.energy.gov/data/10310
    API_driving['avg speed'] = (API_driving['distance']/1000*0.621371) / (API_driving['duration in traffic (min)']/60) #speed in miles per hour
    API_driving['Fuel Consumption (gal)'] = API_driving['duration in traffic (min)']/60* (1/(24.4*(1/API_driving['avg speed']))) #1.025 based on 24.4 MPG car efficienct, 25 mph hour average speed
    API_driving['CO2 Emissions Driving'] = API_driving['Fuel Consumption (gal)'] * 8.887 #kg on CO2 per gallon
    API_bicycling['CO2 Emissions Bicycling'] = API_bicycling['distance']/1000 * 0.0048 #kg CO2 emitted per km  #http://large.stanford.edu/courses/2022/ph240/schutt2/ (old source, see writeup)
    
    #API_df = API_driving.merge(right=API_bicycling, how='inner', on='destinations').fillna(0)
    API_df = pd.concat([API_bicycling, API_driving])
    API_df.set_index('index', drop=True, inplace=True)
    
    #create GeoDataFrame as well
    geometry = [Point(x,y) for x,y in zip(API_df['Lon'], API_df['Lat'])]
    API_gdf = gpd.GeoDataFrame(API_df, crs="EPSG:7131", geometry=geometry)
    
    
    #-------------------------Calc the recharging case-------------------------
    # #set the recharging configuration
    recharge_loc = {'Loc': ['1', '2', '3', '4', '5'], #config for 3.5km radius
                    'Lat': [37.73200, 37.73200, 37.78350, 37.78350, 37.76554], 
                    'Lon': [-122.43000, -122.47750, -122.41200, -122.46000, -122.48750]
                    }
    recharge_loc_df = pd.DataFrame(recharge_loc)
    geometry = [Point(x,y) for x,y in zip(recharge_loc['Lon'], recharge_loc['Lat'])]
    recharge_gdf = gpd.GeoDataFrame(recharge_loc, crs='EPSG:7131', geometry=geometry)
    
   
    #path_ls, shortest_dist_ls = recharging.gen_recharge_trips(sites, depots_df, recharge_loc_df, delivery_radius)
    
    recharge_file = 'C:/Users/Aaron/Documents/GitHub/Hitchhiking/Hitchhiking/recharging results/2024-02-28/recharging_trips.dat'
    with open(recharge_file, 'rb') as temp:
        recharge_trip_df = pickle.load(temp)
        temp.close() 
    
    recharge_trip_df['Shortest Round-Trip Dist (km)'] = recharge_trip_df['Shortest Dist (km)'] * 2
    recharge_trip_df['one-way trip time (min)'] = recharge_trip_df['Shortest Round-Trip Dist (km)'] / 0.00777 / 60 / 2 #divide by 2 to get the one-way trip time
    
    geometry = [Point(x,y) for x,y in zip(recharge_trip_df['Lon'], recharge_trip_df['Lat'])]
    recharge_trip_gdf = gpd.GeoDataFrame(recharge_trip_df, crs="EPSG:7131", geometry=geometry)
    recharge_trip_gdf = recharge_trip_gdf[recharge_trip_gdf['one-way trip time (min)'] < 120] #remove some weird outliers with 4,290 minutes of trip time
    

    #-------------------------Combine the datasets-----------------------------
    combined_gdf = pd.merge(API_gdf, trip_main_feas_gdf, how='inner')
    #combined_gdf = pd.merge(combined_gdf, recharge_trip_gdf, how='left')
    #combined_gdf = combined_gdf[combined_gdf['Total Dist'] <= 7]
    combined_gdf['delta_to_hitch time driving'] = combined_gdf['duration in traffic (min)'] - combined_gdf['full time (min)']
    combined_gdf['delta_to_hitch time driving (perc)'] = (combined_gdf['delta_to_hitch time driving'] / combined_gdf['duration in traffic (min)'])*100
    
    combined_gdf['delta_to_hitch time bicycling'] = combined_gdf['duration no traffic (min)'] - combined_gdf['full time (min)']
    combined_gdf['delta_to_hitch time bicycling (perc)'] = (combined_gdf['delta_to_hitch time bicycling'] / combined_gdf['duration no traffic (min)'])*100
    
    combined_gdf['delta_to_hitch CO2 driving'] = combined_gdf['CO2 Emissions Driving'] - combined_gdf['CO2 Emissions']
    combined_gdf['delta_to_hitch CO2 driving (perc)'] = (combined_gdf['delta_to_hitch CO2 driving'] / combined_gdf['CO2 Emissions Driving'])*100
    
    combined_gdf['delta_to_hitch CO2 bicycling'] = combined_gdf['CO2 Emissions Bicycling'] - combined_gdf['CO2 Emissions']
    combined_gdf['delta_to_hitch CO2 bicycling (perc)'] = (combined_gdf['delta_to_hitch CO2 bicycling'] / combined_gdf['CO2 Emissions Bicycling'])*100
    
    
    combined_gdf_recharging = pd.merge(recharge_trip_gdf, trip_main_gdf, how='left')
    combined_gdf_recharging = combined_gdf_recharging[combined_gdf_recharging['Total Dist'] <= 7]
    combined_gdf_recharging['delta_to_hitch'] = combined_gdf_recharging['one-way trip time (min)'] - combined_gdf_recharging['full time (min)']
    combined_gdf_recharging['delta_to_hitch (perc)'] = (combined_gdf_recharging['delta_to_hitch'] / combined_gdf_recharging['one-way trip time (min)'])*100
    
    combined_gdf_recharging['CO2 Emissions Recharging'] = combined_gdf_recharging['Shortest Round-Trip Dist (km)'] * 0.012  #kg/km  #70 grams/km #https://www.sciencedirect.com/science/article/pii/S2666389922001805#:~:text=Our%20model%20shows%20that%20an,package%20in%20the%20United%20States.
    combined_gdf_recharging['CO2 Emissions Hitchhiking'] = combined_gdf_recharging['Total Dist'] * 0.012  #kg/km  #70 grams/km #https://www.sciencedirect.com/science/article/pii/S2666389922001805#:~:text=Our%20model%20shows%20that%20an,package%20in%20the%20United%20States.
    
    combined_gdf_recharging['delta_to_hitch CO2'] = combined_gdf_recharging['CO2 Emissions Recharging'] - combined_gdf_recharging['CO2 Emissions Hitchhiking']
    combined_gdf_recharging['delta_to_hitch CO2 (perc)'] = (combined_gdf_recharging['delta_to_hitch CO2'] / combined_gdf_recharging['CO2 Emissions Recharging'])*100
    
    
    #-------------------------get the min and max values-----------------------

    
    # vmin_trip = min(min(trip_main_feas_gdf['full time (min)']),
    #            min(recharge_trip_gdf['trip time (min)']),
    #            min(API_gdf['duration (min)'])
    #            )
    # vmax_trip = max(max(trip_main_feas_gdf['full time (min)']),
    #            max(recharge_trip_gdf['trip time (min)']),
    #            max(API_gdf['duration (min)'])
    #            )
    vmin_trip = min(min(combined_gdf['full time (min)']),
                    min(API_bicycling['duration no traffic (min)']),
                    min(API_driving['duration in traffic (min)']),
                    min(recharge_trip_gdf['one-way trip time (min)'])
                    )
    vmax_trip = max(max(combined_gdf['full time (min)']),
                    max(API_bicycling['duration no traffic (min)']),
                    max(API_driving['duration in traffic (min)']),
                    max(recharge_trip_gdf['one-way trip time (min)'])
                    )
    
    # vmax_trip = max(max(trip_main_feas_gdf['full time (min)']),
    #            max(recharge_trip_gdf['trip time (min)']),
    #            max(API_gdf['duration (min)'])
    #            )
    
    #need to subset combined_gdf into driving and bicycling in order to pull out percentages instead of nan values when data is missing
    combined_gdf_driving = combined_gdf[combined_gdf['mode'] == 'driving']
    combined_gdf_bicycling = combined_gdf[combined_gdf['mode'] == 'bicycling']
    
    #get the vmin and vmax values for trip time across modes    
    vmin_trip_perc = min(min(combined_gdf_driving['delta_to_hitch time driving (perc)']),
                         min(combined_gdf_bicycling['delta_to_hitch time bicycling (perc)']),
                         min(combined_gdf_recharging['delta_to_hitch (perc)'])
                         )
    vmax_trip_perc = max(max(combined_gdf_driving['delta_to_hitch time driving (perc)']),
                         max(combined_gdf_bicycling['delta_to_hitch time bicycling (perc)']),
                         max(combined_gdf_recharging['delta_to_hitch (perc)'])
                         )
    
    #add the CO2 min and max value TODO
    
    vmin_CO2 = min(min(combined_gdf['CO2 Emissions']),
                    min(API_bicycling['CO2 Emissions Bicycling']),
                    min(API_driving['CO2 Emissions Driving']),
                    min(combined_gdf_recharging['CO2 Emissions Recharging'])
                    )
    vmax_CO2 = max(max(combined_gdf['CO2 Emissions']),
                    max(API_bicycling['CO2 Emissions Bicycling']),
                    max(API_driving['CO2 Emissions Driving']),
                    max(combined_gdf_recharging['CO2 Emissions Recharging'])
                    )
    
    vmin_CO2_perc = min(min(combined_gdf_driving['delta_to_hitch CO2 driving (perc)']),
                    min(combined_gdf_bicycling['delta_to_hitch CO2 bicycling (perc)']),
                    min(combined_gdf_recharging['delta_to_hitch CO2 (perc)'])
                    )
    vmax_CO2_perc = max(max(combined_gdf_driving['delta_to_hitch CO2 driving (perc)']),
                    max(combined_gdf_bicycling['delta_to_hitch CO2 bicycling (perc)']),
                    max(combined_gdf_recharging['delta_to_hitch CO2 (perc)'])
                    )
    
    
    
    
    
    #-------------------------Making boundaries--------------------------------
    # section is not really needed
    #base = 1
    
    #API_gdf['dur boundary'] = [roundup(API_gdf.iloc[i]['duration (min)'], base=base) for i in range(len(API_gdf))]
    
    API_gdf_driving = API_gdf[API_gdf['mode'] == 'driving']
    API_gdf_bicycling = API_gdf[API_gdf['mode'] == 'bicycling']
    
    # for item in API_gdf_driving['dur boundary'].unique():
    #     print(item)
    #     subset = API_gdf_driving[API_gdf_driving['dur boundary'] == item]
    #     q = subset.unary_union.convex_hull
    #     gpd.GeoSeries(q).plot()
            
    #trip_main_feas_gdf['dur boundary'] = [roundup(trip_main_feas_gdf.iloc[i]['full time (min)'], base=base) for i in range(len(trip_main_feas_gdf))]
    
    
    
    
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
    
    #plot the bus routes
    SF.plot_muni(muni, SF_boundary, SF_zoning)
    SF.plot_muni_freq(bus_freq, SF_boundary, SF_zoning)
    
    
    # #get the vmin and vmax values for trip time across modes
    # trip_main_feas_gdf = trip_main_gdf[trip_main_gdf['Total Dist'] <= 7]
    
    # vmin_trip = min(min(trip_main_feas_gdf['full time (min)']),
    #            min(recharge_trip_gdf['trip time (min)']),
    #            min(API_gdf['duration (min)'])
    #            )
    # vmax_trip = max(max(trip_main_feas_gdf['full time (min)']),
    #            max(recharge_trip_gdf['trip time (min)']),
    #            max(API_gdf['duration (min)'])
    #            )
    
    # #plot the hitchhiking drone travel time graphic
    # fig, ax = init_fig(figsize=(10,10))
    # fig, ax = plot_boundary(SF_boundary)
    # fig, ax = plot_zoning(SF_zoning)
    # fix, ax = plot_bus_routes(fig, ax, bus_freq)
    # fig, ax = plot_travel_time(fig, ax, trips_gdf)
    # fig, ax = plot_depots(fig, ax, depots_gdf)
    # fig, ax = plot_drone_rad(fig, ax, depots_gdf, delivery_radius, figsize=(10,10))
    
    # Test plot
    fig, ax = init_fig(figsize=(8,8))
    fig, ax = plot_boundary(SF_boundary)
    # #fig, ax = plot_zoning(SF_zoning)
    fig, ax = plot_bus_routes(fig, ax, bus_freq)
    fig, ax = plot_depots(fig, ax, depots_gdf)
    fig, ax, radius = plot_drone_rad(fig, ax, depots_gdf, delivery_radius, figsize=(8,8))
    fig, ax = plot_test_points(fig, ax)
    
    #trip_main_feas_gdf.to_crs(7131, inplace=True) #bring trip_main_feas_gdf to the right crs just to be sure
    
    #----Remove bug points from Hitchhiking data
    #trip_main_feas_gdf = trip_main_gdf[trip_main_gdf['Total Dist'] <= 7]
    fig, ax = init_fig(figsize=(10,10))
    fig, ax = plot_boundary(SF_boundary)
    #fig, ax = plot_zoning(SF_zoning)
    fig, ax = plot_bus_routes(fig, ax, bus_freq)
    fig, ax = plot_travel_time(fig, ax, trip_main_feas_gdf, vmin_trip, vmax_trip)
    fig, ax = plot_depots(fig, ax, depots_gdf)
    fig, ax, radius = plot_drone_rad(fig, ax, depots_gdf, delivery_radius, figsize=(10,10))
    
    fig, ax = init_fig(figsize=(10,10))
    fig, ax = plot_boundary(SF_boundary)
    #fig, ax = plot_zoning(SF_zoning)
    fig, ax = plot_bus_routes(fig, ax, bus_freq)
    fig, ax = plot_CO2_emissions(fig, ax, trip_main_feas_gdf, vmin_CO2, vmax_CO2) 
    fig, ax = plot_depots(fig, ax, depots_gdf)
    fig, ax, radius = plot_drone_rad(fig, ax, depots_gdf, delivery_radius, figsize=(10,10))
    
    # #boundaries version
    # fig, ax = init_fig(figsize=(10,10))
    # fig, ax = plot_boundary(SF_boundary)
    # fig, ax = plot_zoning(SF_zoning)
    # fix, ax = plot_bus_routes(fig, ax, bus_freq)
    # fig, ax, gdf_hitchiking_boundaries = plot_hitchhiking_time_boundaries(fig, ax, trip_main_feas_gdf, vmin_trip, vmax_trip)
    # fig, ax = plot_depots(fig, ax, depots_gdf)
    # fig, ax = plot_drone_rad(fig, ax, depots_gdf, delivery_radius, figsize=(10,10))
    
    #plot the recharging locations with data
    fig, ax = init_fig(figsize=(10,10))
    fig, ax = plot_boundary(SF_boundary)
    #fig, ax = plot_zoning(SF_zoning)
    fig, ax = plot_bus_routes(fig, ax, bus_freq)
    fig, ax = plot_depots(fig, ax, depots_gdf)
    fig, ax, radius = plot_drone_rad(fig, ax, depots_gdf, delivery_radius, figsize=(10,10))
    fig, ax = plot_recharging_stations(fig, ax, recharge_gdf, delivery_radius, figsize=(10,10))
    #fig, ax = add_legend_recharging(fig, ax)
    
    #plot the recharging trip
    fig, ax = init_fig(figsize=(10,10))
    fig, ax = plot_boundary(SF_boundary)
    #fig, ax = plot_zoning(SF_zoning)
    fig, ax = plot_bus_routes(fig, ax, bus_freq)
    fig, ax = plot_recharing_trips(fig, ax, combined_gdf_recharging, vmin_trip, vmax_trip)
    fig, ax = plot_depots(fig, ax, depots_gdf)
    fig, ax, radius = plot_drone_rad(fig, ax, depots_gdf, delivery_radius, figsize=(10,10))
    fig, ax = plot_recharging_stations(fig, ax, recharge_gdf, delivery_radius, figsize=(10,10))
    
    fig, ax = init_fig(figsize=(10,10))
    fig, ax = plot_boundary(SF_boundary)
    #fig, ax = plot_zoning(SF_zoning)
    fig, ax = plot_bus_routes(fig, ax, bus_freq)
    fig, ax = plot_recharing_CO2(fig, ax, combined_gdf_recharging, vmin_CO2, vmax_CO2)
    fig, ax = plot_depots(fig, ax, depots_gdf)
    fig, ax, radius = plot_drone_rad(fig, ax, depots_gdf, delivery_radius, figsize=(10,10))
    fig, ax = plot_recharging_stations(fig, ax, recharge_gdf, delivery_radius, figsize=(10,10))

    
    
    #plot the API data for driving
    fig, ax = init_fig(figsize=(10,10))
    fig, ax = plot_boundary(SF_boundary)
    #fig, ax = plot_zoning(SF_zoning)
    fig, ax = plot_bus_routes(fig, ax, bus_freq)
    fig, ax = plot_driving_time(fig, ax, API_gdf, vmin_trip, vmax_trip)
    fig, ax = plot_depots(fig, ax, depots_gdf)
    
    fig, ax = init_fig(figsize=(10,10))
    fig, ax = plot_boundary(SF_boundary)
    #fig, ax = plot_zoning(SF_zoning)
    fig, ax = plot_bus_routes(fig, ax, bus_freq)
    fig, ax = plot_driving_CO2_emissions(fig, ax, API_gdf, vmin_CO2, vmax_CO2) #TODO, update vmin_trip and vmax_trip to CO2 numbers
    fig, ax = plot_depots(fig, ax, depots_gdf)
    
    # #plot the API data for driving
    # fig, ax = init_fig(figsize=(10,10))
    # fig, ax = plot_boundary(SF_boundary)
    # fig, ax = plot_zoning(SF_zoning)
    # fix, ax = plot_bus_routes(fig, ax, bus_freq)
    # fig, ax, gdf_driving_boundaries = plot_driving_time_boundaries(fig, ax, API_gdf, vmin_trip, vmax_trip)
    # fig, ax = plot_depots(fig, ax, depots_gdf)
    
    
    #plot the API data for bicycling
    fig, ax = init_fig(figsize=(10,10))
    fig, ax = plot_boundary(SF_boundary)
    #fig, ax = plot_zoning(SF_zoning)
    fig, ax = plot_bus_routes(fig, ax, bus_freq)
    fig, ax = plot_bicycling_time(fig, ax, API_gdf, vmin_trip, vmax_trip)
    fig, ax = plot_depots(fig, ax, depots_gdf)
    
    fig, ax = init_fig(figsize=(10,10))
    fig, ax = plot_boundary(SF_boundary)
    #fig, ax = plot_zoning(SF_zoning)
    fig, ax = plot_bus_routes(fig, ax, bus_freq)
    fig, ax = plot_bicycling_CO2_emissions(fig, ax, API_gdf, vmin_CO2, vmax_CO2)
    fig, ax = plot_depots(fig, ax, depots_gdf)
    
    # #plot the API data for bicycling
    # fig, ax = init_fig(figsize=(10,10))
    # fig, ax = plot_boundary(SF_boundary)
    # fig, ax = plot_zoning(SF_zoning)
    # fix, ax = plot_bus_routes(fig, ax, bus_freq)
    # fig, ax, gdf_bicycling_boundaries = plot_bicycling_time_boundaries(fig, ax, API_gdf, vmin_trip, vmax_trip)
    # fig, ax = plot_depots(fig, ax, depots_gdf)
    
    
    # #get the vmin and vmax values for trip time across modes    
    # vmin_trip_perc = min(min(combined_gdf['delta_to_hitch (perc)']),
    #                 min(combined_gdf_recharging['delta_to_hitch (perc)'])
    #                 )
    # vmax_trip_perc = max(max(combined_gdf['delta_to_hitch (perc)']),
    #                 max(combined_gdf_recharging['delta_to_hitch (perc)'])
    #                 )
    
    

    
    
    #plot the trip time delta
    fig, ax = init_fig(figsize=(10,10))
    fig, ax = plot_boundary(SF_boundary)
    #fig, ax = plot_zoning(SF_zoning)
    fig, ax = plot_bus_routes(fig, ax, bus_freq)
    fig, ax = plot_delta_driving(fig, ax, combined_gdf, vmin_trip_perc, vmax_trip_perc)
    fig, ax = plot_depots(fig, ax, depots_gdf)
    fig, ax, radius = plot_drone_rad(fig, ax, depots_gdf, delivery_radius, figsize=(10,10))
    
    fig, ax = init_fig(figsize=(10,10))
    fig, ax = plot_boundary(SF_boundary)
    #fig, ax = plot_zoning(SF_zoning)
    fig, ax = plot_bus_routes(fig, ax, bus_freq)
    fig, ax = plot_delta_bicycling(fig, ax, combined_gdf, vmin_trip_perc, vmax_trip_perc)
    fig, ax = plot_depots(fig, ax, depots_gdf)
    fig, ax, radius = plot_drone_rad(fig, ax, depots_gdf, delivery_radius, figsize=(10,10))
    
    fig, ax = init_fig(figsize=(10,10))
    fig, ax = plot_boundary(SF_boundary)
    #fig, ax = plot_zoning(SF_zoning)
    fig, ax = plot_bus_routes(fig, ax, bus_freq)
    fig, ax = plot_delta_recharging(fig, ax, combined_gdf_recharging, vmin_trip_perc, vmax_trip_perc)
    fig, ax = plot_depots(fig, ax, depots_gdf)
    fig, ax = plot_recharging_stations(fig, ax, recharge_gdf, delivery_radius, figsize=(10,10))
    fig, ax, radius = plot_drone_rad(fig, ax, depots_gdf, delivery_radius, figsize=(10,10))
    
    
    #plot the CO2 delta
    fig, ax = init_fig(figsize=(10,10))
    fig, ax = plot_boundary(SF_boundary)
    #fig, ax = plot_zoning(SF_zoning)
    fig, ax = plot_bus_routes(fig, ax, bus_freq)
    fig, ax = plot_delta_CO2_driving(fig, ax, combined_gdf, vmin_CO2_perc, vmax_CO2_perc)
    fig, ax = plot_depots(fig, ax, depots_gdf)
    fig, ax, radius = plot_drone_rad(fig, ax, depots_gdf, delivery_radius, figsize=(10,10))

    fig, ax = init_fig(figsize=(10,10))
    fig, ax = plot_boundary(SF_boundary)
    #fig, ax = plot_zoning(SF_zoning)
    fig, ax = plot_bus_routes(fig, ax, bus_freq)
    fig, ax = plot_delta_CO2_bicycling(fig, ax, combined_gdf, vmin_CO2_perc, vmax_CO2_perc)
    fig, ax = plot_depots(fig, ax, depots_gdf)
    fig, ax, radius = plot_drone_rad(fig, ax, depots_gdf, delivery_radius, figsize=(10,10))
    
    fig, ax = init_fig(figsize=(10,10))
    fig, ax = plot_boundary(SF_boundary)
    #fig, ax = plot_zoning(SF_zoning)
    fig, ax = plot_bus_routes(fig, ax, bus_freq)
    fig, ax = plot_delta_CO2_recharging(fig, ax, combined_gdf_recharging, vmin_CO2_perc, vmax_CO2_perc)
    fig, ax = plot_depots(fig, ax, depots_gdf)
    fig, ax = plot_recharging_stations(fig, ax, recharge_gdf, delivery_radius, figsize=(10,10))
    fig, ax, radius = plot_drone_rad(fig, ax, depots_gdf, delivery_radius, figsize=(10,10))