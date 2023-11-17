# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 09:32:58 2023

@author: Aaron
"""

import geopandas as gpd
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi'] = 500


def load_muni(file_path):
    fp = file_path
    data = gpd.read_file(fp)
    data = data.to_crs(epsg = 4326)
    return data

def load_SF_boundaries(file_path):
    fp = file_path
    data = gpd.read_file(fp)
    data = data.to_crs(epsg = 4326)
    #data_SF = data[data['NAME'] == 'San Francisco']
    return data

def load_SF_zoning(file_path):
    fp = file_path
    data = gpd.read_file(fp)
    data = data.to_crs(epsg = 4326)
    data = data[data['zoning'] == 'PDR-2']
    #data_SF = data[data['NAME'] == 'San Francisco']
    return data

def muni_bus_only(muni_data):
    bus_only_data = muni_data[(muni_data['service_ca'] != 'Cable Car') &
                              (muni_data['service_ca'] != 'Historic') &
                              (muni_data['service_ca'] != 'OWL') &
                              (muni_data['service_ca'] != 'Owl') &
                              (muni_data['direction'] != 'I')
                              ]
    return bus_only_data

def add_bus_freq(bus_only_data):
    #bus_only_data['frequency'] = None
    bus_only_data.loc[bus_only_data['route_name'] == '1', 'frequency'] = 6 #updated based on midday numbers during the week
    bus_only_data.loc[bus_only_data['route_name'] == 'LBUS', 'frequency'] = 8
    bus_only_data.loc[bus_only_data['route_name'] == 'M', 'frequency'] = 10
    bus_only_data.loc[bus_only_data['route_name'] == 'N', 'frequency'] = 10
    bus_only_data.loc[bus_only_data['route_name'] == '2', 'frequency'] = 20
    bus_only_data.loc[bus_only_data['route_name'] == '5', 'frequency'] = 12
    bus_only_data.loc[bus_only_data['route_name'] == '5R', 'frequency'] = 10
    bus_only_data.loc[bus_only_data['route_name'] == '6', 'frequency'] = 20
    bus_only_data.loc[bus_only_data['route_name'] == '7', 'frequency'] = 12
    bus_only_data.loc[bus_only_data['route_name'] == '8', 'frequency'] = 8
    bus_only_data.loc[bus_only_data['route_name'] == '9', 'frequency'] = 10
    bus_only_data.loc[bus_only_data['route_name'] == '9R', 'frequency'] = 12
    bus_only_data.loc[bus_only_data['route_name'] == '12', 'frequency'] = 10
    bus_only_data.loc[bus_only_data['route_name'] == '14', 'frequency'] = 10
    bus_only_data.loc[bus_only_data['route_name'] == '14R', 'frequency'] = 8
    bus_only_data.loc[bus_only_data['route_name'] == '15', 'frequency'] = 10
    bus_only_data.loc[bus_only_data['route_name'] == '18', 'frequency'] = 20
    bus_only_data.loc[bus_only_data['route_name'] == '19', 'frequency'] = 15
    bus_only_data.loc[bus_only_data['route_name'] == '21', 'frequency'] = 20
    bus_only_data.loc[bus_only_data['route_name'] == '22', 'frequency'] = 6
    bus_only_data.loc[bus_only_data['route_name'] == '23', 'frequency'] = 20
    bus_only_data.loc[bus_only_data['route_name'] == '24', 'frequency'] = 10
    bus_only_data.loc[bus_only_data['route_name'] == '25', 'frequency'] = 15
    bus_only_data.loc[bus_only_data['route_name'] == '27', 'frequency'] = 15
    bus_only_data.loc[bus_only_data['route_name'] == '28', 'frequency'] = 10
    bus_only_data.loc[bus_only_data['route_name'] == '28R', 'frequency'] = 12
    bus_only_data.loc[bus_only_data['route_name'] == '29', 'frequency'] = 9
    bus_only_data.loc[bus_only_data['route_name'] == '30', 'frequency'] = 12
    bus_only_data.loc[bus_only_data['route_name'] == '31', 'frequency'] = 20
    bus_only_data.loc[bus_only_data['route_name'] == '33', 'frequency'] = 20
    bus_only_data.loc[bus_only_data['route_name'] == '35', 'frequency'] = 30
    bus_only_data.loc[bus_only_data['route_name'] == '36', 'frequency'] = 30
    bus_only_data.loc[bus_only_data['route_name'] == '37', 'frequency'] = 20
    bus_only_data.loc[bus_only_data['route_name'] == '38', 'frequency'] = 16
    bus_only_data.loc[bus_only_data['route_name'] == '38R', 'frequency'] = 6
    bus_only_data.loc[bus_only_data['route_name'] == '39', 'frequency'] = 20
    bus_only_data.loc[bus_only_data['route_name'] == '43', 'frequency'] = 12
    bus_only_data.loc[bus_only_data['route_name'] == '44', 'frequency'] = 10
    bus_only_data.loc[bus_only_data['route_name'] == '45', 'frequency'] = 10
    bus_only_data.loc[bus_only_data['route_name'] == '48', 'frequency'] = 15
    bus_only_data.loc[bus_only_data['route_name'] == '49', 'frequency'] = 6
    bus_only_data.loc[bus_only_data['route_name'] == '52', 'frequency'] = 20
    bus_only_data.loc[bus_only_data['route_name'] == '54', 'frequency'] = 20
    bus_only_data.loc[bus_only_data['route_name'] == '55', 'frequency'] = 15
    bus_only_data.loc[bus_only_data['route_name'] == '56', 'frequency'] = 20
    bus_only_data.loc[bus_only_data['route_name'] == '57', 'frequency'] = 20
    bus_only_data.loc[bus_only_data['route_name'] == '58', 'frequency'] = 30
    bus_only_data.loc[bus_only_data['route_name'] == '66', 'frequency'] = 20
    bus_only_data.loc[bus_only_data['route_name'] == '67', 'frequency'] = 20
            
    bus_freq_full = bus_only_data
    
    bus_freq = bus_freq_full[bus_freq_full['frequency'].notna()]
    
    return bus_freq, bus_freq_full

def plot_muni(muni_data, SF_boundary, SF_zoning):
    fig, ax = plt.subplots()
    plt.xlim([-122.52480, -122.33907])
    plt.ylim([37.70, 37.84])
    
    SF_boundary.boundary.plot(ax=ax, edgecolor = 'k')
    muni_data.plot(ax=ax, linewidth = 0.5, column = muni_data['service_ca'], legend = True,
                   legend_kwds = {"loc": 'upper left', 'bbox_to_anchor': (1.05, 1.05)})
    SF_zoning.plot(ax=ax, legend = True) #color = 'gray'
    
    plt.xlabel('longitude')
    plt.ylabel('latitude')

    return

def plot_muni_freq(bus_freq, SF_boundary, SF_zoning):
    fig, ax = plt.subplots()
    plt.xlim([-122.52480, -122.33907])
    plt.ylim([37.70, 37.84])
    
    SF_boundary.boundary.plot(ax=ax, edgecolor = 'k')
    SF_zoning.plot(ax=ax, legend = True) #color = 'gray'
    
    bus_freq.plot(ax=ax, linewidth = 0.5, column = 'frequency', legend = True,
                   cmap = 'RdYlGn_r')

    plt.xlabel('longitude')
    plt.ylabel('latitude')
    plt.title('Bus Frequency (minutes)')
    
    return


#execute section

#load the necessary shapefiles
muni_file_path = 'C:/Users/Aaron/Documents/GitHub/Hitchhiking/Hitchhiking/muni_routes/geo_export_6dee9e27-b549-4312-94b4-5ad68950d5fe.shp'
SF_boundary_file_path = 'C:/Users/Aaron/Documents/GitHub/Hitchhiking/Hitchhiking/SF Boundary/s7d02x.shp'
SF_Zoning_file_path = 'C:/Users/Aaron/Documents/GitHub/Hitchhiking/Hitchhiking/SF Zoning/geo_export_486f35c6-390f-4739-8f97-6171581968e5.shp'

muni = load_muni(muni_file_path)
SF_boundary = load_SF_boundaries(SF_boundary_file_path)
SF_zoning = load_SF_zoning(SF_Zoning_file_path)

bus_only = muni_bus_only(muni)
bus_freq, bus_freq_full = add_bus_freq(bus_only)


plot_muni(muni, SF_boundary, SF_zoning)
plot_muni_freq(bus_freq, SF_boundary, SF_zoning)