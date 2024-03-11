# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 10:22:47 2024

@author: Aaron
"""




import pandas as pd
import numpy as np
import geopandas as gpd
from shapely.geometry import Point
from shapely import distance
from shapely import LineString
from geopy.distance import geodesic
#from geopy.distance import great_circle
from scipy.sparse.csgraph import shortest_path
import read_julia_data as read_jl
from datetime import datetime, timedelta, date
import os
import pickle
import warnings


def get_path(Pr, i, j):
    path = [j]
    k = j
    while Pr[i, k] != -9999:
        path.append(Pr[i, k])
        k = Pr[i, k]
    return path[::-1]


if __name__ == '__main__':
    warnings.filterwarnings("ignore")
    
    #----------------------------Read in the Julia Data------------------------
    #set the results parent folder
    file_path = 'C:/Users/Aaron/AppData/Local/Programs/Julia-1.6.7/MultiAgentAllocationTransit.jl/results/2024-02-16 (d_100_s_100_iter_100_2281_sites)' #01-30 (d_100_s_100_iter_20)
    
    #load the depot locations
    d_file_paths = read_jl.data_paths(file_path + '/depots/*.dat')
    depots = read_jl.load_d_and_s(d_file_paths)
    
    #compile the depot locations
    depots_df, depots_gdf = read_jl.compile_depots(depots)
    
    #load the delivery site locations
    s_file_paths = read_jl.data_paths(file_path + '/sites/*.dat')
    sites = read_jl.load_d_and_s(s_file_paths)
    
    
    #-------------------------Calc the recharging case-------------------------
    # #set the recharging configuration
    recharge_loc = {'Loc': ['1', '2', '3', '4', '5'], #config for 3.5km radius
                    'Lat': [37.73200, 37.73200, 37.78350, 37.78350, 37.76554], 
                    'Lon': [-122.43000, -122.47750, -122.41200, -122.46000, -122.48750]
                    }
    recharge_loc_df = pd.DataFrame(recharge_loc)
    geometry = [Point(x,y) for x,y in zip(recharge_loc['Lon'], recharge_loc['Lat'])]
    recharge_gdf = gpd.GeoDataFrame(recharge_loc, crs='EPSG:4326', geometry=geometry)
    
    
    #------------------------Parameters----------------------------------------
    delivery_radius = 3.5 #km
    
    
    #------------------------Run trip estimates--------------------------------
    
    depots_df.drop_duplicates(subset=['geometry'], inplace=True)
    depots_df['Loc'] = 'Depot'
    
    loc_df = pd.concat([depots_df, recharge_loc_df], join='inner', ignore_index=True)
    
    Path_ls = []
    Shortest_dist_ls = []
    
    for s in range(len(sites)):
        #print(s)
        site = sites.iloc[s]
        site['Loc'] = 'Delivery'
        
        df = pd.concat([loc_df, site.to_frame().T], join='inner', ignore_index=True)
        
        
        #calc the distance between each location
        edges_all = []

        for i in range(len(df)):
            edges = []
            start = (df.iloc[i]['Lat'], df.iloc[i]['Lon'])
            for j in range(len(df)):
                end = (df.iloc[j]['Lat'], df.iloc[j]['Lon'])
                dist = geodesic(start, end).km
                
                if (df.iloc[j]['Loc'] == 'Delivery'):
                    if dist > delivery_radius:
                        edges.append(1000)  #make it so that you have to deliver within 3.5km radius  
                    else:
                        edges.append(dist)
                else:
                    if dist > 2*delivery_radius:
                        edges.append(1001)
                    else:            
                        edges.append(dist)


            edges_all.append(edges)
            
            
        M = np.array(edges_all) #https://stackoverflow.com/questions/53074947/examples-for-search-graph-using-scipy

        D, Pr = shortest_path(M, directed=True, method='auto', return_predecessors=True)

        Path = get_path(Pr, 0, 6)
        Shortest_dist = D[0,6]
        
        Path_ls.append(Path)
        Shortest_dist_ls.append(Shortest_dist)
        
        
    #record the final dataframe with delivery location and trip path and distances
    #sites = sites.iloc[0:10]
    recharge_trip_df = pd.DataFrame(sites)
    recharge_trip_df['Shortest Path'] = Path_ls
    recharge_trip_df['Shortest Dist (km)'] = Shortest_dist_ls
        
        
        
    #record the data
    base_path = 'C:/Users/Aaron/Documents/GitHub/Hitchhiking/Hitchhiking/recharging results'
    current_date = date.today().strftime('%Y-%m-%d') 

    folder_path = os.path.join(base_path, current_date)
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)


    save_file = os.path.join(folder_path, 'recharging_trips.dat')

    with open(save_file, 'wb') as file:
        pickle.dump(recharge_trip_df, file)













































# def gen_recharge_trips(sites_df, depots_df, recharge_loc_df, radius=3.5):
    
#     depots_df.drop_duplicates(subset=['geometry'], inplace=True)
#     depots_df['Loc'] = 'Depot'
    
#     loc_df = pd.concat([depots_df, recharge_loc_df], join='inner', ignore_index=True)
    
#     Path_ls = []
#     Shortest_dist_ls = []
    
#     for s in range(10):
#         print(s)
#         site = sites_df.iloc[s]
#         site['Loc'] = 'Delivery'
        
#         df = pd.concat([loc_df, site.to_frame().T], join='inner', ignore_index=True)
        
        
#         #calc the distance between each location
#         edges_all = []

#         for i in range(len(df)):
#             edges = []
#             start = (df.iloc[i]['Lat'], df.iloc[i]['Lon'])
#             for j in range(len(df)):
#                 end = (df.iloc[j]['Lat'], df.iloc[j]['Lon'])
#                 dist = geodesic(start, end).km
                
#                 if (df.iloc[j]['Loc'] == 'Delivery'):
#                     if dist > radius:
#                         edges.append(1000)  #make it so that you have to deliver within 3.5km radius  
#                     else:
#                         edges.append(dist)
#                 else:
#                     if dist > 2*radius:
#                         edges.append(1001)
#                     else:            
#                         edges.append(dist)


#             edges_all.append(edges)
            
            
#         M = np.array(edges_all) #https://stackoverflow.com/questions/53074947/examples-for-search-graph-using-scipy

#         D, Pr = shortest_path(M, directed=True, method='auto', return_predecessors=True)

#         Path = get_path(Pr, 0, 6)
#         Shortest_dist = D[0,6]
        
#         Path_ls.append(Path)
#         Shortest_dist_ls.append(Shortest_dist)
            
    
#     return Path_ls, Shortest_dist_ls


# radius = 3.5 #km
# recharge_loc = {'Recharging Loc': ['Loc_1', 'Loc 2', 'Loc 3', 'Loc 4', 'Loc 5'], #config for 3.5km radius
#                 'Lat': [37.73200, 37.73200, 37.78350, 37.78350, 37.76554], 
#                 'Lon': [-122.43000, -122.47750, -122.41200, -122.46000, -122.48750]
#                 }
# depot_loc=Point(-122.399483, 37.745084)
# #recharge_loc_df = pd.DataFrame(recharge_loc)
# geometry = [Point(x,y) for x,y in zip(recharge_loc['Lon'], recharge_loc['Lat'])]
# recharge_gdf = gpd.GeoDataFrame(recharge_loc, crs='EPSG:4326', geometry=geometry)





# delivery_loc=Point(-122.42673, 37.74518)

# depot_loc = (37.745084, -122.399483)
# delivery_loc = (37.74518, -122.42673)


# delivery_dict = {'Loc': 'Delivery', 'Lat': 37.77872, 'Lon': -122.50594}

# loc_dict = {'Loc': ['Depot', '1', '2', '3', '4', '5'],
#             'Lat': [37.745084, 37.73200, 37.73200, 37.78350, 37.78350, 37.76554],
#             'Lon': [-122.399483, -122.43000, -122.47750, -122.41200, -122.46000, -122.48750]}

# full_dict = loc_dict | delivery_dict

# loc_df = pd.DataFrame(loc_dict)


# edges_all = []

# for i in range(len(loc_df)):
#     edges = []
#     start = (loc_df.iloc[i]['Lat'], loc_df.iloc[i]['Lon'])
#     for j in range(len(loc_df)):
#         end = (loc_df.iloc[j]['Lat'], loc_df.iloc[j]['Lon'])
#         dist = geodesic(start, end).km
        
#         if (loc_df.iloc[j]['Loc'] == 'Delivery'):
#             if dist > radius:
#                 edges.append(1000)  #make it so that you have to deliver within 3.5km radius  
#             else:
#                 edges.append(dist)
#         else:
#             if dist > 2*radius:
#                 edges.append(1001)
#             else:            
#                 edges.append(dist)


#     edges_all.append(edges)


# M = np.array(edges_all) #https://stackoverflow.com/questions/53074947/examples-for-search-graph-using-scipy

# D, Pr = shortest_path(M, directed=True, method='auto', return_predecessors=True)

# get_path(Pr, 0, 6)







# trip_route=['s']

# #update the current location of the drone to the delivery location
# current_loc = delivery_loc

# #Is the delivery loc within range of the depot?
# dist_to_depot = geodesic(depot_loc, current_loc).km
# if dist_to_depot <= radius:
#     trip_route.append('d')
#     shortest_trip = dist_to_depot
    
# else:
#     while 'd' not in trip_route:
#         #how far are each of the recharging spots from the current location?
#         dist_recharge = []
#         for i in range(len(recharge_gdf)):
#             recharge = (recharge_gdf)
            
            
#             print(i)



