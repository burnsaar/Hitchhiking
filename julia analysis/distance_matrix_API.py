

import googlemaps
from datetime import datetime, timedelta, date
import math
import read_julia_data as read_jl
import pandas as pd
import pickle
import os

if __name__ == '__main__':

#--------------Load Julia data---------------------------------------------------------------------------------------------

    #set the results parent folder
    file_path = 'C:/Users/Aaron/AppData/Local/Programs/Julia-1.6.7/MultiAgentAllocationTransit.jl/results/2024-02-16 (d_100_s_100_iter_100_2281_sites)'


    #load the states for each trip, t
    t_file_path = read_jl.data_paths(file_path + '/states/*.dat')
    trips = read_jl.load_trips(t_file_path)

    #compile each of the travel time from the depot to the site
    trip_time_df, trips_gdf = read_jl.compile_trip_time(trips)

    #good code for the batch section
    # destinations_ls = list(zip(trip_time_df['Lat'], trip_time_df['Lon']))

    # destinations_ls = []
    # for i in range(len(trip_time_df)):
    #     destinations_ls.append(str(trip_time_df.iloc[i]['Lat']) + " " + str(trip_time_df.iloc[i]['Lon']))

    # print(destinations_ls)
    # print(len(destinations_ls))
    # print('lol')

    # #break the destination list into chunks to meet the 25 location API requirement
    # chunk_size = 3
    # chunks = [destinations_ls[i:i + chunk_size] for i in range(0, len(destinations_ls), chunk_size)] 



#--------------Connect to the API---------------------------------------------------------------------------------------------

    API_key = googlemaps.Client(key="AIzaSyCKm3v5-ev_n9blU7qt8gwraajN28zJ01A")
    #origins = ["Traverse City, MI"]
    #destinations = ["Old Mission, MI"]

    origins_ls = []
    destinations_ls = []
    lat_ls = []
    lon_ls = []
    response_ls = []
    distance_ls = []
    duration_ls = []
    mode_ls = []
    
    #both manually calculating the time since epoch 1970 and datetime give you the total seconds since Jan 1 1970 which is needed for the API
    # dt = datetime(2024, 5, 15, 10)
    # epoch_time = datetime(1970, 1, 1)
    # datetime_UTC = (dt - epoch_time).total_seconds()
    # print(datetime_UTC)
    
    API_date = datetime(2024, 5, 15, 17) #should be 15 May 2024 at 10am or 5pm as the second case

    for i in range(len(trip_time_df)):
        for mode in ['driving', 'bicycling']:

            origins = ["37.744045, -122.4012"]
            destinations = [str(trip_time_df.iloc[i]['Lat']) + " " + str(trip_time_df.iloc[i]['Lon'])]

            if mode == 'driving': #only the driving mode can include a traffic_model parameter
                response = googlemaps.client.distance_matrix(client=API_key, 
                                                        origins=origins, 
                                                        destinations=destinations,
                                                        mode=mode,
                                                        departure_time=API_date,
                                                        traffic_model="pessimistic", #"best guess"
                                                        units='metric')
            else:
                response = googlemaps.client.distance_matrix(client=API_key, 
                                                        origins=origins, 
                                                        destinations=destinations,
                                                        mode=mode,
                                                        departure_time=API_date,
                                                        units='metric')
            
            origins_ls.append(origins)
            destinations_ls.append(destinations)
            lat_ls.append(trip_time_df.iloc[i]['Lat'])
            lon_ls.append(trip_time_df.iloc[i]['Lon'])
            response_ls.append(response)
            distance_ls.append(response['rows'][0]['elements'][0]['distance']['value']) #distance in meters
            duration_ls.append(response['rows'][0]['elements'][0]['duration']['value']) #time in seconds
            mode_ls.append(mode)

            
            print(i)





    response_dict = {'origins': origins_ls, 'destinations': destinations_ls, 'Lat': lat_ls, 'Lon': lon_ls, 'response': response_ls, 'mode': mode_ls, 'distance': distance_ls, 'duration': duration_ls}
    #print(response_dict)
    df = pd.DataFrame(response_dict)
    #print(df)


    #record the data
    base_path = 'C:/Users/Aaron/Documents/GitHub/Hitchhiking/Hitchhiking/google maps API/results'
    current_date = date.today().strftime('%Y-%m-%d') 

    folder_path = os.path.join(base_path, current_date)
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)


    save_file = os.path.join(folder_path, 'API_results_pessimistic_1700.dat')

    with open(save_file, 'wb') as file:
        pickle.dump(response_dict, file)


# #this batch method works well!  But, the output might be confusing
# total_batches = math.ceil(len(destinations_ls)/25) #API can only consider 25 destination locations per request
# for chunk in chunks:
#     origins = ["37.744045, -122.4012"]
#     destinations = chunk

#     response = googlemaps.client.distance_matrix(client=API_key, 
#                                          origins=origins, 
#                                          destinations=destinations,
#                                          mode='driving',
#                                          units='imperial')
    
#     print(response)



# origins = ["44.756489, -85.619521"]  #44.756489, -85.619521,    
# destinations = ["(44.962540, -85.482259)"]  #44.962540, -85.482259

# test = googlemaps.client.distance_matrix(client=API_key, 
#                                          origins=origins, 
#                                          destinations=destinations,
#                                          mode='driving',
#                                          units='imperial')






