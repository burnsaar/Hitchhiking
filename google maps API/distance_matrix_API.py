

import googlemaps
from datetime import datetime

API_key = googlemaps.Client(key="AIzaSyCKm3v5-ev_n9blU7qt8gwraajN28zJ01A")
#origins = ["Traverse City, MI"]
#destinations = ["Old Mission, MI"]

origins = ["44.756489, -85.619521"]  #44.756489, -85.619521
destinations = ["44.962540, -85.482259"]  #44.962540, -85.482259

test = googlemaps.client.distance_matrix(client=API_key, 
                                         origins=origins, 
                                         destinations=destinations,
                                         mode='driving',
                                         units='imperial')


print(test)
print("hello")


