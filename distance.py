import pandas as pd
from haversine import haversine

# Model prediction data and existing wind farms data
new_sites_df = pd.read_excel('RCP_2.6_top_power_locations.xlsx')  
existing_wind_farms_df = pd.read_csv('updated_validation.csv')  

# Calculate the distance to the closest wind farm
def find_closest_wind_farm(new_site_lat, new_site_lon, wind_farms_df):
    closest_distance = float('inf')
    for index, farm in wind_farms_df.iterrows():
        farm_location = (farm['lat'], farm['lon'])
        new_site_location = (new_site_lat, new_site_lon)
        distance = haversine(farm_location, new_site_location)
        if distance < closest_distance:
            closest_distance = distance
    return closest_distance

# Apply to each new site
new_sites_df['Closest Wind Farm Distance (km)'] = new_sites_df.apply(
    lambda row: find_closest_wind_farm(row['Lat'], row['Lon'], existing_wind_farms_df), axis=1
)

# Save to a new CSV file
new_sites_df.to_csv('new_sites_with_distances_2.6.csv', index=False, float_format='%g')
