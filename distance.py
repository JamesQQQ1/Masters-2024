import pandas as pd
from haversine import haversine

# Load the new sites data and existing wind farms data
new_sites_df = pd.read_excel('RCP_2.6_top_power_locations.xlsx')  # Update the file path accordingly
existing_wind_farms_df = pd.read_csv('updated_validation.csv')  # Update the file path accordingly

# Define a function to calculate the distance to the closest wind farm
def find_closest_wind_farm(new_site_lat, new_site_lon, wind_farms_df):
    closest_distance = float('inf')  # Initialize with an infinitely large number
    for index, farm in wind_farms_df.iterrows():
        farm_location = (farm['lat'], farm['lon'])
        new_site_location = (new_site_lat, new_site_lon)
        distance = haversine(farm_location, new_site_location)
        if distance < closest_distance:
            closest_distance = distance
    return closest_distance

# Apply the function to each new site
new_sites_df['Closest Wind Farm Distance (km)'] = new_sites_df.apply(
    lambda row: find_closest_wind_farm(row['Lat'], row['Lon'], existing_wind_farms_df), axis=1
)

# Save the results to a new CSV file
new_sites_df.to_csv('new_sites_with_distances_2.6.csv', index=False, float_format='%g')
