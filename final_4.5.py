# Section 1: Import and Setup

# Subsection 1.1: Import libaries
import geopandas as gpd
import rasterio
from rasterio.features import rasterize
from rasterio.transform import from_origin
import xarray as xr
import numpy as np
import rioxarray
import os
import netCDF4 as nc
from netCDF4 import Dataset
from geopy.distance import great_circle
import pandas as pd
import simplekml

# Subsection 1.2: Directory Setup
# Base directory and subdirectories
base_directory = '/Users/jamesquessy/Developer/Projects/Masters'
population_directory = os.path.join(base_directory, 'Data/Population')
last_year_avg_directory = os.path.join(base_directory, 'Data/last_year_avg/RCP_4.5')
merged_directory = os.path.join(base_directory, 'RCP_4.5/Code/Merged_Files')
final_files_directory = os.path.join(base_directory, 'RCP_4.5/Code/final_files')
raster_file_directory = os.path.join(base_directory, 'Data/Raster_Data/Raw_Data')


# File paths for orography, land area, and land use data.
orography_file_path = os.path.join(base_directory, 'Data/Raster_Data/Orogrophy/orography_remap.nc')
land_area_file_path = os.path.join(base_directory, 'Data/Raster_Data/Land_Area/land_area_remap.nc')
land_use_file_path = os.path.join(base_directory, 'Data/Raster_Data/land_use/remaped_land.nc')

# File paths for raster files.
airport_mask_file_path = os.path.join(raster_file_directory, 'airport_mask.nc')
spa_mask_file_path = os.path.join(raster_file_directory, 'spa_raster_NetCDF.nc')
nsa_mask_file_path = os.path.join(raster_file_directory, 'nsa_raster_NetCDF.nc')

# Subsection 1.3: Check and Create Directories
# Check all directories exist
os.makedirs(population_directory, exist_ok=True)
os.makedirs(last_year_avg_directory, exist_ok=True)
os.makedirs(merged_directory, exist_ok=True)
os.makedirs(final_files_directory, exist_ok=True)
os.makedirs(raster_file_directory, exist_ok=True)


# Subsection 1.4: Constants
# Years for analysis and variables for climate data
years = ['2020', '2050', '2075', '2099'] 
variables = ['hurs', 'ps', 'sfcWind', 'tas']

# Constants for wind turbine calculations
turbine_area = 2000  # Turbine area in square meters
power_coefficient = 0.35  # Turbine power coefficient
reference_height = 10  # Reference height for wind speed measurement (in meters)
target_height = 80  # Target height for wind speed estimation (in meters)

# Physical constants
Rd = 287.05  # Specific gas constant for dry air (J/kg·K).
Rv = 461.5  # Specific gas constant for water vapor (J/kg·K).
Kelvin = 273.15  # Conversion from Celsius to Kelvin.
power_loss_per_1000km = 0.0035  # Power loss per 1000 km.
days_per_year = 365  # Number of days in a non-leap year.
hours_per_year = 8760  # Number of hours in a non-leap year.
air_density = 1.225  # Air density at sea level (kg/m³).
swept_area = 2000  # Area swept by wind turbine blades (m²).
rated_wind_speed = 14  # Rated wind speed for turbine power calculations (m/s).

# Section 3: Weather Analysis
# Subsection 3.1: Wind speed functions

def calculate_wind_at_80m(wind_speed_10m, friction_coefficient, reference_height, target_height):
    """
    Calculate wind speed at 80 meters
    
    Parameters:
    - wind_speed_10m: Wind speed measured at 10 meters.
    - friction_coefficient: Surface friction coefficient.
    - reference_height: The height at which the reference wind speed is measured.
    - target_height: The height for which the wind speed is to be estimated.

    Returns:
    - Estimated wind speed at 80 meters.
    """
    return wind_speed_10m * (np.log(target_height / friction_coefficient) / np.log(reference_height / friction_coefficient))

def calculate_saturation_vapor_pressure(t):
    """
    Saturation vapor pressure based on temperature.
    
    Parameters:
    - t: Temperature in degrees Celsius.

    Returns:
    - Saturation vapor pressure in Pascals.
    """
    return 6.1094 * np.exp((17.625 * t) / (t + 243.04)) * 100

def calculate_vapor_pressure(t, rh):
    """
    Calculate actual vapor pressure based on temperature and relative humidity.
    
    Parameters:
    - t: Temperature in degrees Celsius.
    - rh: Relative humidity in percentage.

    Returns:
    - Actual vapor pressure in Pascals.
    """
    es = calculate_saturation_vapor_pressure(t)
    return (rh / 100.0) * es

def calculate_air_density(ps, tas, rh, Rd, Rv, Kelvin):
    """
    Calculate air density at surface level.
    
    Parameters:
    - ps: Surface pressure in Pascals.
    - tas: Air temperature in Kelvin.
    - rh: Relative humidity in percentage.
    - Rd: Specific gas constant for dry air (J/kg·K).
    - Rv: Specific gas constant for water vapor (J/kg·K).
    - Kelvin: Conversion constant from Celsius to Kelvin.

    Returns:
    - Air density at the surface level in kg/m³.
    """
    temp_celsius = tas - Kelvin
    e = calculate_vapor_pressure(temp_celsius, rh)
    Pd = ps - e
    return (Pd / (Rd * tas)) + (e / (Rv * tas))

def calculate_power_generation(wind_80m, air_density, turbine_area, power_coefficient):
    """
    Calculate power generation for a single wind turbine.
    
    Parameters:
    - wind_80m: Wind speed at 80 meters.
    - air_density: Air density in kg/m³.
    - turbine_area: Area covered by the wind turbine in square meters.
    - power_coefficient: Power coefficient of the turbine.

    Returns:
    - Power generation in kilowatts.
    """
    wind_power = 0.5 * air_density * turbine_area * (wind_80m ** 3) * power_coefficient
    return wind_power / 1000  # Convert to kW


# Section 4: Data Processing

def merge_datasets(year, nsa_mask_file_path, airport_mask_file_path, spa_mask_file_path):
    """
    Merge climate datasets and apply a mask from a NetCDF file.

    Parameters:
    - year: The year for which the datasets are to be merged.
    - mask_netcdf_path: The file path of the mask NetCDF file.

    Returns:
    - The file path of the merged NetCDF dataset.
    
    """

    # Load datasets
    orography_ds = xr.open_dataset(orography_file_path)
    land_area_ds = xr.open_dataset(land_area_file_path)
    land_use_ds = xr.open_dataset(land_use_file_path)
    datasets = [orography_ds, land_area_ds, land_use_ds]

    # Append climate data
    for variable in variables:
        file_path = os.path.join(last_year_avg_directory, f"{variable}_{year}_yearly_avg.nc")
        if os.path.exists(file_path):
            ds = xr.open_dataset(file_path)
            if 'height' in ds:
                ds = ds.drop_vars('height') 
            datasets.append(ds)

    # Merge datasets 
    merged_ds = xr.merge(datasets)
    if 'sfcWind' in merged_ds and 'friction_coefficient' in merged_ds:
        merged_ds['wind_80m'] = calculate_wind_at_80m(
            merged_ds['sfcWind'], merged_ds['friction_coefficient'], reference_height, target_height
        )
    if 'ps' in merged_ds and 'tas' in merged_ds and 'hurs' in merged_ds:
        merged_ds['air_density'] = calculate_air_density(
            merged_ds['ps'], merged_ds['tas'], merged_ds['hurs'], Rd, Rv, Kelvin
        )
    if 'wind_80m' in merged_ds and 'air_density' in merged_ds:
        merged_ds['power_generation'] = calculate_power_generation(
            merged_ds['wind_80m'], merged_ds['air_density'], turbine_area, power_coefficient
        )

    # Apply the masks
    mask_ds = xr.open_dataset(nsa_mask_file_path)
    mask_aligned = mask_ds['mask'].reindex_like(merged_ds['power_generation'], method='nearest')

    airport_mask_ds = xr.open_dataset(airport_mask_file_path)
    airport_mask_aligned = airport_mask_ds['airport'].reindex_like(merged_ds['power_generation'], method='nearest')
    
    # Load the Special Protection Area mask 
    special_mask_ds = xr.open_dataset(spa_mask_file_path)
    special_mask_aligned = special_mask_ds['mask'].reindex_like(merged_ds['power_generation'], method='nearest')

    # Apply NSA, airport, and SPA mask
    merged_ds['power_generation'] = merged_ds['power_generation'].where(
        (mask_aligned == 0) & (airport_mask_aligned == 0) & (special_mask_aligned == 0), 0)
    
    # Save the dataset
    merged_file_path = os.path.join(merged_directory, f"Merged_{year}.nc")
    merged_ds.to_netcdf(merged_file_path)
    print(f"Merged file for {year} saved at {merged_file_path}")

    return merged_file_path

# Section 5: Data Processing and Analysis

# Loop over each year
for year in years:
    # Merge datasets
    merged_file_path = merge_datasets(year, nsa_mask_file_path, airport_mask_file_path, spa_mask_file_path)
    
    # Drop unnecessary variables
    essential_var_file_path = os.path.join(merged_directory, f"essential_var_{year}.nc")
    if os.path.exists(merged_file_path):
        ds = xr.open_dataset(merged_file_path)
        # Variables not essential for further analysis
        ds = ds.drop_vars([
            "air_density", "change_count", 'friction_coefficient', 'hurs',
            'current_pixel_state', 'observation_count', 'orog', 'processed_flag',
            'ps', 'sfcWind', 'sftlf', 'tas', 'time', 'time_bnds', 'wind_80m'
        ])
        
        # Save dataset
        if not os.path.exists(essential_var_file_path):
            ds.to_netcdf(essential_var_file_path)
            print(f"Essential variables saved for {year}")
        else:
            print(f"Essential variables file already exists for {year}")
        ds.close()
    else:
        print(f"Failed to process file for {year}")

    # Apply land use masks
    if os.path.exists(essential_var_file_path):
        dataset = Dataset(essential_var_file_path, 'r+')
        lccs_class = dataset.variables['lccs_class'][:]
        power_generation = dataset.variables['power_generation'][:]
        urban_mask = lccs_class == 5
        water_mask = lccs_class == 2
        exclusion_mask = np.logical_or(urban_mask, water_mask)
        power_generation_masked = np.ma.array(power_generation, mask=exclusion_mask)
        dataset.variables['power_generation'][:] = power_generation_masked
        dataset.sync()
        dataset.close()
        print(f"Masking applied and saved for {year}")
    else:
        print(f"Failed to apply masks for {year}")

    # Replace NaN values and save
    final_file_path = os.path.join(final_files_directory, f"final_file_{year}.nc")
    if os.path.exists(essential_var_file_path):
        ds = xr.open_dataset(essential_var_file_path)
        for var in ds.variables:
            if ds[var].dtype.kind in 'f':
                ds[var] = ds[var].fillna(0)
        if not os.path.exists(final_file_path):
            ds.to_netcdf(final_file_path)
            print(f"All NaN Values removed and saved in 'final_files' directory for {year}")
        else:
            print(f"Final file already exists in 'final_files' directory for {year}")
        ds.close()
    else:
        print(f"Failed to replace NaN values for {year}")

# Calculate maximum power output at rated wind speed
P_rated = 0.5 * air_density * swept_area * power_coefficient * rated_wind_speed**3
P_rated_kW = P_rated / 1000  # Convert to kilowatts (kW)
max_annual_output = P_rated_kW * hours_per_year 

# Calculate power loss over distance
def calculate_power_loss(power, distance):
    """
    Calculate the power loss over a given distance due to transmission losses.

    Parameters:
    - power: The initial power in kilowatts (kW).
    - distance: The distance the power is transmitted (in meters).

    Returns:
    - The power after accounting for the loss over the distance.
    """
    distance_km = distance / 1000
    loss_fraction = 1 - (power_loss_per_1000km * (distance_km // 1000))
    return power * loss_fraction

# Process data
all_years_top_locations = pd.DataFrame()
all_years_top_locations_no_demand = pd.DataFrame()

for year in years:
    file_path = os.path.join(final_files_directory, f'final_file_{year}.nc')
    dataset = nc.Dataset(file_path)

    # Extract wind power data
    lon = dataset.variables['lon'][:]
    lat = dataset.variables['lat'][:]
    power_generation = dataset.variables['power_generation'][:,:,0].filled(np.nan)

    # Load city energy demand data
    energy_demand_df = pd.read_csv(os.path.join(population_directory, f'city_power_demand_projection_{year}.csv'))

    # DataFrame to store results
    top_locations = pd.DataFrame()

    # Loop over each city
    for index, row in energy_demand_df.iterrows():
        city_name = row['City']
        city_coords = (row['Latitude'], row['Longitude'])
        city_energy_demand_annual = row['Energy Demand (kWh)']

        # Store all locations with power generation
        all_locations = []

        # Loop over each grid point to find potential locations
        for i in range(len(lat)):
            for j in range(len(lon)):
                daily_power_generation = power_generation[i, j]
                if daily_power_generation > 0:
                    wind_farm_coords = (lat[i], lon[j])
                    distance = great_circle(wind_farm_coords, city_coords).kilometers
                    adjusted_daily_power = calculate_power_loss(daily_power_generation, distance)
                    annual_energy_production = adjusted_daily_power * days_per_year
                    
                    # Calculate the annual energy production for the best location
                    annual_energy_production = (adjusted_daily_power * (0.3*24)) * days_per_year
                    
                    # Calculate demand satisfaction percentage
                    demand_satisfaction = (annual_energy_production / city_energy_demand_annual) * 100 if city_energy_demand_annual else 0

                    # Add location and its power generation to the list
                    all_locations.append((adjusted_daily_power, wind_farm_coords, distance, annual_energy_production, demand_satisfaction))

        # Sort locations by power generation and select top 10
        top_10_locations = sorted(all_locations, key=lambda x: x[0], reverse=True)[:10]

        # Add to the DataFrame
        for rank, (power, location, distance, annual_production, satisfaction) in enumerate(top_10_locations, 1):
            new_row = {
                'Year': year,
                'City': city_name,
                'Rank': rank,
                'Lat': location[0],
                'Lon': location[1],
                'Distance_to_City (km)': distance,
                'Adjusted_Daily_Power (kW)': power,
                'Annual_Energy_Production (kWh)': annual_production,
                'City_Energy_Demand (kWh)': city_energy_demand_annual,
                'Demand_Satisfaction (%)': satisfaction,
                'Capacity Factor (%)': (annual_production / max_annual_output) * 100
            }
            top_locations = pd.concat([top_locations, pd.DataFrame([new_row])], ignore_index=True)

    all_years_top_locations = pd.concat([all_years_top_locations, top_locations], ignore_index=True)
    dataset.close()
    print(f"The analysis for {year} has been completed.")
    
    top_power_locations = []

    # Loop over each grid point to find the top locations based on power generation
    for i in range(len(lat)):
        for j in range(len(lon)):
            daily_power_generation = power_generation[i, j]
            if daily_power_generation > 0: 
                wind_farm_coords = (lat[i], lon[j])
                annual_energy_production = daily_power_generation * days_per_year * (0.3 * 24) 
                top_power_locations.append((annual_energy_production, wind_farm_coords))

    # Sort locations by power generation and select top 10
    top_10_power_locations = sorted(top_power_locations, key=lambda x: x[0], reverse=True)[:10]

    # Add to the dataframe
    for rank, (annual_production, location) in enumerate(top_10_power_locations, 1):
        daily_power_generation = annual_production / (days_per_year * 0.3 * 24)

        # Create a new row with required information 
        new_row = {
            'Year': year,
            'Rank': rank,
            'Lat': location[0],
            'Lon': location[1],
            'Daily Power Potential (kW)': daily_power_generation,
            'Annual Energy Production (kWh)': annual_production,
            'Capacity Factor (%)': (annual_production / max_annual_output) * 100
        }

        # Add the new row to the DataFrame
        all_years_top_locations_no_demand = pd.concat([all_years_top_locations_no_demand, pd.DataFrame([new_row])], ignore_index=True)

# Round to one decimal place
all_years_top_locations, all_years_top_locations_no_demand = all_years_top_locations.round(5), all_years_top_locations_no_demand.round(5)

# Save to Excel file
all_years_top_locations.to_excel(os.path.join(base_directory, "RCP_4.5/Code/RCP_4.5_top_locations.xlsx"), index=False)
print("All years processed successfully. Results saved to 'RCP_4.5_top_locations.xlsx'")

# Save new dataframe to a different excel file
all_years_top_locations_no_demand.to_excel(os.path.join(base_directory, "RCP_4.5/Code/RCP_4.5_top_power_locations.xlsx"), index=False)
print("Results for top power generation locations saved to 'RCP_4.5_top_power_locations.xlsx'")

# Section 7: Creating KML Files

def create_kml(df, filename):
    kml = simplekml.Kml()

    for idx, row in df.iterrows():
        pnt = kml.newpoint(name=f"{row['Year']} - Rank {row['Rank']}", 
                           coords=[(row['Lon'], row['Lat'])])
        pnt.description = f"Year: {row['Year']}, Rank: {row['Rank']}"

    kml.save(filename)

# Create and save KML for all_years_top_locations
create_kml(all_years_top_locations, os.path.join(base_directory, "RCP_4.5/Code/top_locations.kml"))

# Create and save KML for all_years_top_locations_no_demand
create_kml(all_years_top_locations_no_demand, os.path.join(base_directory, "RCP_4.5/Code/top_power_locations_no_demand.kml"))
