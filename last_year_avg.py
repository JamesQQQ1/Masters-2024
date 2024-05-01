import os
import xarray as xr

def extract_last_year(file_path, new_folder, var_name, year):

    # Create a folder if it doesn't exist
    if not os.path.exists(new_folder):
        os.makedirs(new_folder)

    # Path to new file
    new_file_path = os.path.join(new_folder, f"{var_name}_{year}_yearly_avg.nc")

    # Load dataset
    data = xr.open_dataset(file_path, engine='netcdf4')

    # Select data for the years
    yearly_data = data.sel(time=data.time.dt.year == int(year))

    # Average for the year
    avg_data = yearly_data[var_name].mean(dim='time')

    # New dataset with the average data
    avg_dataset = xr.Dataset({var_name: avg_data})
    avg_dataset.attrs = data.attrs

    # Save to new file
    avg_dataset.to_netcdf(new_file_path)
    print(f"Saved {new_file_path}")

# Base path
original_file_path = '/Users/jamesquessy/Desktop/Uni Work/Masters/Reasearch Project/Code/Power_Generation /NetCDF_Files'

# New folder path
new_folder = '/Users/jamesquessy/Desktop/Uni Work/Masters/Reasearch Project/Code/Power_Generation/last_year_avg'

# Varialbes used
years = ['2020', '2050', '2075', '2099']
variables = ['tas', 'hurs', 'sfcWind', 'ps']

# Loop through each year and variable
for year in years:
    for var_name in variables:
        file_name = f"{var_name}_{year}_remap.nc" 
        file_path = os.path.join(original_file_path, file_name)
        if os.path.isfile(file_path):
            extract_last_year(file_path, new_folder, var_name, year)
        else:
            print(f"File not found: {file_path}")
