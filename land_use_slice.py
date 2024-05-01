import xarray as xr
import time 

start = time.time()

# Boundaries of the UK
min_lon, max_lon = -10, 2

# Path to files
file_path = '/Users/jamesquessy/Desktop/Uni Work/Masters/Reasearch Project/Code/Power_Generation/land_use/land_use.nc'
output_path = '/Users/jamesquessy/Desktop/Uni Work/Masters/Reasearch Project/Code/Power_Generation/land_use/sliced_land_use_uk.nc'

# Open the file
with xr.open_dataset(file_path) as ds:
    # Check the range of latitude and longitude in the file
    print("Actual latitude range:", ds.lat.min().values, "to", ds.lat.max().values)

    # Select the UK
    sliced_ds1 = ds.sel(lon=slice(min_lon, max_lon))

    # Save the UK to a new file
    sliced_ds1.to_netcdf(output_path)

end = time.time()
print(f'elapsed time is {end - start} seconds')
