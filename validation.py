import pandas as pd
from pyproj import Transformer

try:
    # Try and read the CSV with ISO-8859-1 encoding
    df = pd.read_csv('validation.csv', encoding='ISO-8859-1')
except UnicodeDecodeError:
    # If that fails, try 'utf-16' encoding
    df = pd.read_csv('validation.csv', encoding='utf-16')
except Exception as e:
    # If it fails, print the exception and exit
    print(f"An error occurred: {e}")
    exit()

# Transform coordinates from British National Grid to WGS84
transformer = Transformer.from_crs("EPSG:27700", "EPSG:4326")

# Initialize empty lists for lat and lon
df['lat'] = 0.0
df['lon'] = 0.0

# Loop through the DataFrame and transform coordinates
for index, row in df.iterrows():
    lat, lon = transformer.transform(row['X-coordinate'], row['Y-coordinate'])
    df.at[index, 'lat'] = lat
    df.at[index, 'lon'] = lon

# Save the DataFrame to a new CSV file
df.to_csv('updated_validation.csv', index=False, float_format='%g')

