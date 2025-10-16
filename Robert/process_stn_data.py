import numpy as np
import xarray as xr
import pandas as pd
from glob import glob

# Open the data frame
wx_data_fnames = sorted(glob("/users/rpayne/data/unproc/STN_ECCC/api/hourly/*.csv"))
df = pd.concat(
    (pd.read_csv(f, encoding="utf-8", low_memory=False, index_col=4) for f in wx_data_fnames),
    ignore_index=True
)

# Change the date column to datetime format
df['Date/Time (LST)'] = pd.to_datetime(df['Date/Time (LST)'], format='%Y-%m-%d %H:%M:%S')

# Truncate the data and make sure the index spans
FIRST_YEAR = 1990
LAST_YEAR = 2023
df = df[df['Date/Time (LST)'].dt.year >= FIRST_YEAR]  # Filter for data after FIRST_YEAR
df = df[df['Date/Time (LST)'].dt.year <= LAST_YEAR]  # Filter for data before LAST_YEAR

# Set the datetime column as the index
df.set_index('Date/Time (LST)', inplace=True)

# Expand the index so that they all stns span the same range of time
# Missing data are filled with NaNs
full_index = pd.date_range(f"{FIRST_YEAR}-01-01 00:00:00", f"{LAST_YEAR}-12-31 23:00:00", freq="H")
df = df.reindex(full_index)
df.index.name = 'Date/Time (LST)'

# Remove unnecessary columns
columns_to_keep = ['Longitude (x)', 'Latitude (y)', 'Station Name', 'Temp (deg C)',
       'Dew Point Temp (deg C)', 'Rel Hum (%)', 'Wind Spd (km/h)',  'Stn Press (kPa)',
       'Precip. Amount (mm)',]
df = df[columns_to_keep]

# Save to file
path_out = "/users/rpayne/data/unproc/STN_ECCC/api/hourly/proc/ec_hourly_allstns.csv"
print(f"Saving to {path_out} as csv...")
df.to_csv(path_out)

# Rename wind speed column name to remove slash (required for netcdf save)
# df.rename(columns={"Wind Spd (km/h)": "Wind Spd (kmph)"}, inplace=True)
# path_out = "/users/rpayne/data/unproc/STN_ECCC/api/hourly/proc/ec_hourly_allstns.nc"
# print(f"Saving to {path_out} as netcdf...")
# df.to_xarray().to_netcdf(path_out)
