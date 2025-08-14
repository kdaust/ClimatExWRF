import numpy as np
import xarray as xr
import pandas as pd
from glob import glob

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmcrameri.cm as cm


# Constants

START_DATE = np.datetime64('1999-10-01T00:00:00.000000000')
END_DATE = np.datetime64('2023-09-30T23:00:00.000000000')
N_HOURS = (END_DATE - START_DATE) / np.timedelta64(1, 'h') + 1
COVERAGE_THRESHOLD = 0.9
BAD_MONTHS = [[2005, 4],
             [2005, 7],
             [2005, 10],
             [2006, 3],
             [2006, 9],
             [2008, 7],
             [2015, 10],
             [2017, 8],
             [2017, 11],
             [2019, 9],
             [2021, 7],
             [2022, 2],
             [2023, 6]]


# Functions

def readable_time_coord(ds, time_var='Times'):

    time_vals = ds[time_var].values

    # Extract integer (YYYYMMDD) and fractional day parts
    dates_int = time_vals.astype(int)
    frac_days = time_vals - dates_int

    # Parse integers to datetime
    base_dates = pd.to_datetime(dates_int.astype(str), format="%Y%m%d")

    # Add fractional day as timedelta
    full_times = base_dates + pd.to_timedelta(frac_days, unit="D")
    full_times = full_times.round("1h")

    # Assign back to dataset
    ds = ds.assign_coords({time_var: full_times})
    ds = ds.transpose(time_var, ...)

    ds = ds.rename({time_var: "time"})
    
    return ds

# Find the nearest grid point in winds_ds to the station's lat/lon

def get_nearest_wrf_point(ds, lat, lon):
    abs_lat = np.abs(ds['XLAT'] - lat)
    abs_lon = np.abs(ds['XLONG'] - lon)
    dist = abs_lat + abs_lon
    # Find indices of the minimum distance
    idx = np.unravel_index(np.argmin(dist.values), dist.shape)
    return idx


# WRF DATA

topography = xr.open_dataset("/users/rpayne/data/topography/hr/HGT_orig.nc") # USask WRF Topography for now
winds_ds = readable_time_coord(xr.open_dataset("/users/rpayne/data/unproc/WRF-UBC/WS10/WS10_199910_202309.nc"))
winds_ds = winds_ds.sel(time=slice(START_DATE, END_DATE))

# STN DATA

wx_data_fnames = sorted(glob("/users/rpayne/data/unproc/STN_ECCC/api/hourly/*.csv"))
df = pd.concat(
    (pd.read_csv(f, encoding='latin1', low_memory=False) for f in wx_data_fnames),
    ignore_index=True
)
df['Date/Time (LST)'] = pd.to_datetime(df['Date/Time (LST)'], format='%Y-%m-%d %H:%M:%S')
df = df[df['Date/Time (LST)'] >= START_DATE]  # Filter for data to align w WRF
df = df[df['Date/Time (LST)'] <= END_DATE]
df.set_index('Date/Time (LST)', inplace=True)
df.rename(columns={"Temp (Â°C)": "Temp (°C)"}, inplace=True)
df.rename(columns={"Dew Point Temp (Â°C)": "Dew Point Temp (°C)"}, inplace=True)

# Create a dictionary mapping station name to (lat, lon) from df
stn_latlon = {}
stns = df["Station Name"].unique()
for stn in stns:
    stn_rows = df[df["Station Name"] == stn]
    if not stn_rows.empty:
        lat = stn_rows["Latitude (y)"].iloc[0]
        lon = stn_rows["Longitude (x)"].iloc[0]
        stn_latlon[stn] = (lat, lon)
    else:
        stn_latlon[stn] = (None, None)

def get_stn(stn_name, data=df):
    stn = data[data["Station Name"] == stn_name]
    return stn.sort_index()

# Find all stations that meet the coverage threshold for winds
WS_coverage_dict = {}

for i,stn in enumerate(stns):
    print(f"{i+1}/{len(stns)}")
    
    stn_data = get_stn(stn)

    WS_data = stn_data["Wind Spd (km/h)"]
    WS_coverage = np.count_nonzero(~np.isnan(WS_data)) / N_HOURS
    if WS_coverage >= COVERAGE_THRESHOLD:
        WS_coverage_dict[stn] = WS_coverage

print(f"Number of stations with adequate Wind Speed coverage: {len(WS_coverage_dict)}")

# MAKE THE PLOT

for i, stn in enumerate(WS_coverage_dict.keys()):

    # Get station data
    print(f"Loading in {stn}...")

    data = get_stn(stn)
    lat = stn_latlon[stn][0]
    lon = stn_latlon[stn][1]
    start_time = data.index.min().to_datetime64()
    end_time = data.index.max().to_datetime64()
    stn_winds = data["Wind Spd (km/h)"]
    stn_winds = stn_winds[stn_winds.index < '2023-10-01T00:00:00.000000000']

    for year, month in BAD_MONTHS:
        stn_winds = stn_winds[~((stn_winds.index.year == year) & (stn_winds.index.month == month))]

    # Get model data nearest the station
    print(f"Getting nearest model data...")
    nearest_idx = get_nearest_wrf_point(winds_ds, lat, lon)
    wrf_winds = winds_ds.isel(south_north=nearest_idx[0], west_east=nearest_idx[1])

    for year, month in BAD_MONTHS:
        wrf_winds = wrf_winds.sel(time=~((wrf_winds.time.dt.year == year) & (wrf_winds.time.dt.month == month)))

    # turn into numpy arrays, remove points where stn data is NaN, and sort
    a = stn_winds.values
    b = wrf_winds['WS10'].values
    b = np.sort(b[~np.isnan(a)])
    a = np.sort(a[~np.isnan(a)])

    # plot
    print("Generating plot...")
    fig = plt.figure(dpi=150, figsize=(13, 2.5), facecolor='white')

    lon = topography['XLONG'].squeeze()
    lat = topography['XLAT'].squeeze()
    hgt = topography['HGT'].squeeze()

    ax0 = fig.add_subplot(1, 3, 1, projection=ccrs.PlateCarree())
    mesh = ax0.pcolormesh(lon, lat, hgt, transform=ccrs.PlateCarree(), shading='auto')
    ax0.add_feature(cfeature.BORDERS, linewidth=1)   # Country borders
    provinc_bodr = cartopy.feature.NaturalEarthFeature(category='cultural', 
        name='admin_1_states_provinces_lines', scale='50m', facecolor='none', edgecolor='k')
    ax0.add_feature(provinc_bodr, linestyle='--', linewidth=0.6, edgecolor="k", zorder=10)
    lat, lon = stn_latlon[stn]
    ax0.plot(lon, lat, marker='x', color='red', markersize=4, transform=ccrs.PlateCarree())
    ax0.set_xlim(-134, -114)
    ax0.set_ylim(48, 60.5)
    ax0.coastlines(linewidth=0.5)
    ax0.text(-133.5, 58.8, f"{stn}", fontsize=12, transform=ccrs.PlateCarree(), color='red')

    ax1 = fig.add_subplot(1, 3, 2)
    ax1.scatter(a, b, s=20, facecolors='grey', edgecolors='k', marker='o', linewidth=.5)
    ax1.plot([a.min(), a.max()],
            [a.min(), a.max()],
            'k--', lw=.7)
    ax1.set_xlabel("Station Quantile (km/h)", fontsize=7)
    ax1.set_ylabel("WRF Quantile (km/h)", fontsize=7)
    ax1.grid(alpha=.2, ls='--')
    ax1.set_xticks(np.arange(0,a.max()+5,5))
    ax1.set_yticks(np.arange(0,a.max()+5,5))
    ax1.tick_params(axis='both', which='major', labelsize=5)

    ax2 = fig.add_subplot(1, 3, 3)
    ax2.scatter(a / max(a), b / max(b), s=20, facecolors='grey', edgecolors='k', marker='o', linewidth=.5)
    ax2.plot([0, 1],
            [0, 1],
            'k--', lw=.7)
    ax2.set_xlabel("Normalized Station Quantile", fontsize=7)
    ax2.set_ylabel("Normalized WRF Quantile", fontsize=7)
    ax2.grid(alpha=.2, ls='--')
    ax2.set_xticks(np.arange(0,1.1,0.1))
    ax2.set_yticks(np.arange(0,1.1,0.1))
    ax2.tick_params(axis='both', which='major', labelsize=5)

    plt.subplots_adjust(wspace=0.15,)
    plt.savefig(f"/users/rpayne/ClimatExWRF/Robert/results/winds_qq_plots/{stn}.png")
