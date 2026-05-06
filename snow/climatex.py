"""
Snow Data Extraction and Processing Module

This module provides utilities for extracting and processing snow-related data
from multiple sources (WRF, ERA5-Land, CaSR) for specified station locations
over defined time periods. It supports extraction of snow depth/water equivalent
at target grid points and their 8 neighboring cells.

Data Sources:
-------------
WRF: Daily mean snow data in monthly netCDF files organized by water year
ERA5-Land: Single netCDF file with daily snow depth data
CaSR: Yearly netCDF files with snow water equivalent data

Water Year Convention:
---------------------
Water years run from October (previous calendar year) to September (current calendar year).
Example: Water year 1999 = October 1998 to September 1999
"""
import xarray as xr
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Tuple
import warnings
warnings.filterwarnings('ignore')

# from joblib import Parallel, delayed

def process_all_stations_neighbors(stations_df: pd.DataFrame, base_path: str = '.') -> dict:
    """
    Process all stations and extract SNOW data for the target grid point
    plus its 8 adjacent neighbors.

    Parameters:
    -----------
    stations_df : pd.DataFrame
        DataFrame containing station information
    base_path : str
        Base directory containing the year folders

    Returns:
    --------
    dict
        Dictionary with station_id as keys and DataFrames as values.
        Each DataFrame contains the center cell SNOW data in the first column,
        followed by 8 columns for the neighboring cells.
    """
    results = {}

    for idx, row in stations_df.iterrows():
        station_id = row['station_id']
        print(f"\n{'='*60}")
        print(f"Processing station: {station_id}")
        print(f"{'='*60}")

        df_snow = extract_snow_timeseries_with_neighbors(
            lat=row['lat'],
            lon=row['lon'],
            initial_year=int(row['initial_year']),
            final_year=int(row['final_year']),
            base_path=base_path
        )

        if not df_snow.empty:
            results[station_id] = df_snow
            print(f"  Extracted {len(df_snow)} time steps")
            print(f"  Date range: {df_snow.index.min()} to {df_snow.index.max()}")
            print(f"  Columns: {list(df_snow.columns)}")
        else:
            print(f"  No data extracted for {station_id}")

    return results

# def process_all_stations_neighbors(stations_df: pd.DataFrame, base_path: str = '.', n_jobs=20) -> dict:
#     """
#     Process all stations and extract SNOW data for the target grid point
#     plus its 8 adjacent neighbors - IN PARALLEL.

#     Parameters:
#     -----------
#     stations_df : pd.DataFrame
#         DataFrame containing station information
#     base_path : str
#         Base directory containing the year folders
#     n_jobs : int
#         Number of parallel jobs (-1 means use all processors)

#     Returns:
#     --------
#     dict
#         Dictionary with station_id as keys and DataFrames as values.
#     """

#     def process_single_station(idx, row, base_path):
#         """Helper function to process a single station"""
#         station_id = row['station_id']
#         print(f"\n{'='*60}")
#         print(f"Processing station: {station_id}")
#         print(f"{'='*60}")

#         df_snow = extract_snow_timeseries_with_neighbors(
#             lat=row['lat'],
#             lon=row['lon'],
#             initial_year=int(row['initial_year']),
#             final_year=int(row['final_year']),
#             base_path=base_path
#         )

#         if not df_snow.empty:
#             print(f"  Extracted {len(df_snow)} time steps")
#             print(f"  Date range: {df_snow.index.min()} to {df_snow.index.max()}")
#             print(f"  Columns: {list(df_snow.columns)}")
#             return (station_id, df_snow)
#         else:
#             print(f"  No data extracted for {station_id}")
#             return (station_id, None)

#     # Run in parallel
#     results_list = Parallel(n_jobs=n_jobs)(
#         delayed(process_single_station)(idx, row, base_path)
#         for idx, row in stations_df.iterrows()
#     )

#     # Convert list of tuples to dictionary, filtering out None values
#     results = {station_id: df for station_id, df in results_list if df is not None}

#     return results


def extract_snow_timeseries_with_neighbors(lat: float, lon: float, initial_year: int,
                                           final_year: int, base_path: str = '.') -> pd.DataFrame:
    """
    Extract SNOW time series for a given location and its 8 neighboring grid cells.

    Parameters:
    -----------
    lat : float
        Latitude of the station
    lon : float
        Longitude of the station
    initial_year : int
        Initial water year (e.g., 1999 means water year starting Oct 1999)
    final_year : int
        Final water year (inclusive)
    base_path : str
        Base directory containing the year folders

    Returns:
    --------
    pd.DataFrame
        DataFrame with datetime index and SNOW values from center and 8 neighbors.
        Columns: 'SNOW_center', 'SNOW_n', 'SNOW_ne', 'SNOW_e', 'SNOW_se',
                 'SNOW_s', 'SNOW_sw', 'SNOW_w', 'SNOW_nw'
    """
    base_path = Path(base_path)
    all_data = []
    center_indices = None

    # Define neighbor offsets (relative to center cell)
    # Order: N, NE, E, SE, S, SW, W, NW
    neighbor_offsets = {
        'N': (-1, 0),
        'NE': (-1, 1),
        'E': (0, 1),
        'SE': (1, 1),
        'S': (1, 0),
        'SW': (1, -1),
        'W': (0, -1),
        'NW': (-1, -1)
    }

    # Loop through water years
    for water_year in range(initial_year, final_year + 1):
        print(f"Processing water year {water_year}...")

        # Water year goes from October (year-1) to September (year)
        months_years = []

        # October, November, December of previous year
        for month in [10, 11, 12]:
            months_years.append((water_year - 1, month))

        # January through September of current water year
        for month in range(1, 10):
            months_years.append((water_year, month))

        # Process each month
        for year, month in months_years:
            # Construct filename
            filename = f"SNOW_d03_metgrid_{year}_{month:02d}_daily_mean.nc"
            filepath = base_path / str(water_year) / filename

            if not filepath.exists():
                print(f"  Warning: File not found: {filepath}")
                continue

            try:
                # Open dataset
                ds = xr.open_dataset(filepath)

                # Find nearest grid point (only need to do this once)
                if center_indices is None:
                    sn_idx, we_idx = find_nearest_grid_point(ds, lat, lon)
                    center_indices = (sn_idx, we_idx)
                    actual_lat = float(ds.XLAT.isel(south_north=sn_idx, west_east=we_idx).values)
                    actual_lon = float(ds.XLONG.isel(south_north=sn_idx, west_east=we_idx).values)
                    print(f"  Target coordinates: ({lat:.4f}, {lon:.4f})")
                    print(f"  Nearest grid point: ({actual_lat:.4f}, {actual_lon:.4f})")
                    print(f"  Grid indices: south_north={sn_idx}, west_east={we_idx}")

                sn_idx, we_idx = center_indices

                # Extract SNOW data for center cell
                snow_center = ds.SNOW.isel(south_north=sn_idx, west_east=we_idx)
                df_month = snow_center.to_dataframe().reset_index()
                df_month = df_month[['Times', 'SNOW']]
                df_month.rename(columns={'Times': 'datetime', 'SNOW': 'SNOW_center'}, inplace=True)

                # Extract SNOW data for each neighbor
                for neighbor_name, (offset_sn, offset_we) in neighbor_offsets.items():
                    neighbor_sn = sn_idx + offset_sn
                    neighbor_we = we_idx + offset_we

                    # Check if neighbor indices are within bounds
                    if (0 <= neighbor_sn < ds.dims['south_north'] and
                        0 <= neighbor_we < ds.dims['west_east']):
                        snow_neighbor = ds.SNOW.isel(south_north=neighbor_sn, west_east=neighbor_we)
                        df_month[f'SNOW_{neighbor_name}'] = snow_neighbor.values
                    else:
                        print(f"  Warning: Neighbor {neighbor_name} at indices ({neighbor_sn}, {neighbor_we}) out of bounds")
                        df_month[f'SNOW_{neighbor_name}'] = np.nan

                all_data.append(df_month)

                ds.close()

            except Exception as e:
                print(f"  Error processing {filepath}: {e}")
                continue

    # Combine all data
    if not all_data:
        print("No data was extracted!")
        return pd.DataFrame()

    df_final = pd.concat(all_data, ignore_index=True)
    df_final = df_final.sort_values('datetime').reset_index(drop=True)
    df_final.set_index('datetime', inplace=True)

    return df_final


def extract_snow_timeseries_era5(lat: float, lon: float, initial_year: int,
                                  final_year: int, filepath: str) -> pd.DataFrame:
    """
    Extract snow depth time series from ERA5-Land netCDF file.

    Parameters:
    -----------
    lat : float
        Latitude of the location
    lon : float
        Longitude of the location
    initial_year : int
        Initial year to extract
    final_year : int
        Final year to extract (inclusive)
    filepath : str
        Path to ERA5-Land netCDF file

    Returns:
    --------
    pd.DataFrame
        DataFrame with datetime index and 'sd' (snow depth) column
    """
    # Open dataset
    ds = xr.open_dataset(filepath)

    # Find nearest grid point using 2D coordinates
    # Calculate distances to all grid points
    lat_diff = (ds.XLAT - lat) ** 2
    lon_diff = (ds.XLONG - lon) ** 2
    distance = np.sqrt(lat_diff + lon_diff)

    # Find indices of minimum distance
    min_idx = distance.argmin()
    south_north_idx, west_east_idx = np.unravel_index(
        min_idx.values, distance.shape
    )

    actual_lat = float(ds.XLAT.isel(south_north=south_north_idx,
                                     west_east=west_east_idx).values)
    actual_lon = float(ds.XLONG.isel(south_north=south_north_idx,
                                      west_east=west_east_idx).values)

    print(f"Target coordinates: ({lat:.4f}, {lon:.4f})")
    print(f"Nearest grid point: ({actual_lat:.4f}, {actual_lon:.4f})")

    # Extract snow data for grid point
    snow_data = ds.sd.isel(south_north=south_north_idx, west_east=west_east_idx)
    df = snow_data.to_dataframe().reset_index()
    df = df[['valid_time', 'sd']].rename(columns={'valid_time': 'datetime'})
    df['sd'] = df['sd'] * 1000  # m to mm

    # Filter by year range
    df['year'] = df['datetime'].dt.year
    df = df[(df['year'] >= initial_year) & (df['year'] <= final_year)]
    df = df.drop('year', axis=1)

    # Set datetime as index and sort
    df.set_index('datetime', inplace=True)
    df = df.sort_index()

    ds.close()

    return df


def find_nearest_grid_point(ds: xr.Dataset, target_lat: float, target_lon: float) -> Tuple[int, int]:
    """
    Find the nearest grid point to the target lat/lon coordinates.

    Parameters:
    -----------
    ds : xarray.Dataset
        Dataset containing XLAT and XLONG coordinates
    target_lat : float
        Target latitude
    target_lon : float
        Target longitude

    Returns:
    --------
    Tuple[int, int]
        Indices (south_north, west_east) of nearest grid point
    """
    # Calculate distance from target point
    distance = np.sqrt((ds.XLAT - target_lat)**2 + (ds.XLONG - target_lon)**2)

    # Find indices of minimum distance
    min_idx = distance.argmin(dim=['south_north', 'west_east'])
    south_north_idx = int(min_idx['south_north'].values)
    west_east_idx = int(min_idx['west_east'].values)

    return south_north_idx, west_east_idx

def extract_snow_timeseries(lat: float, lon: float, initial_year: int,
                            final_year: int, base_path: str = '.') -> pd.DataFrame:
    """
    Extract SNOW time series for a given location and time period.

    Parameters:
    -----------
    lat : float
        Latitude of the station
    lon : float
        Longitude of the station
    initial_year : int
        Initial water year (e.g., 1999 means water year starting Oct 1999)
    final_year : int
        Final water year (inclusive)
    base_path : str
        Base directory containing the year folders

    Returns:
    --------
    pd.DataFrame
        DataFrame with datetime index and SNOW values
    """
    base_path = Path(base_path)
    all_data = []

    # Loop through water years
    for water_year in range(initial_year, final_year + 1):
        print(f"Processing water year {water_year}...")

        # Water year goes from October (year-1) to September (year)
        # Months for this water year
        months_years = []

        # October, November, December of previous year
        for month in [10, 11, 12]:
            months_years.append((water_year - 1, month))

        # January through September of current water year
        for month in range(1, 10):
            months_years.append((water_year, month))

        # Process each month
        for year, month in months_years:
            # Construct filename
            filename = f"SNOW_d03_metgrid_{year}_{month:02d}_daily_mean.nc"
            filepath = base_path / str(water_year) / filename

            if not filepath.exists():
                print(f"  Warning: File not found: {filepath}")
                continue

            try:
                # Open dataset
                ds = xr.open_dataset(filepath)

                # Find nearest grid point (only need to do this once)
                if not all_data:  # First file
                    sn_idx, we_idx = find_nearest_grid_point(ds, lat, lon)
                    actual_lat = float(ds.XLAT.isel(south_north=sn_idx, west_east=we_idx).values)
                    actual_lon = float(ds.XLONG.isel(south_north=sn_idx, west_east=we_idx).values)
                    print(f"  Target coordinates: ({lat:.4f}, {lon:.4f})")
                    print(f"  Nearest grid point: ({actual_lat:.4f}, {actual_lon:.4f})")
                    print(f"  Grid indices: south_north={sn_idx}, west_east={we_idx}")

                # Extract SNOW data for this grid point
                snow_data = ds.SNOW.isel(south_north=sn_idx, west_east=we_idx)

                # Convert to dataframe
                df_month = snow_data.to_dataframe().reset_index()
                df_month = df_month[['Times', 'SNOW']]
                df_month.rename(columns={'Times': 'datetime'}, inplace=True)

                all_data.append(df_month)

                ds.close()

            except Exception as e:
                print(f"  Error processing {filepath}: {e}")
                continue

    # Combine all data
    if not all_data:
        print("No data was extracted!")
        return pd.DataFrame()

    df_final = pd.concat(all_data, ignore_index=True)
    df_final = df_final.sort_values('datetime').reset_index(drop=True)
    df_final.set_index('datetime', inplace=True)

    return df_final


def extract_snow_timeseries_casr(lat: float, lon: float, initial_year: int,
                                  final_year: int, base_path: str = '.') -> pd.DataFrame:
    """
    Extract snow water equivalent time series from CaSR netCDF files (yearly format).
    """
    base_path = Path(base_path)
    all_data = []
    lat_idx = lon_idx = None
    actual_lat = actual_lon = None

    for year in range(initial_year, final_year + 1):
        filename = f"swe_{year}.nc"
        filepath = base_path / filename

        if not filepath.exists():
            print(f"File not found: {filepath}")
            continue

        try:
            ds = xr.open_dataset(filepath)

            # Find nearest grid point only on first successful file
            if lat_idx is None:
                # Calculate distances to all grid points
                lat_diff = (ds['XLAT'] - lat) ** 2
                lon_diff = (ds['XLONG'] - lon) ** 2
                distance = lat_diff + lon_diff

                # Find minimum distance indices
                min_idx = distance.argmin()
                lat_idx, lon_idx = np.unravel_index(min_idx, distance.shape)

                actual_lat = float(ds['XLAT'].isel(south_north=lat_idx, west_east=lon_idx).values)
                actual_lon = float(ds['XLONG'].isel(south_north=lat_idx, west_east=lon_idx).values)
                print(f"Target coordinates: ({lat:.4f}, {lon:.4f})")
                print(f"Nearest grid point: ({actual_lat:.4f}, {actual_lon:.4f})")

            # Extract SWE data for all times in this year
            swe_values = ds['CaSR_v3.1_P_SWE_LAND'].isel(south_north=lat_idx, west_east=lon_idx).values
            times = ds['time'].values

            # Append data for each timestep
            for time, swe_value in zip(times, swe_values):
                all_data.append({
                    'datetime': pd.Timestamp(time),
                    'swe': float(swe_value)
                })

            ds.close()

        except Exception as e:
            print(f"Error processing {filepath}: {e}")
            continue

    if not all_data:
        return pd.DataFrame()

    # Create dataframe once at the end
    df_final = pd.DataFrame(all_data)
    df_final.set_index('datetime', inplace=True)
    df_final.sort_index(inplace=True)

    return df_final


def find_nearest_grid_point_casr(ds: xr.Dataset, target_lat: float,
                                  target_lon: float) -> tuple:
    """
    Find nearest grid point in CaSR rotated grid.

    Parameters:
    -----------
    ds : xarray.Dataset
        Dataset with 2D lat/lon coordinates
    target_lat : float
        Target latitude
    target_lon : float
        Target longitude

    Returns:
    --------
    tuple
        Indices (rlat, rlon) of nearest grid point
    """
    # Convert target_lon from -180/180 to 0/360 format
    target_lon = target_lon if target_lon >= 0 else target_lon + 360

    # Calculate distance
    distance = np.sqrt((ds.lat.values - target_lat)**2 +
                       (ds.lon.values - target_lon)**2)

    # Find minimum
    min_idx = np.unravel_index(np.argmin(distance), distance.shape)

    return int(min_idx[0]), int(min_idx[1])


def process_all_stations(stations_df: pd.DataFrame, source: str = 'wrf',
                         base_path: str = '.', filepath: str = None) -> dict:
    """
    Process all stations from WRF, ERA5-Land, or CaSR data.

    Parameters:
    -----------
    stations_df : pd.DataFrame
        DataFrame with columns: station_id, lat, lon, initial_year, final_year
    source : str
        Data source: 'wrf', 'era5', or 'casr'
    base_path : str
        Base directory for WRF/CaSR (ignored if source='era5')
    filepath : str
        Path to ERA5-Land netCDF file (required if source='era5')

    Returns:
    --------
    dict
        Dictionary with station_id as keys and DataFrames as values
    """
    sources = {'wrf': extract_snow_timeseries,
               'era5': extract_snow_timeseries_era5,
               'casr': extract_snow_timeseries_casr}

    if source not in sources:
        raise ValueError(f"source must be one of {list(sources.keys())}")
    if source == 'era5' and filepath is None:
        raise ValueError("filepath required for ERA5-Land")

    extract_func = sources[source]
    results = {}

    for idx, row in stations_df.iterrows():
        station_id = row['station_id']
        print(f"Processing {station_id}...")

        kwargs = {
            'lat': row['lat'],
            'lon': row['lon'],
            'initial_year': int(row['initial_year']),
            'final_year': int(row['final_year'])
        }

        if source == 'era5':
            kwargs['filepath'] = filepath
        else:
            kwargs['base_path'] = base_path

        df_snow = extract_func(**kwargs)

        if not df_snow.empty:
            results[station_id] = df_snow
            print(f"  {len(df_snow)} timesteps: {df_snow.index.min()} to {df_snow.index.max()}")
        else:
            print(f"  No data")

    return results
