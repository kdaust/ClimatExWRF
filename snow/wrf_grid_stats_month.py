"""
Compute monthly climatological statistics for snow water equivalent (SWE) data.

This module provides functions to calculate monthly statistics (mean, max, min, standard
deviation) from multiple snow datasets and save them in a standardized NetCDF format.
It supports three data sources:

1. CaSR - swe_YYYY.nc files
2. WRF (Weather Research and Forecasting) - SNOW_d03_metgrid daily files organized by water year
3. ERA5-Land - Reanalysis snow depth data

Functions:
-----------
compute_casr_monthly(casr_path, initial_year, final_year, months, output_path)
    Processes CaSR SWE data and computes monthly climatological statistics.

compute_monthly_statistics(base_path, initial_year, final_year, months, output_path)
    Memory-efficient processor for WRF SNOW data organized by water year.

compute_era5land_monthly_statistics(era5_path, initial_year, final_year, months, output_path)
    Processes ERA5-Land snow depth data and computes monthly statistics.

Output Format:
--------------
All functions produce standardized NetCDF files containing:
- Monthly climatological statistics (mean, max, min, std)
- Spatial coordinates (XLAT, XLONG for CaSR/ERA5; indices for WRF)
- Month dimension (Times) for temporal organization
- Variable naming convention: SNOW_[stat] (mean, max, min, std)

Notes:
------
- Default month range is Nov-Apr (hydrologic winter season)
- Functions skip missing files with warning messages
- WRF data uses water year convention (Oct-Sep)
"""
import xarray as xr
import numpy as np
from pathlib import Path

def compute_casr_monthly(casr_path, initial_year, final_year,
                                    months=[11, 12, 1, 2, 3, 4],
                                    output_path="casr_monthly_stats.nc"):
    """
    Compute monthly statistics for CaSR SWE data.

    Parameters:
    -----------
    casr_path : str or Path
        Path to folder containing CaSR swe_YYYY.nc files
    initial_year : int
        Starting year
    final_year : int
        Ending year
    months : list
        Months to process (default: Nov-Apr)
    output_path : str
        Output netcdf file path
    """

    casr_path = Path(casr_path)
    results = {}
    spatial_coords = None
    swe_var_name = 'SNOW'  # Use 'SNOW' to match WRF output

    # Loop through months
    for month in months:
        print(f"Processing month {month}...")
        monthly_stats = []

        for year in range(initial_year, final_year + 1):
            filename = f"swe_{year}.nc"
            filepath = casr_path / filename

            if not filepath.exists():
                print(f"  Warning: File not found: {filepath}")
                continue

            try:
                ds = xr.open_dataset(filepath)

                # Get the SWE variable name
                casr_swe_var = 'CaSR_v3.1_P_SWE_LAND' if 'CaSR_v3.1_P_SWE_LAND' in ds.data_vars else list(ds.data_vars)[0]

                # Filter data for the specific month
                ds['month'] = ds['time'].dt.month
                monthly_data = ds.where(ds['month'] == month, drop=True)

                if len(monthly_data['time']) == 0:
                    print(f"  No data for month {month} in year {year}")
                    ds.close()
                    continue

                # Get spatial coordinates from first file
                if spatial_coords is None:
                    spatial_coords = {
                        'XLAT': ds['XLAT'].values,
                        'XLONG': ds['XLONG'].values,
                        'south_north': ds.dims['south_north'],
                        'west_east': ds.dims['west_east'],
                    }

                # Compute stats for this month in this year
                stats = {
                    'mean': monthly_data[casr_swe_var].mean(dim='time').values,
                    'max': monthly_data[casr_swe_var].max(dim='time').values,
                    'min': monthly_data[casr_swe_var].min(dim='time').values,
                    'std': monthly_data[casr_swe_var].std(dim='time').values,
                }
                monthly_stats.append(stats)
                ds.close()

            except Exception as e:
                print(f"  Error reading {filepath}: {e}")
                continue

        # Average statistics across years
        if monthly_stats:
            results[month] = {
                'mean': np.mean([s['mean'] for s in monthly_stats], axis=0),
                'max': np.mean([s['max'] for s in monthly_stats], axis=0),
                'min': np.mean([s['min'] for s in monthly_stats], axis=0),
                'std': np.mean([s['std'] for s in monthly_stats], axis=0),
            }

    # Create output dataset with same format as WRF version
    print("\nCreating output NetCDF file...")

    sorted_months = sorted(results.keys())

    output_ds = xr.Dataset(
        data_vars={
            f'{swe_var_name}_mean': (['Times', 'south_north', 'west_east'],
                                      np.array([results[m]['mean'] for m in sorted_months])),
            f'{swe_var_name}_max': (['Times', 'south_north', 'west_east'],
                                     np.array([results[m]['max'] for m in sorted_months])),
            f'{swe_var_name}_min': (['Times', 'south_north', 'west_east'],
                                     np.array([results[m]['min'] for m in sorted_months])),
            f'{swe_var_name}_std': (['Times', 'south_north', 'west_east'],
                                     np.array([results[m]['std'] for m in sorted_months])),
        },
        coords={
            'Times': (['Times'], sorted_months),
            'south_north': (['south_north'], np.arange(spatial_coords['south_north'])),
            'west_east': (['west_east'], np.arange(spatial_coords['west_east'])),
            'XLAT': (['south_north', 'west_east'], spatial_coords['XLAT']),
            'XLONG': (['south_north', 'west_east'], spatial_coords['XLONG']),
        }
    )

    # Save to file
    output_ds.to_netcdf(output_path)
    print(f"Saved results to {output_path}")

    # Display the result
    print("\nResult dataset:")
    print(output_ds)

    return output_ds


def compute_monthly_statistics(base_path, initial_year, final_year,
                                         months=[11, 12, 1, 2, 3, 4],
                                         output_path="swe_monthly_stats.nc"):
    """Memory-efficient version that computes and saves monthly statistics"""

    base_path = Path(base_path)
    results = {}
    spatial_dims = None
    snow_var_name = None

    for month in months:
        print(f"Processing month {month}...")
        monthly_stats = []

        for water_year in range(initial_year, final_year + 1):
            # Determine which year for this month
            year = water_year - 1 if month in [10, 11, 12] else water_year

            filename = f"SNOW_d03_metgrid_{year}_{month:02d}_daily_mean.nc"
            filepath = base_path / str(water_year) / filename

            if not filepath.exists():
                continue

            try:
                ds = xr.open_dataset(filepath)
                snow_var = 'SNOW' if 'SNOW' in ds.data_vars else list(ds.data_vars)[0]
                snow_var_name = snow_var

                # Get spatial dimensions from first file
                if spatial_dims is None:
                    spatial_dims = {dim: ds.dims[dim] for dim in ds[snow_var].dims if dim != 'Time'}

                # Compute stats for each year
                stats = {
                    'mean': ds[snow_var].mean(dim='Times').values,
                    'max': ds[snow_var].max(dim='Times').values,
                    'min': ds[snow_var].min(dim='Times').values,
                    'std': ds[snow_var].std(dim='Times').values,
                }
                monthly_stats.append(stats)
                ds.close()

            except Exception as e:
                print(f"  Error with {filepath}: {e}")

        # Average across years
        if monthly_stats:
            results[month] = {
                'mean': np.mean([s['mean'] for s in monthly_stats], axis=0),
                'max': np.mean([s['max'] for s in monthly_stats], axis=0),
                'min': np.mean([s['min'] for s in monthly_stats], axis=0),
                'std': np.mean([s['std'] for s in monthly_stats], axis=0),
            }

    # Get a reference dataset for coordinates
    print("\nCreating output NetCDF file...")

    # Get coordinate information from first file
    ref_filepath = None
    for water_year in range(initial_year, final_year + 1):
        year = water_year - 1 if 11 in months else water_year
        for m in months:
            year = water_year - 1 if m in [10, 11, 12] else water_year
            filename = f"SNOW_d03_metgrid_{year}_{m:02d}_daily_mean.nc"
            filepath = base_path / str(water_year) / filename
            if filepath.exists():
                ref_filepath = filepath
                break
        if ref_filepath:
            break

    if ref_filepath:
        ref_ds = xr.open_dataset(ref_filepath)
        south_north = ref_ds.dims.get('south_north', ref_ds.dims.get('y', None))
        west_east = ref_ds.dims.get('west_east', ref_ds.dims.get('x', None))

        # Get coordinate values if they exist
        south_north_coords = ref_ds.coords.get('south_north', np.arange(south_north)).values
        west_east_coords = ref_ds.coords.get('west_east', np.arange(west_east)).values

        ref_ds.close()
    else:
        raise FileNotFoundError("Could not find any reference file to extract coordinates")

    # Stack all statistics for each month
    sorted_months = sorted(results.keys())
    mean_data = np.array([results[m]['mean'] for m in sorted_months])
    max_data = np.array([results[m]['max'] for m in sorted_months])
    min_data = np.array([results[m]['min'] for m in sorted_months])
    std_data = np.array([results[m]['std'] for m in sorted_months])

    # Create dataset with proper structure - use .data to extract numpy arrays
    output_ds = xr.Dataset(
        data_vars={
            f'{snow_var_name}_mean': (['Times', 'south_north', 'west_east'], mean_data),
            f'{snow_var_name}_max': (['Times', 'south_north', 'west_east'], max_data),
            f'{snow_var_name}_min': (['Times', 'south_north', 'west_east'], min_data),
            f'{snow_var_name}_std': (['Times', 'south_north', 'west_east'], std_data),
        },
        coords={
            'Times': (['Times'], sorted_months),
            'south_north': (['south_north'], south_north_coords),
            'west_east': (['west_east'], west_east_coords),
        }
    )

    # Save to file
    output_ds.to_netcdf(output_path)
    print(f"Saved results to {output_path}")

    # Display the result
    print("\nResult dataset:")
    print(output_ds)

    return output_ds


def compute_era5land_monthly_statistics(era5_path, initial_year, final_year,
                                        months=[11, 12, 1, 2, 3, 4],
                                        output_path="era5land_monthly_stats.nc"):
    """
    Compute monthly statistics for ERA5Land snow depth data.

    Parameters:
    -----------
    era5_path : str or Path
        Path to ERA5Land netcdf file (e.g., era5land_data.nc)
    initial_year : int
        Starting year
    final_year : int
        Ending year
    months : list
        Months to process (default: Nov-Apr)
    output_path : str
        Output netcdf file path
    """

    era5_path = Path(era5_path)
    results = {}
    spatial_coords = None
    swe_var_name = 'SNOW'  # Use 'SNOW' to match other outputs

    # Load the entire ERA5Land dataset once
    print(f"Loading ERA5Land data from {era5_path}...")
    try:
        ds = xr.open_dataset(era5_path)
    except Exception as e:
        print(f"Error loading {era5_path}: {e}")
        return None

    # Get the snow depth variable name
    era5_swe_var = 'sd' if 'sd' in ds.data_vars else list(ds.data_vars)[0]

    # Store spatial coordinates
    spatial_coords = {
        'XLAT': ds['XLAT'].values,
        'XLONG': ds['XLONG'].values,
        'south_north': ds.dims['south_north'],
        'west_east': ds.dims['west_east'],
    }

    # Loop through months
    for month in months:
        print(f"Processing month {month}...")
        monthly_stats = []

        # Filter data for the specific month and year range
        ds['month'] = ds['valid_time'].dt.month
        ds['year'] = ds['valid_time'].dt.year

        # Get data for this month across all years in range
        monthly_data = ds.where(
            (ds['month'] == month) &
            (ds['year'] >= initial_year) &
            (ds['year'] <= final_year),
            drop=True
        )

        if len(monthly_data['valid_time']) == 0:
            print(f"  No data available for month {month}")
            continue

        # Group by year and compute statistics
        for year in range(initial_year, final_year + 1):
            year_month_data = monthly_data.where(monthly_data['year'] == year, drop=True)*1000

            if len(year_month_data['valid_time']) == 0:
                continue

            # Compute stats for this month in this year
            stats = {
                'mean': year_month_data[era5_swe_var].mean(dim='valid_time').values,
                'max': year_month_data[era5_swe_var].max(dim='valid_time').values,
                'min': year_month_data[era5_swe_var].min(dim='valid_time').values,
                'std': year_month_data[era5_swe_var].std(dim='valid_time').values,
            }
            monthly_stats.append(stats)

        # Average statistics across years
        if monthly_stats:
            results[month] = {
                'mean': np.mean([s['mean'] for s in monthly_stats], axis=0),
                'max': np.mean([s['max'] for s in monthly_stats], axis=0),
                'min': np.mean([s['min'] for s in monthly_stats], axis=0),
                'std': np.mean([s['std'] for s in monthly_stats], axis=0),
            }
            print(f"  Processed {len(monthly_stats)} years for month {month}")
        else:
            print(f"  No valid data found for month {month}")

    ds.close()

    # Create output dataset with same format as other functions
    print("\nCreating output NetCDF file...")

    sorted_months = sorted(results.keys())

    if not sorted_months:
        print("Error: No data was processed!")
        return None

    output_ds = xr.Dataset(
        data_vars={
            f'{swe_var_name}_mean': (['Times', 'south_north', 'west_east'],
                                      np.array([results[m]['mean'] for m in sorted_months])),
            f'{swe_var_name}_max': (['Times', 'south_north', 'west_east'],
                                     np.array([results[m]['max'] for m in sorted_months])),
            f'{swe_var_name}_min': (['Times', 'south_north', 'west_east'],
                                     np.array([results[m]['min'] for m in sorted_months])),
            f'{swe_var_name}_std': (['Times', 'south_north', 'west_east'],
                                     np.array([results[m]['std'] for m in sorted_months])),
        },
        coords={
            'Times': (['Times'], sorted_months),
            'south_north': (['south_north'], np.arange(spatial_coords['south_north'])),
            'west_east': (['west_east'], np.arange(spatial_coords['west_east'])),
            'XLAT': (['south_north', 'west_east'], spatial_coords['XLAT']),
            'XLONG': (['south_north', 'west_east'], spatial_coords['XLONG']),
        }
    )

    # Save to file
    output_ds.to_netcdf(output_path)
    print(f"Saved results to {output_path}")

    # Display the result
    print("\nResult dataset:")
    print(output_ds)

    return output_ds


# Usage example
if __name__ == "__main__":
    # Your parameters
    # base_path = "/path/to/wrf/data"
    initial_year = 2000
    final_year = 2023
    months = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]

    # # Compute statistics
    # result = compute_monthly_statistics(
    #     base_path='/climatex/',
    #     initial_year=initial_year,
    #     final_year=final_year,
    #     months=months,
    #     output_path='/climatex/wrf_data/snow_stats_month.nc'
    # )

    result = compute_casr_monthly(
        casr_path='/climatex/casr_dataset/year/',
        initial_year=initial_year,
        final_year=final_year,
        months=months,
        output_path='climatex/casr_dataset/snow_stats_month.nc'
    )


    # result = compute_era5land_monthly_statistics(
    #     era5_path='climatex/era5_dataset/swe_daily_1999-2023.nc',
    #     initial_year=initial_year,
    #     final_year=final_year,
    #     months=months,
    #     output_path='climatex/era5_dataset/snow_stats_month.nc'
    # )

    # Display the result
    print("\nResult dataset:")
    print(result)
