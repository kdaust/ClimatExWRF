"""
Parallel Daily Mean NetCDF Processor using Dask

This script processes NetCDF files to generate daily mean values across
multiple years using distributed Dask processing.

Usage:
    python script.py /path/to/base/directory start_year end_year [num_workers]
"""

import os
import sys
import xarray as xr
import pandas as pd
from datetime import datetime
import dask
from dask.distributed import Client, LocalCluster

def process_daily_mean(input_file, output_file):
    """
    Process a single NetCDF file to generate daily mean values using Dask.

    Args:
        input_file (str): Path to the input NetCDF file
        output_file (str): Path to save the daily mean NetCDF file

    Returns:
        bool: True if processing successful, False otherwise
    """
    try:
        # Read the NetCDF file with Dask
        ds = xr.open_dataset(input_file, chunks={'Times': 'auto'})

        # Parse the start date
        start_date_str = ds.attrs['START_DATE']
        start_date = datetime.strptime(start_date_str, '%Y-%m-%d_%H:%M:%S')

        # Create a time coordinate based on hourly frequency and start date
        time_coords = pd.date_range(start=start_date, periods=len(ds.Times), freq='H')

        # Add time coordinate to the dataset
        ds = ds.assign_coords(Times=time_coords)

        # Calculate daily mean for each variable
        ds_daily = ds.resample(Times='D').mean()

        # Compute and save the daily mean dataset
        ds_daily.compute()
        ds_daily.to_netcdf(output_file)

        print(f"Daily mean file created: {output_file}")

        # Close the datasets
        ds.close()
        ds_daily.close()

        return True
    except Exception as e:
        print(f"Error processing {input_file}: {e}")
        return False

def process_yearly_data(base_path, start_year, end_year):
    """
    Process NetCDF files across multiple years in a base directory using Dask.

    Args:
        base_path (str): Base directory containing year folders
        start_year (int): Starting year to process
        end_year (int): Ending year to process
    """
    # Convert years to integers
    start_year = int(start_year)
    end_year = int(end_year)

    # Process each year
    for year in range(start_year, end_year + 1):
        year_path = os.path.join(base_path, str(year))

        # Check if year directory exists
        if not os.path.isdir(year_path):
            print(f"Directory not found: {year_path}")
            continue

        # Process files in the year directory
        for filename in os.listdir(year_path):
            # Skip compressed files
            if filename.startswith('COMPRESSED_'):
                continue

            # Process only NetCDF files starting with SNOW
            if filename.startswith('SNOW_') and filename.endswith('.nc'):
                # Full input file path
                input_file = os.path.join(year_path, filename)

                # Create output filename
                output_filename = filename.replace('.nc', '_daily_mean.nc')
                output_file = os.path.join(year_path, output_filename)

                # Process the file
                print(f"Processing: {input_file}")
                process_daily_mean(input_file, output_file)

def main():
    """
    Main function to handle command-line arguments and initiate processing.
    Sets up Dask client and initiates distributed processing.
    """
    # Check if correct number of arguments is provided
    if len(sys.argv) < 4 or len(sys.argv) > 5:
        print("Usage: python script.py /path/to/base/directory start_year end_year [num_workers]")
        print("Example: python script.py /path/to/snow 1999 2010 4")
        sys.exit(1)

    # Parse arguments
    base_path = sys.argv[1]
    start_year = sys.argv[2]
    end_year = sys.argv[3]

    # Set up Dask client
    if len(sys.argv) == 5:
        num_workers = int(sys.argv[4])
        cluster = LocalCluster(n_workers=num_workers)
        client = Client(cluster)
        print(f"Dask Client created with {num_workers} workers")
    else:
        client = Client()
        print("Dask Client created with default settings")

    # Validate base path
    if not os.path.isdir(base_path):
        print(f"Error: {base_path} is not a valid directory")
        sys.exit(1)

    try:
        # Process the data
        process_yearly_data(base_path, start_year, end_year)
    finally:
        # Ensure Dask client is closed
        client.close()

if __name__ == '__main__':
    main()
