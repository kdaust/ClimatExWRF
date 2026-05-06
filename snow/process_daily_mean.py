"""
Daily Mean NetCDF Processor

This script processes NetCDF files to generate daily mean values across
multiple years. It handles files in a specified directory structure,
skipping compressed files and creating new daily mean files.

Usage:
    python script.py /path/to/base/directory start_year end_year
"""

import os
import sys
import xarray as xr
import pandas as pd
from datetime import datetime

def process_daily_mean(input_file, output_file):
    """
    Process a single NetCDF file to generate daily mean values.

    Args:
        input_file (str): Path to the input NetCDF file
        output_file (str): Path to save the daily mean NetCDF file

    Returns:
        bool: True if processing successful, False otherwise
    """
    try:
        # Read the NetCDF file
        ds = xr.open_dataset(input_file)

        # Parse the start date
        start_date_str = ds.attrs['START_DATE']
        start_date = datetime.strptime(start_date_str, '%Y-%m-%d_%H:%M:%S')

        # Create a time coordinate based on hourly frequency and start date
        time_coords = pd.date_range(start=start_date, periods=len(ds.Times), freq='H')

        # Add time coordinate to the dataset
        ds = ds.assign_coords(Times=time_coords)

        # Calculate daily mean for each variable
        ds_daily = ds.resample(Times='D').mean()

        # Save the daily mean dataset
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
    Process NetCDF files across multiple years in a base directory.

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
    Validates input and calls the processing function.
    """
    # Check if correct number of arguments is provided
    if len(sys.argv) != 4:
        print("Usage: python script.py /path/to/base/directory start_year end_year")
        print("Example: python script.py /path/to/snow 1999 2010")
        sys.exit(1)

    # Parse arguments
    base_path = sys.argv[1]
    start_year = sys.argv[2]
    end_year = sys.argv[3]

    # Validate base path
    if not os.path.isdir(base_path):
        print(f"Error: {base_path} is not a valid directory")
        sys.exit(1)

    # Process the data
    process_yearly_data(base_path, start_year, end_year)

if __name__ == '__main__':
    main()
