"""
Time Points Consistency Checker for NetCDF Files

This script analyzes time points across multiple years and months for SNOW NetCDF files.
It checks the number of time points for each month and reports any inconsistencies.

Usage:
    python time_points_check.py /path/to/base/directory start_year end_year

"""
import os
import sys
import netCDF4 as nc
from collections import defaultdict

def count_time_points(file_path):
    """Count the number of time points in a NetCDF file."""
    try:
        with nc.Dataset(file_path, 'r') as ds:
            # Try to find the time dimension
            time_vars = ['time', 'Time', 'XTIME', 'Times']
            for var in time_vars:
                if var in ds.dimensions:
                    return len(ds.dimensions[var])

            # If no standard time dimension found, look for 4D variable
            for var_name in ds.variables:
                var = ds.variables[var_name]
                if len(var.dimensions) == 4:
                    return var.shape[0]

            # If no time dimension found
            return None
    except Exception as e:
        print(f"Error processing {file_path}: {e}")
        return None

def check_time_points(base_dir, start_year, end_year):
    """Check time points for files across specified years for each month."""
    # Dictionary to store time points for each month
    monthly_time_points = defaultdict(list)

    # Walk through directories for specified years
    for year in range(int(start_year), int(end_year) + 1):
        year_dir = str(year)
        year_path = os.path.join(base_dir, year_dir)

        # Ensure the directory exists
        if not os.path.isdir(year_path):
            print(f"Warning: Directory for year {year} not found")
            continue

        # Find all SNOW files in this directory
        for filename in os.listdir(year_path):
            # Skip files that are still compressed
            if filename.startswith('COMPRESSED_SNOW_') or filename.endswith('_daily_mean.nc'):
                continue

            # Look for SNOW files
            if filename.startswith('SNOW_') and filename.endswith('.nc'):
                # Extract month from filename
                try:
                    month = filename.split('_')[-1].split('.')[0]
                except:
                    continue

                # Full file path
                file_path = os.path.join(year_path, filename)

                # Count time points
                time_points = count_time_points(file_path)

                # Store results
                if time_points is not None:
                    monthly_time_points[month].append((year_dir, time_points))

    # Check and report results
    print("Time Points Analysis:")
    for month, time_point_list in sorted(monthly_time_points.items()):
        print(f"\nMonth {month}:")
        # Group by time points
        grouped = defaultdict(list)
        for year, time_points in time_point_list:
            grouped[time_points].append(year)

        # Print results
        for time_points, years in grouped.items():
            print(f"  {time_points} time points in years: {', '.join(years)}")

        # Check if all files for this month have the same time points
        if len(set(tp for _, tp in time_point_list)) > 1:
            print(f"  WARNING: Inconsistent time points for month {month}")

            # Detailed breakdown of inconsistencies
            for time_points, years in grouped.items():
                print(f"    Subset with {time_points} time points: {', '.join(years)}")

def main():
    # Check if correct number of arguments is provided
    if len(sys.argv) != 4:
        print("Usage: python script.py /path/to/base/directory start_year end_year")
        print("Example: python script.py /path/to/data 2005 2010")
        sys.exit(1)

    base_dir = sys.argv[1]
    start_year = sys.argv[2]
    end_year = sys.argv[3]

    # Validate base directory
    if not os.path.isdir(base_dir):
        print(f"Error: {base_dir} is not a valid directory")
        sys.exit(1)

    # Validate years
    try:
        int(start_year)
        int(end_year)
        if int(start_year) > int(end_year):
            raise ValueError("Start year must be less than or equal to end year")
    except ValueError as e:
        print(f"Error: Invalid year input - {e}")
        sys.exit(1)

    # Run the analysis
    check_time_points(base_dir, start_year, end_year)

if __name__ == '__main__':
    main()
