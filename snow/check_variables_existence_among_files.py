"""
Script to analyze and report variable consistency across NetCDF SNOW files.

This script examines NetCDF files in a directory structure organized by year,
extracts variables from SNOW_*.nc files, and identifies inconsistencies in
variable presence across different years for the same month.

The script processes files organized as:
    base_dir/
    ├── year1/
    │   ├── SNOW_*.nc
    │   └── ...
    ├── year2/
    │   ├── SNOW_*.nc
    │   └── ...
    └── ...

For each month, it reports which variables are present/missing across years,
helping identify data quality issues or structural inconsistencies.

Usage:
    python script.py /path/to/base/directory start_year end_year

"""
import os
import sys
import netCDF4 as nc
from collections import defaultdict

def get_variables(file_path):
    """
    Extract variables from a NetCDF file.

    Args:
        file_path (str): Path to the NetCDF file

    Returns:
        set: Set of variable names in the file
    """
    try:
        with nc.Dataset(file_path, 'r') as ds:
            # Get all variable names, excluding dimension variables
            variables = set(ds.variables.keys()) - set(ds.dimensions.keys())
            return variables
    except Exception as e:
        print(f"Error processing {file_path}: {e}")
        return set()

def check_variables(base_dir, start_year, end_year):
    """
    Check variable presence across files for specified years and months.

    Args:
        base_dir (str): Base directory containing year folders
        start_year (int): Starting year
        end_year (int): Ending year
    """
    # Dictionary to store variables for each month
    monthly_variables = defaultdict(list)

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
            if filename.startswith('COMPRESSED_SNOW_'):
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

                # Get variables
                variables = get_variables(file_path)

                # Store results
                if variables:
                    monthly_variables[month].append((year_dir, variables))

    # Check and report results
    print("Variable Presence Analysis:")
    for month, var_list in sorted(monthly_variables.items()):
        print(f"\nMonth {month}:")

        # Collect all unique variables across years
        all_variables = set().union(*[vars_ for _, vars_ in var_list])

        # Track variables missing in some years
        missing_vars = defaultdict(list)

        # Check each variable's presence across years
        for var in sorted(all_variables):
            years_with_var = [year for year, vars_ in var_list if var in vars_]
            years_without_var = [year for year, vars_ in var_list if var not in vars_]

            if years_without_var:
                print(f"  WARNING: Variable '{var}':")
                print(f"    Present in years: {', '.join(years_with_var)}")
                print(f"    Missing in years: {', '.join(years_without_var)}")

                # Track missing variables
                missing_vars[var] = years_without_var

        # Summary of missing variables
        if missing_vars:
            print("\n  Summary of Missing Variables:")
            for var, years in missing_vars.items():
                print(f"    {var}: Missing in {len(years)} year(s)")
        else:
            print("  All variables consistent across years.")

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
    check_variables(base_dir, start_year, end_year)

if __name__ == '__main__':
    main()
