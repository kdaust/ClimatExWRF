#!/usr/bin/env bash
#!/usr/bin/env bash
################################################################################
# casr_to_wrf_regrid.sh
#
# DESCRIPTION:
#   Regridds daily snow water equivalent (SWE) data from the CaSR dataset to
#   match a WRF (Weather Research and Forecasting)
#   model grid for a specified date range. The script performs bilinear
#   remapping and spatial subsetting to extract the region of interest.
#
# USAGE:
#   bash casr_to_wrf_regrid.sh <start_year> <end_year>
#
# ARGUMENTS:
#   start_year      : Starting year (YYYY format) for data processing
#   end_year        : Ending year (YYYY format) for data processing (inclusive)
#
# PREREQUISITES:
#   - CDO (Climate Data Operators) must be installed and in PATH
#   - WRF grid description file at: /casr_dataset/wrf_grid.txt
#   - CaSR SWE data at: CaSR_v3.1/swe/ with naming convention swe_YYYYMMDD12.nc
#
# PROCESSING STEPS:
#   1. Validates input arguments and required files
#   2. Iterates through each day from start_year-01-01 to end_year-12-31
#   3. For each available CaSR file:
#      a. Remaps data to WRF grid using bilinear interpolation
#      b. Subsets to geographic boundaries (lon: -126.9 to -124.8,
#         lat: 43.10 to 66.6)
#      c. Saves final output to /casr_dataset/swe_YYYYMMDD12.nc
#
# OUTPUT:
#   Regridded and subset NetCDF files saved to: /casr_dataset/
#
# EXIT CODES:
#   0 : Successful completion
#   1 : Incorrect number of arguments or missing required files
#
################################################################################
set -euo pipefail

WRF_FILE="/climatex/static_d03.nc"
CASR_PATH="/CaSR_v3.1/swe/"
OUTPUT_DIR="climatex/casr_dataset/"
WRF_GRID="climatex/casr_dataset/wrf_grid.txt"

#Usage: bash casr_to_wrf_regrid.sh <wrf_grid.txt> <casr_path> <output_dir> <start_year> <end_year>
if [ "$#" -ne 2 ]; then
  echo "Usage: $0  <start_year> <end_year>"
  exit 1
fi

echo "Started"

START_YEAR=$1
END_YEAR=$2

if [ ! -f "$WRF_GRID" ]; then
  echo "Error: WRF grid description file not found: $WRF_GRID"
  exit 1
fi
if [ ! -d "$CASR_PATH" ]; then
  echo "Error: Casr data path not found: $CASR_PATH"
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_DIR="${OUTPUT_DIR%/}"

echo "Parsing"
# Parse 4 numeric boundaries from wrf_grid.txt (order: lon_min, lon_max, lat_min, lat_max)
#nums=($(grep -Eo '[-]?[0-9]+\.?[0-9]*' "$WRF_GRID"))
# Use grep and awk to extract specific values
lon_min=-126.9
lon_max=-124.8
lat_min=43.10
lat_max=66.6

start_date=$(date -u -d "${START_YEAR}-01-01" +%Y-%m-%d)
end_date=$(date -u -d "${END_YEAR}-12-31" +%Y-%m-%d)

current_date="$start_date"
echo "Loop!"
while :; do
  if [ "$(date -u -d "$current_date" +%Y%m%d)" -gt "$(date -u -d "$end_date" +%Y%m%d)" ]; then
    break
  fi

  ymd=$(date -u -d "$current_date" +%Y%m%d)
  casr_file="${CASR_PATH}/swe_${ymd}12.nc"
  echo "Processing: ${casr_file}"

  if [ -f "$casr_file" ]; then
    subset_file="${SCRIPT_DIR}/swe_${ymd}_subset.nc"
    remapped_file="${SCRIPT_DIR}/swe_${ymd}_remap.nc"
    final_file="${OUTPUT_DIR}/swe_${ymd}12.nc"

    cdo -L remapbil,"${WRF_GRID}" "$casr_file" "$remapped_file"

    echo "${lon_min}","${lon_max}","${lat_min}","${lat_max}"

    cdo -L sellonlatbox,"${lon_min}","${lon_max}","${lat_min}","${lat_max}" "$remapped_file" "$subset_file"

    mv "$subset_file" "$final_file"
    rm -f "$remapped_file"
  fi

  current_date=$(date -u -d "$current_date + 1 day" +%Y-%m-%d)
done

echo "Done. Regridded data saved to: $OUTPUT_DIR"
