import os
import re
import glob
import subprocess
from pathlib import Path

# ----- Configuration ----- #
root_dir = "/blue/data/WRF/kikridaust/Share_Data/SUBSETTED/" 
output_dir = os.path.join(root_dir, "processed_daily_precip")
os.makedirs(output_dir, exist_ok=True)

pattern = "**/COMPRESSED_RAIN_d03_metgrid_*.nc"
varname = "HPRECIPNC"
# -------------------------- #

file_re = re.compile(r"COMPRESSED_RAIN_d03_metgrid_(\d{4})_(\d{2})\.nc$")

def process_file(filepath):
    filename = os.path.basename(filepath)
    match = file_re.search(filename)
    if not match:
        print(f"Skipping {filename}: couldn't parse date.")
        return

    year, month = match.groups()


    start_date = f"{year}-{month}-01"
    print(f"Processing {filename} → start date: {start_date}")

    # Temp and output paths
    outfile = os.path.join(output_dir, f"daily_precip_{year}_{month}.nc")

    # Set time axis (hourly timestep)
    # Chain: settaxis → selname → daysum
    subprocess.run([
        "cdo",
        "daysum",
        f"-selname,{varname}",
        f"-settaxis,{start_date},00:01:00,1hour",
        filepath,
        outfile
    ], check=True)

    print(f"Saved: {outfile}")

# Process all matching files
filepaths = glob.glob(os.path.join(root_dir, pattern), recursive=True)
for fp in sorted(filepaths):
    process_file(fp)

