##calculate wind speed and statistics

import os
import re
import glob
import subprocess
from pathlib import Path

# ----- Configuration ----- #
root_dir = "/blue/data/WRF/kikridaust/Share_Data/SUBSETTED/"
output_dir = os.path.join(root_dir, "processed_daily_wind")
os.makedirs(output_dir, exist_ok=True)

pattern = "**/COMPRESSED_SUBSETTED_d03_metgrid_*.nc"
varname = "ws"
# -------------------------- #

file_re = re.compile(r"COMPRESSED_SUBSETTED_d03_metgrid_(\d{4})_(\d{2})\.nc$")

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
    tempfile_ws = os.path.join(output_dir, f"hourly_ws_{year}_{month}.nc")
    outfile_wsmax = os.path.join(output_dir, f"daily_wsmax_{year}_{month}.nc")
    outfile_wsmean = os.path.join(output_dir, f"daily_wsmean_{year}_{month}.nc")
    outfile_wsmin = os.path.join(output_dir, f"daily_wsmin_{year}_{month}.nc")
    outfile_wssd = os.path.join(output_dir, f"daily_wssd_{year}_{month}.nc")

    #calculate wind speed
    subprocess.run([
        "cdo",
        f"-selname,{varname}",
        f"-expr,ws=sqrt(U10*U10+V10*V10)",
        f"-settaxis,{start_date},00:00:00,1hour",
        filepath,
        tempfile_ws
    ], check=True)


    subprocess.run([
        "cdo",
        "daymax",
        tempfile_ws,
        outfile_wsmax
    ], check=True)

    subprocess.run([
        "cdo",
        "daymean",
        tempfile_ws,
        outfile_wsmean
    ], check=True)

    subprocess.run([
        "cdo",
        "daymin",
        tempfile_ws,
        outfile_wsmin
    ], check=True)

    subprocess.run([
        "cdo",
        "daystd",
        tempfile_ws,
        outfile_wssd
    ], check=True)

    print(f"Saved: {outfile_wssd} and others")

# Process all matching files
filepaths = glob.glob(os.path.join(root_dir, pattern), recursive=True)
for fp in sorted(filepaths):
    process_file(fp)
