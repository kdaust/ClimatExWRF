### Original script by Kiri Daust
### modified by Robert Payne

import os
import re
import glob
import subprocess

# ----- Configuration ----- #
root_dir = "/blue/data/WRF/kikridaust/Share_Data/SUBSETTED/"
output_dir = "/users/rpayne/data/unproc/WRF-UBC/"
os.makedirs(output_dir, exist_ok=True)

pattern = "**/COMPRESSED_SUBSETTED_d03_metgrid_*.nc"
expression = "WS10=sqrt(U10*U10 + V10*V10)"
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

    # output paths
    outfile_ws10 = os.path.join(output_dir, "WS10", f"ws10_{year}_{month}.nc")
    outfile_q2 = os.path.join(output_dir, "Q2", f"q2_{year}_{month}.nc")

    # Set time axis (hourly timestep)
    # and calculate WS10 from components
    try:
        subprocess.run([
            "cdo",
            f"-settaxis,{start_date},00:00:00,1hour",
            f"expr,{expression}",
            filepath,
            outfile_ws10
        ], check=True)
        print(f"Saved: {outfile_ws10}")

    except subprocess.CalledProcessError as e:
        print(f"Error processing {filepath} for WS10: {e}")

    try:
        subprocess.run([
            "cdo",
            f"-settaxis,{start_date},00:00:00,1hour",
            "Q2",
            filepath,
            outfile_q2
        ], check=True)
        print(f"Saved: {outfile_q2}")

    except subprocess.CalledProcessError as e:
        print(f"Error processing {filepath} for Q2: {e}")


# 🔁 Process all matching files
filepaths = glob.glob(os.path.join(root_dir, pattern), recursive=True)
for fp in sorted(filepaths):
    process_file(fp)
