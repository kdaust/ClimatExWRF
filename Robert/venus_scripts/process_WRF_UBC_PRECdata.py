# Original script by Kiri Daust
# Modified by Robert Payne

# This script extracts precIP from UBC WRF data

import os
import re
import glob
import subprocess
import sys
import logging

logging.basicConfig(level=logging.INFO)
logging.info("Processing Script started")
logging.info("Python executable: %s", sys.executable)

# ----- Configuration ----- #
root_dir = "/net/venus/blue/data/WRF/kikridaust/Share_Data/SUBSETTED/"
output_dir = "/home/rpayne/data-rpayne/unproc/WRF-UBC/"
os.makedirs(os.path.join(output_dir, "PREC/unmerged"), exist_ok=True)

pattern = "**/COMPRESSED_RAIN_d03_metgrid_*.nc"
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
    logging.info("Processing %s → start date: %s", filename, start_date)

    # output paths
    outfile_prec = os.path.join(output_dir, "PREC/unmerged/", f"{year}{month}.nc")

    # PREC
    try:
        subprocess.run([
            "cdo",
            f"-settaxis,{start_date},01:00:00,1hour",
            "-selvar,HPRECIPNC",
            filepath,
            outfile_prec,
        ], check=True)
        logging.info("Saved: %s", outfile_prec)

    except Exception as e:
        logging.error("Error processing %s for PREC: %s", filepath, e)


# 🔁 Process all matching files
filepaths = glob.glob(os.path.join(root_dir, pattern), recursive=True)
logging.info("Found %d files matching pattern.", len(filepaths))

for fp in sorted(filepaths):
    process_file(fp)
