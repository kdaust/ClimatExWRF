# Original script by Kiri Daust
# Modified by Robert Payne

# This script extracts hrly wind components, Q2, PSFC from WRF data
# and calculates the wind speed (WS10)

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
os.makedirs(os.path.join(output_dir, "WS10/unmerged"), exist_ok=True)
os.makedirs(os.path.join(output_dir, "Q2/unmerged"), exist_ok=True)
os.makedirs(os.path.join(output_dir, "T2/unmerged"), exist_ok=True)
os.makedirs(os.path.join(output_dir, "PSFC/unmerged"), exist_ok=True)
os.makedirs(os.path.join(output_dir, "RH/unmerged"), exist_ok=True)

pattern = "**/COMPRESSED_SUBSETTED_d03_metgrid_*.nc"

L = 2.457e6 # J/kg
m_v = 18.01528e-3 # kg/mol
m_d = 28.9647e-3 # kg/mol
emm = m_d/m_v - 1
e_0 = 6.112e2
T_0 = 273.15
R = 8.314462618 # J/(mol·K)
expr_RH = f"RH=100*Q2*PSFC*{m_d}/((1+Q2*{emm})*{m_v}*{e_0}*exp({m_v}*{L}/{R}*(1/{T_0}-1/T2)))"
expr_WS10 = "WS10=sqrt(U10*U10+V10*V10)"
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
    logging.info("Processing %s → start date: %s", filename, start_date)

    # output paths
    outfile_ws = os.path.join(output_dir, "WS10/unmerged/", f"{year}{month}.nc")
    outfile_t2 = os.path.join(output_dir, "T2/unmerged/", f"{year}{month}.nc")
    outfile_q2 = os.path.join(output_dir, "Q2/unmerged/", f"{year}{month}.nc")
    outfile_psfc = os.path.join(output_dir, "PSFC/unmerged/", f"{year}{month}.nc")
    outfile_rh = os.path.join(output_dir, "RH/unmerged/", f"{year}{month}.nc")

    # Set time axis (hourly timestep)
    # and calculate WS10 from components
    # try:
    #     subprocess.run([
    #         "cdo",
    #         f"-expr,{expr_WS10}",
    #         "-delete,timestep=-1",
    #         f"-settaxis,{start_date},00:00:00,1hour",
    #         filepath,
    #         outfile_ws,
    #     ], check=True)
    #     logging.info("Saved: %s", outfile_ws)

    # except Exception as e:
    #     logging.error("Error processing %s for WS: %s", filepath, e)

    # Q2
    try:
        subprocess.run([
            "cdo",
            "-delete,timestep=-1",
            f"-settaxis,{start_date},00:00:00,1hour",
            "-selvar,Q2",
            filepath,
            outfile_q2,
        ], check=True)
        logging.info("Saved: %s", outfile_q2)

    except Exception as e:
        logging.error("Error processing %s for Q2: %s", filepath, e)

    # T2
    try:
        subprocess.run([
            "cdo",
            "-delete,timestep=-1",
            f"-settaxis,{start_date},00:00:00,1hour",
            "-selvar,T2",
            filepath,
            outfile_t2,
        ], check=True)
        logging.info("Saved: %s", outfile_t2)

    except Exception as e:
        logging.error("Error processing %s for T2: %s", filepath, e)

    # PSFC
    try:
        subprocess.run([
            "cdo",
            "-delete,timestep=-1",
            f"-settaxis,{start_date},00:00:00,1hour",
            "-selvar,PSFC",
            filepath,
            outfile_psfc,
        ], check=True)
        logging.info("Saved: %s", outfile_psfc)

    except Exception as e:
        logging.error("Error processing %s for PSFC: %s", filepath, e)

    # RH
    try:
        subprocess.run([
            "cdo",
            f"-expr,{expr_RH}",
            "-delete,timestep=-1",
            f"-settaxis,{start_date},00:00:00,1hour",
            filepath,
            outfile_rh,
        ], check=True)
        logging.info("Saved: %s", outfile_rh)

    except Exception as e:
        logging.error("Error processing %s for RH: %s", filepath, e)


# 🔁 Process all matching files
filepaths = glob.glob(os.path.join(root_dir, pattern), recursive=True)
logging.info("Found %d files matching pattern.", len(filepaths))

for fp in sorted(filepaths):
    process_file(fp)
