import os
import gc
import pandas as pd
import xarray as xr
import logging
from glob import glob


def readable_time_coord(ds):

    # Determine the time coordinate name
    time_coord = "Times" if "Times" in ds.coords else "XTIME" if "XTIME" in ds.coords else "time" if "time" in ds.coords else None
    if time_coord is None:
        raise ValueError("No time coordinate named 'Times' or 'XTIME' or 'time' found in dataset.")

    time_vals = ds[time_coord].values

    # Extract integer (YYYYMMDD) and fractional day parts
    dates_int = time_vals.astype(int)
    frac_days = time_vals - dates_int

    # Parse integers to datetime
    base_dates = pd.to_datetime(dates_int.astype(str), format="%Y%m%d")

    # Add fractional day as timedelta
    full_times = base_dates + pd.to_timedelta(frac_days, unit="D")
    full_times = full_times.round("1h")

    # Assign back to dataset
    ds = ds.assign_coords({time_coord: full_times})
    ds = ds.transpose(time_coord, ...)

    return ds

logging.basicConfig(level=logging.INFO)
# vars = ['WS10', 'T2', 'Q2', 'PSFC', 'PREC', 'RH']
vars = ['T2', 'Q2', 'PSFC', 'RH']
for var in vars:
    logging.info(f"Processing variable: {var}")
    root = f"/home/rpayne/data-rpayne/unproc/WRF-UBC/{var}/unmerged/"
    filepaths = sorted(glob(os.path.join(root, "*.nc")))
    for filepath in filepaths:
        try:
            ds = readable_time_coord(xr.open_dataset(filepath))
            ds.to_netcdf(filepath.replace(".nc", "_r.nc"))
            os.remove(filepath)
            os.rename(filepath.replace(".nc", "_r.nc"), filepath)
        except Exception as e:
            logging.error(f"Error processing {filepath}: {e}")
    gc.collect()


logging.info("Done.")
