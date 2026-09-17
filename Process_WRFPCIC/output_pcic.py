"""Process one-record hourly WRF files, with independent adjacent-pair workers."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timedelta
import multiprocessing
import os
from pathlib import Path
import re
import tempfile

# Set before importing numerical libraries, including in spawned workers.
# Explicit environment settings still take precedence.
for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_name, "1")

from netCDF4 import Dataset
import numpy as np
import wrf
from compute_ipw_ivt import compute_ipw_ivt

TIME_FORMAT = "%Y-%m-%d_%H:%M:%S"
BUCKET_J = 1.0e9


def discover_files(input_dir, domain):
    """Accept native names and _compressed suffixes, never duplicate an hour."""
    pattern = re.compile(
        rf"wrfout_{domain}_(\d{{4}}-\d{{2}}-\d{{2}}_\d{{2}}:\d{{2}}:\d{{2}})(?:_compressed)?$"
    )
    by_time = {}
    for path in Path(input_dir).glob(f"wrfout_{domain}_*"):
        match = pattern.fullmatch(path.name)
        if not match or not path.is_file():
            continue
        timestamp = datetime.strptime(match[1], TIME_FORMAT)
        if timestamp in by_time:
            raise ValueError(f"Duplicate input hour: {by_time[timestamp]} and {path}")
        by_time[timestamp] = path
    ordered = sorted(by_time)
    if len(ordered) < 2:
        raise ValueError("At least two hourly WRF files are required")
    for previous, current in zip(ordered, ordered[1:]):
        if current - previous != timedelta(hours=1):
            raise ValueError(f"Inputs are not one hour apart: {previous} -> {current}")
    return [by_time[t] for t in ordered]


def read_timestamp(dataset):
    if len(dataset.dimensions["Time"]) != 1:
        raise ValueError("Each input must contain exactly one Time record")
    raw = dataset["Times"][0]
    value = b"".join(raw.tolist()).decode("ascii").strip("\x00 ")
    return datetime.strptime(value, TIME_FORMAT)


def copy_time_coordinates(source, target):
    """Keep timestamps available to cdo mergetime and downstream consumers."""
    for name in ("Times", "XTIME"):
        if name not in source.variables:
            continue
        original = source[name]
        for dim in original.dimensions:
            if dim not in target.dimensions:
                target.createDimension(dim, len(source.dimensions[dim]))
        attrs = {key: original.getncattr(key) for key in original.ncattrs()}
        kwargs = {}
        if "_FillValue" in attrs:
            kwargs["fill_value"] = attrs.pop("_FillValue")
        copied = target.createVariable(name, original.datatype, original.dimensions, **kwargs)
        copied.setncatts(attrs)
        copied[:] = original[:]


def write_to_ncfile(ncfile, in_array, var_string, var_dims, var_units,
                    var_coordinates, var_description):
    newvar = ncfile.createVariable(var_string, np.float32, var_dims)
    newvar.units = var_units
    newvar.coordinates = var_coordinates
    newvar.description = var_description
    values = np.asanyarray(in_array[:])
    if values.ndim == len(var_dims):
        values = values[0]
    newvar[0, ...] = values


def uncompressed_name(name):
    """Remove the optional suffix using operations supported by Python 3.8."""
    suffix = "_compressed"
    return name[:-len(suffix)] if name.endswith(suffix) else name


def process_pair(previous_file, current_file, output_dir):
    """Only filenames cross process boundaries; each worker owns its handles."""
    current_file = Path(current_file)
    stem = uncompressed_name(current_file.name)
    output_path = Path(output_dir) / ("wrfpcic_" + stem[len("wrfout_"):])
    # Hidden temporary names cannot enter the driver's wrfpcic* merge glob.
    fd, temporary = tempfile.mkstemp(prefix=".pcic-", suffix=".tmp", dir=output_dir)
    os.close(fd)
    try:
        with Dataset(previous_file, "r") as pf, Dataset(current_file, "r") as cf:
            previous_time, current_time = read_timestamp(pf), read_timestamp(cf)
            if current_time - previous_time != timedelta(hours=1):
                raise ValueError(f"Input Times are not one hour apart: {previous_file}, {current_file}")
            for path, timestamp in ((Path(previous_file), previous_time), (current_file, current_time)):
                expected = uncompressed_name(path.name).split("_", 2)[2]
                if timestamp.strftime(TIME_FORMAT) != expected:
                    raise ValueError(f"Filename and Times disagree: {path}")
            with Dataset(temporary, "w", format="NETCDF4_CLASSIC") as ncfile:
                _compute_pair(cf, pf, ncfile)
        os.replace(temporary, output_path)
        return str(output_path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def _compute_pair(cf, pf, ncfile):
    # Spatial dimensions
    south_north = ncfile.createDimension("south_north", cf["ACSWUPB"].shape[1])
    west_east = ncfile.createDimension("west_east", cf["ACSWUPB"].shape[2])
    soil_layers_stag = ncfile.createDimension("soil_layers_stag", cf["TSLB"].shape[1])
    time = ncfile.createDimension("Time", None)

    copy_time_coordinates(cf, ncfile)

    lat = ncfile.createVariable('XLAT', np.float32, ("Time","south_north","west_east",))
    lat.units = 'degree_north'
    lat.long_name = 'latitude'
    lat.coordinates = 'XLONG XLAT'
    lon = ncfile.createVariable('XLONG', np.float32, ("Time","south_north","west_east",))
    lon.units = 'degree_east'
    lon.long_name = 'longitude'
    lon.coordinates = 'XLONG XLAT' 
   
    lat[:] = np.array(cf["XLAT"])
    lon[:] = np.array(cf["XLONG"])
    
    ################### HOURLY RADIATION AND SNOWMELT ######################
    #ACSWUPB,ACSWUPBC,ACSWDNB,ACSWDNBC,ACLWUPB,ACLWUPBC,ACLWDNB,ACLWDNBC

    radvars = ["ACSWUPB","ACSWUPBC","ACSWDNB","ACSWDNBC","ACLWUPB","ACLWUPBC","ACLWDNB","ACLWDNBC"]
    dims = ("Time","south_north","west_east",)
    for radvar in radvars:
        cfrad = cf[radvar]
        hourly_radvar = (
            np.array(cf["I_"+radvar])*BUCKET_J + np.array(cfrad)
            - (np.array(pf["I_"+radvar])*BUCKET_J + np.array(pf[radvar]))
        )
        
        write_to_ncfile(ncfile,hourly_radvar,"h"+radvar.lower(),dims,
                        cfrad.units,cfrad.coordinates,"HOURLY "+cfrad.description)

    
    othervars = ["ACGRDFLX","ACSNOM"]
    for othervar in othervars:
        cfother = cf[othervar]
        hourly_othervar = np.array(cfother) - np.array(pf[othervar])

        write_to_ncfile(ncfile,hourly_othervar,"h"+othervar.lower(),dims,
                        cfother.units,cfother.coordinates,"HOURLY "+cfother.description)

    #######################################
    
    wrfout_variables = ["SFROFF","UDROFF","CANWAT","SR","TSLB","SMOIS","SH2O"]
    
    for wrfout_var_string in wrfout_variables:
        wrfout_var_in = cf[wrfout_var_string]
        if wrfout_var_in.ndim == 3:
            dims = ("Time","south_north","west_east",)
        else:
            dims = ("Time","soil_layers_stag","south_north","west_east",)
        
        write_to_ncfile(ncfile,wrfout_var_in,wrfout_var_string,dims,
            wrfout_var_in.units,wrfout_var_in.coordinates,wrfout_var_in.description)
    
    
    getvar_variables = ["ua","va","wa","temp","height","height_agl","pressure","rh","QVAPOR","slp"]

    winds = {name: wrf.getvar(cf, name) for name in ("ua", "va")}
    for wrfout_var_string in getvar_variables:
        if wrfout_var_string == "slp":
            units = "Pa"
            wrfout_var_in = wrf.getvar(cf,wrfout_var_string,units="Pa")[:,:] 
        elif wrfout_var_string in winds:
            wrfout_var_in = winds[wrfout_var_string][0,:,:]
        else:
            wrfout_var_in = wrf.getvar(cf,wrfout_var_string)[0,:,:] # take lowest level

        
        dims = ("Time","south_north","west_east",)
        
       
        if wrfout_var_string == "slp":
            write_to_ncfile(ncfile,wrfout_var_in,wrfout_var_string,dims,
                wrfout_var_in.units,wrfout_var_in.coordinates,wrfout_var_in.description)

        else:             
            write_to_ncfile(ncfile,wrfout_var_in,wrfout_var_string+"_b",dims,
                wrfout_var_in.units,wrfout_var_in.coordinates,"lowest level "+wrfout_var_in.description)

    

    # Extract required 3D fields (taking first time step)
    qv = cf.variables["QVAPOR"][0, :, :, :]  # Water vapor mixing ratio (kg/kg)
    p = cf.variables["P"][0, :, :, :]  # Perturbation pressure (Pa)
    pb = cf.variables["PB"][0, :, :, :]  # Base state pressure (Pa)
    psfc = cf.variables["PSFC"][0, :, :]  # Surface pressure (Pa)

    ua = np.asarray(winds["ua"])[:,:,:]
    va = np.asarray(winds["va"])[:,:,:]

    ipw_50kPa,ivt_50kPa,ivtx_50kPa,ivty_50kPa = compute_ipw_ivt(qv, p, pb, psfc, ua, va)
    ipw_total = wrf.getvar(cf,"pw")[:,:]
 
    write_to_ncfile(ncfile,ipw_50kPa,"IPW_50kPA",dims,
        "mm","XLONG XLAT XTIME","Integrated precipitable water to 50 kPa")
    write_to_ncfile(ncfile,ipw_total,"IPW_total",dims,
        "mm","XLONG XLAT XTIME","Integrated precipitable water to model top")
    write_to_ncfile(ncfile,ivt_50kPa,"IVT_50kPa",dims,
        "kg m-1 s-1","XLONG XLAT XTIME","Integrated vapour transport to 50 kPa")
    write_to_ncfile(ncfile,ivtx_50kPa,"IVTX_50kPa",dims,
        "kg m-1 s-1","XLONG XLAT XTIME","Integrated x-direction vapour transport to 50 kPa")
    write_to_ncfile(ncfile,ivty_50kPa,"IVTY_50kPa",dims,
        "kg m-1 s-1","XLONG XLAT XTIME","Integrated y-direction vapour transport to 50 kPa")

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("domain", choices=("d01", "d02", "d03", "d04"))
    parser.add_argument("--workers", type=int, default=1,
                        help="Number of independent processes (default: 1)")
    parser.add_argument("--input-dir", type=Path, default=Path("."))
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    args = parser.parse_args(argv)
    if args.workers < 1:
        parser.error("--workers must be at least 1")
    try:
        files = discover_files(args.input_dir, args.domain)
    except ValueError as error:
        parser.error(str(error))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    jobs = [(previous, current, args.output_dir)
            for previous, current in zip(files, files[1:])]
    print(f"Processing {len(jobs)} hours with {args.workers} worker(s)", flush=True)
    if args.workers == 1:
        for job in jobs:
            print(process_pair(*job), flush=True)
    else:
        # Spawn avoids inheriting NetCDF/HDF5 state from the parent process.
        with ProcessPoolExecutor(max_workers=args.workers,
                                 mp_context=multiprocessing.get_context("spawn")) as pool:
            futures = [pool.submit(process_pair, *job) for job in jobs]
            try:
                for future in as_completed(futures):
                    print(future.result(), flush=True)
            except BaseException:
                for future in futures:
                    future.cancel()
                raise


if __name__ == "__main__":
    main()
