# WRF PCIC processing

Requires Python 3.8+, NumPy, netCDF4, and WRF-Python versions compatible with your
Python installation. The shell driver also
requires the existing rclone configuration, CDO, and nccopy on the processing host.

## Run

```bash
# Full download/process/merge workflow; start with 2 workers if RAM is limited.
PCIC_WORKERS=4 bash process_pcic.sh

# Store all data on a different drive; the scripts can live anywhere.
PCIC_WORKERS=4 bash /path/to/scripts/process_pcic.sh "/mnt/data/WRF/2021-11"

# Process files already downloaded, without the rest of the driver.
python output_pcic.py d03 --workers 4 --input-dir . --output-dir ./pcic-output

# Serial baseline with the same calculations.
python output_pcic.py d03 --workers 1 --output-dir ./pcic-serial
```

The driver accepts an optional data-directory argument, creates it if necessary,
and stores downloads, hourly/intermediate/final outputs, and success markers there.
Relative paths are resolved from the directory where you invoke the command;
omitting the argument keeps the original current-directory behaviour. The Python
script is located alongside the Bash driver, independently of the data directory.

The Python command defaults to one worker; the driver defaults to four.
Each worker opens the previous and current input and produces only the current
hour's output. Inputs may be native `wrfout_d03_YYYY-MM-DD_HH:MM:SS` files or
have an `_compressed` suffix. NetCDF decompresses data as it is read. Keep
exactly one input per timestamp: duplicates, missing hours, filename/Times
disagreements, and files containing multiple Time records cause an error.
Include the preceding hour when processing a month or other time range.

Workers write separate temporary files and publish each output after closing it.
A failed worker makes the command fail; the driver stops before merging or
creating its success marker. Completed outputs may remain after a failure.
Every driver run uses a fresh `pcic-hours.*` directory, retained for inspection,
so older hourly outputs cannot enter its merge. Direct Python runs overwrite
matching outputs; use a fresh output directory when changing the input range.

The existing `WRFOUT.OK` marker still skips processing. Move it aside before
regenerating outputs with this version. The hourly radiation, ground heat flux,
and snowmelt differences now correctly subtract the previous hour; these results
intentionally differ from the original script's cumulative outputs.

## Performance and validation

The horizontal interpolation is vectorised, while the original interpolation
weights and vertical integration method are preserved. Full `ua` and `va`
diagnostics are computed once per hour and reused for IVT. This optimisation
does not independently validate the original scientific integration method.

Start with 2–4 workers and compare elapsed time and peak memory with one worker.
Each worker holds several full 3D fields, so more workers can exhaust RAM or
saturate storage. OpenMP/OpenBLAS/MKL thread counts default to one unless already
set in the environment; avoid combining many processes with many library threads.

```bash
python -m unittest discover -s tests -v
```

Tests compare all four integrated outputs to a frozen copy of the original
algorithm using float32/float64 and pressure-boundary cases. They also check
hourly subtraction (including radiation bucket rollover), wind reuse, input
pairing, and failure cleanup using mocked NetCDF/WRF interfaces. On the WRF host,
compare one-worker and multi-worker runs on the same small real-file subset,
including a month boundary, before running a full month. Check timestamps and
the final CDO merge there as well.

## Pressure-level output

Each hourly file now includes `U_850`, `V_850`, `T_850`, `Q_850`, `Z_850`,
and the corresponding variables at 700, 500 and 250 hPa (20 new variables).
The existing CDO merges include these automatically in the final monthly file.

- U/V: destaggered, grid-relative wind components in m s-1, matching the
  existing `ua_b`/`va_b` convention (not rotated to east/north).
- T: actual air temperature in K, not WRF perturbation potential temperature.
- Q: QVAPOR water-vapour mixing ratio in kg kg-1, not specific humidity.
- Z: geopotential height above mean sea level in decametres (`dm`).

Interpolation uses `wrf.interplevel` with pressure in hPa and all four levels
in one call per field. Full pressure, winds, temperature and height diagnostics
are reused from the lowest-level calculations. Levels outside the resolved
vertical range, including below terrain, are missing (`_FillValue`); there is
no extrapolation. Each variable also records `pressure_level` and its units.

Move aside `WRFOUT.OK` to regenerate hourly and monthly outputs with the new
variables. Existing output files are not updated merely by replacing the script.

Pressure-output unit tests (mocked WRF/NetCDF interfaces):
`python -m unittest discover -s tests -p 'test_pressure_levels.py' -v`.
A real WRF/CDO integration run is still needed on the processing host.
