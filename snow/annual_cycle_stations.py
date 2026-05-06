"""
Annual Snow Water Equivalent (SWE) Cycle Comparison Tool

This script processes and compares seasonal snow water equivalent data from multiple
climate models (WRF, ERA5Land, CaSR) against observed measurements from the
CanSWE-CanEEN dataset. It calculates mean annual cycles for the water year (August-July)
and generates a comparative visualization.

The script performs the following operations:
1. Loads observation data from a NetCDF file (CanSWE-CanEEN dataset)
2. Processes model output CSV files for matching weather stations
3. Aligns model and observation data by common dates
4. Calculates mean monthly SWE values across the water year cycle
5. Generates a multi-model comparison plot with observations

Input Requirements:
    - Observation file: NetCDF dataset with 'snw' variable and 'station_id' dimension
    - Model data: CSV files named '{station_id}_snow_timeseries.csv' in model-specific folders
    - File structure: {folder}/{station_id}_snow_timeseries.csv for each model

Output:
    - Console output: Mean SWE values for each model and observations
    - PNG plot: Comparative annual SWE cycle visualization saved to 'annual_swe_cycle_stations.png'

Configuration:
    - MODELS: Dictionary defining model folders, plot colors, and markers
    - OBS_FILE: Path to observation NetCDF file
    - OBS_VAR: Variable name for SWE in observation file
    - WATER_YEAR_MONTHS: Month indices defining water year (Aug-Jul)
"""
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')


MODELS = {
    'WRF': {'folder': '/wrf_data/', 'color': 'blue', 'marker': 'o'},
    'ERA5Land': {'folder': '/era5_data/', 'color': 'red', 'marker': 's'},
    'CaSR': {'folder': '/casr_data/', 'color': 'green', 'marker': '^'}
}
OBS_FILE = '/CanSWE-CanEEN_1928-2023_v6.nc'
OBS_VAR = 'snw'
WATER_YEAR_MONTHS = [7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6]
MONTHS = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul']

def load_csv(filepath):
    """Load CSV and return aligned SWE series."""
    df = pd.read_csv(filepath, parse_dates=[0], index_col=0)
    df.index = pd.to_datetime(df.index.strftime('%Y-%m-%d')).normalize()
    return df.iloc[:, 0].dropna()

def get_water_year_month(date):
    m = date.month
    return m - 8 if m >= 8 else m + 4

def calc_annual_cycle(swe_series):
    if len(swe_series) == 0:
        return None
    df = pd.DataFrame({'swe': swe_series}).assign(wym=swe_series.index.map(get_water_year_month))
    return df.groupby('wym')['swe'].mean().reindex(range(12)).values

def process_folder(folder_path, obs_ds):
    """Process CSVs and align with obs for same stations only."""
    all_cycles = []

    for csv_file in Path(folder_path).glob('*_snow_timeseries.csv'):
        station_id = csv_file.name.split('_')[0]
        model_df = load_csv(csv_file)
        if len(model_df) == 0 or obs_ds is None:
            continue

        try:
            obs_series = obs_ds[OBS_VAR].sel(station_id=station_id).to_series().dropna()
            obs_series.index = pd.to_datetime(obs_series.index.strftime('%Y-%m-%d')).normalize()

            common_dates = model_df.index.intersection(obs_series.index)
            if len(common_dates) > 0:
                aligned = model_df.loc[common_dates]
                mean = calc_annual_cycle(aligned)
                if mean is not None:
                    all_cycles.append(mean)
        except:
            continue

    return np.nanmean(np.array(all_cycles), axis=0) if all_cycles else None

def process_obs(folder_path, obs_ds):
    """Calculate obs annual cycle for stations in folder only."""
    all_cycles = []

    for csv_file in Path(folder_path).glob('*_snow_timeseries.csv'):
        station_id = csv_file.name.split('_')[0]
        try:
            obs_series = obs_ds[OBS_VAR].sel(station_id=station_id).to_series().dropna()
            if len(obs_series) > 0:
                mean = calc_annual_cycle(obs_series)
                if mean is not None:
                    all_cycles.append(mean)
        except:
            continue

    return np.nanmean(np.array(all_cycles), axis=0) if all_cycles else None

def plot_results(results, obs_mean, output='climatex/annual_swe_cycle_stations.png'):
    fig, ax = plt.subplots(figsize=(14, 8))
    months = np.arange(12)

    if obs_mean is not None:
        ax.plot(months, obs_mean, marker='D', linewidth=2.5, markersize=8,
               label='Observations', color='black')

    for model, (mean, color, marker) in results.items():
        ax.plot(months, mean, marker=marker, linewidth=2.5, markersize=8,
               label=model, color=color)

    ax.set(xlabel='Month', ylabel='SWE (mm)', xticks=range(12), xticklabels=MONTHS)
    ax.set_title('Annual Snow Water Equivalent Cycle (Aug-Jul)', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.legend(fontsize=11)
    plt.tight_layout()
    plt.savefig(output, dpi=300, bbox_inches='tight')
    print(f"Plot saved to {output}")
    plt.show()

# Main
obs_ds = xr.open_dataset(OBS_FILE) if Path(OBS_FILE).exists() else None
results = {}

first_model_folder = list(MODELS.values())[0]['folder']
obs_mean = process_obs(first_model_folder, obs_ds) if obs_ds else None

for model, config in MODELS.items():
    mean = process_folder(config['folder'], obs_ds)
    if mean is not None:
        results[model] = (mean, config['color'], config['marker'])
        print(f"✓ {model}: {np.nanmean(mean):.1f} mm")
    else:
        print(f"✗ {model}")

if obs_mean is not None:
    print(f"✓ Observations: {np.nanmean(obs_mean):.1f} mm")

if results:
    plot_results(results, obs_mean)
