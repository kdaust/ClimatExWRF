"""
Analyze and compare annual snow water equivalent (SWE) cycles across multiple climate models.

This script loads monthly snow statistics from three climate models (WRF, ERA5Land, and CaSR),
applies an ice mask to exclude pixels with persistent snow cover, computes area-averaged SWE,
and generates a comparative visualization of the annual cycle aligned to the water year
(August-July).

Key Features:
    - Loads NetCDF datasets from three different climate models
    - Applies ice masking based on minimum snow depth to exclude glaciated regions
    - Computes spatially-averaged snow water equivalent (SWE) for each timestep
    - Reorders monthly data to water year convention (Aug-Jul)
    - Calculates long-term climatology for each model
    - Generates a multi-model comparison plot with distinct colors and markers

Data Processing Steps:
    1. Load monthly snow statistics from each model's NetCDF file
    2. Create ice mask where minimum snow depth > 0 throughout the study period
    3. Mask out ice-covered pixels and compute spatial mean SWE
    4. Extract water year (Aug-Jul) climatology
    5. Plot all models on a single figure for visual comparison

Output:
    - Console messages indicating successful data loading for each model
    - PNG figure showing annual SWE cycle for all models (default: /annual_swe_cycle.png)

Author notes:
    - Models are processed independently; failure of one does not halt others
    - Ice mask based on SNOW_min variable ensures analysis focuses on seasonal snow
    - Spatial averaging removes orographic complexity for regional-scale comparison
"""
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

MODELS = {
    'WRF': {'file': '/wrf_data/snow_stats_month.nc', 'color': 'blue', 'marker': 'o'},
    'ERA5Land': {'file': '/era5_dataset/snow_stats_month.nc', 'color': 'red', 'marker': 's'},
    'CaSR': {'file': '/casr_dataset/snow_stats_month.nc', 'color': 'green', 'marker': '^'}
}
VAR_NAME = 'SNOW_mean'
WATER_YEAR_ORDER = [7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 6]
MONTHS = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul']

def load_and_process(filepath, var_name, ice_var='SNOW_min'):
    """Load data, apply ice mask, and return area-averaged SWE."""
    ds = xr.open_dataset(filepath)
    ice_mask = ds[ice_var].min(dim='Times') > 0
    masked_data = ds[var_name].where(~ice_mask, np.nan)
    return masked_data.mean(dim=['south_north', 'west_east'])

def water_year_climatology(data):
    """Reorder to water year (Aug-Jul) and compute climatology."""
    reshaped = data.values.reshape(-1, 12)[:, WATER_YEAR_ORDER]
    return reshaped.mean(axis=0)

def plot_annual_cycle(results, output_file='climatex/annual_swe_cycle.png'):
    """Plot annual SWE cycle for all models."""
    fig, ax = plt.subplots(figsize=(12, 7))

    for model, (data, color, marker) in results.items():
        ax.plot(data, marker=marker, linewidth=2.5, markersize=8,
                label=model, color=color)

    ax.set(xlabel='Month', ylabel='Area-Averaged SWE (mm)',
           xticks=range(12), xticklabels=MONTHS)
    ax.set_title('Annual Snow Water Equivalent (Aug-Jul)', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.legend(fontsize=11, loc='upper right')
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()

# Process and plot
results = {}
for model, config in MODELS.items():
    try:
        data = load_and_process(config['file'], VAR_NAME)
        climatology = water_year_climatology(data)
        results[model] = (climatology, config['color'], config['marker'])
        print(f"✓ {model}")
    except Exception as e:
        print(f"✗ {model}: {e}")

if results:
    plot_annual_cycle(results)
