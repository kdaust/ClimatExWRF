"""
WRF Snow Water Equivalent (SWE) Visualization Script

This script loads monthly snow water equivalent (SWE) data from a WRF (Weather
Research and Forecasting) model NetCDF output file and generates a 2x3 subplot
visualization covering the winter season (November through April).

The script:
- Loads SWE data from a specified NetCDF file
- Extracts data for six consecutive months spanning the snow season
- Creates a multi-panel figure with individual maps for each month
- Uses a progressive 'Reds' colormap to show SWE accumulation
- Caps displayed values at 1000 mm for consistent visualization
- Saves the output figure at high resolution (300 DPI)

Configuration:
    WRF_FILE: Path to the input NetCDF file containing monthly SWE statistics
    SWE_VAR: Name of the SWE variable in the NetCDF file
    MONTHS: List of month numbers to extract (11=Nov, 12=Dec, 1=Jan, etc.)
    MONTH_NAMES: Corresponding month names for plot titles
    OUTPUT_DIR: Directory where the output PNG will be saved
    CMAP_SWE: Matplotlib colormap name ('Reds' for progressive display)
    SWE_MAX: Maximum SWE value for colorbar scaling (mm)
    DPI: Output figure resolution (dots per inch)

Output:
    Generates a PNG file showing 6-panel grid of monthly SWE maps with
    individual colorbars, saved to OUTPUT_DIR.
"""
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# Configuration
WRF_FILE = 'climatex/wrf_data/snow_stats_month.nc'
SWE_VAR = 'SNOW_mean'

MONTHS = [11, 12, 1, 2, 3, 4]
MONTH_NAMES = ['November', 'December', 'January', 'February', 'March', 'April']

OUTPUT_DIR = Path('datasets/climatex/')
OUTPUT_DIR.mkdir(exist_ok=True)

CMAP_SWE = 'Reds'  # Progressive colormap for SWE
SWE_MAX = 1000  # mm
DPI = 300

def load_wrf_data():
    """Load WRF SWE data."""
    ds = xr.open_dataset(WRF_FILE)
    return ds[SWE_VAR].sel(Times=MONTHS)

def plot_wrf_swe(data, months, month_names, output_dir):
    """Plot SWE data from WRF model with progressive colormap."""
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    vmin, vmax = 0, SWE_MAX

    for idx, (month, month_name) in enumerate(zip(months, month_names)):
        data_slice = np.ma.masked_invalid(data.sel(Times=month).values)
        # Cap values at SWE_MAX
        data_slice = np.ma.masked_greater(data_slice, SWE_MAX)

        im = axes[idx].imshow(data_slice, cmap=CMAP_SWE, origin='lower', vmin=vmin, vmax=vmax)
        axes[idx].set_title(f'Month: {month_name}', fontsize=12, fontweight='bold')
        axes[idx].set_xlabel('West-East')
        axes[idx].set_ylabel('South-North')

        cbar = plt.colorbar(im, ax=axes[idx])
        cbar.set_label(f'{SWE_VAR} (mm)')
        cbar.set_ticks([0, 250, 500, 750, 1000])
        cbar.set_ticklabels(['0', '250', '500', '750', f'>{SWE_MAX}'])

    title = f'{SWE_VAR} - WRF Model'
    fig.suptitle(title, fontsize=16, fontweight='bold', y=0.995)
    plt.tight_layout()

    filename = output_dir / f"WRF_SWE_{SWE_VAR.lower().replace('_', '')}.png"
    plt.savefig(filename, dpi=DPI, bbox_inches='tight')
    print(f"Saved: {filename}")
    plt.show()

def main():
    """Main execution."""
    print(f"Loading WRF data from {WRF_FILE}...")
    data = load_wrf_data()

    print("Plotting WRF SWE data...")
    plot_wrf_swe(data, MONTHS, MONTH_NAMES, OUTPUT_DIR)

    print("Done!")

if __name__ == '__main__':
    main()
