"""
Compare snow statistics between two climate models by computing and visualizing their spatial ratio.

This script loads monthly snow data from two climate models (WRF and ERA5Land),
computes the pixel-wise ratio between them, and generates a 6-panel comparison plot
(one for each month from November to April). The ratio is visualized using a
diverging colormap centered at 1.0 (100%), where red indicates higher values in the
numerator model and blue indicates lower values.

The script performs the following operations:
1. Loads snow statistics from NetCDF files for both models
2. Extracts the specified snow variable for winter months (Nov-Apr)
3. Computes the ratio: model2 / model1
4. Creates a 2x3 subplot figure with one panel per month
5. Applies a TwoSlopeNorm colormap centered at 100% (ratio = 1.0)
6. Saves the comparison plot as a high-resolution PNG file

Configuration:
    - MODELS: Dictionary specifying file paths and variable names for each model
    - MONTHS: List of months to include (11, 12, 1, 2, 3, 4)
    - VMIN/VMAX: Ratio limits (0.5 to 2.0, representing 50% to 200%)
    - CMAP: Diverging colormap (RdBu_r) for intuitive interpretation
    - OUTPUT_DIR: Directory where output plots are saved

Output:
    - PNG file containing 6 subplots with colorbars showing model ratio comparisons
"""
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from matplotlib.colors import TwoSlopeNorm

# Dictionary with model information
MODELS = {
    'WRF': {
        'file': 'climatex/wrf_data/snow_stats_month.nc',
        'var': 'SNOW_mean',
    },
    # 'CaSR': {
    #     'file': 'climatex/casr_dataset/snow_stats_month.nc',
    #     'var': 'SNOW_mean',
    # }
    'ERA5Land': {
        'file': 'climatex/era5_dataset/snow_stats_month.nc',
        'var': 'SNOW_max',
    }
}

OUTPUT_DIR = Path('/climatex/')

MONTHS = [11, 12, 1, 2, 3, 4]
MONTH_NAMES = ['November', 'December', 'January', 'February', 'March', 'April']

# Ratio limits: 50% (0.5) to 200% (2.0), centered at 100% (1.0)
VMIN, VMAX = 0.5, 2.0
VCENTER = 1.0
CMAP = 'RdBu_r'  # RdBu reversed: red for >100%, blue for <100%
DPI = 300

def load_datasets(models_dict):
    return {name: xr.open_dataset(info['file']) for name, info in models_dict.items()}

def extract_model_data(datasets, models_dict, months):
    return {name: datasets[name][models_dict[name]['var']].sel(Times=months)
            for name in models_dict.keys()}

def plot_ratio(model_names, data_dict, months, month_names, var_name, output_dir):
    """Plot ratio between two models (model1/model2)."""
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    if len(model_names) != 2:
        print(f"Error: This function requires exactly 2 models, got {len(model_names)}")
        return

    # Calculate ratio: CaSR/WRF (or model_names[1]/model_names[0])
    # Avoid division by zero
    data_plot = np.divide(
        data_dict[model_names[1]].values,
        data_dict[model_names[0]].values,
        where=data_dict[model_names[0]].values != 0,
        out=np.full_like(data_dict[model_names[0]].values, np.nan, dtype=float)
    )

    # Create TwoSlopeNorm to center the colormap at 1.0 (100%)
    norm = TwoSlopeNorm(vmin=VMIN, vcenter=VCENTER, vmax=VMAX)

    for idx, (month, month_name) in enumerate(zip(months, month_names)):
        data_slice = np.ma.masked_invalid(data_plot[idx])
        im = axes[idx].imshow(data_slice, cmap=CMAP, origin='lower', norm=norm)
        axes[idx].set_title(f'Month: {month_name}', fontsize=12, fontweight='bold')
        axes[idx].set_xlabel('West-East')
        axes[idx].set_ylabel('South-North')

        cbar = plt.colorbar(im, ax=axes[idx])
        cbar.set_label('Ratio')
        cbar.set_ticks([0.5, 0.75, 1.0, 1.5, 2.0])
        cbar.set_ticklabels(['50%', '75%', '100%', '150%', '200%'])

    title = f'{var_name} Ratio ({model_names[1]} / {model_names[0]})\nRed: Higher in {model_names[1]}, Blue: Lower in {model_names[1]}'
    fig.suptitle(title, fontsize=16, fontweight='bold', y=0.995)
    plt.tight_layout()

    filename = output_dir / f"{'_'.join(model_names)}_ratio_{var_name.lower().replace('_', '')}.png"
    plt.savefig(filename, dpi=DPI, bbox_inches='tight')
    print(f"Saved: {filename}")
    plt.show()

def main():
    datasets = load_datasets(MODELS)
    model_names = list(MODELS.keys())
    var_name = MODELS[model_names[0]]['var']

    if len(model_names) != 2:
        print(f"Error: Expected 2 models for ratio plot, got {len(model_names)}")
        return

    model_data = extract_model_data(datasets, MODELS, MONTHS)
    plot_ratio(model_names, model_data, MONTHS, MONTH_NAMES, var_name, OUTPUT_DIR)

    for ds in datasets.values():
        ds.close()
    print("Done!")

if __name__ == '__main__':
    main()
