"""
Compare snow statistics between two climate/weather models.

This script loads snow data from multiple model datasets, computes the difference
between the first two models, and generates a 6-panel visualization showing monthly
differences in maximum snow depth from November through April.

The script is configured to work with WRF and ERA5-Land datasets, though other
models can be added to the MODELS dictionary. The output is a PNG figure with
symmetric color scaling centered at zero, where red indicates positive differences
(first model higher) and blue indicates negative differences (second model lower).

Configuration:
    MODELS: Dictionary mapping model names to their dataset file paths and variable names
    MONTHS: List of month numbers to process (11, 12, 1-4 for Nov-Apr winter season)
    MONTH_NAMES: Display names corresponding to each month
    OUTPUT_DIR: Directory where output figures will be saved
    VMIN, VMAX: Symmetric limits for difference colormap (-200 to 200 mm)
    CMAP: Colormap name ('RdBu_r' reversed so red=positive, blue=negative)
    DPI: Output figure resolution in dots per inch

Output:
    - PNG image with 6 subplots showing monthly snow depth differences
"""
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# Dictionary with model information
MODELS = {
    'WRF': {
        'file': 'climatex/wrf_data/snow_stats_month.nc',
        'var': 'SNOW_max',
    },
    # 'CaSR': {
    #     'file': '/climatex/casr_dataset/snow_stats_month.nc',
    #     'var': 'SNOW_max',
    # }
    'ERA5Land': {
        'file': 'climatex/era5_dataset/snow_stats_month.nc',
        'var': 'SNOW_max',
    }
}

MONTHS = [11, 12, 1, 2, 3, 4]
MONTH_NAMES = ['November', 'December', 'January', 'February', 'March', 'April']

OUTPUT_DIR = Path('/datasets/climatex/')

VMIN, VMAX = -200, 200  # Difference limits (symmetric)
CMAP = 'RdBu_r'  # RdBu reversed so red is positive, blue is negative
DPI = 300

def load_datasets(models_dict):
    return {name: xr.open_dataset(info['file']) for name, info in models_dict.items()}

def extract_model_data(datasets, models_dict, months):
    return {name: datasets[name][models_dict[name]['var']].sel(Times=months)
            for name in models_dict.keys()}

def plot_comparison(model_names, data_dict, months, month_names, var_name, output_dir):
    """Plot difference between two models."""
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    if len(model_names) != 2:
        print(f"Error: This function requires exactly 2 models, got {len(model_names)}")
        return

    data_plot = data_dict[model_names[0]] - data_dict[model_names[1]]
    vmin, vmax = VMIN, VMAX

    for idx, (month, month_name) in enumerate(zip(months, month_names)):
        data_slice = np.ma.masked_invalid(data_plot.sel(Times=month).values)
        im = axes[idx].imshow(data_slice, cmap=CMAP, origin='lower', vmin=vmin, vmax=vmax)
        axes[idx].set_title(f'Month: {month_name}', fontsize=12, fontweight='bold')
        axes[idx].set_xlabel('West-East')
        axes[idx].set_ylabel('South-North')

        cbar = plt.colorbar(im, ax=axes[idx])
        cbar.set_label('Difference (mm)')
        cbar.set_ticks([VMIN, -100, 0, 100, VMAX])
        cbar.set_ticklabels([f'{VMIN}', '-100', '0', '100', f'{VMAX}'])

    title = f'{var_name} Difference ({model_names[0]} - {model_names[1]})'
    fig.suptitle(title, fontsize=16, fontweight='bold', y=0.995)
    plt.tight_layout()

    filename = output_dir / f"{'_'.join(model_names)}_difference_{var_name.lower().replace('_', '')}.png"
    plt.savefig(filename, dpi=DPI, bbox_inches='tight')
    print(f"Saved: {filename}")
    plt.show()

def main():
    datasets = load_datasets(MODELS)
    model_names = list(MODELS.keys())
    var_name = MODELS[model_names[0]]['var']

    if len(model_names) != 2:
        print(f"Error: Expected 2 models for difference plot, got {len(model_names)}")
        return

    model_data = extract_model_data(datasets, MODELS, MONTHS)
    plot_comparison(model_names, model_data, MONTHS, MONTH_NAMES, var_name, OUTPUT_DIR)

    for ds in datasets.values():
        ds.close()
    print("Done!")

if __name__ == '__main__':
    main()
