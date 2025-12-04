"""
Snow Water Equivalent (SWE) Model Validation and Spatial Analysis Script

This script processes snow model predictions against observed SWE data from the
CanSWE dataset across multiple stations in British Columbia. It generates
comprehensive validation metrics and visualizations including:

1. Spatial maps showing the geographic distribution of model performance metrics
   (RMSE, R², MAE) across all stations
2. Time series comparisons of best and worst performing stations for each metric
3. Summary statistics and CSV export of validation results

The analysis filters observations based on a SWE threshold (default: 10 mm) to
focus on periods with meaningful snow cover. Output includes high-resolution
PNG maps, time series plots, and a metrics summary CSV file.

Workflow:
---------
1. Load observed SWE data from NetCDF file
2. Process all model output CSV files, matching them to observation stations
3. Calculate validation metrics (RMSE, R², MAE) for each station
4. Generate spatial scatter plots showing metric distribution across BC region
5. Create time series plots highlighting best/worst case scenarios
6. Export summary statistics and metrics to CSV

Configuration:
---------------
MODEL_DIR : Path
    Directory containing model output files (named as *_snow_timeseries.csv)
OBS_FILE : str
    Path to CanSWE NetCDF observation file (1928-2023)
OUTPUT_DIR : Path
    Directory where output plots and CSV files are saved
SWE_THRESHOLD : float
    Minimum SWE value (mm) for metric calculation (filters noise at low SWE)

Output Files:
-------------
- rmse_spatial_map.png : Geographic scatter plot of RMSE values
- r2_spatial_map.png : Geographic scatter plot of R² scores
- mae_spatial_map.png : Geographic scatter plot of MAE values
- rmse_best_worst_timeseries.png : Time series comparison for RMSE metric
- r2_best_worst_timeseries.png : Time series comparison for R² metric
- mae_best_worst_timeseries.png : Time series comparison for MAE metric
- metrics_summary.csv : Tabular summary of all station metrics
"""
import pandas as pd
import xarray as xr
import numpy as np
from pathlib import Path
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Configuration
MODEL_DIR = Path("climatex/snow/stations_data/")  # Update this path
OBS_FILE = "/CanSWE-CanEEN_1928-2023_v6.nc"  # Update this path
OUTPUT_DIR = Path("outputs")
SWE_THRESHOLD = 10  # mm

OUTPUT_DIR.mkdir(exist_ok=True)

# Load observations
obs_ds = xr.open_dataset(OBS_FILE)

def calculate_metrics(obs, model, threshold=SWE_THRESHOLD):
    """Calculate RMSE, R2, MAE for observations above threshold."""
    df = pd.DataFrame({'obs': obs, 'model': model}).dropna()

    if threshold > 0:
        df = df[df['obs'] > threshold]

    if len(df) < 10:  # Minimum samples required
        return None, None, None

    rmse = np.sqrt(mean_squared_error(df['obs'], df['model']))
    mae = mean_absolute_error(df['obs'], df['model'])
    r2 = r2_score(df['obs'], df['model'])

    return rmse, r2, mae

def process_station(model_file, obs_ds):
    """Process single station: read model data and calculate metrics."""
    station_id = model_file.stem.split('_')[0]

    # Check if station exists in observations
    if station_id not in obs_ds.station_id.values:
        return None

    # Read model data
    model_df = pd.read_csv(model_file, parse_dates=[0], index_col=0)
    model_df.columns = ['model']

    # Extract observations for this station
    obs_data = obs_ds.sel(station_id=station_id)
    obs_df = obs_data['snw'].to_dataframe()[['snw']].rename(columns={'snw': 'obs'})

    # Merge on time
    merged = model_df.join(obs_df, how='inner')

    if len(merged) == 0:
        return None

    rmse, r2, mae = calculate_metrics(merged['obs'], merged['model'])

    if rmse is None:
        return None

    return {
        'station_id': station_id,
        'lat': float(obs_data['lat'].values),
        'lon': float(obs_data['lon'].values),
        'rmse': rmse,
        'r2': r2,
        'mae': mae,
        'timeseries': merged
    }

# Process all stations
print("Processing stations...")
results = []
for model_file in MODEL_DIR.glob("*_snow_timeseries.csv"):
    result = process_station(model_file, obs_ds)
    if result:
        results.append(result)
        print(f"  {result['station_id']}: RMSE={result['rmse']:.2f}, R2={result['r2']:.3f}, MAE={result['mae']:.2f}")

results_df = pd.DataFrame([{k: v for k, v in r.items() if k != 'timeseries'} for r in results])

print(f"\nProcessed {len(results)} stations successfully")

# Plotting function for spatial maps
def plot_metric_map(results_df, metric, title, cmap='RdYlGn_r', vmin=None, vmax=None):
    """Create spatial map of a metric for BC region."""
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

    # BC boundaries (approximate)
    ax.set_extent([-139, -114, 48, 60], crs=ccrs.PlateCarree())

    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    ax.add_feature(cfeature.STATES, linewidth=0.3, edgecolor='gray')
    ax.gridlines(draw_labels=True, alpha=0.3)

    # Reverse colormap for R2 (higher is better)
    if metric == 'r2':
        cmap = 'RdYlGn'
        vmin = vmin or 0
        vmax = vmax or 1

    scatter = ax.scatter(results_df['lon'], results_df['lat'],
                        c=results_df[metric], cmap=cmap,
                        s=100, alpha=0.7, edgecolors='black', linewidth=0.5,
                        vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree())

    plt.colorbar(scatter, ax=ax, label=metric.upper(), shrink=0.7)
    ax.set_title(title, fontsize=14, fontweight='bold')

    return fig

# Create spatial plots
print("\nGenerating spatial plots...")
metrics_config = [
    ('rmse', 'Root Mean Square Error (mm)', 'RdYlGn_r', None, None),
    ('r2', 'R² Score', 'RdYlGn', 0, 1),
    ('mae', 'Mean Absolute Error (mm)', 'RdYlGn_r', None, None)
]

for metric, title, cmap, vmin, vmax in metrics_config:
    fig = plot_metric_map(results_df, metric, title, cmap, vmin, vmax)
    plt.savefig(OUTPUT_DIR / f'{metric}_spatial_map.png', dpi=300, bbox_inches='tight')
    plt.show()
    print(f"  Saved {metric}_spatial_map.png")

# Time series plots: best and worst cases
def plot_best_worst_timeseries(results, metric):
    """Plot time series for best and worst stations for a given metric."""
    results_sorted = sorted(results, key=lambda x: x[metric], reverse=(metric != 'r2'))

    # For R2, higher is better; for RMSE/MAE, lower is better
    if metric == 'r2':
        best = results_sorted[-1]
        worst = results_sorted[0]
    else:
        best = results_sorted[0]
        worst = results_sorted[-1]

    fig, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=True)

    for ax, station, label in zip(axes, [best, worst], ['Best', 'Worst']):
        ts = station['timeseries']
        ax.plot(ts.index, ts['obs'], label='Observed', alpha=0.7, linewidth=1.5)
        ax.plot(ts.index, ts['model'], label='Model', alpha=0.7, linewidth=1.5)
        ax.set_ylabel('SWE (mm)', fontsize=11)
        ax.legend(loc='upper right')
        ax.grid(alpha=0.3)

        title_text = (f"{label} {metric.upper()}: {station['station_id']} "
                     f"(RMSE={station['rmse']:.2f}, R²={station['r2']:.3f}, MAE={station['mae']:.2f})")
        ax.set_title(title_text, fontsize=11, fontweight='bold')

    axes[-1].set_xlabel('Date', fontsize=11)
    plt.tight_layout()

    return fig

print("\nGenerating time series plots...")
for metric in ['rmse', 'r2', 'mae']:
    fig = plot_best_worst_timeseries(results, metric)
    plt.savefig(OUTPUT_DIR / f'{metric}_best_worst_timeseries.png', dpi=300, bbox_inches='tight')
    plt.show()
    print(f"  Saved {metric}_best_worst_timeseries.png")

# Save summary statistics
results_df.to_csv(OUTPUT_DIR / 'metrics_summary.csv', index=False)
print(f"\nSaved metrics_summary.csv")

print("\n=== Summary Statistics ===")
print(results_df[['rmse', 'r2', 'mae']].describe())
