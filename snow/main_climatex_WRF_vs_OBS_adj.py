"""
Evaluate snow water equivalent (SWE) model predictions against observational data.

This script compares ERA5 model snow predictions with CanSWE-CanEEN observational
data across multiple weather stations. It calculates validation metrics (RMSE, R², MAE),
generates time series and spatial comparison visualizations, and produces a summary
statistics report.

Workflow:
    1. Load observational SWE data from NetCDF file
    2. Process each station's model output CSV file:
       - Read model predictions and extract corresponding observations
       - Calculate performance metrics for SWE above a specified threshold
       - Interpolate observation gaps and compute extended SWE metrics
    3. Generate visualizations:
       - Time series plots for each station
       - Spatial maps of RMSE, R², and MAE across stations
       - Best/worst performing station comparisons
    4. Export metrics summary to CSV

Configuration:
    MODEL_DIR (Path): Directory containing model snow timeseries CSV files
    OBS_FILE (str): Path to NetCDF file with observational SWE data
    OUTPUT_DIR (Path): Directory for output plots and summary files
    SWE_THRESHOLD (float): Minimum SWE (mm) for metric calculation; filters out
                          low-snow periods to focus on meaningful accumulation events
    height_filter (str): Optional CSV file to filter stations by elevation differences

Output:
    - Individual station time series plots (PNG)
    - Spatial metric maps for RMSE, R², MAE (PNG)
    - Best/worst station comparison plots (PNG)
    - metrics_summary.csv: Aggregated performance statistics across all stations
"""
import pandas as pd
import xarray as xr
import numpy as np
from pathlib import Path
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
import matplotlib.pyplot as plt
from plots_WRF_vs_OBS import (plot_timeseries_all_gridcells, plot_metric_map,
                              plot_best_worst_timeseries)
from utils.evaluation import compute_metrics_swe


def calculate_metrics(obs, model, threshold):
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


def process_station(model_file, obs_ds, threshold):
    """Process single station: read model data and calculate metrics."""
    station_id = model_file.stem.split('_')[0]

    # Check if station exists in observations
    if station_id not in obs_ds.station_id.values:
        return None

    # Read model data
    model_df = pd.read_csv(model_file, parse_dates=[0], index_col=0)
    model_df.index = pd.to_datetime(model_df.index.strftime('%Y-%m-%d'))
    model_df.index = model_df.index.normalize()

    # Extract observations for this station
    obs_data = obs_ds.sel(station_id=station_id)
    obs_df = obs_data['snw'].to_dataframe()[['snw']].rename(columns={'snw': 'obs'})
    # Restrict obs_df to the same date range
    obs_df = obs_df.loc[(obs_df.index >= model_df.index.min()) & (obs_df.index <= model_df.index.max())]

    column_names = ["SNOW_center", "sd",'swe']
    found_col = [col for col in column_names if col in model_df.columns]

    # Merge on time - use found_col for metrics calculation
    merged = model_df[found_col].join(obs_df, how='inner')
    merged.rename(columns={found_col[0]: 'model'}, inplace=True)

    if len(merged) == 0:
        return None

    rmse, r2, mae = calculate_metrics(merged['obs'], merged['model'], threshold)

    if rmse is None:
        return None

    # For time series plotting, keep all grid cell data
    full_merged = model_df.join(obs_df, how='inner')

    # Interpolate observations linearly between first and last valid observation
    obs_interpolated = interpolate_observations(full_merged['obs'].copy())

    # Calculate extended metrics using interpolated observations
    swe_metrics = compute_metrics_swe(
        obs_interpolated,
        full_merged[found_col[0]],
        full_merged.index
    )

    return {
        'station_id': station_id,
        'lat': float(obs_data['lat'].values),
        'lon': float(obs_data['lon'].values),
        'rmse': rmse,
        'r2': r2,
        'mae': mae,
        'timeseries': full_merged,
        'metrics_timeseries': merged,
        **swe_metrics
    }


def interpolate_observations(obs_series):
    """
    Linearly interpolate observations between first and last valid observation.

    Args:
        obs_series (pd.Series): Series with observation values (may contain NaNs).

    Returns:
        pd.Series: Series with linearly interpolated values between first and last valid obs.
    """
    # Find first and last valid (non-NaN) indices
    valid_mask = obs_series.notna()
    valid_indices = valid_mask[valid_mask].index

    if len(valid_indices) < 2:
        # Not enough data to interpolate
        return obs_series

    first_valid_idx = valid_indices[0]
    last_valid_idx = valid_indices[-1]

    # Create a copy and interpolate only between first and last valid observations
    obs_interp = obs_series.copy()
    obs_interp = obs_interp.loc[first_valid_idx:last_valid_idx]
    obs_interp = obs_interp.interpolate(method='linear')

    # Pad back to original index if needed
    obs_result = obs_series.copy()
    obs_result.loc[first_valid_idx:last_valid_idx] = obs_interp

    return obs_result

if __name__ == '__main__':

    # Configuration
    MODEL_DIR = Path("/climatex/era5_data/")
    OBS_FILE = "/CanSWE-CanEEN_1928-2023_v6.nc"
    OUTPUT_DIR = Path("climatex/comparison_era5")
    SWE_THRESHOLD = 10  # mm
    height_filter = '/comparison_wrf/metrics_with_elevations.csv'

    OUTPUT_DIR.mkdir(exist_ok=True)

    # Load observations
    obs_ds = xr.open_dataset(OBS_FILE)

    # Process all stations
    print("Processing stations...")

    # if the difference between station and WRF height would be consider
    height_df = pd.read_csv(height_filter, index_col=0) if height_filter else None
    model_files = [f for f in MODEL_DIR.glob("*_snow_timeseries.csv")
                   if height_df is None or f.stem.split('_')[0] in height_df.index]
    results = []
    for model_file in model_files:
        result = process_station(model_file, obs_ds, SWE_THRESHOLD)
        if result:
            results.append(result)
            print(f"  {result['station_id']}: RMSE={result['rmse']:.2f}, R2={result['r2']:.3f}, MAE={result['mae']:.2f}")

    results_df = pd.DataFrame([{
        k: v for k, v in r.items()
        if k not in ['timeseries', 'metrics_timeseries', 'swe_metrics']
    } for r in results])

    print(f"\nProcessed {len(results)} stations successfully")

    # Generate time series plots for all stations
    print("\nGenerating time series plots for all stations...")
    for station_result in results:
        fig = plot_timeseries_all_gridcells(station_result)
        filename = f"timeseries_{station_result['station_id']}.png"
        plt.savefig(OUTPUT_DIR / filename, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Saved {filename}")

    # Create spatial plots
    print("\nGenerating spatial plots...")
    metrics_config = [
        ('rmse', 'Root Mean Square Error (mm)', 'viridis', None, None),
        #('r2', 'R² Score', 'RdYlGn', 0, 1),
        ('r2', 'R² Score', 'viridis', None, None),
        ('mae', 'Mean Absolute Error (mm)', 'viridis', None, None)
    ]

    for metric, title, cmap, vmin, vmax in metrics_config:
        fig = plot_metric_map(results_df, metric, title, cmap, vmin, vmax)
        plt.savefig(OUTPUT_DIR / f'{metric}_spatial_map.png', dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Saved {metric}_spatial_map.png")

    print("\nGenerating best/worst time series plots...")
    for metric in ['rmse', 'r2', 'mae']:
        fig = plot_best_worst_timeseries(results, metric)
        plt.savefig(OUTPUT_DIR / f'{metric}_best_worst_timeseries.png', dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Saved {metric}_best_worst_timeseries.png")

    # Save summary statistics
    results_df.to_csv(OUTPUT_DIR / 'metrics_summary.csv', index=False)
    print(f"\nSaved metrics_summary.csv")

    print("\n=== Summary Statistics ===")
    print(results_df[['rmse', 'r2', 'mae']].describe())
