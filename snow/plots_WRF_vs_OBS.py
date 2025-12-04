"""
Snow Water Equivalent (SWE) Model Visualization and Analysis Script

This module provides visualization tools for analyzing snow water equivalent model
performance across the British Columbia region. It includes functions to create
spatial maps of model metrics and time series plots comparing model predictions
against observations.

Module Contents:
    - ADJACENT_CELLS: List of adjacent grid cell identifiers used in spatial analysis
    - COLUMN_NAMES: Standard column names for snow metrics (center, standard deviation, SWE)
    - plot_metric_map(): Create spatial scatter plots of metrics over BC region
    - plot_timeseries_all_gridcells(): Time series visualization for a single station
    - plot_best_worst_timeseries(): Comparative time series of best and worst performing stations

Key Features:
    - Geographic visualization using Cartopy with BC regional boundaries
    - Time series comparison of observations vs. model predictions
    - Support for center grid cell and 8 adjacent grid cells (N, NE, E, SE, S, SW, W, NW)
    - Model performance metrics display (RMSE, R², MAE)
    - Best/worst case analysis based on user-selected metrics
"""
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

ADJACENT_CELLS = ['SNOW_N', 'SNOW_NE', 'SNOW_E', 'SNOW_SE',
                  'SNOW_S', 'SNOW_SW', 'SNOW_W', 'SNOW_NW']

COLUMN_NAMES = ["SNOW_center", "sd",'swe']

# Plotting function for spatial maps
def plot_metric_map(results_df, metric, title, cmap='viridis', vmin=None, vmax=None):
    """Create spatial map of a metric for BC region."""
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

    # BC boundaries (approximate)
    ax.set_extent([-139, -114, 48, 60], crs=ccrs.PlateCarree())

    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    ax.add_feature(cfeature.STATES, linewidth=0.3, edgecolor='gray')
    ax.gridlines(draw_labels=True, alpha=0.3)

    scatter = ax.scatter(results_df['lon'], results_df['lat'],
                        c=results_df[metric], cmap=cmap,
                        s=100, alpha=0.7, edgecolors='black', linewidth=0.5,
                        vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree())

    plt.colorbar(scatter, ax=ax, label=metric.upper(), shrink=0.7)
    ax.set_title(title, fontsize=14, fontweight='bold')

    return fig


# Time series plots with all grid cells
def plot_timeseries_all_gridcells(station_result):
    """
    Plot time series for a station with observations and all grid cells.
    - Solid lines: Observations and SNOW_center
    - Dashed lines: Adjacent grid cells
    """
    ts = station_result['timeseries']

    fig, ax = plt.subplots(figsize=(14, 6))

    # Plot observations (solid line)
    if 'obs' in ts.columns:
        ax.plot(ts.index, ts['obs'], label='Observed', linewidth=2,
                linestyle='-', color='black', alpha=0.8)

    found_col = [col for col in COLUMN_NAMES if col in ts.columns]

    # Plot SNOW_center (solid line)
    if found_col[0] in ts.columns:
        ax.plot(ts.index, ts[found_col[0]], label='Model (Center)',
                linewidth=2, linestyle='-', alpha=0.7)

    for cell in ADJACENT_CELLS:
        if cell in ts.columns:
            ax.plot(ts.index, ts[cell], label=f'Model ({cell.replace("SNOW_", "")})',
                   linewidth=1.5, linestyle='--', alpha=0.6)

    ax.set_ylabel('SWE (mm)', fontsize=12)
    ax.set_xlabel('Date', fontsize=12)
    ax.legend(loc='upper right', fontsize=9, ncol=2)
    ax.grid(alpha=0.3)

    title_text = (f"Station {station_result['station_id']}: "
                 f"RMSE={station_result['rmse']:.2f}, "
                 f"R²={station_result['r2']:.3f}, MAE={station_result['mae']:.2f}")
    ax.set_title(title_text, fontsize=12, fontweight='bold')

    plt.tight_layout()
    return fig


# Time series plots: best and worst cases
def plot_best_worst_timeseries(results, metric):
    """Plot time series for best and worst stations for a given metric."""
    results_sorted = sorted(results, key=lambda x: x[metric], reverse=(metric == 'r2'))

    # For R2, higher is better; for RMSE/MAE, lower is better
    if metric == 'r2':
        best = results_sorted[-1]
        worst = results_sorted[0]
    else:
        best = results_sorted[0]
        worst = results_sorted[-1]

    fig, axes = plt.subplots(2, 1, figsize=(14, 10), sharex=True)

    for ax, station, label in zip(axes, [best, worst], ['Best', 'Worst']):
        ts = station['timeseries']

        # Plot observations (solid line)
        if 'obs' in ts.columns:
            ax.plot(ts.index, ts['obs'], label='Observed', linewidth=2.5,
                   linestyle='-', color='black', alpha=0.8)

        found_col = [col for col in COLUMN_NAMES if col in ts.columns]
        if found_col[0] in ts.columns:
            ax.plot(ts.index, ts[found_col[0]], label='Model (Center)',
                   linewidth=2.5, linestyle='-', alpha=0.7)

        for cell in ADJACENT_CELLS:
            if cell in ts.columns:
                ax.plot(ts.index, ts[cell], label=f'Model ({cell.replace("SNOW_", "")})',
                       linewidth=1.5, linestyle='--', alpha=0.6)

        ax.set_ylabel('SWE (mm)', fontsize=11)
        ax.legend(loc='upper right', fontsize=9, ncol=2)
        ax.grid(alpha=0.3)

        title_text = (f"{label} {metric.upper()}: {station['station_id']} "
                     f"(RMSE={station['rmse']:.2f}, R²={station['r2']:.3f}, MAE={station['mae']:.2f})")
        ax.set_title(title_text, fontsize=11, fontweight='bold')

    axes[-1].set_xlabel('Date', fontsize=11)
    plt.tight_layout()

    return fig
