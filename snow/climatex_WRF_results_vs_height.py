"""
Module for comparing observed and WRF model elevations and analyzing their impact on performance metrics.

This script provides utilities to:
- Extract elevation data from observational datasets and WRF model outputs.
- Compute elevation differences between observed stations and WRF grid points using nearest-neighbor matching.
- Integrate elevation data into a metrics DataFrame and optionally filter stations based on a maximum elevation difference.
- Generate visualizations to explore relationships between elevation and performance metrics, including:
    * Scatter plots of metrics vs. elevation (observed, WRF, and elevation difference).
    * Histograms and scatter plots of elevation differences.
    * Combined comparison plots for multiple metrics against multiple elevation types.

Usage:
------
This module is intended for analyzing how elevation discrepancies between observed stations and WRF model grids
affect performance metrics such as precipitation bias, RMSE, or correlation.
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from pathlib import Path
from numpy.linalg import LinAlgError

# Configuration
METRICS_FILE = "/comparison_wrf/metrics_summary.csv"
OBS_FILE = "/CanSWE-CanEEN_1928-2023_v6.nc"
WRF_FILE = "/climatex/static_d03.nc"
OUTPUT_DIR = Path("/climatex/comparison_wrf/")
MAX_ELEVATION_DIFF = 200  # Maximum allowed absolute elevation difference (m), set to None to disable

OUTPUT_DIR.mkdir(exist_ok=True)


# ============================================================================
# DATA LOADING AND PROCESSING
# ============================================================================

def load_datasets(metrics_file, obs_file, wrf_file):
    """Load all required datasets."""
    print("Loading data...")
    metrics_df = pd.read_csv(metrics_file)
    obs_ds = xr.open_dataset(obs_file)
    wrf_ds = xr.open_dataset(wrf_file)
    return metrics_df, obs_ds, wrf_ds


def extract_obs_elevations(metrics_df, obs_ds):
    """Extract observed elevations for each station."""
    print("Extracting observation elevations...")
    elevations = []
    for station_id in metrics_df['station_id']:
        elev = float(obs_ds.sel(station_id=station_id)['elevation'].values)
        elevations.append(elev)
    return elevations


def extract_wrf_elevations(metrics_df, obs_ds, wrf_ds):
    """Extract WRF elevations for each station using nearest neighbor."""
    print("Extracting WRF elevations for each station...")

    # Get WRF height data - it's time-independent, so take first time step
    if 'HGT' in wrf_ds:
        if 'Times' in wrf_ds['HGT'].dims or 'XTIME' in wrf_ds['HGT'].dims:
            wrf_hgt = wrf_ds['HGT'].isel({list(wrf_ds['HGT'].dims)[0]: 0})
        else:
            wrf_hgt = wrf_ds['HGT']
    else:
        raise ValueError("HGT variable not found in WRF file. Please check your WRF dataset.")

    # Get lat/lon coordinates (as numpy arrays for easier calculation)
    wrf_lat = wrf_ds['XLAT'].values
    wrf_lon = wrf_ds['XLONG'].values
    wrf_hgt_values = wrf_hgt.values

    wrf_elevations = []
    for station_id in metrics_df['station_id']:
        station_data = obs_ds.sel(station_id=station_id)
        station_lat = float(station_data['lat'].values)
        station_lon = float(station_data['lon'].values)

        # Calculate distance to all grid points
        lat_diff = wrf_lat - station_lat
        lon_diff = wrf_lon - station_lon
        distance = np.sqrt(lat_diff**2 + lon_diff**2)

        # Find the indices of the minimum distance
        min_idx = np.unravel_index(distance.argmin(), distance.shape)

        # Get the corresponding height value
        wrf_elev = float(wrf_hgt_values[min_idx])
        wrf_elevations.append(wrf_elev)

    return wrf_elevations


def prepare_elevation_data(metrics_df, obs_ds, wrf_ds, max_diff=None):
    """Add all elevation data to metrics dataframe and filter by elevation difference."""
    metrics_df['obs_elevation'] = extract_obs_elevations(metrics_df, obs_ds)
    metrics_df['wrf_elevation'] = extract_wrf_elevations(metrics_df, obs_ds, wrf_ds)
    metrics_df['elevation_diff'] = metrics_df['wrf_elevation'] - metrics_df['obs_elevation']
    metrics_df['elevation_abs_diff'] = metrics_df['elevation_diff'].abs()

    print(f"\nBefore filtering:")
    print(f"  Total stations: {len(metrics_df)}")
    print(f"  Elevation difference: mean={metrics_df['elevation_diff'].mean():.2f}m, "
          f"std={metrics_df['elevation_diff'].std():.2f}m")
    print(f"  Absolute difference: mean={metrics_df['elevation_abs_diff'].mean():.2f}m, "
          f"max={metrics_df['elevation_abs_diff'].max():.2f}m")

    # Filter by elevation difference if specified
    if max_diff is not None:
        mask = metrics_df['elevation_abs_diff'] <= max_diff
        n_filtered = (~mask).sum()
        metrics_df = metrics_df[mask].copy()

        print(f"\nAfter filtering (max abs diff = {max_diff}m):")
        print(f"  Stations removed: {n_filtered}")
        print(f"  Stations remaining: {len(metrics_df)}")
        print(f"  Elevation difference: mean={metrics_df['elevation_diff'].mean():.2f}m, "
              f"std={metrics_df['elevation_diff'].std():.2f}m")
        print(f"  Absolute difference: mean={metrics_df['elevation_abs_diff'].mean():.2f}m, "
              f"max={metrics_df['elevation_abs_diff'].max():.2f}m")

    return metrics_df


# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def safe_polyfit(x, y, deg=1):
    """Safely fit polynomial, return None if fails."""
    try:
        return np.polyfit(x, y, deg)
    except (LinAlgError, np.linalg.LinAlgError, ValueError) as e:
        print(f"    Warning: Polyfit failed - {str(e)}")
        return None


def safe_correlation(df, col1, col2):
    """Safely calculate correlation, return NaN if fails."""
    try:
        return df[[col1, col2]].corr().iloc[0, 1]
    except Exception as e:
        print(f"    Warning: Correlation calculation failed - {str(e)}")
        return np.nan


def clean_data(df, columns):
    """Remove NaN and infinite values from specified columns."""
    return df[columns].replace([np.inf, -np.inf], np.nan).dropna()


# ============================================================================
# PLOTTING FUNCTIONS
# ============================================================================

def add_trend_line(ax, x, y, color='red'):
    """Add trend line to plot if possible."""
    z = safe_polyfit(x, y, 1)
    if z is not None:
        p = np.poly1d(z)
        x_trend = np.linspace(x.min(), x.max(), 100)
        ax.plot(x_trend, p(x_trend), "--", color=color, alpha=0.8, linewidth=2,
                label=f'Trend: y={z[0]:.2e}x+{z[1]:.2f}')
    return z is not None


def add_stats_box(ax, n, corr=None, position=(0.05, 0.95)):
    """Add statistics text box to plot."""
    textstr = f'n = {n}'
    if corr is not None and not np.isnan(corr):
        textstr += f'\nCorr = {corr:.3f}'

    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(position[0], position[1], textstr, transform=ax.transAxes,
            fontsize=10, verticalalignment='top', bbox=props)


def plot_metric_vs_elevation(df, metric, elevation_col, xlabel, title, color):
    """Create scatter plot of metric vs elevation."""
    fig, ax = plt.subplots(figsize=(10, 6))

    # Clean data
    df_clean = clean_data(df, [elevation_col, metric])

    if len(df_clean) < 3:
        print(f"    Warning: Not enough valid data points ({len(df_clean)}) for {metric} vs {elevation_col}")
        plt.close(fig)
        return None

    # Scatter plot
    ax.scatter(df_clean[elevation_col], df_clean[metric],
               c=color, s=80, alpha=0.6,
               edgecolors='black', linewidth=0.5)

    # Add trend line
    has_trend = add_trend_line(ax, df_clean[elevation_col], df_clean[metric])

    # Calculate correlation
    corr = safe_correlation(df_clean, elevation_col, metric)

    # Labels and title
    ax.set_xlabel(xlabel, fontsize=12, fontweight='bold')
    ylabel = metric.upper() + (' (mm)' if metric != 'r2' else '')
    ax.set_ylabel(ylabel, fontsize=12, fontweight='bold')

    title_str = title
    if not np.isnan(corr):
        title_str += f' (r={corr:.3f})'
    ax.set_title(title_str, fontsize=14, fontweight='bold')

    # Formatting
    ax.grid(alpha=0.3, linestyle='--')
    if has_trend:
        ax.legend(fontsize=10)

    add_stats_box(ax, len(df_clean), corr)

    plt.tight_layout()
    return fig


def plot_elevation_difference(df, max_diff=None):
    """Plot histogram and scatter of elevation differences."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Clean data
    df_clean = clean_data(df, ['obs_elevation', 'wrf_elevation', 'elevation_diff'])

    # Histogram
    axes[0].hist(df_clean['elevation_diff'], bins=30, color='steelblue',
                 alpha=0.7, edgecolor='black')
    axes[0].axvline(0, color='red', linestyle='--', linewidth=2, label='Zero difference')

    # Add filter threshold lines if applicable
    if max_diff is not None:
        axes[0].axvline(max_diff, color='orange', linestyle=':', linewidth=2,
                       label=f'Filter threshold: ±{max_diff}m')
        axes[0].axvline(-max_diff, color='orange', linestyle=':', linewidth=2)

    axes[0].set_xlabel('Elevation Difference (WRF - Obs) [m]',
                       fontsize=11, fontweight='bold')
    axes[0].set_ylabel('Frequency', fontsize=11, fontweight='bold')

    title_str = 'Distribution of Elevation Differences'
    if max_diff is not None:
        title_str += f' (filtered: |diff| ≤ {max_diff}m)'
    axes[0].set_title(title_str, fontsize=12, fontweight='bold')
    axes[0].legend()
    axes[0].grid(alpha=0.3)

    # Statistics
    stats_text = (f'n = {len(df_clean)}\n'
                  f'Mean = {df_clean["elevation_diff"].mean():.2f} m\n'
                  f'Std = {df_clean["elevation_diff"].std():.2f} m\n'
                  f'Median = {df_clean["elevation_diff"].median():.2f} m')
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    axes[0].text(0.72, 0.95, stats_text, transform=axes[0].transAxes,
                fontsize=10, verticalalignment='top', bbox=props)

    # Scatter: WRF vs Obs elevation
    axes[1].scatter(df_clean['obs_elevation'], df_clean['wrf_elevation'],
                   c='steelblue', s=80, alpha=0.6,
                   edgecolors='black', linewidth=0.5)

    # 1:1 line
    min_elev = min(df_clean['obs_elevation'].min(), df_clean['wrf_elevation'].min())
    max_elev = max(df_clean['obs_elevation'].max(), df_clean['wrf_elevation'].max())
    axes[1].plot([min_elev, max_elev], [min_elev, max_elev],
                'r--', linewidth=2, label='1:1 line')

    # Add tolerance band if filter is active
    if max_diff is not None:
        x_range = np.linspace(min_elev, max_elev, 100)
        axes[1].fill_between(x_range, x_range - max_diff, x_range + max_diff,
                            alpha=0.2, color='orange', label=f'±{max_diff}m band')

    # Correlation
    corr = safe_correlation(df_clean, 'obs_elevation', 'wrf_elevation')

    axes[1].set_xlabel('Observed Elevation (m)', fontsize=11, fontweight='bold')
    axes[1].set_ylabel('WRF Elevation (m)', fontsize=11, fontweight='bold')

    title_str = 'WRF vs Observed Elevation'
    if not np.isnan(corr):
        title_str += f' (r={corr:.3f})'
    axes[1].set_title(title_str, fontsize=12, fontweight='bold')

    axes[1].legend()
    axes[1].grid(alpha=0.3)

    plt.tight_layout()
    return fig


def plot_combined_comparison(df, metrics_config):
    """Create combined plot with all metrics vs all elevation types."""
    fig, axes = plt.subplots(3, 3, figsize=(18, 15))

    elevation_cols = ['obs_elevation', 'wrf_elevation', 'elevation_diff']
    elevation_labels = ['Observed Elevation (m)', 'WRF Elevation (m)', 'Elevation Diff (WRF - Obs) [m]']

    for i, (metric, metric_name, color) in enumerate(metrics_config):
        # Clean data once per metric
        cols_needed = [metric] + elevation_cols
        df_clean = clean_data(df, cols_needed)

        for j, (elev_col, elev_label) in enumerate(zip(elevation_cols, elevation_labels)):
            axes[i, j].scatter(df_clean[elev_col], df_clean[metric],
                              c=color, s=60, alpha=0.6,
                              edgecolors='black', linewidth=0.5)

            corr = safe_correlation(df_clean, elev_col, metric)

            axes[i, j].set_xlabel(elev_label, fontsize=10)
            axes[i, j].set_ylabel(metric.upper(), fontsize=10)

            title_str = f'{metric_name} vs {elev_label.split("(")[0].strip()}'
            if not np.isnan(corr):
                title_str += f' (r={corr:.3f})'
            axes[i, j].set_title(title_str, fontsize=11)
            axes[i, j].grid(alpha=0.3)

    plt.tight_layout()
    return fig


# ============================================================================
# PLOTTING WORKFLOW
# ============================================================================

def save_and_show(fig, filename, output_dir):
    """Save figure and display it."""
    if fig is not None:
        filepath = output_dir / filename
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        plt.show()
        print(f"  Saved {filename}")
        return True
    return False


def create_all_plots(metrics_df, output_dir, max_diff=None):
    """Create and save all plots."""
    metrics_config = [
        ('rmse', 'RMSE', 'steelblue'),
        ('r2', 'R²', 'forestgreen'),
        ('mae', 'MAE', 'coral')
    ]

    elevation_types = [
        ('wrf_elevation', 'WRF Elevation (m)', 'WRF Elevation'),
        ('obs_elevation', 'Observed Elevation (m)', 'Observed Elevation'),
        ('elevation_diff', 'Elevation Difference (WRF - Obs) [m]', 'Elevation Difference')
    ]

    # 1. Elevation difference plots
    print("\nGenerating elevation difference plots...")
    try:
        fig = plot_elevation_difference(metrics_df, max_diff)
        save_and_show(fig, 'elevation_differences.png', output_dir)
    except Exception as e:
        print(f"  Error: {str(e)}")

    # 2. Metrics vs elevation plots
    for elev_col, elev_label, elev_name in elevation_types:
        print(f"\nGenerating metrics vs {elev_name.lower()} plots...")
        for metric, metric_name, color in metrics_config:
            try:
                fig = plot_metric_vs_elevation(
                    metrics_df, metric, elev_col,
                    elev_label, f'{metric_name} vs {elev_name}', color
                )
                filename = f'{metric}_vs_{elev_col}.png'
                save_and_show(fig, filename, output_dir)
            except Exception as e:
                print(f"  Error creating {metric} vs {elev_name} plot: {str(e)}")

    # 3. Combined comparison plot
    print("\nGenerating combined comparison plot...")
    try:
        fig = plot_combined_comparison(metrics_df, metrics_config)
        save_and_show(fig, 'all_elevation_comparisons.png', output_dir)
    except Exception as e:
        print(f"  Error: {str(e)}")


def print_correlation_summary(metrics_df):
    """Print correlation statistics."""
    print("\n=== Elevation Correlation Summary ===")

    try:
        print("\nCorrelations with Observed Elevation:")
        print(metrics_df[['obs_elevation', 'rmse', 'r2', 'mae']].corr()['obs_elevation'])

        print("\nCorrelations with WRF Elevation:")
        print(metrics_df[['wrf_elevation', 'rmse', 'r2', 'mae']].corr()['wrf_elevation'])

        print("\nCorrelations with Elevation Difference:")
        print(metrics_df[['elevation_diff', 'rmse', 'r2', 'mae']].corr()['elevation_diff'])
    except Exception as e:
        print(f"Error calculating correlations: {str(e)}")


# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    """Main execution function."""
    # Load data
    metrics_df, obs_ds, wrf_ds = load_datasets(METRICS_FILE, OBS_FILE, WRF_FILE)

    # Process elevation data with filtering
    metrics_df = prepare_elevation_data(metrics_df, obs_ds, wrf_ds, max_diff=MAX_ELEVATION_DIFF)

    # Create all plots
    create_all_plots(metrics_df, OUTPUT_DIR, max_diff=MAX_ELEVATION_DIFF)

    # Save updated metrics
    try:
        output_file = OUTPUT_DIR / 'metrics_with_elevations.csv'
        metrics_df.to_csv(output_file, index=False)
        print(f"\nSaved metrics_with_elevations.csv")
    except Exception as e:
        print(f"Error saving CSV: {str(e)}")

    # Print summary
    print_correlation_summary(metrics_df)


if __name__ == "__main__":
    main()
