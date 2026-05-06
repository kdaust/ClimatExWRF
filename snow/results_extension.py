"""
Analyze and visualize WRF vs ERA5 climate metrics across weather stations.

This script compares Weather Research and Forecasting (WRF) model outputs against
ERA5 reanalysis data by analyzing performance metrics at multiple weather stations.
It filters stations based on data coverage thresholds, generates boxplots for metric
distribution, and creates geographic maps showing spatial patterns of model performance.

Configuration:
    FILE1 (str): Path to CSV file containing model performance metrics summary
    FILE2 (str): Path to CSV file containing station metadata and coverage statistics
    OUTPUT_DIR (str): Directory path where output visualizations will be saved
    MIN_COVERAGE (float): Minimum seasonal data coverage threshold for station inclusion (0.0-1.0)
    METRICS (list): List of performance metric names to analyze and visualize

Output:
    - boxplots.png: Box and swarm plots showing distribution of all metrics
    - map_*.png: Geographic distribution maps for each metric
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from plots_WRF_vs_OBS import plot_metric_map

# Configuration
FILE1 = '/comparison_era5/metrics_summary.csv'
FILE2 = '/climatex/stations_summary.csv'
OUTPUT_DIR = '/climatex/comparison_era5'
MIN_COVERAGE = 0.94
METRICS = ['rmse', 'r2', 'mae', 'MAPE Max', 'MAPE Snow Duration', 'Bias Daily']


# Create output directory
Path(OUTPUT_DIR).mkdir(exist_ok=True)

# Load and filter data
file1 = pd.read_csv(FILE1)
file2 = pd.read_csv(FILE2)

filtered_stations = file2[file2['seasonal_coverage'] > MIN_COVERAGE]['station_id']
results_df = file1[file1['station_id'].isin(filtered_stations)]
print(filtered_stations)

results_df = results_df.loc[~np.isinf(results_df['MAPE Snow Duration'])]

print(f"Filtered: {len(results_df)} / {len(file1)} stations")
print(results_df[['station_id','MAPE Max', 'MAPE Snow Duration']])
fig, axes = plt.subplots(1, len(METRICS), figsize=(5*len(METRICS), 5))
for ax, metric in zip(axes, METRICS):
    # Create boxplot
    sns.boxplot(data=results_df, y=metric, ax=ax, width=0.5)
    # Add swarm plot on top
    sns.swarmplot(data=results_df, y=metric, ax=ax, color='black', alpha=0.5, size=6)
    ax.set_title(metric)

plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/boxplots.png', dpi=300, bbox_inches='tight')
plt.close()


for metric in METRICS:
    fig = plot_metric_map(results_df, metric, f'{metric.upper()} Distribution')
    plt.savefig(f'{OUTPUT_DIR}/map_{metric.lower().replace(" ", "_")}.png',
                dpi=300, bbox_inches='tight')
    plt.close()

print(f"Results saved to {OUTPUT_DIR}/")
