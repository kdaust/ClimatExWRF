"""
Skill Score Comparison Analysis for Model Performance Evaluation

This script compares the performance of two models (CaSR and ERA5) against a benchmark
using skill score metrics. Skill scores quantify the relative improvement or degradation
of model predictions compared to a baseline benchmark across multiple error metrics.

Methodology:
    - For error metrics (RMSE, MAE): SS = 1 - (error_model / error_benchmark)
      * SS > 0 indicates better performance than benchmark
      * SS < 0 indicates worse performance than benchmark

    - For R² (coefficient of determination): SS = (R²_model - R²_benchmark) / (1 - R²_benchmark)
      * SS > 0 indicates better predictive skill than benchmark
      * Scaled by benchmark uncertainty (1 - R²_benchmark)

Workflow:
    1. Load three CSV files containing metrics (benchmark, CaSR model, ERA5 model)
    2. Remove invalid data (rows with infinite MAPE Snow Duration values)
    3. Calculate skill scores for RMSE, MAE, and R² metrics
    4. Generate boxplot visualizations comparing model performance
    5. Print detailed summary statistics for each metric and model

Outputs:
    - Boxplot comparison figure saved as 'skill_score_comparison.png'
    - Console output with summary statistics including mean, median, std dev, min/max,
      and percentage of cases where each model outperforms the benchmark

Returns:
    - ss_data (dict): Nested dictionary containing skill score Series for each metric
      and model combination
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load the three CSV files
file1 = pd.read_csv('climatex/comparison_wrf/metrics_summary.csv')  # Benchmark
file2 = pd.read_csv('climatex/comparison_casr/metrics_summary.csv')
file3 = pd.read_csv('climatex/comparison_era5/metrics_summary.csv')

# Calculate Skill Scores (SS) for each error metric
# SS = 1 - (error_model / error_benchmark)
# For R2, since higher is better, we use: SS = (R2_model - R2_benchmark) / (1 - R2_benchmark)

def calculate_skill_scores(benchmark, model, metric):
    """Calculate skill score for a given metric"""
    if metric == 'r2':
        # For R2: positive SS means better than benchmark
        ss = (model[metric] - benchmark[metric]) / (1 - benchmark[metric])
    else:  # rmse and mae: lower is better
        ss = 1 - (model[metric] / benchmark[metric])
    return ss

# Calculate skill scores for all metrics
metrics = ['rmse', 'mae', 'r2']
#metrics = ['MAPE Max','MAPE Snow Duration','Bias Daily']
ss_data = {}


file2 = file2.loc[~np.isinf(file1['MAPE Snow Duration'])]
file3 = file3.loc[~np.isinf(file1['MAPE Snow Duration'])]
file1 = file1.loc[~np.isinf(file1['MAPE Snow Duration'])]

for metric in metrics:
    ss_data[metric] = {
        'CaSR': calculate_skill_scores(file1, file2, metric),
        'ERA5': calculate_skill_scores(file1, file3, metric)
    }

# Create boxplots
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

for idx, metric in enumerate(metrics):
    data_to_plot = [ss_data[metric]['CaSR'], ss_data[metric]['ERA5']]

    bp = axes[idx].boxplot(data_to_plot, labels=['CaSR', 'ERA5'], patch_artist=True)

    # Color the boxes
    colors = ['lightblue', 'lightcoral']
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)

    axes[idx].set_ylabel('Skill Score', fontsize=11)
    axes[idx].set_title(f'Skill Score for {metric.upper()}', fontsize=12, fontweight='bold')
    axes[idx].grid(axis='y', alpha=0.3)
    axes[idx].axhline(y=0, color='k', linestyle='--', linewidth=1, alpha=0.5)
    axes[idx].set_ylim(-3, 3)

plt.tight_layout()
plt.savefig('/climatex/skill_score_comparison.png', dpi=300, bbox_inches='tight')
plt.show()

# Print summary statistics
print("=" * 70)
print("SKILL SCORE SUMMARY STATISTICS")
print("=" * 70)
for metric in metrics:
    print(f"\n{metric.upper()}:")
    print("-" * 70)
    for model in ['CaSR', 'ERA5']:
        ss = ss_data[metric][model]
        print(f"\n{model}:")
        print(f"  Mean SS:     {ss.mean():.4f}")
        print(f"  Median SS:   {ss.median():.4f}")
        print(f"  Std Dev:     {ss.std():.4f}")
        print(f"  Min SS:      {ss.min():.4f}")
        print(f"  Max SS:      {ss.max():.4f}")
        positive = (ss > 0).sum()
        print(f"  % Better than benchmark: {(positive/len(ss)*100):.1f}%")
