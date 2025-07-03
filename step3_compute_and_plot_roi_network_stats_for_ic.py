# -*- coding: utf-8 -*-
"""
This script reads regional statistics from an independent component map (z > 3 threshold),
computes network-level statistics (total voxel count, percentage, and weighted average signal),
and creates a horizontal bar plot of voxel contributions by network using the Schaefer 400-parcel
and Tian 32-parcel atlas. It processes regions with significant voxels (Num_Voxels > 0).
@author: Huan Huang
email: lexie_hh@163.com
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Define relative file paths
data_csv = os.path.join("results", "ic_schaefer400_tian32_stats.csv")  # Input CSV with regional statistics
output_csv = os.path.join("results", "ic_schaefer400_tian32_network_stats.csv")  # Output CSV with network statistics
output_figure_path = os.path.join("results", "ic_schaefer400_tian32_network_plot.png")  # Output plot

# Create output directory if it doesn't exist
os.makedirs("results", exist_ok=True)

# Load the data from CSV file
data = pd.read_csv(data_csv)

# Filter out rows where Num_Voxels is zero
data = data[data["Num_Voxels"] > 0]

# Get unique network names
network_names = data["Network"].unique()
num_networks = len(network_names)

# Initialize arrays for results
network_voxel_count = np.zeros(num_networks)
network_voxel_percentage = np.zeros(num_networks)
network_weighted_avg_signal = np.zeros(num_networks)

# Calculate total voxel count across all regions
total_num_voxels = data["Num_Voxels"].sum()

# Process each network
for i, current_network in enumerate(network_names):
    # Filter rows for the current network
    network_data = data[data["Network"] == current_network]
    
    # Calculate total voxel count for this network
    network_voxel_count[i] = network_data["Num_Voxels"].sum()
    
    # Calculate percentage of total Num_Voxels for this network
    network_voxel_percentage[i] = (network_voxel_count[i] / total_num_voxels) * 100
    
    # Calculate weighted average of Average_Signal, weighted by Num_Voxels
    network_weighted_avg_signal[i] = (network_data["Average_Signal"] * network_data["Num_Voxels"]).sum() / network_voxel_count[i]

# Create results table
results = pd.DataFrame({
    "Network": network_names,
    "Total_Num_Voxels": network_voxel_count,
    "Voxel_Percentage": network_voxel_percentage,
    "Weighted_Average_Signal": network_weighted_avg_signal
})

# Save results to CSV
results.to_csv(output_csv, index=False)
print(f"Results saved to: {output_csv}")

# Specify the order of the networks for plotting
network_order = ['Vis', 'SomMot', 'DorsAttn', 'SalVentAttn', 'Limbic', 'Cont', 'Default', 'SubCor']

# Ensure the data is in the specified order
results['Network'] = pd.Categorical(results['Network'], categories=network_order, ordered=True)
results = results.sort_values('Network')

# Filter out networks not in network_order or with NaN values
results = results[results['Network'].notna()]

# Create a horizontal bar plot
plt.figure(figsize=(10, 6))
bar_colors = '#90c3d4'  # Lighter blue color

# Plot each network bar
barh = plt.barh(results['Network'], results['Total_Num_Voxels'], color=bar_colors, height=0.8)

# Add Voxel_Percentage annotations outside the bars
for i, (x, y) in enumerate(zip(results['Total_Num_Voxels'], results['Voxel_Percentage'])):
    plt.text(x + max(results['Total_Num_Voxels']) * 0.02, i, f'{y:.2f}%', 
             va='center', ha='left', fontsize=14, color='black', weight='bold')

# Customize plot appearance
plt.xlabel('Number of Voxels (x$10^4$)', fontsize=20, labelpad=10)
plt.suptitle('Functional Network Contributions (z>3)', fontsize=22, y=0.95, weight='bold')
plt.xticks(ticks=[0, 0.2e4, 0.4e4, 0.6e4, 0.8e4, 1.0e4, 1.2e4, 1.4e4, 1.6e4],
           labels=['0', '0.2', '0.4', '0.6', '0.8', '1.0', '1.2', '1.4', '1.6'], fontsize=16)
plt.yticks(fontsize=18)
plt.grid(axis='x', linestyle='--', alpha=0.5)
plt.gca().invert_yaxis()  # Match MATLAB's default behavior
plt.xlim([0, 1.6e4])  # Set x-axis maximum to 1.6 × 10^4

# Adjust plot layout
plt.tight_layout(pad=2)

# Save the figure
plt.savefig(output_figure_path, dpi=300)
plt.close()

print(f"Figure saved to: {output_figure_path}")