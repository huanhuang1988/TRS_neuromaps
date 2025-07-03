#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script computes Spearman correlations between parcellated receptor data and mean values from an
independent component map, using the Schaefer 400-parcel and Tian 32-parcel atlas. It performs spin tests
to assess significance, applies FDR corrections, and generates permutation test plots for each receptor.
Inputs include receptor data and mean values from a previous step. Outputs include correlation results and plots.
@author: Huan Huang
email: lexie_hh@163.com
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, zscore
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import neuromaps
import neuromaps.nulls
from netneurotools import datasets, utils, plotting

# Define constants
NNODES = 432  # Total number of nodes (cortical + subcortical)
NCORTICALNODES = 400  # Number of cortical nodes
NSUBCORTNODES = 32  # Number of subcortical nodes
N_PERM = 10000  # Number of permutations for spin tests
ATLAS_NAME = "schaefer400parcel7networks_subcor32"  # Atlas name for file paths

# Define relative file paths
output_dir = os.path.join("results", "receptor_correlations")  # Directory for output files
receptor_data_file = os.path.join("data", f"{ATLAS_NAME}_receptor_data.csv")  # Parcellated receptor data
receptor_names_file = os.path.join("data", "receptor_names_pet.npy")  # Receptor names
mean_values_file = os.path.join("results", f"ic_{ATLAS_NAME}_values.csv")  # Mean values from step 1
zscore_output_file = os.path.join("data", f"{ATLAS_NAME}_receptor_data_zscore.csv")  # Z-scored receptor data
results_file = os.path.join(output_dir, "receptor_ic_correlations_fdr05.csv")  # Correlation results

# Create output directory
os.makedirs(output_dir, exist_ok=True)

# Load Schaefer atlas for cortical parcellation
schaefer_fs = datasets.fetch_schaefer2018(version='fsaverage')['400Parcels7Networks']
schaefer_fs = neuromaps.images.annot_to_gifti(schaefer_fs)

# Load receptor data
parcellated_receptors_data_array = np.genfromtxt(receptor_data_file, delimiter=',')
receptor_names = np.load(receptor_names_file)

# Z-score receptor data
parcellated_receptors_data_array_zscore = zscore(parcellated_receptors_data_array, axis=0)

# Save z-scored receptor data
np.savetxt(zscore_output_file, parcellated_receptors_data_array_zscore, delimiter=',')
print(f"Z-scored receptor data saved to: {zscore_output_file}")

# Load mean values from fMRI component
data = pd.read_csv(mean_values_file)

neuroimage_data = data["Mean_Value"].values

# Initialize results dictionary
rotated_results = {}

# Compute correlations and spin tests for each receptor
for idx, receptor_name in enumerate(receptor_names):
    receptor_vector = parcellated_receptors_data_array_zscore[:, idx]
    
    # Compute Spearman correlation
    observed_corr, observed_pval = spearmanr(neuroimage_data, receptor_vector)
    
    # Perform spin test for cortical nodes
    rotated_cort = neuromaps.nulls.vasa(
        neuroimage_data[:NCORTICALNODES],
        atlas='fsaverage',
        density='164k',
        parcellation=schaefer_fs,
        n_perm=N_PERM,
        seed=1234
    )
    
    # Random permutation for subcortical nodes
    rotated_sub = np.zeros((NSUBCORTNODES, N_PERM))
    np.random.seed(1234)
    for p in range(N_PERM):
        shuffled_indices = np.random.permutation(np.arange(NCORTICALNODES, NNODES))
        rotated_sub[:, p] = neuroimage_data[shuffled_indices[:NSUBCORTNODES]]
    
    # Combine cortical and subcortical rotations
    rotated = np.vstack((rotated_cort, rotated_sub))
    null_corrs = np.array([spearmanr(rotated[:, p], receptor_vector)[0] for p in range(N_PERM)])
    spin_pval = (np.sum(np.abs(null_corrs) >= np.abs(observed_corr)) + 1) / (N_PERM + 1)
    
    print(f'{receptor_name}: r = {observed_corr:.3f}, p_spin = {spin_pval:.5f}')
    
    # Store results
    rotated_results[receptor_name] = {
        'observed_corr': observed_corr,
        'original_pval': observed_pval,
        'spin_pval': spin_pval,
        'null_distributions': null_corrs
    }
    
    # Create permutation test plot
    plt.figure(figsize=(8, 5))
    plt.hist(null_corrs, bins=50, color='lightgray', edgecolor='black', alpha=0.7)
    plt.axvline(observed_corr, color='red', linestyle='--', linewidth=2,
                label=f'Observed r = {observed_corr:.3f}, p_spin = {spin_pval:.5f}')
    plt.title(f'Permutation Test for {receptor_name}', fontsize=16)
    plt.xlabel('Spearman Correlation', fontsize=14)
    plt.ylabel('Frequency', fontsize=14)
    plt.legend(fontsize=12)
    plt.tight_layout()
    plot_filename = os.path.join(output_dir, f'{receptor_name}_permutation_test_plot.png')
    plt.savefig(plot_filename, dpi=300)
    plt.close()
    
    print(f'Plot saved for {receptor_name}: {plot_filename}')

# Create results table
results_list = []
for receptor_name, result in rotated_results.items():
    results_list.append([
        receptor_name,
        result['observed_corr'],
        result['original_pval'],
        result['spin_pval']
    ])

results_df = pd.DataFrame(results_list, columns=['Receptor', 'Observed_Corr', 'Original_P_Value', 'Spin_P_Value'])

# Apply FDR and FWE corrections
spin_p_values = results_df['Spin_P_Value'].dropna().values
if len(spin_p_values) > 0:
    fdr_corrected = multipletests(spin_p_values, alpha=0.05, method='fdr_bh')[1]
    fwe_corrected = multipletests(spin_p_values, alpha=0.05, method='bonferroni')[1]
    results_df['FDR_Corrected_Spin_P'] = np.nan
    results_df['FWE_Corrected_Spin_P'] = np.nan
    results_df.loc[results_df['Spin_P_Value'].notna(), 'FDR_Corrected_Spin_P'] = fdr_corrected
    results_df.loc[results_df['Spin_P_Value'].notna(), 'FWE_Corrected_Spin_P'] = fwe_corrected
else:
    results_df['FDR_Corrected_Spin_P'] = np.nan
    results_df['FWE_Corrected_Spin_P'] = np.nan

# Save correlation results
results_df.to_csv(results_file, index=False)
print(f"Correlation results saved to: {results_file}")