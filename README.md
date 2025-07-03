This repository contains Python scripts for processing independent component (IC) maps using the Schaefer 400-parcel and Tian 32-parcel atlas, computing regional and network-level statistics, and correlating these with receptor data. 

**Overview**
The workflow consists of four scripts that process an IC map (ic.nii) through various analyses:

step1_extract_roi_mean_values_for_ic.py: Extracts mean values for each atlas region from the IC map.
step2_compute_roi_stats_for_ic.py: Computes regional statistics with a z > 3 threshold.
step3_compute_and_plot_roi_network_stats_for_ic.py: Aggregates regional statistics into network-level statistics and generates a bar plot of voxel contributions.
step4_compute_spatial_correlations_between_ic_and_receptors.py: Computes Spearman correlations between IC mean values and receptor data, performs spin tests, and generates permutation plots.

**Usage**
Run the scripts in the following order, as each depends on the output of the previous step:

python step1_extract_roi_mean_values_for_ic.py
python step2_compute_roi_stats_for_ic.py
python step3_compute_and_plot_roi_network_stats_for_ic.py
python step4_compute_spatial_correlations_between_ic_and_receptors.py

**Directory Structure**

Inputs (data/):
atlas-schaefer400parcel7networks_subcor32_1mm.csv: Atlas labels with ROILabel and ROIName columns.
atlas-schaefer400parcel7networks_subcor32_1mm.nii.gz: Schaefer 400-parcel and Tian 32-parcel atlas (NIfTI format).
ic.nii: fMRI independent component map derived from SBM analysis.
receptor_names_pet.npy: Array of receptor names.
schaefer400parcel7networks_subcor32_receptor_data.csv: Parcellated receptor data (432 nodes).

Outputs (results/):
ic_schaefer400_tian32_values.csv: Regional mean values from step 1.
ic_schaefer400_tian32_stats.csv: Regional statistics with z > 3 from step 2.
ic_schaefer400_tian32_network_stats.csv: Network-level statistics from step 3.
ic_schaefer400_tian32_network_plot.png: Bar plot of network voxel contributions from step 3.

receptor_correlations/:
schaefer400parcel7networks_subcor32_receptor_data_zscore.csv: Z-scored receptor data from step 4.
receptor_ic_correlations_fdr05.csv: Correlation results with FDR/FWE corrections from step 4.
<receptor_name>_permutation_test_plot.png: Per-receptor permutation test plots from step 4.

Input Requirements
data/ic.nii: NIfTI file containing the fMRI independent component map.
data/atlas-schaefer400parcel7networks_subcor32_1mm.nii.gz: Atlas file with 432 regions (400 cortical + 32 subcortical).
data/atlas-schaefer400parcel7networks_subcor32_1mm.csv: CSV file with columns ROILabel and ROIName (e.g., 7Networks_LH_Vis_1).
data/schaefer400parcel7networks_subcor32_receptor_data.csv: CSV file with 432 rows of receptor data.
data/receptor_names_pet.npy: NumPy array of receptor names.
Note: The atlas and fMRI data must have matching dimensions (432 nodes for receptor and mean value data).

Output Details
results/ic_schaefer400_tian32_values.csv: Contains Label_ID, Region_Name, and Mean_Value columns.
results/ic_schaefer400_tian32_stats.csv: Contains Label_ID, Region_Name, Hemisphere, Network, Num_Voxels, Total_Voxels, Percentage_NonZero_Voxels, and Average_Signal columns.
results/ic_schaefer400_tian32_network_stats.csv: Contains Network, Total_Num_Voxels, Voxel_Percentage, and Weighted_Average_Signal columns.
results/ic_schaefer400_tian32_network_plot.png: Horizontal bar plot of voxel contributions by network (ordered: Vis, SomMot, DorsAttn, SalVentAttn, Limbic, Cont, Default, SubCor).
results/receptor_correlations/schaefer400parcel7networks_subcor32_receptor_data_zscore.csv: Z-scored receptor data (432 nodes).
results/receptor_correlations/receptor_ic_correlations_fdr05.csv: Contains Receptor, Observed_Corr, Original_P_Value, Spin_P_Value, FDR_Corrected_Spin_P, and FWE_Corrected_Spin_P columns.
results/receptor_correlations/<receptor_name>_permutation_test_plot.png: Histogram of null correlations with observed correlation line for each receptor.
