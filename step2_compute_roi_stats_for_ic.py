# -*- coding: utf-8 -*-
"""
This script processes an independent component map derived from SBM analysis, applies a z > 3 threshold,
and computes regional statistics (number of significant voxels, percentage, average signal, hemisphere, and network)
using the Schaefer 400-parcel and Tian 32-parcel atlas. Hemisphere and Network are extracted from Region_Name.
Results are saved to a CSV file.
@author: Huan Huang
email: lexie_hh@163.com
"""

import os
import numpy as np
import pandas as pd
import nibabel as nib

# Define constants
ATLAS_FILE = "atlas-schaefer400parcel7networks_subcor32_1mm.nii"  # Schaefer 400-parcel and Tian 32-parcel atlas
LABEL_FILE = "atlas-schaefer400parcel7networks_subcor32_1mm.csv"   # Atlas labels CSV
THRESHOLD = 3.0  # Z-score threshold for significant voxels

# Define relative file paths
nii_file = os.path.join("data", "ic.nii")  # Input an independent component map derived from SBM analysis
atlas_file = os.path.join("data", ATLAS_FILE)  # Atlas file
label_file = os.path.join("data", LABEL_FILE)  # Atlas labels CSV
output_csv = os.path.join("results", "ic_schaefer400_tian32_stats.csv")  # Output CSV with z > 3 statistics

# Create output directory if it doesn't exist
os.makedirs("results", exist_ok=True)

# Load NIfTI volumes
nii_img = nib.load(nii_file)
atlas_img = nib.load(atlas_file)

# Load volume data
nii_data = nii_img.get_fdata()
atlas_data = atlas_img.get_fdata()

# Load atlas labels
atlas_labels = pd.read_csv(label_file)

# Initialize results table with additional Hemisphere and Network columns
results = pd.DataFrame(columns=["Label_ID", "Region_Name", "Hemisphere", "Network", "Num_Voxels",
                               "Total_Voxels", "Percentage_NonZero_Voxels", "Average_Signal"])

# Get unique labels from atlas (excluding background)
template_labels = np.unique(atlas_data[atlas_data > 0])

# Count total voxels per label
total_voxels_per_label = {lbl: np.sum(atlas_data == lbl) for lbl in template_labels}

# Process each label
for label_id in template_labels:
    # Create region mask
    region_mask = atlas_data == label_id
    region_voxels = nii_data[region_mask]
    
    # Apply z-score threshold
    non_zero_voxels = region_voxels[region_voxels > THRESHOLD]
    non_zero_count = len(non_zero_voxels)
    total_voxels = total_voxels_per_label[label_id]
    
    # Calculate statistics
    percentage_nonzero = non_zero_count / total_voxels if total_voxels > 0 else 0
    average_signal = np.mean(non_zero_voxels) if non_zero_count > 0 else np.nan
    
    # Get region name
    region_name = atlas_labels[atlas_labels["ROILabel"] == label_id]["ROIName"].iloc[0]
    
    # Extract Hemisphere and Network from Region_Name (e.g., 7Networks_LH_Vis_1)
    name_parts = region_name.split("_")
    if len(name_parts) >= 3 and name_parts[0] == "7Networks":
        hemisphere = name_parts[1]  # LH or RH
        network = name_parts[2]     # e.g., Vis, SomMot, etc.
    else:
        hemisphere = "Unknown"
        network = "Unknown"
    
    # Append results
    new_row = pd.DataFrame([{
        "Label_ID": label_id,
        "Region_Name": region_name,
        "Hemisphere": hemisphere,
        "Network": network,
        "Num_Voxels": non_zero_count,
        "Total_Voxels": total_voxels,
        "Percentage_NonZero_Voxels": percentage_nonzero,
        "Average_Signal": average_signal
    }])
    results = pd.concat([results, new_row], ignore_index=True)

# Save results to CSV
results.to_csv(output_csv, index=False)
print(f"Results saved to: {output_csv}")