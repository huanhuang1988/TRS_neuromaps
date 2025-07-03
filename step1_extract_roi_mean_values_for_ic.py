# -*- coding: utf-8 -*-
"""

This script extracts mean values for each region in the Schaefer 400-parcel and Tian 32-parcel atlas 
from an independent component based on SBM analysis and saves them to a CSV file, without applying a threshold.
@author: Huan Huang
email：lexie_hh@163.com

"""

import os
import nibabel as nib
import numpy as np
import pandas as pd

# Define constants
ATLAS_FILE = "atlas-schaefer400parcel7networks_subcor32_1mm.nii"  # Schaefer 400-parcel and Tian 32-parcel atlas
LABEL_FILE = "atlas-schaefer400parcel7networks_subcor32_1mm.csv"   # Atlas labels CSV

# Define relative file paths
nii_file = os.path.join("data", "ic.nii")  # Input an independent component map derived from SBM analysis
atlas_file = os.path.join("data", ATLAS_FILE)  # Atlas file
label_file = os.path.join("data", LABEL_FILE)  # Atlas labels CSV
output_csv = os.path.join("results", "ic_schaefer400_tian32_values.csv")  # Output CSV

# Create output directory if it doesn't exist
os.makedirs("results", exist_ok=True)

# Load NIfTI files
nii_img = nib.load(nii_file)
atlas_img = nib.load(atlas_file)

# Load volume data
nii_data = nii_img.get_fdata()
atlas_data = atlas_img.get_fdata()

# Round atlas values to nearest integer
atlas_data = np.round(atlas_data)

# Load atlas labels
atlas_labels = pd.read_csv(label_file)

# Initialize results table
results = pd.DataFrame(columns=["Label_ID", "Region_Name", "Mean_Value"])

# Extract unique labels (excluding 0)
roi_labels = np.unique(atlas_data[atlas_data > 0])

# Process each ROI
for label_id in roi_labels:
    # Create region mask
    region_mask = atlas_data == label_id
    region_values = nii_data[region_mask]
    
    # Calculate mean value, ignoring NaNs
    mean_value = np.nanmean(region_values)
    
    # Get region name
    region_name = atlas_labels[atlas_labels["ROILabel"] == label_id]["ROIName"].iloc[0]
    
    # Append to results
    new_row = pd.DataFrame([{
        "Label_ID": label_id,
        "Region_Name": region_name,
        "Mean_Value": mean_value
    }])
    results = pd.concat([results, new_row], ignore_index=True)

# Save results to CSV
results.to_csv(output_csv, index=False)
print(f"Results saved to: {output_csv}")