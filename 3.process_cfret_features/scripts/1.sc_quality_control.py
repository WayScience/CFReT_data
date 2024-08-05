#!/usr/bin/env python
# coding: utf-8

# # Perform single-cell quality control
# 
# In this notebook, we perform single-cell quality control using coSMicQC. To filter the single-cells, we use z-score to find outliers using the values from only one feature at a time. We use features from the AreaShape and Intensity modules to assess the quality of the segmented single-cells:
# 
# ### Assessing poor nuclei segmentation
# 
# Due to high confluence, sometimes nuclei overlap on top of each other, creating highly intense clusters within the Hoechst channel. To identify these nuclei, we use:
# 
# - **Nuclei Area:** This metric quantifies the number of pixels in a nucleus segmentation. 
# We detect nuclei that are abnormally large, which likely indicates poor nucleus segmentation where overlapping nuclei are merged into one segmentation. 
# - **Nuclei Intensity:** This metric quantifies the total intensity of all pixels in a nucleus segmentation. 
# In combination with abnormally large nuclei, we detect nuclei that are also highly intense, likely indicating that this a group of overlapped nuclei.
# 
# As well, there are times where nuclei are under or over-segmented, leading to non-circular nuclei shapes. To identify these nuclei, we use:
# 
# - **Nuclei FormFactor:** This metric quantifies how circular an object is, with 1 meaning perfect circle and 0 meaning non-circular. 
# We are detecting nuclei that are not circular and have rough edges or shapes that look like budding yeast.
# 
# ### Assessing poor cell segmentation
# 
# Also due to high confluence, images with large, intense clusters of cells leads to errors in the segmentation algorithm that causes cells around the cluster to segmented incorrectly. 
# When this happens, a cell is segmented around the same segmentation as the nucleus, giving it the same area which is very small for a normal cardiac fibroblast cell. To detect poorly segmented cells, we use:
# 
# - **Cells area in conjunction with nuclei area:** The Cells Area metric quantifies the number of pixels in a cell segmentation.
# We detect cells that are abnormally small using nuclei area as the threshold to try and avoid removing too many properly segmented cells. 
# We find cells that are with 2 standard deviations above the nuclei mean are more likely  to be poor cell segmentation due to high confluence clusters.

# In[1]:


import pathlib
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from cosmicqc import find_outliers


# ## Set paths and variables

# In[2]:


# Directory with data
data_dir = pathlib.Path("./data/converted_profiles/")

# Directory to save cleaned data
cleaned_dir = pathlib.Path("./data/cleaned_profiles/")
cleaned_dir.mkdir(exist_ok=True)

# Directory to save qc figures
qc_fig_dir = pathlib.Path("./qc_figures")
qc_fig_dir.mkdir(exist_ok=True)

# Create an empty dictionary to store data frames for each plate
all_qc_data_frames = {}

# metadata columns to include in output data frame
metadata_columns = [
    "Image_Metadata_Plate",
    "Image_Metadata_Well",
    "Image_Metadata_Site",
    "Metadata_Nuclei_Location_Center_X",
    "Metadata_Nuclei_Location_Center_Y",
]

# Set outlier threshold to use for all metrics (# of std deviations away from the mean)
outlier_threshold = 2


# ## Load in plate to perform QC on

# In[3]:


# Set plate as variable to load in
plate = "localhost231120090001"

# Load in converted plate data
plate_df = pd.read_parquet(f"{data_dir}/{plate}_converted.parquet")

# Add plate values if plate_4 is loaded in as it is absent due to error
if plate == 'localhost231120090001':
    plate_df['Image_Metadata_Plate'] = 'localhost231120090001'

print(plate_df.shape)
plate_df.head()


# ## Oversegmented nuclei with very high intensity (clusters)

# In[4]:


# find large nuclei and high intensity
feature_thresholds = {
    "Nuclei_AreaShape_Area": outlier_threshold,
    "Nuclei_Intensity_IntegratedIntensity_Hoechst": outlier_threshold,
}

large_nuclei_high_int_outliers = find_outliers(
    df=plate_df,
    metadata_columns=metadata_columns,
    feature_thresholds=feature_thresholds
)

print(large_nuclei_high_int_outliers.shape)
large_nuclei_high_int_outliers.head()


# ### Scatterplot of area to integrated intensity

# In[5]:


# Set the default value to 'inlier'
plate_df['Outlier_Status'] = 'Single-cell passed QC'

# Update the 'Outlier_Status' column based on the outliers DataFrame using index
plate_df.loc[plate_df.index.isin(large_nuclei_high_int_outliers.index), 'Outlier_Status'] = 'Single-cell failed QC'

# Create scatter plot
plt.figure(figsize=(10, 6))
plot = sns.scatterplot(
    data=plate_df,
    x="Nuclei_AreaShape_Area",
    y="Nuclei_Intensity_IntegratedIntensity_Hoechst",
    hue="Outlier_Status",
    palette={'Single-cell passed QC': '#006400', 'Single-cell failed QC': '#990090'},  # Specify colors
    alpha=0.6
)

# Add threshold lines
plt.axvline(
    x=large_nuclei_high_int_outliers['Nuclei_AreaShape_Area'].min(),
    color="r",
    linestyle="--",
    label='Min. threshold for Nuclei Area'
)
plt.axhline(
    y=large_nuclei_high_int_outliers['Nuclei_Intensity_IntegratedIntensity_Hoechst'].min(),
    color="b",
    linestyle="--",
    label='Min. threshold for Nuclei Intensity'
)

plt.title(f"Nuclei Area vs. Nuclei Integrated Intensity for {plate}")
plt.xlabel("Nuclei Area")
plt.ylabel("Nuclei Integrated Intensity (Hoechst)")
plt.tight_layout()

# Show the legend
plt.legend(loc='lower right', bbox_to_anchor=(1.0, 0), prop={'size': 10})

# Save figure
plt.savefig(pathlib.Path(f"{qc_fig_dir}/{plate}_nuclei_outliers.png"), dpi=500)

plt.show()


# ## Identify all non-circular (poorly-segmented) nuclei

# In[6]:


# find non-circular nuclei (low formfactor means below mean so must add a negative to the outlier threshold)
feature_thresholds = {
    "Nuclei_AreaShape_FormFactor": -outlier_threshold,
}

low_formfactor_outliers = find_outliers(
    df=plate_df,
    metadata_columns=metadata_columns,
    feature_thresholds=feature_thresholds
)

print(low_formfactor_outliers.shape)
low_formfactor_outliers.head()


# ### Plot FormFactor distribution labelling outlier status

# In[7]:


# Reset the default value to 'inlier'
plate_df["Outlier_Status"] = "Single-cell passed QC"

# Update the 'Outlier_Status' column based on the outliers DataFrame using index
plate_df.loc[plate_df.index.isin(low_formfactor_outliers.index), "Outlier_Status"] = (
    "Single-cell failed QC"
)

# Create histogram
plt.figure(figsize=(10, 6))
sns.histplot(
    data=plate_df,
    x="Nuclei_AreaShape_FormFactor",
    hue="Outlier_Status",
    multiple="stack",
    bins=50,  # Adjust the number of bins as needed
    palette={"Single-cell passed QC": "#006400", "Single-cell failed QC": "#990090"},
    legend=True,
)

plt.title(f"Distribution of Nuclei FormFactor all outlier statuses for plate {plate}")
plt.xlabel("Nuclei FormFactor")
plt.ylabel("Count")
plt.tight_layout()

# Show the legend
plt.legend(
    title="QC Status",
    loc="upper left",
    prop={"size": 10},
    labels=["Failed QC", "Passed QC"],
)

# save figure
plt.savefig(pathlib.Path(f"{qc_fig_dir}/{plate}_nuclei_formfactor_outliers.png"), dpi=500)

plt.show()


# ## Mis-segmented cells due to high confluence (segmentation for cells around the nuclei)

# In[8]:


# find large nuclei and high intensity
feature_thresholds = {
    "Cells_AreaShape_Area": -1.08,
}

small_cells_outliers = find_outliers(
    df=plate_df,
    metadata_columns=metadata_columns,
    feature_thresholds=feature_thresholds
)

print(small_cells_outliers.shape)
small_cells_outliers.sort_values(by="Cells_AreaShape_Area", ascending=False).head()


# In[9]:


# Create a density plot
plt.figure(figsize=(10, 6))
sns.kdeplot(x='Cells_AreaShape_Area', data=plate_df, fill=True)

# Add threshold line
plt.axvline(
    x=small_cells_outliers["Cells_AreaShape_Area"].max(),
    color="r",
    linestyle="--",
    label=f'Threshold for Outliers: < {small_cells_outliers["Cells_AreaShape_Area"].max()}',
)

# Set labels and title
plt.ylabel('Count')
plt.xlabel('Cells Area')
plt.title(f'Distribution of cells area for {plate}')
plt.legend()
plt.tight_layout()

# save figure
plt.savefig(pathlib.Path(f"{qc_fig_dir}/{plate}_cells_outliers.png"), dpi=500)

# Show the plot
plt.show()


# ## Remove all outliers and save cleaned data frame

# In[10]:


# Assuming nuclei_outliers_df and cells_outliers_df have the same index
outlier_indices = pd.concat([large_nuclei_high_int_outliers, small_cells_outliers, low_formfactor_outliers]).index

# Remove rows with outlier indices from plate_4_df
plate_4_df_cleaned = plate_df.drop(outlier_indices)

# Save cleaned data for this plate
plate_name = plate_df['Image_Metadata_Plate'].iloc[0]
plate_4_df_cleaned.to_parquet(f"{cleaned_dir}/{plate_name}_cleaned.parquet")

# Verify the result
print(plate_4_df_cleaned.shape)
plate_4_df_cleaned.head()

