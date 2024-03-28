#!/usr/bin/env python
# coding: utf-8

# # Whole image quality control metric evaluation
# 
# In this notebook, we will use the outputted QC metrics to start working on developing thresholds using z-score to flag and skip images during CellProfiler illumination correction. These are currently not proven to be generalizable and will only be applied to one plate.
# 
# **Note**: We will be using Plate 4 measurements only. 
# 
# **Blur metric to detect out of focus images** -> PowerLogLogSlope
# 
# **Saturation metric to detect large smudges** -> PercentMaximal

# ## Import libraries

# In[1]:


import pathlib
import pandas as pd

from scipy.stats import zscore
import matplotlib.pyplot as plt
import seaborn as sns


# ## Set paths and load in data frame

# In[2]:


# Directory for figures to be outputted
figure_dir = pathlib.Path("./qc_figures")
figure_dir.mkdir(exist_ok=True)

# Focus on Plate 4
plate = "localhost231120090001"

# Read in CSV with all image quality metrics per image
qc_df = pd.read_csv(pathlib.Path(f"./qc_results/{plate}/Image.csv"))

# Adding 'Metadata_Well' and 'Metadata_Site' columns using filename from FileName_OrigActin (same for all channels)
qc_df['Metadata_Well'] = qc_df['FileName_OrigActin'].str.extract(r'_(\w{3})\w*\.TIF')
qc_df['Metadata_Site'] = qc_df['FileName_OrigActin'].str.extract(r'_\w{3}(\w{3})\w*\.TIF')

print(qc_df.shape)
qc_df.head()


# ## Create concat dataframe combining blur and saturation metrics from all channels

# In[3]:


# List of channels
channels = ["Actin", "DNA", "ER", "PM", "Mito"]

# Create DataFrames for each channel with all Metadata columns
Actin_df = qc_df.filter(like="Metadata_").drop(columns=["Metadata_Series", "Metadata_Frame"]).copy()
DNA_df = Actin_df.copy()
ER_df = Actin_df.copy()
PM_df = Actin_df.copy()
Mito_df = Actin_df.copy()

# Iterate through each channel and add the specified columns
for channel in channels:
    # Add PowerLogLogSlope column
    globals()[f"{channel}_df"][f"ImageQuality_PowerLogLogSlope_Orig{channel}"] = qc_df[
        f"ImageQuality_PowerLogLogSlope_Orig{channel}"
    ]

    # Add PercentMaximal column
    globals()[f"{channel}_df"][f"ImageQuality_PercentMaximal_Orig{channel}"] = qc_df[
        f"ImageQuality_PercentMaximal_Orig{channel}"
    ]

    # Rename columns for each channel to concat
    globals()[f"{channel}_df"] = globals()[f"{channel}_df"].rename(
        columns={
            f"ImageQuality_PowerLogLogSlope_Orig{channel}": "ImageQuality_PowerLogLogSlope",
            f"ImageQuality_PercentMaximal_Orig{channel}": f"ImageQuality_PercentMaximal"
        }
    )

    # Add "Channel" column
    globals()[f"{channel}_df"]["Channel"] = channel

# Concat the channel data frames together for plotting
df = pd.concat([Actin_df, DNA_df, ER_df, PM_df, Mito_df], ignore_index=True)

print(df.shape)
df.head()


# ## Visualize blur metric
# 
# Based on the plot below, we can see the `actin` channel distribution is more different than the rest of the channels. We found that this difference is due to this channel on average being more dim, but not blurry. We decided that the DNA and PM channels were best at identifying blurry images and minimized removing good quality images. **We used the z-scoring method to identify max and min thresholds using 2 standard deviations.**

# In[4]:


sns.set_style('whitegrid')
sns.kdeplot(data=df, x='ImageQuality_PowerLogLogSlope', hue='Channel', palette=['r', 'b', 'g', 'y', 'magenta'],fill=True, common_norm=False)
plt.title(f'Density plots per channel for {plate}')
plt.xlabel('ImageQuality_PowerLogLogSlope')
plt.ylabel('Density')

plt.tight_layout()
plt.savefig(pathlib.Path(f"{figure_dir}/{plate}_channels_blur_density.png"), dpi=500)
plt.show()


# In[5]:


summary_statistics = df["ImageQuality_PowerLogLogSlope"].describe()
print(summary_statistics)


# In[6]:


# Calculate Z-scores for the column
z_scores = zscore(df['ImageQuality_PowerLogLogSlope'])

# Set a threshold for Z-scores (adjust as needed for number of standard deviations away from the mean)
threshold_z = 2

# Identify outlier rows based on Z-scores above and below the mean since we are using absolute values of the z-scores
blur_outliers = df[abs(z_scores) > threshold_z]

# Remove any blur outliers detected in for all channels (currently we have concluded DNA and PM channels are best able to detect blurry images)
blur_outliers = blur_outliers[(blur_outliers['Channel'] == 'DNA') | (blur_outliers['Channel'] == 'PM')]

print(blur_outliers.shape)
print(blur_outliers['Channel'].value_counts())
blur_outliers.sort_values(by='ImageQuality_PowerLogLogSlope', ascending=False).head()


# In[7]:


# Calculate the mean and standard deviation
mean_value = df["ImageQuality_PowerLogLogSlope"].mean()
std_dev = df["ImageQuality_PowerLogLogSlope"].std()

# Set the threshold multiplier for above and below the mean
threshold = 2

# Calculate the threshold values
threshold_value_above_mean = mean_value + threshold * std_dev
threshold_value_below_mean = mean_value - threshold * std_dev

# Print the calculated threshold values
print("Threshold for outliers above the mean:", threshold_value_above_mean)
print("Threshold for outliers below the mean:", threshold_value_below_mean)


# In[8]:


sns.set_style('whitegrid')
sns.kdeplot(data=df[df['Channel'].isin(['DNA', 'PM'])], x='ImageQuality_PowerLogLogSlope', hue='Channel', palette=['b', 'magenta'], fill=True, common_norm=False)
plt.title(f'DNA and PM channels are best for identifying blur\n for {plate}')
plt.xlabel('ImageQuality_PowerLogLogSlope')
plt.ylabel('Density')

plt.axvline(x=-1.2827898535807258, color='k', linestyle='--')
plt.axvline(x=-2.3840811649380393, color='k', linestyle='--')

plt.tight_layout()
plt.savefig(pathlib.Path(f"{figure_dir}/{plate}_DNA_PM_blur_density.png"), dpi=500)
plt.show()


# ## Saturation metric
# 
# For saturation metrics, we are looking for:
# 
# - Smudged images (usually seen in `DNA` channel regardless of stain)
# - Specifically for the CFReT project, some cells will start clustering where they are overlapping, causing issues during segmentation. Based on the below threshold method, we can see that intensely clustered cells are easily identified in the `PM` channel.
# 
# This means that we will only be setting threshold for saturation outliers that are detected from the DNA and PM channels in the CellProfiler IC pipeline.

# In[9]:


summary_statistics = df["ImageQuality_PercentMaximal"].describe()
print(summary_statistics)


# In[10]:


# Calculate Z-scores for the column
z_scores = zscore(df['ImageQuality_PercentMaximal'])

# Set a threshold for Z-scores (adjust as needed for number of standard deviations away from the mean)
threshold_z = 2

# Identify outlier rows based on Z-scores greater than as to identify whole images with abnormally high saturated pixels
sat_outliers = df[abs(z_scores) > threshold_z]

# Remove any outliers detected in other channels (currently we have concluded that those other channels don't detect artifacts)
sat_outliers = sat_outliers[(sat_outliers['Channel'] == 'DNA') | (sat_outliers['Channel'] == 'PM')]

# Append saturation outliers to blur outliers, drop any duplicate FOVs and save outliers for QC report
outliers = pd.concat([sat_outliers, blur_outliers], ignore_index=True).drop_duplicates()

# Save outliers data frame to view in report
outliers.to_csv("./qc_results/qc_outliers.csv")

print(sat_outliers.shape)
print(sat_outliers['Channel'].value_counts())
sat_outliers.sort_values(by='ImageQuality_PercentMaximal', ascending=True).head()


# In[11]:


# Calculate the mean and standard deviation
mean_value = df["ImageQuality_PercentMaximal"].mean()
std_dev = df["ImageQuality_PercentMaximal"].std()

# Set the threshold multiplier for above and below the mean
threshold = 2

# Calculate the threshold values
threshold_value_above_mean = mean_value + threshold * std_dev

# Print the calculated threshold values
print("Threshold for outliers above the mean:", threshold_value_above_mean)


# In[12]:


# Create a histogram plot
plt.figure(figsize=(10, 6))
plt.hist(df['ImageQuality_PercentMaximal'], bins=80, color='skyblue', edgecolor='black', alpha=0.7)

# Add mean and threshold lines
plt.axvline(
    x=df["ImageQuality_PercentMaximal"].mean(),
    color="b",
    linestyle="--",
    label=f'Mean Percent Maximal: {df["ImageQuality_PercentMaximal"].mean():.2f}',
)
plt.axvline(
    x=0.26,
    color="r",
    linestyle="--",
    label='Threshold for Outliers: > 0.26',
)

# Set labels and title
plt.ylabel('Count')
plt.xlabel('Percent Maximal')
plt.title(f'Histogram plot of percent maximal for {plate}')
plt.legend()
plt.tight_layout()

plt.savefig(pathlib.Path(f"{figure_dir}/{plate}_percent_maximal.png"), dpi=500)

# Show the plot
plt.show()


# ### Zoomed-in plot

# In[13]:


# Create a histogram plot
plt.figure(figsize=(10, 6))
plt.hist(df['ImageQuality_PercentMaximal'], bins=50, color='skyblue', edgecolor='black', alpha=0.7)

# Add mean and threshold lines
plt.axvline(
    x=df["ImageQuality_PercentMaximal"].mean(),
    color="b",
    linestyle="--",
    label=f'Mean Percent Maximal: {df["ImageQuality_PercentMaximal"].mean():.2f}',
)
plt.axvline(
    x=0.26,
    color="r",
    linestyle="--",
    label='Threshold for Outliers: > 0.26',
)

# Set labels
plt.ylabel('Count')
plt.xlabel('Percent Maximal')

# Set the zoomed-in axis ranges
plt.xlim(0, 1)
plt.ylim(0, 100)

plt.tight_layout()

plt.savefig(pathlib.Path(f"{figure_dir}/{plate}_percent_maximal_zoomed.png"), dpi=500)

# Show the plot
plt.show()

