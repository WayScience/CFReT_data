#!/usr/bin/env python
# coding: utf-8

# # Find one top representative single-cell crop visualize per top healthy and failing coefficients
# 
# The features with the highest coefficients per channel per cell type are identified *manually*. We do not use correlation features due to having to use two features and those are harder to visualize.
# 
# We first filter the plate 4 data frame to only include isolated cells (0 cell neighbors adjacent) and then filter out any single-cell that is too close to an edge (based on crop_size).
# 
# ### **Top healthy features per channel**
# 
# | Organelle | Feature | Coefficient |
# |-----------------|-----------------|-----------------|
# | Nucleus | Nuclei_Intensity_MeanIntensityEdge_Hoechst | 1.41883228668943 |
# | ER | Nuclei_Texture_AngularSecondMoment_ER_3_01_256 | 0.247086864563417 |
# | Golgi/Plasma membrane | Nuclei_Intensity_IntegratedIntensity_PM | 0.551412935417769 |
# | Mitochondria | Cells_Intensity_MinIntensityEdge_Mitochondria | 0.247078061220504 |
# | F-actin | Nuclei_Intensity_MinIntensity_Actin | 0.623372719535748 |
# 
# ### **Top failing features per channel**
# 
# | Organelle | Feature | Coefficient |
# |-----------------|-----------------|-----------------|
# | Nucleus | Cells_Intensity_StdIntensityEdge_Hoechst | 0.736873911951805 |
# | ER | Cytoplasm_Texture_InverseDifferenceMoment_ER_3_01_256 | 0.357607115585286 |
# | Golgi/Plasma membrane | Cytoplasm_Texture_InfoMeas1_PM_3_00_256 | 0.518741372506253 |
# | Mitochondria | Nuclei_Intensity_IntegratedIntensity_Mitochondria | 0.26804586648328 |
# | F-actin | Cells_Intensity_IntegratedIntensityEdge_Actin | 1.35720412209369 |
# 

# In[1]:


import pathlib
from pprint import pprint

import cv2
import pandas as pd
from typing import List


# In[2]:


# Function for formatting min/max row data frames into dictionaries
def create_sc_dict(dfs: List[pd.DataFrame], names: List[str]) -> dict:
    """Format lists of data frames and names into a dictionary with all relevant metadata to find single-cell images.

    Args:
        dfs (List[pd.DataFrame]): List of data frames each containing a single cell and relevant metadata.
        names (List[str]): List of names corresponding to the data frames.

    Returns:
        dict: Dictionary containing info relevant for finding single-cell crops.
    """
    sc_dict = {}
    for df, name in zip(dfs, names):
        for i, (_, row) in enumerate(df.iterrows()):
            key = f"{name}"
            sc_dict[key] = {
                "plate": row["Metadata_Plate"],
                "well": row["Metadata_Well"],
                "site": row["Metadata_Site"],
                "location_center_x": row["Metadata_Nuclei_Location_Center_X"],
                "location_center_y": row["Metadata_Nuclei_Location_Center_Y"],
            }
    return sc_dict


# In[3]:


# Function for generating and saving single-cell crops per channel as PNGs
def generate_sc_crops(
    sc_dict: dict,
    images_dir: pathlib.Path,
    output_img_dir: pathlib.Path,
    crop_size: int,
) -> None:
    """Using a dictionary with single-cell metadata info per image set, single-cell crops per channel are generated 
    and saved as PNGs in an image set folder.

    Args:
        sc_dict (dict): Dictionary containing info relevant for finding single-cell crops.
        images_dir (pathlib.Path): Directory where illumination corrected images are found.
        output_img_dir (pathlib.Path): Main directory to save each image set single-cell crops
        crop_size (int): Size of the box in pixels (example: setting crop_size as 250 will make a 250x250 pixel crop around the single-cell center coordinates)
    """
    for key, info in sc_dict.items():
        # Initialize a list to store file paths for every image set
        file_paths = []

        # Create file paths with well, site, and channel
        for i in range(5):  # Update the range to start from 0 and end at 4
            filename = f"{images_dir}/{info['plate']}_{info['well']}{info['site']}d{i}_illumcorrect.tiff"
            file_paths.append(filename)

            # Read the image
            channel_image = cv2.imread(filename, cv2.IMREAD_UNCHANGED)

            # Use the location_center_x and location_center_y to create a crop
            center_x = info.get("location_center_x")
            center_y = info.get("location_center_y")

            # Crop dimensions (including crop_size)
            half_crop = crop_size // 2

            # Ensure the center coordinates are valid
            if center_x is not None and center_y is not None:
                # Calculate crop boundaries
                top_left_x = max(int(center_x - half_crop), 0)
                top_left_y = max(int(center_y - half_crop), 0)
                bottom_right_x = min(int(center_x + half_crop), channel_image.shape[1])
                bottom_right_y = min(int(center_y + half_crop), channel_image.shape[0])

                # Perform cropping
                cropped_channel = channel_image[
                    top_left_y:bottom_right_y, top_left_x:bottom_right_x
                ]

                # Ensure the cropped image is of size 250x250
                cropped_channel = cv2.resize(cropped_channel, (crop_size, crop_size))

                # Make directory for the key to keep all channels for an image in one folder
                key_dir = pathlib.Path(f"{output_img_dir}/{key}")
                key_dir.mkdir(exist_ok=True, parents=True)

                # Save the cropped image with single_cell and channel information
                output_filename = str(pathlib.Path(f"{key_dir}/{key}_d{i}_cropped.png"))
                cv2.imwrite(output_filename, cropped_channel)


# In[4]:


# Images directory for plate 4
images_dir = pathlib.Path(
    "../1.preprocessing_data/Corrected_Images/localhost231120090001/"
).resolve(strict=True)

# Output dir for cropped images
output_img_dir = pathlib.Path("./sc_crops")
output_img_dir.mkdir(exist_ok=True)

# Define the size of the cropping box (250x250 pixels)
crop_size = 250

# Create open list for one row data frames for each top feature per channel per cell type
list_of_dfs = []

# Create open list of names to assign each data frame in a list relating to the feature, channel, and cell type
list_of_names = []


# In[5]:


# Load in normalized + feature selected data as data frame
plate4_df = pd.read_parquet(
    pathlib.Path(
        "../3.process_cfret_features/data/single_cell_profiles/localhost231120090001_sc_feature_selected.parquet"
    )
)

# Load in annotated dataframe to extract neighbors
annot_df = pd.read_parquet(
    pathlib.Path(
        "../3.process_cfret_features/data/single_cell_profiles/localhost231120090001_sc_annotated.parquet"
    ),
    columns=[
        "Metadata_Well",
        "Metadata_Site",
        "Metadata_Nuclei_Number_Object_Number",
        "Cells_Neighbors_NumberOfNeighbors_Adjacent",
    ],
)

plate4_df = plate4_df.merge(
    annot_df,
    on=["Metadata_Well", "Metadata_Site", "Metadata_Nuclei_Number_Object_Number"],
    how="inner",
)

plate4_df.rename(
    columns={
        "Cells_Neighbors_NumberOfNeighbors_Adjacent": "Metadata_Number_of_Cells_Neighbors_Adjacent"
    },
    inplace=True,
)

print(plate4_df.shape)
plate4_df.head()


# In[6]:


# Filter the DataFrame directly
filtered_plate4_df = plate4_df[
    (plate4_df['Metadata_Number_of_Cells_Neighbors_Adjacent'].isin([0])) &
    (plate4_df['Metadata_Nuclei_Location_Center_X'] > crop_size // 2) &
    (plate4_df['Metadata_Nuclei_Location_Center_X'] < (plate4_df['Metadata_Nuclei_Location_Center_X'].max() - crop_size // 2)) &
    (plate4_df['Metadata_Nuclei_Location_Center_Y'] > crop_size // 2) &
    (plate4_df['Metadata_Nuclei_Location_Center_Y'] < (plate4_df['Metadata_Nuclei_Location_Center_Y'].max() - crop_size // 2))
]

print(filtered_plate4_df.shape)
filtered_plate4_df.head()


# ## Healthy features

# In[7]:


# Get data frame with the top single-cell from the top healthy nuclei coefficient
top_healthy_nuclei = filtered_plate4_df[filtered_plate4_df["Metadata_cell_type"] == "Healthy"].nlargest(
    1, "Nuclei_Intensity_MeanIntensityEdge_Hoechst"
)[
    [
        "Nuclei_Intensity_MeanIntensityEdge_Hoechst",
        "Metadata_Well",
        "Metadata_Plate",
        "Metadata_Site",
        "Metadata_Number_of_Cells_Neighbors_Adjacent",
        "Metadata_Nuclei_Location_Center_X",
        "Metadata_Nuclei_Location_Center_Y",
        "Metadata_heart_number",
        "Metadata_cell_type",
    ]
]

# Append the DataFrame and its name to the lists
list_of_dfs.append(top_healthy_nuclei)
list_of_names.append("top_healthy_nuclei")

print(top_healthy_nuclei.shape)
top_healthy_nuclei


# In[8]:


# Get data frame with the top single-cell from the top healthy ER coefficient
top_healthy_ER = filtered_plate4_df[filtered_plate4_df["Metadata_cell_type"] == "Healthy"].nlargest(
    1, "Nuclei_Texture_AngularSecondMoment_ER_3_01_256"
)[
    [
        "Nuclei_Texture_AngularSecondMoment_ER_3_01_256",
        "Metadata_cell_type",
        "Metadata_Well",
        "Metadata_Plate",
        "Metadata_Site",
        "Metadata_Nuclei_Location_Center_X",
        "Metadata_Nuclei_Location_Center_Y",
        "Metadata_heart_number",
    ]
]

# Append the DataFrame and its name to the lists
list_of_dfs.append(top_healthy_ER)
list_of_names.append("top_healthy_ER")

print(top_healthy_ER.shape)
top_healthy_ER


# In[9]:


# Get data frame with the top single-cell from the top healthy PM coefficient
top_healthy_PM = filtered_plate4_df[filtered_plate4_df["Metadata_cell_type"] == "Healthy"].nlargest(
    1, "Nuclei_Intensity_IntegratedIntensity_PM"
)[
    [
        "Nuclei_Intensity_IntegratedIntensity_PM",
        "Metadata_cell_type",
        "Metadata_Well",
        "Metadata_Plate",
        "Metadata_Site",
        "Metadata_Nuclei_Location_Center_X",
        "Metadata_Nuclei_Location_Center_Y",
        "Metadata_heart_number",
    ]
]

# Append the DataFrame and its name to the lists
list_of_dfs.append(top_healthy_PM)
list_of_names.append("top_healthy_PM")

print(top_healthy_PM.shape)
top_healthy_PM


# In[10]:


# Get data frame with the top single-cell from the top healthy Mitochondria coefficient
top_healthy_mito = filtered_plate4_df[filtered_plate4_df["Metadata_cell_type"] == "Healthy"].nlargest(
    1, "Cells_Intensity_MinIntensityEdge_Mitochondria"
)[
    [
        "Cells_Intensity_MinIntensityEdge_Mitochondria",
        "Metadata_cell_type",
        "Metadata_Well",
        "Metadata_Plate",
        "Metadata_Site",
        "Metadata_Nuclei_Location_Center_X",
        "Metadata_Nuclei_Location_Center_Y",
        "Metadata_heart_number",
    ]
]

# Append the DataFrame and its name to the lists
list_of_dfs.append(top_healthy_mito)
list_of_names.append("top_healthy_mito")

print(top_healthy_mito.shape)
top_healthy_mito


# In[11]:


# Get data frame with the top single-cell from the top healthy Actin coefficient
top_healthy_actin = filtered_plate4_df[filtered_plate4_df["Metadata_cell_type"] == "Healthy"].nlargest(
    1, "Nuclei_Intensity_MinIntensity_Actin"
)[
    [
        "Nuclei_Intensity_MinIntensity_Actin",
        "Metadata_cell_type",
        "Metadata_Well",
        "Metadata_Plate",
        "Metadata_Site",
        "Metadata_Nuclei_Location_Center_X",
        "Metadata_Nuclei_Location_Center_Y",
        "Metadata_heart_number",
    ]
]

# Append the DataFrame and its name to the lists
list_of_dfs.append(top_healthy_actin)
list_of_names.append("top_healthy_actin")

print(top_healthy_actin.shape)
top_healthy_actin


# ## Failing features

# In[12]:


# Get data frame with the top single-cell from the top failing nuclei coefficient
top_failing_nuclei = filtered_plate4_df[filtered_plate4_df["Metadata_cell_type"] == "Failing"].nlargest(
    1, "Cells_Intensity_StdIntensityEdge_Hoechst"
)[
    [
        "Cells_Intensity_StdIntensityEdge_Hoechst",
        "Metadata_cell_type",
        "Metadata_Well",
        "Metadata_Plate",
        "Metadata_Site",
        "Metadata_Nuclei_Location_Center_X",
        "Metadata_Nuclei_Location_Center_Y",
        "Metadata_heart_number",
    ]
]

# Append the DataFrame and its name to the lists
list_of_dfs.append(top_failing_nuclei)
list_of_names.append("top_failing_nuclei")

print(top_failing_nuclei.shape)
top_failing_nuclei


# In[13]:


# Get data frame with the top single-cell from the top failing ER coefficient
top_failing_ER = filtered_plate4_df[filtered_plate4_df["Metadata_cell_type"] == "Failing"].nlargest(
    1, "Cytoplasm_Texture_InverseDifferenceMoment_ER_3_01_256"
)[
    [
        "Cytoplasm_Texture_InverseDifferenceMoment_ER_3_01_256",
        "Metadata_cell_type",
        "Metadata_Well",
        "Metadata_Plate",
        "Metadata_Site",
        "Metadata_Nuclei_Location_Center_X",
        "Metadata_Nuclei_Location_Center_Y",
        "Metadata_heart_number",
    ]
]

# Append the DataFrame and its name to the lists
list_of_dfs.append(top_failing_ER)
list_of_names.append("top_failing_ER")

print(top_failing_ER.shape)
top_failing_ER


# In[14]:


# Get data frame with the top single-cell from the top failing PM coefficient
top_failing_PM = filtered_plate4_df[filtered_plate4_df["Metadata_cell_type"] == "Failing"].nlargest(
    1, "Cytoplasm_Texture_InfoMeas1_PM_3_00_256"
)[
    [
        "Cytoplasm_Texture_InfoMeas1_PM_3_00_256",
        "Metadata_cell_type",
        "Metadata_Well",
        "Metadata_Plate",
        "Metadata_Site",
        "Metadata_Nuclei_Location_Center_X",
        "Metadata_Nuclei_Location_Center_Y",
        "Metadata_heart_number",
    ]
]

# Append the DataFrame and its name to the lists
list_of_dfs.append(top_failing_PM)
list_of_names.append("top_failing_PM")

print(top_failing_PM.shape)
top_failing_PM


# In[15]:


# Get data frame with the top single-cell from the top failing Mitochondria coefficient
top_failing_mito = filtered_plate4_df[filtered_plate4_df["Metadata_cell_type"] == "Failing"].nlargest(
    1, "Nuclei_Intensity_IntegratedIntensity_Mitochondria"
)[
    [
        "Nuclei_Intensity_IntegratedIntensity_Mitochondria",
        "Metadata_cell_type",
        "Metadata_Well",
        "Metadata_Plate",
        "Metadata_Site",
        "Metadata_Nuclei_Location_Center_X",
        "Metadata_Nuclei_Location_Center_Y",
        "Metadata_heart_number",
    ]
]

# Append the DataFrame and its name to the lists
list_of_dfs.append(top_failing_mito)
list_of_names.append("top_failing_mito")

print(top_failing_mito.shape)
top_failing_mito


# In[16]:


# Get data frame with the top single-cell from the top failing Actin coefficient
top_failing_actin = filtered_plate4_df[filtered_plate4_df["Metadata_cell_type"] == "Failing"].nlargest(
    1, "Cells_Intensity_IntegratedIntensityEdge_Actin"
)[
    [
        "Cells_Intensity_IntegratedIntensityEdge_Actin",
        "Metadata_cell_type",
        "Metadata_Well",
        "Metadata_Plate",
        "Metadata_Site",
        "Metadata_Nuclei_Location_Center_X",
        "Metadata_Nuclei_Location_Center_Y",
        "Metadata_heart_number",
    ]
]

# Append the DataFrame and its name to the lists
list_of_dfs.append(top_failing_actin)
list_of_names.append("top_failing_actin")

print(top_failing_actin.shape)
top_failing_actin


# In[17]:


sc_dict = create_sc_dict(dfs=list_of_dfs, names=list_of_names)

# Check the created dictionary for the first two items
pprint(list(sc_dict.items())[:2], indent=4)


# In[18]:


generate_sc_crops(sc_dict=sc_dict, images_dir=images_dir, output_img_dir=output_img_dir, crop_size=crop_size)

