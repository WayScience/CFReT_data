#!/usr/bin/env python
# coding: utf-8

# ## Extract UMAP embeddings for CFReT CP Features

# In[1]:


import glob
import pathlib
import pandas as pd
import umap
import numpy as np

from pycytominer import feature_select
from pycytominer.cyto_utils import infer_cp_features


# ## Generate Embeddings for Whole Plates

# ### Set constant for whole plates
# 
# Note: All plates (1-4) without filtering had a random seed of 1234. For plates with filtering, we use a random seed of 0 which is a standard for the Way lab.

# In[2]:


# Set constants (previously set prior, normally use 0 but the change in coordinates will impact already generated single-cell crops)
umap_random_seed = 1234
umap_n_components = 2

output_dir = pathlib.Path("results")
output_dir.mkdir(parents=True, exist_ok=True)


# ### Set paths to all plates

# In[3]:


# Set input paths
data_dir = pathlib.Path("..", "..", "..", "3.process_cfret_features", "data", "single_cell_profiles")

# Select only the feature selected files
file_suffix = "*sc_feature_selected.parquet"

# Obtain file paths for all feature selected plates
fs_files = glob.glob(f"{data_dir}/{file_suffix}")
fs_files


# ### Generate dictionary with plate and data

# In[4]:


# Load feature data into a dictionary, keyed on plate name
cp_dfs = {x.split("/")[-1]: pd.read_parquet(x) for x in fs_files}

# Print out useful information about each dataset
print(cp_dfs.keys())
[cp_dfs[x].shape for x in cp_dfs]


# In[5]:


cp_dfs


# ### Fit UMAP for whole plates

# In[6]:


# Fit UMAP features per dataset and save
for plate in cp_dfs:
    # Set plate name
    plate_name = pathlib.Path(plate).stem
    # Set output file for the UMAP
    output_umap_file = pathlib.Path(output_dir, f"UMAP_{plate_name}.tsv.gz")

    # # Check if the output file already exists
    # if output_umap_file.exists():
    #     print(f"Skipping {output_umap_file.stem} as it already exists.")
    #     continue

    # Make sure to reinitialize UMAP instance per plate
    umap_fit = umap.UMAP(
        random_state=umap_random_seed,
        n_components=umap_n_components
    )
    
    # Remove NA columns
    cp_df = cp_dfs[plate]
    cp_df = feature_select(
        cp_df,
        operation="drop_na_columns",
        na_cutoff=0
    )
    
    # Process cp_df to separate features and metadata
    cp_features = infer_cp_features(cp_df)
    meta_features = infer_cp_features(cp_df, metadata=True)
    
    # Fit UMAP and convert to pandas DataFrame
    embeddings = pd.DataFrame(
        umap_fit.fit_transform(cp_df.loc[:, cp_features]),
        columns=[f"UMAP{x}" for x in range(0, umap_n_components)]
    )
    print(embeddings.shape)
    
    # Combine with metadata
    cp_umap_with_metadata_df = pd.concat([
        cp_df.loc[:, meta_features],
        embeddings
    ], axis=1)
    
    # Generate output file, drop unnamed column, and save 
    cp_umap_with_metadata_df.to_csv(output_umap_file, index=False, sep="\t")

    # Print an example output file
    cp_umap_with_metadata_df.head()


# ## Generate embeddings for filtered data
# 
# Note: We are filtering out single-cells from plates 3 and 4 where there is more than 1 single-cell adjacent. We are looking to see the impact on the UMAP when only including "isolated" single-cells.

# In[7]:


# Set random seed as 0 for filtered datasets
filtered_umap_random_seed = 0

# Select only the feature selected files
file_suffix = "*sc_annotated.parquet"

# Obtain file paths for all annotated plates (contains neighbors data)
annot_files = glob.glob(f"{data_dir}/{file_suffix}")

plate_names = []

for file_path in pathlib.Path("../../../0.download_data/Images").iterdir():
    plate_names.append(str(file_path.stem))

print(plate_names)


# In[8]:


# create plate info dictionary
plate_info_dictionary = {
    name: {
        "fs_data": pd.read_parquet(
            pathlib.Path(
                list(data_dir.rglob(f"{name}_sc_feature_selected.parquet"))[0]
            ).resolve(strict=True)
        ),
        "annot_data": pd.read_parquet(
            pathlib.Path(
                list(data_dir.rglob(f"{name}_sc_annotated.parquet"))[0]
            ).resolve(strict=True)
        ),
    }
    for name in plate_names
    if name == "localhost230405150001" or name == "localhost231120090001"
}

# view the dictionary info to assess that all info is added correctly
print(plate_info_dictionary.keys())
print(
    "The shapes of the feature selected data frames are:",
    [plate_info_dictionary[x]["fs_data"].shape for x in plate_info_dictionary],
)
print(
    "The shapes of the annotated data frames are:",
    [plate_info_dictionary[x]["annot_data"].shape for x in plate_info_dictionary],
)


# In[9]:


for plate, info in plate_info_dictionary.items():
    # Set output file for the UMAP
    output_umap_file = pathlib.Path(output_dir, f"UMAP_{plate}_fs_filtered.tsv.gz")

    # # Check if the output file already exists
    # if output_umap_file.exists():
    #     print(f"Skipping {output_umap_file.stem} as it already exists.")
    #     continue

    # Give variable names to data frames
    fs_df = info["fs_data"]
    annot_df = info["annot_data"]

    # Merging neighbor column onto fs_df from annot_df
    fs_df = fs_df.merge(
        annot_df[[
            "Metadata_Well", "Metadata_Site", "Metadata_Nuclei_Number_Object_Number", "Cells_Neighbors_NumberOfNeighbors_Adjacent"
        ]],
        on=["Metadata_Well", "Metadata_Site", "Metadata_Nuclei_Number_Object_Number"],
        how="inner",
    )

    # Rename neighbors column to include as metadata
    fs_df = fs_df.rename(columns={"Cells_Neighbors_NumberOfNeighbors_Adjacent": "Metadata_Neighbors_Adjacent"})

    # Only including rows where Metadata_Neighbors_Adjacent is less than or equal to 1 neighbor
    filtered_fs_df = fs_df[fs_df['Metadata_Neighbors_Adjacent'] <= 1]

    # Reset index to avoid any issues with concat
    filtered_fs_df.reset_index(drop=True, inplace=True)

    # Make sure to reinitialize UMAP instance per plate (uses random seed 0 and same umap components as above)
    umap_fit = umap.UMAP(
        random_state=filtered_umap_random_seed,
        n_components=umap_n_components
    )

    # Remove NA columns
    filtered_fs_df = feature_select(
        filtered_fs_df,
        operation="drop_na_columns",
        na_cutoff=0
    )
    
    # Process filtered_fs_df to separate features and metadata
    cp_features = infer_cp_features(filtered_fs_df)
    meta_features = infer_cp_features(filtered_fs_df, metadata=True)
    
    # Fit UMAP and convert to pandas DataFrame
    embeddings = pd.DataFrame(
        umap_fit.fit_transform(filtered_fs_df.loc[:, cp_features]),
        columns=[f"UMAP{x}" for x in range(0, umap_n_components)]
    )
    print(embeddings.shape)
    
    # Combine with metadata
    filtered_umap_with_metadata_df = pd.concat([
        filtered_fs_df.loc[:, meta_features],
        embeddings
    ], axis=1)
    
    # Generate output file, drop unnamed column, and save 
    filtered_umap_with_metadata_df.to_csv(output_umap_file, index=False, sep="\t")

    # Print an example output file
    filtered_umap_with_metadata_df.head()

