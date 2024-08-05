#!/usr/bin/env python
# coding: utf-8

# # Convert SQLite outputs to parquet files with cytotable

# ## Import libraries

# In[1]:


import pathlib
import pandas as pd

# cytotable will merge objects from SQLite file into single cells and save as parquet file
from cytotable import convert, presets

import logging

# Set the logging level to a higher level to avoid outputting unnecessary errors from config file in convert function
logging.getLogger().setLevel(logging.ERROR)


# ## Set paths and variables

# In[2]:


# preset configurations based on typical CellProfiler outputs
preset = "cellprofiler_sqlite_pycytominer"

# update preset to include site metadata and cell counts
joins = presets.config["cellprofiler_sqlite_pycytominer"]["CONFIG_JOINS"].replace(
    "Image_Metadata_Well,",
    "Image_Metadata_Well, Image_Metadata_Site, Image_Count_Cells,",
)

# type of file output from cytotable (currently only parquet)
dest_datatype = "parquet"

# set path to directory with SQLite files
sqlite_dir = pathlib.Path("../2.cellprofiler_processing/cp_output")

# directory for processed data
output_dir = pathlib.Path("data")
output_dir.mkdir(parents=True, exist_ok=True)

plate_names = []

for file_path in sqlite_dir.iterdir():
    plate_names.append(file_path.stem)

# print the plate names and how many plates there are (confirmation)
print(f"There are {len(plate_names)} plates in this dataset. Below are the names:")
for name in plate_names:
    print(name)


# ## Convert SQLite to parquet files

# In[3]:


for file_path in sqlite_dir.iterdir():
    output_path = pathlib.Path(
        f"{output_dir}/converted_profiles/{file_path.stem}_converted.parquet"
    )
    print("Starting conversion with cytotable for plate:", file_path.stem)
    # Merge single cells and output as parquet file
    convert(
        source_path=str(file_path),
        dest_path=str(output_path),
        dest_datatype=dest_datatype,
        preset=preset,
        joins=joins,
        chunk_size=5000,
    )

print("All plates have been converted with cytotable!")


# # Load in converted profiles to update

# In[4]:


# Directory with converted profiles
converted_dir = pathlib.Path(f"{output_dir}/converted_profiles")

for file_path in converted_dir.iterdir():
    # Load the DataFrame from the Parquet file
    df = pd.read_parquet(file_path)

    # If any, drop rows where "Metadata_ImageNumber" is NaN (artifact of cytotable)
    df = df.dropna(subset=["Metadata_ImageNumber"])

    # Columns to move to the front
    columns_to_move = [
        "Nuclei_Location_Center_X",
        "Nuclei_Location_Center_Y",
        "Cells_Location_Center_X",
        "Cells_Location_Center_Y",
        "Image_Count_Cells",
    ]

    # Rearrange columns and add "Metadata" prefix in one line
    df = df[
        columns_to_move + [col for col in df.columns if col not in columns_to_move]
    ].rename(
        columns=lambda col: "Metadata_" + col if col in columns_to_move else col
    )

    # Save the processed DataFrame as Parquet in the same path
    df.to_parquet(file_path, index=False)


# ## Check output to confirm process worked
# 
# To confirm the number of single cells is correct, please use any database browser software to see if the number of rows in the "Per_Cells" compartment matches the number of rows in the data frame.

# In[5]:


converted_df = pd.read_parquet(
    "./data/converted_profiles/localhost220512140003_KK22-05-198_converted.parquet"
)

print(converted_df.shape)
converted_df.head()

