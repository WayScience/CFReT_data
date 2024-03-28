#!/usr/bin/env python
# coding: utf-8

# # CellProfiler Segmentation and Feature Extraction of CFReT data

# ## Import libraries

# In[1]:


import pathlib
import pprint
import sys

sys.path.append("../utils/")
import cp_parallel as cp


# ## Set paths and variables

# In[2]:


# set the run type for the parallelization
run_name = "analysis"

# path to analysis pipeline
path_to_pipeline = pathlib.Path("./pipeline/CFReT_project_CL.cppipe").resolve(strict=True)

# path to output for SQLite database files per plate folder (create if does not already exist)
output_dir = pathlib.Path("./cp_output/")
output_dir.mkdir(exist_ok=True)

# Directory where all images are separated by folder per plate
images_dir = pathlib.Path("../1.preprocessing_data/Corrected_Images").resolve(strict=True)

# list for plate names based on folders to use to create dictionary
plate_names = []

# iterate through 0.download_data and append plate names from folder names that contain image data from that plate
for file_path in images_dir.iterdir():
    plate_names.append(str(file_path.stem))

print("There are a total of", len(plate_names), "plates. The names of the plates are:")
for plate in plate_names:
    print(plate)


# ## Create dictionary with all of the necessary paths to run CellProfiler analysis

# In[3]:


# create plate info dictionary with all parts of the CellProfiler CLI command to run in parallel
plate_info_dictionary = {
    name: {
        "path_to_images": pathlib.Path(list(images_dir.rglob(name))[0]).resolve(
            strict=True
        ),
        "path_to_output": pathlib.Path(f"{output_dir}/{name}"),
        "path_to_pipeline": path_to_pipeline,

    }
    for name in plate_names if name == 'localhost231120090001'
}

# view the dictionary to assess that all info is added correctly
pprint.pprint(plate_info_dictionary, indent=4)


# ## Run CellProfiler analysis on all plates
# 
# **Note:** This code cell will not be run in this notebook due to the instability of jupyter notebooks compared to running as a python script. All CellProfiler SQLite outputs will have the same name but outputted into their respective plate folder (due to parallelization).

# In[ ]:


cp.run_cellprofiler_parallel(
    plate_info_dictionary=plate_info_dictionary, run_name=run_name
)

