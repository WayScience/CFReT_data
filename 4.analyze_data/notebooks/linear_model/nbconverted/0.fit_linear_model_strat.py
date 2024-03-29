#!/usr/bin/env python
# coding: utf-8

# # Stratify to perform linear modeling on certain data

# In[1]:


import pathlib
import pandas as pd

from sklearn.linear_model import LinearRegression

from pycytominer import feature_select
from pycytominer.cyto_utils import infer_cp_features


# In[2]:


# Define inputs and outputs
plate = "localhost231120090001"  # Focusing on plate 4
file_suffix = "_sc_feature_selected.parquet"

data_dir = pathlib.Path("../../../3.process_cfret_features/data/single_cell_profiles")

data_df = pd.read_parquet(pathlib.Path(data_dir, f"{plate}{file_suffix}"))

output_dir = pathlib.Path("results")
output_cp_file = pathlib.Path(output_dir, f"{plate}_linear_model_failing_healthy_no_treatment_neighbors.tsv")

# Replace NA values with "None"
data_df['Metadata_treatment'].fillna('None', inplace=True)

print(data_df.shape)
data_df.head()


# ## Stratify data

# In[3]:


# Filter by cell type and only cells without DMSO treatment
specific_type = ["None"]
specific_cell_types = ["Failing", "Healthy"]

filtered_df = data_df[
    (data_df["Metadata_treatment"].isin(specific_type))
    & (data_df["Metadata_cell_type"].isin(specific_cell_types))
]

# Drop NA columns
cp_df = feature_select(filtered_df, operation="drop_na_columns", na_cutoff=0)

# Count number of cells per well and add to data frame as metadata
cell_count_df = pd.DataFrame(
    cp_df.groupby("Metadata_Well").count()["Metadata_treatment"]
).reset_index()
cell_count_df.columns = ["Metadata_Well", "Metadata_cell_count_per_well"]
cp_df = cell_count_df.merge(cp_df, on=["Metadata_Well"])

# Add cell neighbors per single-cell to use as co-variate
neighbors_df = pd.read_parquet(
    pathlib.Path(
        "../../../3.process_cfret_features/data/single_cell_profiles/localhost231120090001_sc_annotated.parquet"
    )
)

# Selecting the only the joining columns and cell neighbors column
neighbors_subset = neighbors_df[
    [
        "Cells_Neighbors_NumberOfNeighbors_Adjacent",
        "Metadata_Well",
        "Metadata_Site",
        "Metadata_Nuclei_Number_Object_Number",
    ]
]

# Merging the cell neighbors column to the cp_df
cp_df = pd.merge(
    cp_df,
    neighbors_subset,
    on=["Metadata_Well", "Metadata_Site", "Metadata_Nuclei_Number_Object_Number"],
    how="inner",
)

# Rename the neighbors column so it isn't considered a CellProfiler feature
cp_df = cp_df.rename(
    columns={
        "Cells_Neighbors_NumberOfNeighbors_Adjacent": "Metadata_Cells_Neighbors_NumberOfNeighbors_Adjacent"
    }
)

# Define CellProfiler features
cp_features = infer_cp_features(cp_df)

print(f"We are testing {len(cp_features)} CellProfiler features")
print(cp_df.shape)
cp_df.head()


# ## Fit linear model

# In[4]:


# Setup linear modeling framework -> in plate 4 we are looking at the treatments or cell type
variables = ["Metadata_Cells_Neighbors_NumberOfNeighbors_Adjacent", "Metadata_cell_type"]
X = cp_df.loc[:, variables]

print(X.shape)
X.head()


# In[5]:


# Set the variables and treatments used for LM
variables = ["Metadata_Cells_Neighbors_NumberOfNeighbors_Adjacent", "Metadata_cell_type"]
treatments_to_select = ["Failing", "Healthy"]

# Select rows with specific treatment values
selected_rows = X[X["Metadata_cell_type"].isin(treatments_to_select)]

# Create dummy variables
dummies = pd.get_dummies(selected_rows["Metadata_cell_type"])

# Concatenate dummies with the selected rows DataFrame
X = pd.concat([selected_rows, dummies], axis=1)

# Drop the original treatment column
X.drop("Metadata_cell_type", axis=1, inplace=True)

print(X.shape)
X.head()


# In[6]:


# Fit linear model for each feature
lm_results = []
for cp_feature in cp_features:
    # Create a boolean mask to filter rows with the specified treatments
    mask = cp_df["Metadata_cell_type"].isin(treatments_to_select)

    # Apply the mask to Subset CP data to each individual feature (univariate test)
    cp_subset_df = cp_df.loc[mask, cp_feature]

    # Fit linear model
    lm = LinearRegression()
    lm_result = lm.fit(X=X, y=cp_subset_df)
    
    # Extract Beta coefficients
    # (contribution of feature to X co-variates)
    coef = lm_result.coef_
    
    # Estimate fit (R^2)
    r2_score = lm.score(X=X, y=cp_subset_df)
    
    # Add results to a growing list
    lm_results.append([cp_feature, r2_score] + list(coef))

# Convert results to a pandas DataFrame
lm_results = pd.DataFrame(
    lm_results,
    columns=["feature", "r2_score", "cell_neighbors_coef", "failing_coef", "healthy_coef"]
)

# Output file
lm_results.to_csv(output_cp_file, sep="\t", index=False)

print(lm_results.shape)
lm_results.head()

