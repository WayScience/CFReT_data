#!/usr/bin/env python
# coding: utf-8

# ## Fit a linear model on cell morphology features
# 
# We aim to determine which features are significantly impacted by the drug treatment, adjusted by cell count.
# 
# We use normalized and feature selected data.

# In[1]:


import pathlib
import pandas as pd

from sklearn.linear_model import LinearRegression

from pycytominer import feature_select
from pycytominer.cyto_utils import infer_cp_features


# In[2]:


# Define inputs and outputs
plate = "localhost230405150001"  # Focusing on plate 3
file_suffix = "_sc_norm_fs_cellprofiler.csv.gz"

data_dir = pathlib.Path("../../../3.process_cfret_features/data/")

cp_file = pathlib.Path(data_dir, f"{plate}{file_suffix}")

output_dir = pathlib.Path("results")
output_cp_file = pathlib.Path(output_dir, f"{plate}_linear_model_TGFRi_drug.tsv")


# In[3]:


# Load data
cp_df = pd.read_csv(cp_file)

# Drop NA columns
cp_df = feature_select(
    cp_df,
    operation="drop_na_columns",
    na_cutoff=0
)

# Count number of cells per well and add to dataframe as metadata
cell_count_df = pd.DataFrame(
    cp_df.groupby("Metadata_Well").count()["Metadata_treatment"]
).reset_index()
cell_count_df.columns = ["Metadata_Well", "Metadata_cell_count_per_well"]
cp_df = cell_count_df.merge(cp_df, on=["Metadata_Well"])

# # Only for plates 1, 2, and 3: Clean the dose column to extract numeric value
# cp_df = cp_df.assign(Metadata_dose_numeric=cp_df.Metadata_dose.str.strip("uM").astype(float))

# Define CellProfiler features
cp_features = infer_cp_features(cp_df)

print(f"We are testing {len(cp_features)} CellProfiler features")
print(cp_df.shape)
cp_df.head()


# ## Fit linear model

# In[4]:


# Setup linear modeling framework -> in plate 3 we are looking at the treatments
variables = ["Metadata_cell_count_per_well", "Metadata_treatment"]
X = cp_df.loc[:, variables]

print(X.shape)
X.head()


# In[5]:


# Assuming cp_df is your DataFrame
variables = ["Metadata_cell_count_per_well", "Metadata_treatment"]
treatments_to_select = ["drug_x", "TGFRi"]

# Select rows with specific treatment values
selected_rows = X[X["Metadata_treatment"].isin(treatments_to_select)]

# Create dummy variables
dummies = pd.get_dummies(selected_rows["Metadata_treatment"])

# Concatenate dummies with the selected rows DataFrame
X = pd.concat([selected_rows, dummies], axis=1)

# Drop the original treatment column
X.drop("Metadata_treatment", axis=1, inplace=True)

print(X.shape)
X.head()


# In[6]:


# Fit linear model for each feature
lm_results = []
for cp_feature in cp_features:
    # Create a boolean mask to filter rows with the specified treatments
    mask = cp_df["Metadata_treatment"].isin(treatments_to_select)

    # Apply the mask to Subset CP data to each individual feature (univariate test)
    cp_subset_df = cp_df.loc[mask, cp_feature]

    # Fit linear model
    lm = LinearRegression(fit_intercept=True)
    lm_result = lm.fit(X=X, y=cp_subset_df)
    
    # Extract Beta coefficients
    # (contribution of feature to X covariates)
    coef = lm_result.coef_
    
    # Estimate fit (R^2)
    r2_score = lm.score(X=X, y=cp_subset_df)
    
    # Add results to a growing list
    lm_results.append([cp_feature, r2_score] + list(coef))

# Convert results to a pandas DataFrame
lm_results = pd.DataFrame(
    lm_results,
    columns=["feature", "r2_score", "cell_count_coef", "TGFRi_coef", "drug_x_coef"]
)

# Output file
lm_results.to_csv(output_cp_file, sep="\t", index=False)

print(lm_results.shape)
lm_results.head()

