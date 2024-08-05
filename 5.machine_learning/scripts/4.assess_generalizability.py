#!/usr/bin/env python
# coding: utf-8

# # Assess generalizability of the model by applying to different plates
# 
# **NOTE:** We currently will only be applying the model to Plate 3, split by treatments.

# ## Import libraries

# In[1]:


import pathlib
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import seaborn as sns
from joblib import load
from sklearn.metrics import precision_recall_curve

sys.path.append("../utils")
from eval_utils import generate_confusion_matrix_df, generate_f1_score_df
from training_utils import get_X_y_data


# ## Set paths and variables

# In[2]:


# Directory with plate datasets 
data_dir = pathlib.Path("../3.process_cfret_features/data/single_cell_profiles")

# Directory with models
models_dir = pathlib.Path("./models")

# Directory for model figures output
fig_dir = pathlib.Path("./figures")
fig_dir.mkdir(exist_ok=True)

# Directory for probability data to be saved
prob_dir = pathlib.Path("./prob_data")
prob_dir.mkdir(exist_ok=True)

# Load in each model individually
final_model = load(
    pathlib.Path(f"{models_dir}/log_reg_fs_plate_4_final_downsample.joblib")
)
shuffled_model = load(
    pathlib.Path(f"{models_dir}/log_reg_fs_plate_4_shuffled_downsample.joblib")
)


# ## Load in Plate 4 fs data to extract column names to filter from the other plates

# In[3]:


# Load in Plate 4 normalized feature selected data metadata (used with model) to get the feature columns to filter the plate data
parquet_metadata = pq.read_metadata(pathlib.Path(f"{data_dir}/localhost231120090001_sc_feature_selected.parquet"))

# Get the column names from the metadata
all_column_names = parquet_metadata.schema.names

# Filter out the column names that start with "Metadata_"
model_column_names = [col for col in all_column_names if not col.startswith("Metadata_")]

print(len(model_column_names))
print(model_column_names)


# ## Load in Plate 3 normalized data and drop any rows that have NaNs in the model feature columns

# In[4]:


# Load in Plate 3 data
plate_3_df = pd.read_parquet(pathlib.Path(f"{data_dir}/localhost230405150001_sc_normalized.parquet"))

# Drop rows with NaN values in feature columns that the model uses
plate_3_df = plate_3_df.dropna(subset=model_column_names)

# Capitalize the cell type values to match the model
plate_3_df["Metadata_cell_type"] = plate_3_df["Metadata_cell_type"].str.capitalize()

print(plate_3_df["Metadata_treatment"].unique())

print(plate_3_df.shape)
plate_3_df.head()


# ## Filter the Plate 3 data to only include metadata and the model feature columns

# In[5]:


#Extract metadata columns from the plate
metadata_columns = [col for col in plate_3_df.columns if col.startswith("Metadata_")]

# Extract feature columns that don't start with "Metadata_"
feature_columns = [col for col in plate_3_df.columns if not col.startswith("Metadata_")]

# Filter columns in data frame to only include those in the model
filtered_feature_columns = [col for col in plate_3_df.columns if col in model_column_names]

# Filter the DataFrame to keep only the desired columns
plate_3_df = plate_3_df[metadata_columns + filtered_feature_columns]

plate_3_df


# ## Create dictionary for Plate 3 data to split the data by the treatment

# In[6]:


# Split the plate data into different data frames with different data for applying the model to

# Define a dictionary
plate_3_dfs_dict = {}

# Filter the DataFrame to a data frame per treatment
DMSO_df = plate_3_df[plate_3_df['Metadata_treatment'] == 'DMSO']
drug_x_df = plate_3_df[plate_3_df['Metadata_treatment'] == 'drug_x']
TGFRi_df = plate_3_df[plate_3_df['Metadata_treatment'] == 'TGFRi']

# Add each DataFrame to the dictionary with a corresponding key
plate_3_dfs_dict['DMSO'] = {'data_df': DMSO_df}
plate_3_dfs_dict['drug_x'] = {'data_df': drug_x_df}
plate_3_dfs_dict['TGFRi'] = {'data_df': TGFRi_df}


# In[7]:


DMSO_df['Metadata_cell_type'].value_counts()


# ## Create a data frame with precision recall data

# In[8]:


# Initialize empty lists to store data for each iteration
precision_list = []
recall_list = []
threshold_list = []
model_type_list = []
data_type_list = []

for model_path in models_dir.iterdir():
    print("Evaluating", model_path.stem.split("_")[5], "model...")
    model_type = model_path.stem.split("_")[5]  # Extract model type
    
    # Initialize empty lists to store data for each model
    model_precision_list = []
    model_recall_list = []
    model_threshold_list = []
    
    for data, info in plate_3_dfs_dict.items():
        print("Applying model to", data, "...")
        # Load in model to apply to datasets
        model = load(model_path)

        # Load in label encoder
        le = load(
            pathlib.Path("./encoder_results/label_encoder_log_reg_fs_plate_4.joblib")
        )

        # Load in data frame associated with the data split
        data_df = info["data_df"]

        # Load in X and y data from dataset
        X, y = get_X_y_data(df=data_df, label="Metadata_cell_type")

        # Assign y classes to correct binary using label encoder results
        y_binary = le.transform(y)

        # Predict class probabilities for morphology feature data
        predicted_probs = model.predict_proba(X)

        # Calculate the precision, recall data
        precision, recall, threshold = precision_recall_curve(
            y_binary, predicted_probs[:, -1]
        )
        threshold = np.append(threshold, np.nan)

        # Append data to lists for the current model
        model_precision_list.extend(precision.tolist())
        model_recall_list.extend(recall.tolist())
        model_threshold_list.extend(threshold.tolist())

        # Append the corresponding data type to data_type_list
        data_type_list.extend([data] * len(precision))

    # Extend model_type_list once per model
    model_type_list.extend([model_type] * len(model_precision_list))

    # Extend precision_list, recall_list, and threshold_list with the data for the current model
    precision_list.extend(model_precision_list)
    recall_list.extend(model_recall_list)
    threshold_list.extend(model_threshold_list)

# Create a DataFrame from the accumulated data
pr_df = pd.DataFrame(
    {
        "Precision": precision_list,
        "Recall": recall_list,
        "Threshold": threshold_list,
        "Model_Type": model_type_list,
        "Data_Type": data_type_list,
    }
)

# Drop any NA data
pr_df.dropna(inplace=True)

# Show output of all data
print(pr_df.shape)
pr_df.head()


# ## Create PR curve with only DMSO (control) data to assess performance

# In[9]:


# PR curves with only testing and training data
plt.figure(figsize=(12, 10))
sns.set_style("whitegrid")

# Combine model and data type as one column for plotting
pr_df["data_split"] = pr_df["Model_Type"] + " (" + pr_df["Data_Type"] + ")"

# Filter data frame to only show one data_type
filtered_df = pr_df[pr_df["Data_Type"].isin(["DMSO"])]

sns.lineplot(
    x="Recall",
    y="Precision",
    hue="data_split",
    style="Model_Type",
    dashes={"final": (1, 0), "shuffled": (2, 2)},
    drawstyle='steps',
    data=filtered_df,
)

plt.legend(loc="lower right", fontsize=15)
plt.ylim(bottom=0.0, top=1.02)
plt.xlabel("Recall", fontsize=18)
plt.ylabel("Precision", fontsize=18)
plt.title("Precision vs. Recall Plate 3 DMSO Treatment Cell Type Classification", fontsize=18)

# Adjust x-axis ticks font size
plt.xticks(fontsize=14)

# Adjust y-axis ticks font size and labels
plt.yticks(fontsize=14)

plt.tight_layout()
plt.savefig(f"{fig_dir}/precision_recall_plate3_DMSO_only.png", dpi=500)

plt.show()


# ## Extract final model predicted probabilities for each treatment

# In[10]:


# Create an empty DataFrame to store the results
combined_prob_df = pd.DataFrame()

for model_path in models_dir.iterdir():
    model_type = model_path.stem.split("_")[5]  # Get the model type
    
    for data, info in plate_3_dfs_dict.items():
        # Ensure that the file is named the correct data split
        data_split = data
        print(f"Extracting {model_type} probabilities from {data} data...")

        # Load in model to apply to datasets
        model = load(model_path)

        # Load in label encoder
        le = load(
            pathlib.Path("./encoder_results/label_encoder_log_reg_fs_plate_4.joblib")
        )

        # Get unique cell types and their corresponding encoded values
        unique_labels = le.classes_
        encoded_values = le.transform(unique_labels)

        # Create a dictionary mapping encoded values to original labels
        label_dict = dict(zip(encoded_values, unique_labels))

        # Load in data frame associated with the data split
        data_df = info["data_df"].reset_index(drop=True)

        # Load in X data to get predicted probabilities
        X, _ = get_X_y_data(df=data_df, label="Metadata_cell_type")

        # Predict class probabilities for morphology feature data
        predicted_probs = model.predict_proba(X)

        # Storing probabilities in a pandas DataFrame
        prob_df = pd.DataFrame(predicted_probs, columns=model.classes_)

        # Update column names in prob_df using the dictionary and add suffix "_probas"
        prob_df.columns = [label_dict[col] + '_probas' for col in prob_df.columns]

        # Add a new column called predicted_label for each row
        prob_df['predicted_label'] = prob_df.apply(lambda row: row.idxmax()[:-7], axis=1)

        # Select metadata columns from the data
        metadata_columns = data_df.filter(like='Metadata_')

        # Combine metadata columns with predicted probabilities DataFrame based on index
        prob_df = prob_df.join(metadata_columns)
        
        # Add a new column for model_type
        prob_df['model_type'] = model_type
        
        # Append the probability DataFrame to the combined DataFrame
        combined_prob_df = pd.concat([combined_prob_df, prob_df], ignore_index=True)

# Save combined prob data
combined_prob_df.to_csv(f"{prob_dir}/combined_plate_3_predicted_proba.csv", index=False)

