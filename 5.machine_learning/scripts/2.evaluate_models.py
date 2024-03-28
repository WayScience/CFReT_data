#!/usr/bin/env python
# coding: utf-8

# # Evaluate model performance with PR curves, confusion matrixes, and F1 scores

# ## Import libraries

# In[1]:


import pathlib
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from joblib import load
from sklearn.metrics import precision_recall_curve

sys.path.append("../utils")
from eval_utils import generate_confusion_matrix_df, generate_f1_score_df
from training_utils import load_data


# ## Set paths to different datasets and models

# In[2]:


# Directory with datasets
data_dir = pathlib.Path("./data")

# Directory with models
models_dir = pathlib.Path("./models")

# Directory for model figures output
fig_dir = pathlib.Path("./figures")
fig_dir.mkdir(exist_ok=True)

# Load in each model individually
final_model = load(
    pathlib.Path(f"{models_dir}/log_reg_fs_plate_4_final_downsample.joblib")
)
shuffled_model = load(
    pathlib.Path(f"{models_dir}/log_reg_fs_plate_4_shuffled_downsample.joblib")
)


# ## Precision-Recall Curves

# ### Output PR data per dataset and combine into one df

# In[3]:


# Initialize empty lists to store data for each iteration
precision_list = []
recall_list = []
threshold_list = []
model_type_list = []
data_type_list = []

for model_path in models_dir.iterdir():
    print("Evaluating", model_path.stem.split("_")[5], "model...")
    for data_path in data_dir.iterdir():
        print("Applying model to", data_path.stem, "...")
        # load in model to apply to datasets
        model = load(model_path)

        # load in label encoder
        le = load(
            pathlib.Path("./encoder_results/label_encoder_log_reg_fs_plate_4.joblib")
        )

        # Load in X and y data from dataset
        X, y = load_data(path_to_data=data_path, label="Metadata_cell_type")

        # Assign y classes to correct binary using label encoder results
        y_binary = le.transform(y)

        # predict class probabilities for morphology feature data
        predicted_probs = model.predict_proba(X)

        # Calculate the precision, recall data
        precision, recall, threshold = precision_recall_curve(
            y_binary, predicted_probs[:, -1]
        )
        threshold = np.append(threshold, np.nan)

        # Append data to lists
        precision_list.extend(precision.tolist())
        recall_list.extend(recall.tolist())
        threshold_list.extend(threshold.tolist())
        model_type_list.extend([model_path.stem.split("_")[5]] * len(precision))
        data_type_list.extend([data_path.stem.split("_")[0]] * len(precision))

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


# ### Plot precision-recall curves for testing and training data per model type

# In[4]:


# PR curves with only testing and training data
plt.figure(figsize=(12, 10))
sns.set_style("whitegrid")

# Combine model and data type as one column for plotting
pr_df["data_split"] = pr_df["Model_Type"] + " (" + pr_df["Data_Type"] + ")"

# Filter data frame to only show training versus testing
filtered_df = pr_df[pr_df["Data_Type"].isin(["training", "testing"])]

sns.lineplot(
    x="Recall",
    y="Precision",
    hue="data_split",
    style="Model_Type",
    dashes={"final": (1, 0), "shuffled": (2, 2)},
    data=filtered_df,
)

plt.legend(loc="lower right", fontsize=15)
plt.ylim(bottom=0.0, top=1.02)
plt.xlabel("Recall", fontsize=18)
plt.ylabel("Precision", fontsize=18)
plt.title("Precision vs. Recall Plate 4 Cell Type Classification", fontsize=18)

# Adjust x-axis ticks font size
plt.xticks(fontsize=14)

# Adjust y-axis ticks font size and labels
plt.yticks(fontsize=14)

plt.tight_layout()
plt.savefig(f"{fig_dir}/precision_recall_plate4_downsample.png", dpi=500)

plt.show()


# ### Plot precision-recall curves for holdout data per model type

# In[5]:


# PR curves with hold-out data
plt.figure(figsize=(12, 10))
sns.set_style("whitegrid")

# Combine model and data type as one column for plotting
pr_df["data_split"] = pr_df["Model_Type"] + " (" + pr_df["Data_Type"] + ")"

# Filter data frame to only show holdout datasets
filtered_df = pr_df[pr_df["Data_Type"].isin(["holdout1", "holdout2"])]

sns.lineplot(
    x="Recall",
    y="Precision",
    hue="data_split",
    style="Model_Type",
    dashes={"final": (1, 0), "shuffled": (2, 2)},
    data=filtered_df,
)

plt.legend(loc="lower right", fontsize=13)
plt.ylim(bottom=0.0, top=1.02)
plt.xlabel("Recall", fontsize=18)
plt.ylabel("Precision", fontsize=18)
plt.title("Precision vs. Recall Plate 4 Cell Type Classification", fontsize=18)

# Adjust x-axis ticks font size
plt.xticks(fontsize=14)

# Adjust y-axis ticks font size and labels
plt.yticks(fontsize=14)

plt.tight_layout()
plt.savefig(f"{fig_dir}/precision_recall_plate4_holdout_data_downsample.png", dpi=500)

plt.show()


# ## Confusion Matrices

# ### Generate confusion matrix plots for each data set and model type

# In[6]:


# List of paths that contains each model
path_list = list(models_dir.glob("*"))

# List of data sets to apply models to
data_set_list = ["training", "testing", "holdout1", "holdout2"]

# Path to encoder file
encoder_path = pathlib.Path("./encoder_results/label_encoder_log_reg_fs_plate_4.joblib")

# Metadata column used for classifying
label = "Metadata_cell_type"

# Iterate over each model path
for model_path in path_list:
    # Iterate over each dataset
    for data_set in data_set_list:
        print(
            "Generating plot for",
            data_set.capitalize(),
            model_path.stem.split("_")[5].capitalize(),
        )
        # Generate confusion matrix data frame
        df = generate_confusion_matrix_df(
            model_path=model_path,
            data_dir=data_dir,
            encoder_path=encoder_path,
            label=label,
            data_set=data_set,
        )

        # Rename binary labels to failing versus healthy
        df["True_Label"] = df["True_Label"].replace({0: "Failing", 1: "Healthy"})
        df["Predicted_Label"] = df["Predicted_Label"].replace(
            {0: "Failing", 1: "Healthy"}
        )

        # Plot confusion matrix
        plt.figure(figsize=(12, 10))
        sns.set_style("whitegrid")

        sns.heatmap(
            data=df.pivot(
                index="True_Label", columns="Predicted_Label", values="Recall"
            ),
            annot=df.pivot(
                index="True_Label", columns="Predicted_Label", values="Count"
            ),
            fmt=".0f",
            cmap="Reds",
            square=True,
            vmin=0,  # Set vmin here
            vmax=1,  # Set vmax here
            cbar_kws={"label": "Recall"},
            xticklabels=["Failing", "Healthy"],
            yticklabels=["Failing", "Healthy"],
            linewidths=0.5,
            annot_kws={"size": 14},
            mask=(
                df.pivot(
                    index="True_Label", columns="Predicted_Label", values="Count"
                ).isna()
            ),
        )

        plt.title(
            f"Confusion Matrix for {data_set.capitalize()} Data - Model: {model_path.stem.split('_')[5].capitalize()}"
        )
        plt.tight_layout()
        plt.savefig(
            f"{fig_dir}/plate4_confusion_matrix_{data_set}_{model_path.stem.split('_')[5]}_downsample.png",
            dpi=500,
        )
        plt.show()


# ## F1 scores

# ### Concat F1 scores from all datasets and model types

# In[7]:


# List of paths that contains each model
path_list = list(models_dir.glob("*"))

# Path to encoder file
encoder_path = pathlib.Path("./encoder_results/label_encoder_log_reg_fs_plate_4.joblib")

# Metadata column used for classifying
label = "Metadata_cell_type"

dfs = []

# Iterate over each model path
for model_path in path_list:
    # Iterate over each dataset
    for data_set in data_set_list:
        print(data_set.capitalize(), model_path.stem.split("_")[5].capitalize())
        # Generate confusion matrix data frame
        df = generate_f1_score_df(
            model_path=model_path,
            data_dir=data_dir,
            encoder_path=encoder_path,
            label=label,
            data_set=data_set,
        )

        # Rename binary labels to failing versus healthy
        df = df.rename(columns={0: "Failing", 1: "Healthy"})

        # Add columns for model and data_set
        df["Model"] = model_path.stem.split("_")[5].capitalize()
        df["Data_Set"] = data_set.capitalize()

        # Append the DataFrame to the list
        dfs.append(df)

concat_f1_scores = pd.concat(dfs, ignore_index=True)


# ### Check if output looks correct

# In[8]:


print(concat_f1_scores.shape)
concat_f1_scores.head()


# ### Plot F1 scores for only testing and training data

# In[9]:


plt.figure(figsize=(10, 6))

# Set y-axis limits
plt.ylim(0, 1)

# Filter the data to include only training and testing sets
filtered_data = concat_f1_scores[
    concat_f1_scores["Data_Set"].isin(["Training", "Testing"])
]

ax = sns.barplot(x="Model", y="Weighted", hue="Data_Set", data=filtered_data)

# Add the weighted values above the bars
for p in ax.patches:
    height = p.get_height()
    if height > 0:
        ax.annotate(
            f"{height:.2f}",
            (p.get_x() + p.get_width() / 2.0, height),
            ha="center",
            va="center",
            xytext=(0, 10),
            textcoords="offset points",
        )

plt.title("Weighted F1 Score by Model Type")
plt.xlabel("Model Type")
plt.ylabel("Weighted F1 Score")
plt.legend(title="Data Set")
plt.tight_layout()

plt.savefig(f"{fig_dir}/f1_scores_plate4_downsample.png", dpi=500)

plt.show()


# ### Plot F1 scores for holdout data

# In[10]:


plt.figure(figsize=(10, 6))

# Set y-axis limits
plt.ylim(0, 1)

# Filter the data to include only hold-out sets
filtered_data = concat_f1_scores[
    concat_f1_scores["Data_Set"].isin(["Holdout1", "Holdout2"])
]

ax = sns.barplot(x="Model", y="Weighted", hue="Data_Set", data=filtered_data)

# Add the weighted values above the bars
for p in ax.patches:
    height = p.get_height()
    if height > 0:
        ax.annotate(
            f"{height:.2f}",
            (p.get_x() + p.get_width() / 2.0, height),
            ha="center",
            va="center",
            xytext=(0, 10),
            textcoords="offset points",
        )

plt.title("Weighted F1 Score by Model Type")
plt.xlabel("Model Type")
plt.ylabel("Weighted F1 Score")
plt.legend(title="Data Set")
plt.tight_layout()

plt.savefig(f"{fig_dir}/f1_scores_plate4_holdout_data_downsample.png", dpi=500)

plt.show()

