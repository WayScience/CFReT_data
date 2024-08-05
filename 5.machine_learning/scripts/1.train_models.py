#!/usr/bin/env python
# coding: utf-8

# ## Train machine learning models to predict failing or healthy cell status

# ## Import libraries

# In[1]:


import pathlib
import sys
import warnings

import numpy as np
import pandas as pd
from joblib import dump
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import RandomizedSearchCV, StratifiedKFold
from sklearn.preprocessing import LabelEncoder
from sklearn.utils import parallel_backend

sys.path.append("../utils")
from training_utils import downsample_data, get_X_y_data


# ## Set paths and variables

# In[2]:


# set numpy seed to make sure any random operations performs are reproducible
np.random.seed(0)

# path to training data set
training_data_path = pathlib.Path("./data/training_data.csv")

# Metadata column used for prediction class
label = "Metadata_cell_type"

# Directory for models to be outputted
model_dir = pathlib.Path("./models")
model_dir.mkdir(exist_ok=True, parents=True)

# Directory for label encoder
encoder_dir = pathlib.Path("./encoder_results")
encoder_dir.mkdir(exist_ok=True, parents=True)


# ## Load in training data

# In[3]:


df = pd.read_csv(training_data_path)

print(df.shape)
df.head()


# ## Perform downsampling on training data and output as data frame

# In[4]:


# load in training plate 4 data as downsampled to lowest class
downsample_df = downsample_data(path_to_data=training_data_path, label=label)

print(downsample_df.shape)
print(downsample_df["Metadata_cell_type"].value_counts())
downsample_df.head()


# ## Get X and y data for both final and shuffled models

# In[5]:


# Get not shuffled training data from downsampled df (e.g., "final")
X_train, y_train = get_X_y_data(df=downsample_df, label=label, shuffle=False)

# Get shuffled training data from downsampled df(e.g., "shuffled_baseline")
X_shuffled_train, y_shuffled_train = get_X_y_data(
    df=downsample_df, label=label, shuffle=True
)


# ## Encode labels in both shuffled and non-shuffled
# 
# **Note:** Failing will be considered as 0 and Healthy will be 1.

# In[6]:


# Encode classes
le = LabelEncoder()
le.fit(y_train)
# Fit the labels onto the shuffled and non-shuffled data
y_train = le.transform(y_train)
y_shuffled_train = le.transform(y_shuffled_train)

# Print the original classes and their corresponding encoded values
class_mapping = dict(zip(le.classes_, le.transform(le.classes_)))
print("Class Mapping:")
print(class_mapping)


# ## Train the models
# 
# **Note:** We will be using RandomizedSearchCV to hyperparameterize the model since we have a larger dataset and it will be easier to try random combinations than all combinations.

# ### Set up the model and hyper parameter method

# In[7]:


# Set folds for k-fold cross validation (default is 5)
straified_k_folds = StratifiedKFold(n_splits=10, shuffle=False)

# Set Logistic Regression model parameters (use default for max_iter)
logreg_params = {
    "penalty": "elasticnet",
    "solver": "saga",
    "max_iter": 1000,
    "n_jobs": -1,
    "random_state": 0,
    "class_weight": "balanced",
}

# Define the hyperparameter search space for RandomizedSearchCV
param_dist = {
    "C": np.logspace(-3, 3, 7),
    "l1_ratio": np.linspace(0, 1, 11),
}

# Set the random search hyperparameterization method parameters (used default for "cv" and "n_iter" parameter)
random_search_params = {
    "param_distributions": param_dist,
    "scoring": "f1_weighted",
    "random_state": 0,
    "n_jobs": -1,
    "cv": straified_k_folds,
}


# ### Train final model

# In[8]:


# Check if the "models" folder contains a file with "final" in its name
if any(model_dir.glob("*final*")):
    print("Model training skipped as a 'final' model already exists.")
else:
    # Generate logistic regression model for non-shuffled training data
    final_logreg = LogisticRegression(**logreg_params)

    # Initialize the RandomizedSearchCV
    final_random_search = RandomizedSearchCV(final_logreg, **random_search_params)

    # Prevent the convergence warning in sklearn, it does not impact the result
    with parallel_backend("multiprocessing"):
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", category=ConvergenceWarning, module="sklearn"
            )
            # Perform the random hyperparameter search
            final_random_search.fit(X_train, y_train)

# Print the best parameters
print("Optimal parameters:", final_random_search.best_params_)


# ### Train shuffled baseline model

# In[9]:


# Check if the "models" folder contains a file with "final" in its name
if any(model_dir.glob("*shuffled*")):
    print("Model training skipped as a 'shuffled' model already exists.")
else:
    # Generate logistic regression model for shuffled training data
    shuffled_logreg = LogisticRegression(**logreg_params)

    # Initialize the RandomizedSearchCV
    shuffled_random_search = RandomizedSearchCV(shuffled_logreg, **random_search_params)

    # Prevent the convergence warning in sklearn, it does not impact the result
    with parallel_backend("multiprocessing"):
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", category=ConvergenceWarning, module="sklearn"
            )
            # Perform the random hyperparameter search
            shuffled_random_search.fit(X_shuffled_train, y_shuffled_train)


# ## Save the models and label encoder

# In[10]:


data_prefix = "log_reg_fs_plate_4"

# Check if there are models with "final" or "shuffled" in its name that exists in the models folder
if any(model_dir.glob("*final*")) or any(model_dir.glob("*shuffled*")):
    print(
        "No models were generated or saved because 'final' and/or 'shuffled' files already exist."
    )
else:
    # Save the models
    dump(
        final_random_search.best_estimator_,
        f"{model_dir}/{data_prefix}_final_downsample.joblib",
    )
    dump(
        shuffled_random_search.best_estimator_,
        f"{model_dir}/{data_prefix}_shuffled_downsample.joblib",
    )

    # Save label encoder
    dump(le, f"{encoder_dir}/label_encoder_{data_prefix}.joblib")

print("Models and label encoder have been saved!")

