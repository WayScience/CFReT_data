#!/usr/bin/env python
# coding: utf-8

# # Extract logistic regression model coefficients per CellProfiler feature
# 
# The coefficients from the machine learning model will either be positive (x > 0), where if the feature value increases, the more likely the feature is related to the Healthy cell type, negative (x < 0), where if the feature value increases, the more likely the feature is the Failing cell type, and zero (x = 0), where that feature has no impact on predicting if a cell is Healthy or Failing.

# ## Import libraries

# In[1]:


from joblib import load
import pathlib
import pandas as pd


# ## Load in the training data to collect the CellProfiler feature columns

# In[2]:


# Path to data dir
data_dir = pathlib.Path("./data/")

# Output dir for coefficients
coeff_dir = pathlib.Path("./coeff_data/")
coeff_dir.mkdir(exist_ok=True)

# path to training data to access the feature columns
path_to_training_data = pathlib.Path(f"{data_dir}/training_data.csv")

# Get all the column names from the training data
all_columns = pd.read_csv(path_to_training_data, nrows=0).columns

# Filter columns that start with 'Metadata_' to only get feature columns from CellProfiler
feature_columns = [col for col in all_columns if not col.startswith('Metadata_')]

print(len(feature_columns))
feature_columns
print(f"Examples of first five feature columns:")
for example_column in feature_columns[:5]:
    print(example_column)


# ## Set paths and load in the final model

# In[3]:


# path to the final model joblib file
path_to_final_model = pathlib.Path("./models/log_reg_fs_plate_4_final_downsample.joblib").resolve(strict=True)

# load in final model
final_model = load(path_to_final_model)


# ## Collect coefficients from the model and match with the correct feature in a dataframe

# In[4]:


# Get the coefficients
coefficients = final_model.coef_

# Print the coefficients shape and confirm it is the same number as feature columns from training data
print(coefficients.shape)
# Confirm it is the same number as feature columns from training data
if coefficients.shape[1] == len(feature_columns):
    print("The number of coefficients matches the number of feature columns.")
else:
    print("Warning: The number of coefficients does not match the number of feature columns.")

# Create a DataFrame with the coefficients and features
coefficients_df = pd.DataFrame({'Feature': feature_columns, 'Coefficient': coefficients.flatten()})

# Print the DataFrame
coefficients_df.head()


# ## Split the data frame by negative, positive. and zero coefficients which relate to different class importance

# In[5]:


# Split into DataFrames with positive, negative, and zero coefficients
positive_coeffs_df = coefficients_df[coefficients_df['Coefficient'] > 0].copy()
negative_coeffs_df = coefficients_df[coefficients_df['Coefficient'] < 0].copy()
zero_coeffs_df = coefficients_df[coefficients_df['Coefficient'] == 0].copy()

# Make the values in negative_coeffs_df positive to prevent issues during plotting
negative_coeffs_df['Coefficient'] = abs(negative_coeffs_df['Coefficient'])

# Rename the columns
positive_coeffs_df.columns = ['Feature', 'Healthy_Coeffs']
negative_coeffs_df.columns = ['Feature', 'Failing_Coeffs']
zero_coeffs_df.columns = ['Feature', 'Zero_Coeffs']

# Save the coef data into the "/coeff_data" folder
positive_coeffs_df.to_csv(f'{coeff_dir}/healthy_coeffs.csv', index=False)
negative_coeffs_df.to_csv(f'{coeff_dir}/failing_coeffs.csv', index=False)
zero_coeffs_df.to_csv(f'{coeff_dir}/zero_coeffs.csv', index=False)


# Print or use the resulting DataFrames
print("Positive Coefficients:", positive_coeffs_df.shape[0])
print("\nNegative Coefficients:", negative_coeffs_df.shape[0])
print("\nZero Coefficients:", zero_coeffs_df.shape[0])
negative_coeffs_df.head()


# ## Explore the coefficients

# In[6]:


# Find the row with the highest coefficient value
max_row = coefficients_df.loc[coefficients_df['Coefficient'].idxmax()]

# Extract the feature and coefficient values
max_feature = max_row['Feature']
max_coefficient_value = max_row['Coefficient']

# Print or use the result
print("Feature with the highest coefficient:", max_feature)
print("Coefficient value:", max_coefficient_value)


# In[7]:


# Sort the DataFrame based on the coefficient values (from most positive to most negative)
coefficients_healthy_df = coefficients_df.sort_values(by='Coefficient', ascending=False)

# Show the top ten ranking features for predicting "Healthy" class
coefficients_healthy_df.head(10)


# In[8]:


# Find the row with the most negative coefficient value
min_row = coefficients_df.loc[coefficients_df['Coefficient'].idxmin()]

# Extract the feature and coefficient values
min_feature = min_row['Feature']
min_coefficient_value = min_row['Coefficient']

# Print or use the result
print("Feature with the most negative coefficient:", min_feature)
print("Coefficient value:", min_coefficient_value)


# In[9]:


# Sort the DataFrame based on the coefficient values (from most negative to most positive)
coefficients_failing_df = coefficients_df.sort_values(by='Coefficient', ascending=True)

# Show the top ten ranking features for predicting "Failing" class
coefficients_failing_df.head(10)


# ## Add ranking column with sorted descending values and save the CSV for visualization
# 
# Rank is based on the highest positive coefficient which will have rank one and then descending from there. We expect to see that the model will take into account many different features (positive and negative which relate to different classes) and there will be many features at zero meaning they are redundant to the model.

# In[10]:


# Sort coefficients_df by descending order
coefficients_df = coefficients_df.sort_values(by='Coefficient', ascending=False)

# Add a new column 'Rank'
coefficients_df['Rank'] = range(1, len(coefficients_df) + 1)

# Save the ranked df
coefficients_df.to_csv(f'{coeff_dir}/ranked_coeffs.csv', index=False)

# Show df to assess if the ranking was performed correctly
print(coefficients_df.shape)
coefficients_df.head()

