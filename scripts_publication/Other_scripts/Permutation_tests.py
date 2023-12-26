#!/usr/bin/env python
# coding: utf-8

# In[3]:


import pandas as pd
import numpy as np
import random


# In[4]:


random.seed(18671107)


# In[5]:


# Load your data from the CSV file
data = pd.read_csv("#path_umap.txt", sep='\t')


# In[6]:
print(data)


# In[7]:


# Define the observed test statistic
observed_statistic = data.groupby("group")[["x", "y"]].mean().diff().abs().sum().sum()


# In[8]:


# Number of permutations (adjust if needed)
num_permutations = 1000


# In[6]:


# Initialize an array to store permuted statistics
permuted_statistics = []

# Perform the permutation test
for _ in range(num_permutations):
    # Shuffle the group labels
    shuffled_data = data.copy()
    shuffled_data["group"] = np.random.permutation(shuffled_data["group"])

    # Calculate the statistic on the shuffled data
    shuffled_statistic = shuffled_data.groupby("group")[["x", "y"]].mean().diff().abs().sum().sum()

    permuted_statistics.append(shuffled_statistic)

# Calculate the p-value
p_value = (np.sum(np.array(permuted_statistics) >= observed_statistic) + 1) / (num_permutations + 1)

print(f"Observed Statistic: {observed_statistic}")
print(f"P-Value: {p_value}")
