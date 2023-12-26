#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Load libraries
import pandas as pd
import numpy as np
import os
from sklearn.impute import KNNImputer
from sklearn.impute import SimpleImputer
# from natsort import natsorted, ns
from datetime import datetime
import matplotlib.pyplot as plt
from pybedtools import BedTool
import multiprocess as mp
from functools import partial
import random
from tqdm import tqdm
# Set random seed for reproducibility
random.seed(18671107)
np.random.seed(18671107)


# In[2]:


# Define basic function for loading and organising data
def binning_bed(cellname, bedfileDF, windows):
    bedfile = bedfileDF[bedfileDF['Cell'] == cellname] # get mini bed (= of each Cell)
    bedfile = BedTool.from_dataframe(bedfile).sort() # sort bed
    mapped = BedTool.map(BedTool(), o='median', f=0.5, c=4, a=windows, b=bedfile, bed = True) # overlap CN with windows
    mapped = BedTool.to_dataframe(mapped) # back to pandas dataframe
    mapped['Cell'] = cellname # regive Cell column
    return(mapped)

def import_file(input_file, matrix, sep):
    if matrix == True:
        print('Importing file as Matrix. . .')
        BED = pd.read_csv(input_file, sep=sep, na_values='.', header=0, index_col = 0)
        BED.index.names = ['Cell']
        BED = BED.stack().reset_index(name='copy_number').rename(columns={'level_1':'pos'})
        BED[['chr', 'start']] = BED.pos.str.split(":", expand = True)
        BED[['start', 'end']] = BED.start.str.split("-", expand = True)
        BED = BED[['chr','start', 'end','copy_number','Cell']]
    else:
        print('Importing BED file. . .')
        BED = pd.read_csv(input_file, sep=sep, na_values='.', header=0)[['chr','start', 'end','copy_number','Cell']]
    return(BED)


# In[3]:


# BED scCNV file to load (MCF-7 Normal population from Gnan et al. 2022)
file_to_load = '/Volumes/SC_Data/Kronos_data_analyses/CNV/20kb_2023/MCF7_Normal_cnv_calls.bed'
# Define output variables
output_dir = 'data_tmp/'
data_name = 'imputation-test_dta_2023'
file_name = data_name+'_scMatrix'+'.tsv.gz'
file_path = output_dir+'/'+file_name
file_path = file_path.replace("//", "/")
Reference_genome = 'hg38'
window_size = 100000
sep = '\t'
# Make output directories
os.makedirs(output_dir,exist_ok=True)


# In[4]:


# Make windows from reference genome and sort
print('Binning the reference genome ('+Reference_genome+'). . .')
start=datetime.now()
if len (Reference_genome) > 5:
    windows = BedTool.window_maker(BedTool(), g=Reference_genome, w=window_size)
else:
    windows = BedTool.window_maker(BedTool(), genome=Reference_genome, w=window_size)
windows = windows.sort()
print(datetime.now()-start)


# In[5]:


# Load BED File
start=datetime.now()
bedfileDF = import_file(file_to_load, False, sep)
bedfile = BedTool.from_dataframe(bedfileDF).sort()
print(datetime.now()-start)


# In[6]:


# Find unique Cell names
Cells = bedfileDF["Cell"].unique()


# In[7]:


# Left join Windows + BED file copy-numbers by median per Cell.
print('Binning BED file using '+str(6)+' cores . . .')
start=datetime.now()
pool = mp.Pool(6)
binning_bed_1=partial(binning_bed, bedfileDF=bedfileDF)
binning_bed_2=partial(binning_bed_1, windows=windows)
CN_data = pd.DataFrame()
with mp.Pool() as pool:
    results = pool.map(binning_bed_2, Cells)
CN_data = pd.concat(results)
pool.close()

# Merged data from all cells : replace for correct Nan values.
CN_data.rename(columns={"name": "CN"}, inplace = True)
CN_data.replace({'CN': "."}, np.nan, inplace = True)
print(datetime.now()-start)


# In[8]:


# Create unique position column (chromosome + start + end)
CN_data["pos"] = (CN_data["chrom"] + ':' + CN_data["start"].astype('str') + '-' + CN_data["end"].astype('str'))
all_regions = CN_data["pos"]


# In[9]:


# Convert data to matrix
print('Creating matrix. . .')
start=datetime.now()
CN_data["Cell"] = CN_data["Cell"].astype("category")
CN_data["pos"] = CN_data["pos"].astype("category")
CN_data = CN_data.pivot(index="Cell", columns="pos", values="CN")
CN_data = CN_data.dropna(axis=1, how='all')
print(datetime.now()-start)


# In[10]:


# Save cell names for later use
Cells = CN_data.index


# In[11]:


# Check matrix size
CN_data.shape


# In[12]:


# Drop any regions that contain NaNs
CN_data = CN_data.dropna(axis=1, how='any')


# In[13]:


# Check matrix size after removing those regions
CN_data.shape


# In[14]:


# Verify that there are no more NaN values
print(CN_data.isna().sum().sum())


# In[15]:


# Create a backup of the data
bckup = CN_data


# In[16]:


#Set the missing value percentages
missing_values_percentage = range(5,56,5)


# In[17]:


# Define random imputation function
def random_imputation(df):
    non_nan_values = df.stack().dropna().values

    total_nan_count = df.isna().sum().sum()
    progress_bar = tqdm(total=total_nan_count, desc="Imputing NaN values")

    while total_nan_count > 0:
        nan_rows, nan_cols = np.where(df.isna())

        for i in range(len(nan_rows)):
            random_value = np.random.choice(non_nan_values)
            df.iloc[nan_rows[i], nan_cols[i]] = random_value
            progress_bar.update(1)

        total_nan_count = df.isna().sum().sum()

    progress_bar.close()
    return df


# In[18]:


# Define accuracy and invariance function
def imp_accuracy_inv (Original_matrix, Imputed_matrix, removed_indices):
    accuracy = np.sum(Original_matrix.values.flatten()[random_indices] == Imputed_matrix.values.flatten()[random_indices])
    accuracy = (accuracy / num_elements) * 100
    invariance = np.sum(Original_matrix.values == Imputed_matrix.values)
    invariance = (invariance / total_elements) * 100
    return (accuracy, invariance)


# In[19]:


# Create accuracy lists
KNN_accuracy = []
KNN_accuracy_1 = []
RANDOM_accuracy = []
MEDIAN_accuracy = []
# Create invariance lists
KNN_invariance = []
KNN_1_invariance = []
RANDOM_invariance = []
MEDIAN_invariance = []

# Calculate the total number of elements in the matrix
total_elements = bckup.size
#Perform imputations and save accuracy and invariance rates
for i in missing_values_percentage:
    start=datetime.now()
    print('Imputation after removing '+str(i)+'% of values.')
    CN_data = bckup.copy()
    # Determine the number of elements to replace
    num_elements = int(CN_data.size * i/100)
    # Randomly select x% of the values and replace them with NaNs
    random_indices = np.random.choice(CN_data.size, size=num_elements, replace=False)
    matrix_flattened = CN_data.values.flatten()
    matrix_flattened[random_indices] = np.nan
    matrix_with_nans = pd.DataFrame(matrix_flattened.reshape(CN_data.shape), columns=CN_data.columns)
    # Force integer values
    CN_data = CN_data.astype(int)
    # KNN Imputation
    imputer = KNNImputer(n_neighbors=10, weights='distance')
    imp_matrix = pd.DataFrame(imputer.fit_transform(matrix_with_nans.copy()), columns = CN_data.columns, index = CN_data.index)
    imp_matrix = imp_matrix.round(0).astype(int)
    # Median imputation
    imp = SimpleImputer(strategy='median')
    median_imp_matrix = pd.DataFrame(imp.fit_transform(matrix_with_nans.copy()), columns = CN_data.columns, index = CN_data.index)
    median_imp_matrix = median_imp_matrix.round(0).astype(int)
    # Random imputation
    random_matrix = random_imputation(matrix_with_nans.copy())  # Perform random imputation
    random_matrix = random_matrix.round(0).astype(int)
    # Calculate the number of exact matches with the KNN imputation
    KNN_accuracy_df, KNN_invariance_df = imp_accuracy_inv(CN_data, imp_matrix, random_indices)
    # Calculate the number of equivalent values (with a CN difference of 1 or less)
    KNN_1_accuracy_df = np.sum(np.abs(CN_data.values.flatten()[random_indices] - imp_matrix.values.flatten()[random_indices]) <= 1)
    KNN_1_accuracy_df = (KNN_1_accuracy_df / num_elements) * 100
    KNN_1_invariance_df = np.sum(np.abs(CN_data.values - imp_matrix.values) <= 1)
    KNN_1_invariance_df = (KNN_1_invariance_df / total_elements) * 100
    # Calculate the number of matches with the random imputation
    random_accuracy_df, random_invariance_df = imp_accuracy_inv(CN_data, random_matrix, random_indices)
    # Calculate the number of matches wit the median imputation
    median_accuracy_df, median_invariance_df = imp_accuracy_inv(CN_data, median_imp_matrix, random_indices)
    # Add accuracy to list
    KNN_accuracy.append(KNN_accuracy_df)
    KNN_accuracy_1.append(KNN_1_accuracy_df)
    RANDOM_accuracy.append(random_accuracy_df)
    MEDIAN_accuracy.append(median_accuracy_df)
    # Add invariance rate to list
    KNN_invariance.append(KNN_invariance_df)
    KNN_1_invariance.append(KNN_1_invariance_df)
    RANDOM_invariance.append(random_invariance_df)
    MEDIAN_invariance.append(median_invariance_df)
    # Print accuracy
    print(f"KNN accuracy: {KNN_accuracy_df}%")
    print(f"KNN +/- 1 accuracy: {KNN_1_accuracy_df}%")
    print(f"Random accuracy: {random_accuracy_df}%")
    print(f"Median accuracy: {median_accuracy_df}%")
    # Print runtime
    print(datetime.now()-start)


# In[20]:


# Averages
print('KNN mean: ', np.mean(KNN_accuracy) )
print('Median mean: ', np.mean(MEDIAN_accuracy) )
print('Random mean: ', np.mean(RANDOM_accuracy) )


# In[21]:


# Print Results
print(KNN_accuracy,'\n',KNN_accuracy_1,'\n',RANDOM_accuracy,'\n',MEDIAN_accuracy)
# Create invariance lists
print(' ')
print(KNN_invariance,'\n',KNN_1_invariance,'\n',RANDOM_invariance,'\n',MEDIAN_invariance)


# In[22]:


from scipy import stats
# Perform paired t-test KNN vs Random
t_statistic, p_value = stats.ttest_rel(KNN_accuracy, RANDOM_accuracy)
print("Paired t-test results Knn vs Random:")
print(f"T-statistic: {t_statistic}")
print(f"P-value: {p_value}")


# In[23]:


from scipy import stats
# Perform paired t-test KNN vs Median
t_statistic, p_value = stats.ttest_rel(KNN_accuracy, MEDIAN_accuracy)
print("Paired t-test results Knn vs Median:")
print(f"T-statistic: {t_statistic}")
print(f"P-value: {p_value}")


# In[24]:


from scipy import stats
# Perform paired t-test Random vs Median
t_statistic, p_value = stats.ttest_rel(RANDOM_accuracy, MEDIAN_accuracy)
print("Paired t-test results Random vs Median:")
print(f"T-statistic: {t_statistic}")
print(f"P-value: {p_value}")


# In[25]:


# Import colorblind-friendly colors
from matplotlib.cm import cividis
# Create the plot
plt.figure(figsize=(5, 4), dpi=200)
plt.plot(list(missing_values_percentage), KNN_accuracy, marker='o', color=cividis(0.99), label='KNN accuracy')
plt.plot(list(missing_values_percentage), KNN_accuracy_1, marker='o', color=cividis(0.7), label='KNN ±1 similarity')
plt.plot(list(missing_values_percentage), MEDIAN_accuracy, marker='o', color=cividis(0.5), label='Median accuracy')
plt.plot(list(missing_values_percentage), RANDOM_accuracy, marker='o', color=cividis(0.0), label='Random accuracy')
plt.xlabel('Percentage of removed values')
plt.ylabel('Accuracy rate of imputation method')
plt.legend()
plt.grid(True)
plt.show()


# In[26]:


# Import colorblind-friendly colors
from matplotlib.cm import cividis
# Create the plot
plt.figure(figsize=(5, 4), dpi=200)
plt.plot(list(missing_values_percentage), KNN_invariance, marker='o', color=cividis(0.99), label='KNN invariance')
plt.plot(list(missing_values_percentage), KNN_1_invariance, marker='o', color=cividis(0.7), label='KNN ±1 invariance')
plt.plot(list(missing_values_percentage), MEDIAN_invariance, marker='o', color=cividis(0.5), label='Median invariance')
plt.plot(list(missing_values_percentage), RANDOM_invariance, marker='o', color=cividis(0.0), label='Random invariance')
plt.xlabel('Percentage of removed values')
plt.ylabel('Percentage of matching copy-numbers')
plt.legend()
plt.grid(True)
plt.show()


# In[27]:


# Set different random seeds for reproducibility
random_seeds = [2023, 1911, 190312]

# Make a list to include all KNN accuracies, with different seeds
ALL_KNN_acc = []
# Add the 1st KNN with the original random seed
ALL_KNN_acc.append(KNN_accuracy)

for seed in random_seeds:
    start=datetime.now()
    print(seed)
    # Set random seed
    np.random.seed(seed)  
    random.seed(seed)
    # Create accuracy lists
    KNN_accuracy = []
    # Create invariance lists
    KNN_invariance = []
    # Calculate the total number of elements in the matrix
    total_elements = bckup.size
    #Perform imputations and save accuracy and invariance rates
    for i in missing_values_percentage:
        print('Imputation after removing '+str(i)+'% of values.')
        CN_data = bckup.copy()
        # Determine the number of elements to replace
        num_elements = int(CN_data.size * i/100)
        # Randomly select x% of the values and replace them with NaNs
        random_indices = np.random.choice(CN_data.size, size=num_elements, replace=False)
        matrix_flattened = CN_data.values.flatten()
        matrix_flattened[random_indices] = np.nan
        matrix_with_nans = pd.DataFrame(matrix_flattened.reshape(CN_data.shape), columns=CN_data.columns)
        # Force integer values
        CN_data = CN_data.astype(int)
        # KNN Imputation
        imputer = KNNImputer(n_neighbors=10, weights='distance')
        imp_matrix = pd.DataFrame(imputer.fit_transform(matrix_with_nans.copy()), columns = CN_data.columns, index = CN_data.index)
        imp_matrix = imp_matrix.round(0).astype(int)
        # Calculate the number of exact matches with the KNN imputation
        KNN_accuracy_df, KNN_invariance_df = imp_accuracy_inv(CN_data, imp_matrix, random_indices)
        # Add accuracy to list
        KNN_accuracy.append(KNN_accuracy_df)
        # Add invariance rate to list
        KNN_invariance.append(KNN_invariance_df)
    # Print runtime
    print(datetime.now()-start)
    ALL_KNN_acc.append(KNN_accuracy)


# In[28]:


# Calculate standard deviation between corresponding positions in lists across seeds
std_dev_between_seeds = np.std(ALL_KNN_acc, axis=0)

# Print standard deviation between lists across seeds
print(f"Standard deviation between lists across seeds:\n{std_dev_between_seeds}")


# In[29]:


np.mean(std_dev_between_seeds)


# In[30]:


min(std_dev_between_seeds)


# In[31]:


max(std_dev_between_seeds)


# In[32]:


np.median(std_dev_between_seeds)


# In[ ]:




