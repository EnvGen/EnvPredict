#!/usr/bin/env python
# coding: utf-8

# ## TabPFN Predictions for the EnvPredict project
# 
# Let's load the packages first

import pandas as pd
import numpy as np
from sklearn.model_selection import KFold
from tabpfn import TabPFNRegressor
import os
from datetime import datetime
import math
from scipy.stats import pearsonr, spearmanr
from sklearn.datasets import fetch_openml
from sklearn.model_selection import train_test_split
from joblib import Parallel, delayed

# Define functions

# Function to use alternative sample names
def use_alternative_sample_names(matrix, alt_ids):
    ns = alt_ids.set_index('sample_id').loc[matrix.columns, 'station_id_date'].values
    matrix.columns = ns
    return matrix

# Function to extract shared samples
def extract_shared_samples(features_matrix_full, responses_matrix_full):
    shared_samples = responses_matrix_full.columns.intersection(features_matrix_full.columns)
    if any('NA_' in sample for sample in shared_samples):
        shared_samples = shared_samples[~shared_samples.str.contains('NA_')]
    features_matrix = features_matrix_full.loc[:, shared_samples]
    responses_matrix = responses_matrix_full.loc[:, shared_samples]
    responses_matrix_bin = responses_matrix.applymap(lambda x: 1 if x > 0 else 0)
    responses_matrix = responses_matrix.loc[responses_matrix_bin.sum(axis=1) > 0, :]
    return features_matrix, responses_matrix

# Function to perform feature selection
def do_feature_selection(features_matrix, occupancy):
    binary_features_matrix = features_matrix.map(lambda x: 1 if x > 0 else 0)
    ix = binary_features_matrix.sum(axis=1) > occupancy * features_matrix.shape[1]
    return features_matrix.loc[ix, :]

def process_fold(train_index, test_index, df, min_samples):
    train_data, test_data = df.iloc[train_index, :], df.iloc[test_index, :]
    train_data = train_data.dropna(subset=[0])
    if train_data.shape[0] >= min_samples:
        model = TabPFNRegressor(device='auto', ignore_pretraining_limits=True)
        model.fit(train_data.iloc[:, 1:], train_data.iloc[:, 0])
        return test_index, model.predict(test_data.iloc[:, 1:])
    return test_index, None

def run_tabpfn(features_matrix, responses_matrix, numfolds, min_samples):
    predicted_responses_matrix = responses_matrix.copy()
    predicted_responses_matrix[:] = np.nan
    kf = KFold(n_splits=numfolds)
    
    # Fetch the number of available cores
    n_jobs = min(os.cpu_count(), numfolds)
    
    for i in range(responses_matrix.shape[0]):
        response = responses_matrix.iloc[i, :]
        if response.dropna().shape[0] < min_samples:
            continue
        df = pd.DataFrame(np.column_stack((response, features_matrix.T)))
        
        results = Parallel(n_jobs=n_jobs, backend="loky")(delayed(process_fold)(train_index, test_index, df, min_samples) 
                                                          for train_index, test_index in kf.split(df))
        
        for test_index, prediction in results:
            if prediction is not None:
                predicted_responses_matrix.iloc[i, test_index] = prediction
    
    return predicted_responses_matrix

def process_file(infile, alt_ids, RF_output_files_path, output_files_path, num_folds=10, min_samples=10, occupancy=0.1):
    outfile_actual = os.path.basename(infile).replace(".tsv", "_RF10fold_Actual.tsv")
    outfile_predicted = os.path.basename(infile).replace(".tsv", "_RF10fold_Predictions.tsv")

    # Read the responses matrix from the actual file
    responses_matrix_full = pd.read_csv(os.path.join(RF_output_files_path, outfile_actual), delimiter="\t", index_col=0)
    
    # Read the features matrix
    features_matrix_full = pd.read_csv(infile, delimiter="\t", index_col=0)

    # Use alternative sample names
    features_matrix_full = use_alternative_sample_names(features_matrix_full, alt_ids)

    # Extract shared samples
    features_matrix, responses_matrix = extract_shared_samples(features_matrix_full, responses_matrix_full)

    # Perform feature selection
    features_matrix = do_feature_selection(features_matrix, occupancy)

    # Run TabPFN
    predicted_responses_matrix = run_tabpfn(features_matrix, responses_matrix, num_folds, min_samples)
    
    # Write the predicted responses matrix to a file
    predicted_responses_matrix.to_csv(os.path.join(output_files_path, outfile_predicted), sep="\t")

# Define the output directory and create it if it doesn't exist
output_files_path = "../output/TabPFN"
os.makedirs(output_files_path, exist_ok=True)

## Define the RF output file folder to get the response variables cut and in proper format
RF_output_files_path = "../output/DifferentTaxonomicLevels"

# Define the input directories
features_files_path_16S = "../seq_files/16S"
features_files_path_18S = "../seq_files/18S"

# Get the list of input files
infiles = sorted(
    [os.path.join(features_files_path_16S, f) for f in os.listdir(features_files_path_16S) if f.startswith("norm_seqtab") and f.endswith(".tsv")] +
    [os.path.join(features_files_path_18S, f) for f in os.listdir(features_files_path_18S) if f.startswith("norm_seqtab") and f.endswith(".tsv")]
)
# infiles = [f for f in infiles if not f.endswith("_1.tsv")]

# Load the phys_chem data
phys_chem = pd.read_csv("../env_files/physical_chemical_processed_translation.tsv", delimiter="\t")

# Load the alternative sample IDs
alt_ids = pd.DataFrame({'sample_id': phys_chem['sample_id'], 'station_id_date': phys_chem['station_id_date']})

## Call the model on irrelevant to secure checkopoints (common bug fix)

df = fetch_openml(data_id=531, as_frame=True)  # Boston Housing dataset
X = df.data
y = df.target.astype(float)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.5, random_state=42)
regressor = TabPFNRegressor()
regressor.fit(X_train, y_train)

## Parallelize the processing of input files
Parallel(n_jobs=2)(delayed(process_file)(infile, alt_ids, RF_output_files_path, output_files_path) for infile in infiles)