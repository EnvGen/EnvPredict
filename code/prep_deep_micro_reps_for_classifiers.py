import pandas as pd
from pathlib import Path
import numpy as np
import random

# step 1: read files
base_path = Path('/proj/berzelius-2023-48/envpredict')
path_to_deep_micro_files = base_path / 'DeepMicro' / 'results'

# loop through all .csv files
all_files = []
file_names = []
for file in path_to_deep_micro_files.glob('*.csv'):
    all_files.append(pd.read_csv(file, delimiter=',',header=None))
    file_names.append(file.name)

print("Random DataFrame Shape:", random.choice(all_files).shape)


# create folder
path_to_save = base_path / 'RepresentationsFromDeepMicro'
path_to_save.mkdir(exist_ok=True)

# read sample names
names = pd.read_csv(base_path / 'code' / 'SampleNames16S.csv', sep = ',')
names = names.iloc[:,1]

# Save as comma-separated files in new folder
for ind in np.arange(len(all_files)):
    file = all_files[ind]
    file.index = names
    
    file_name = file_names[ind]
    file_name = file_name.replace("[", "_")
    file_name = file_name.replace("]", "_")
    file_name = file_name.replace("__", "_")

    file_name = file_name.replace("_pivoted_rep", "")

    colnames = ['feature_' + str(i+1) for i in np.arange(file.shape[1])]
    file.columns = colnames

    file.to_csv(path_to_save / file_name, index=True, index_label='Sample ID', header=True, sep = '\t')
    print(f'{file_name} file saved')
