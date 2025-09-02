import pandas as pd
from pathlib import Path
import numpy as np
import random

print('🟢 Script started')

# Step 1: Define paths
base_path = Path('/cfs/klemming/projects/supr/snic2020-6-126/projects/envpredict/from_berzelius_nov_2024')
path_to_deep_micro_files = base_path / 'DeepMicro' / 'results'
print(f'📂 Looking for .csv files in: {path_to_deep_micro_files}')

# Step 2: Read all .csv files
all_files = []
file_names = []

csv_files = list(path_to_deep_micro_files.glob('*.csv'))
print(f'🔍 Found {len(csv_files)} CSV files.')

for file in csv_files:
    try:
        df = pd.read_csv(file, delimiter=',', header=None)
        all_files.append(df)
        file_names.append(file.name)
        print(f'✅ Loaded: {file.name} (shape: {df.shape})')
    except Exception as e:
        print(f'❌ Error reading {file.name}: {e}')

# Show a random file shape
if all_files:
    print("🎲 Example DataFrame Shape:", random.choice(all_files).shape)
else:
    print("⚠️ No files loaded. Exiting script.")
    exit()

# Step 3: Create output folder
path_to_save = base_path / 'RepresentationsFromDeepMicro_2025_04_01'
path_to_save.mkdir(exist_ok=True)
print(f'📁 Output directory: {path_to_save}')

# Step 4: Read sample names
print('📄 Reading sample names...')
names_path = base_path / 'code' / 'SampleNames16S.csv'

try:
    names = pd.read_csv(names_path, sep=',').iloc[:, 1]
    print(f'✅ Loaded sample names: {len(names)} entries')
except Exception as e:
    print(f'❌ Error reading sample names: {e}')
    exit()

# Step 5: Save processed files
print('💾 Saving processed files...')

for ind in range(len(all_files)):
    file = all_files[ind]
    file.index = names

    file_name = file_names[ind]
    file_name = file_name.replace("[", "_").replace("]", "_").replace("__", "_")
    file_name = file_name.replace("_pivoted_rep", "")

    colnames = [f'feature_{i+1}' for i in range(file.shape[1])]
    file.columns = colnames

    output_path = path_to_save / file_name
    file.to_csv(output_path, index=True, index_label='Sample ID', header=True, sep='\t')
    print(f'✅ Saved: {output_path.name} (shape: {file.shape})')

print('🏁 All done!')
