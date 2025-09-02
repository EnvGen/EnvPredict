import pandas as pd

# Load 16S data
with open('/cfs/klemming/projects/supr/snic2020-6-126/projects/envpredict/EnvPredict/seq_data_march_2025/combined/16S/norm_seqtab_16S.tsv', 'r') as tsvfile:
    seqtab16S = pd.read_csv(tsvfile, delimiter='\t')  # Assuming first column is ASV IDs

# Load 18S data
with open('/cfs/klemming/projects/supr/snic2020-6-126/projects/envpredict/EnvPredict/seq_data_march_2025/combined/18S/norm_seqtab_18S.tsv', 'r') as tsvfile:
    seqtab18S = pd.read_csv(tsvfile, delimiter='\t')  # Assuming first column is ASV IDs

print('Files read')

print("16S DataFrame Shape:", seqtab16S.shape)
print("18S DataFrame Shape:", seqtab18S.shape)

# 16S filtering
# Convert to presence/absence (1 if present, 0 otherwise)
presence16S = (seqtab16S > 0).astype(int)

# Calculate the proportion of columns (samples) where each ASV (row) is present
asv_presence_ratio16S = presence16S.sum(axis=1) / seqtab16S.shape[1]

# Filter rows where the presence is at least 10% of columns
filtered_seqtab16S = seqtab16S.loc[asv_presence_ratio16S >= 0.1, :]

print("16S DataFrame Shape after filtering:", filtered_seqtab16S.shape)

# 18S filtering
presence18S = (seqtab18S > 0).astype(int)
asv_presence_ratio18S = presence18S.sum(axis=1) / seqtab18S.shape[1]
filtered_seqtab18S = seqtab18S.loc[asv_presence_ratio18S >= 0.1, :]

print("18S DataFrame Shape after filtering:", filtered_seqtab18S.shape)

# Save filtered data
filtered_seqtab16S.to_csv('/cfs/klemming/projects/supr/snic2020-6-126/projects/envpredict/EnvPredict/seq_data_march_2025/combined/16S/filtered_10_percent_norm_seqtab_16S.tsv', sep='\t')
print('Filtered 16S file saved')

filtered_seqtab18S.to_csv('/cfs/klemming/projects/supr/snic2020-6-126/projects/envpredict/EnvPredict/seq_data_march_2025/combined/18S/filtered_10_percent_norm_seqtab_18S.tsv', sep='\t')
print('Filtered 18S file saved')

# Pivot the DataFrames if required
pivot_filtered_seqtab16S = filtered_seqtab16S.T
pivot_filtered_seqtab18S = filtered_seqtab18S.T


# pivoted dimensions

print("16S DataFrame Shape after pivoting:", pivot_filtered_seqtab16S.shape)
print("18S DataFrame Shape after pivoting:", pivot_filtered_seqtab18S.shape)

# Save the pivoted files
pivot_filtered_seqtab16S.to_csv('filtered_norm_seqtab_16S_pivoted_march_2025.csv', index=False, header=False)
print('Pivoted 16S file saved')

pivot_filtered_seqtab18S.to_csv('filtered_norm_seqtab_18S_pivoted_march_2025.csv', index=False, header=False)
print('Pivoted 18S file saved')
