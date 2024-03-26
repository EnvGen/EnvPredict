import pandas as pd

with open('/proj/berzelius-2023-48/envpredict/data/norm_seqtab_16S.tsv', 'r') as tsvfile:
    # Create a pandas DataFrame from the TSV file
    seqtab16S = pd.read_csv(tsvfile, delimiter='\t')

with open('/proj/berzelius-2023-48/envpredict/data/norm_seqtab_18S.tsv', 'r') as tsvfile:
    # Create a pandas DataFrame from the TSV file
    seqtab18S = pd.read_csv(tsvfile, delimiter='\t')

print('files read')

print("16S DataFrame Shape:", seqtab16S.shape)
print(type(seqtab16S))
print("18S DataFrame Shape:", seqtab18S.shape)
print(type(seqtab18S))


pd.DataFrame(seqtab16S.columns).to_csv("SampleNames16S.csv")
pd.DataFrame(seqtab18S.columns).to_csv("SampleNames18S.csv")


# Pivot the DataFrames 
pivot_seqtab16S = seqtab16S.T
pivot_seqtab18S = seqtab18S.T


print('after pivoting:')
print("16S DataFrame Shape:", pivot_seqtab16S.shape)

print("18S DataFrame Shape:", pivot_seqtab18S.shape)

# Save as comma-separated file
pivot_seqtab16S.to_csv('norm_seqtab_16S_filt_pivoted.csv', index=False, header=False)
print('16S file saved')

pivot_seqtab18S.to_csv('norm_seqtab_18S_filt_pivoted.csv', index=False, header=False)
print('18S file saved')


