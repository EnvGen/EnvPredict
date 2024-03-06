#Set wd
setwd("C:/Users/johan/OneDrive/R/Merged metabarcoding data")

#Load smhi data
smhi_data=read.table("C:/Users/johan/OneDrive/R/SMHI data/Merged_smhi_data_2015_2024/zooplankton_2015_2020_2024-01-25_utf8.txt", header = TRUE, sep = "\t")

#Load packages
library(dplyr)
library(tidyr)

#Create column for sample name
smhi_data$sample=paste0(smhi_data$station_id, "_", smhi_data$sample_date)

#COUNT TABLE FOR MIXED TAXA LEVELS
unique_taxa=unique(smhi_data$scientific_name)
unique_samples=unique(smhi_data$sample)

num_columns = length(unique_samples)
num_rows = length(unique_taxa)

count_df = data.frame(matrix(0, nrow = num_rows, ncol = num_columns))
colnames(count_df) = unique_samples
rownames(count_df) = unique_taxa

#Loop through smhi data and add ocurrances to the count_df
for (i in seq(nrow(smhi_data))) {
  sample_name = smhi_data$sample[i]
  taxa = smhi_data$scientific_name[i]
  
  if (sample_name %in% colnames(count_df) && taxa %in% rownames(count_df)) {
    count_df[taxa, sample_name] <- count_df[taxa, sample_name] + 1
  }
}

colnames(count_df) = gsub(colnames(count_df), pattern = '^X', replacement = '')
colnames(count_df) = gsub(colnames(count_df), pattern = '\\.', replacement = '-')

write.table(count_df, "C:/Users/johan/OneDrive/R/ML/count_table_zooplankton")

#COUNT TABLE FOR GENUS LEVEL 

#Create count table 
unique_taxa=unique(smhi_data$taxon_genus)
unique_samples=unique(smhi_data$sample)

num_columns = length(unique_samples)
num_rows = length(unique_taxa)

count_df = data.frame(matrix(0, nrow = num_rows, ncol = num_columns))
colnames(count_df) = unique_samples
rownames(count_df) = unique_taxa

#Loop through smhi data and add ocurrances to the count_df
for (i in seq(nrow(smhi_data))) {
  sample_name = smhi_data$sample[i]
  taxa = smhi_data$taxon_genus[i]
  
  if (sample_name %in% colnames(count_df) && taxa %in% rownames(count_df)) {
    count_df[taxa, sample_name] <- count_df[taxa, sample_name] + 1
  }
}

colnames(count_df) = gsub(colnames(count_df), pattern = '^X', replacement = '')
colnames(count_df) = gsub(colnames(count_df), pattern = '\\.', replacement = '-')
count_df=count_df[-4, , drop = FALSE]

write.table(count_df, "C:/Users/johan/OneDrive/R/ML/count_table_zooplankton_genus")

#COUNT TABLE FOR SPECIES LEVEL 

#Create count table 
unique_taxa=unique(smhi_data$taxon_species)
unique_samples=unique(smhi_data$sample)

num_columns = length(unique_samples)
num_rows = length(unique_taxa)

count_df = data.frame(matrix(0, nrow = num_rows, ncol = num_columns))
colnames(count_df) = unique_samples
rownames(count_df) = unique_taxa

#Loop through smhi data and add ocurrances to the count_df
for (i in seq(nrow(smhi_data))) {
  sample_name = smhi_data$sample[i]
  taxa = smhi_data$taxon_species[i]
  
  if (sample_name %in% colnames(count_df) && taxa %in% rownames(count_df)) {
    count_df[taxa, sample_name] <- count_df[taxa, sample_name] + 1
  }
}

colnames(count_df) = gsub(colnames(count_df), pattern = '^X', replacement = '')
colnames(count_df) = gsub(colnames(count_df), pattern = '\\.', replacement = '-')
count_df=count_df[-1, , drop = FALSE]


write.table(count_df, "C:/Users/johan/OneDrive/R/ML/count_table_zooplankton_species")
