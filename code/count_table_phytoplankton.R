#Set wd
setwd("C:/Users/johan/OneDrive/R/Merged metabarcoding data")

#Load smhi data
phyto_data=read.csv("C:/Users/johan/Downloads/phytoplankton_2015_2020_windows.txt", sep = "\t")

#Create column for sample name
phyto_data$sample=paste0(phyto_data$Nationellt.övervakningsstations.ID, "_", phyto_data$Provtagningsdatum..start.)

#Create count table 
unique_taxa_phyto=unique(phyto_data$Vetenskapligt.namn)
unique_samples_phyto=unique(phyto_data$sample)

num_columns = length(unique_samples_phyto)
num_rows = length(unique_taxa_phyto)

count_df_phyto = data.frame(matrix(0, nrow = num_rows, ncol = num_columns))
colnames(count_df_phyto) = unique_samples_phyto
rownames(count_df_phyto) = unique_taxa_phyto

#Loop through smhi data and add ocurrances to the count_df
for (i in seq(nrow(phyto_data))) {
  sample_name = phyto_data$sample[i]
  taxa = phyto_data$Vetenskapligt.namn[i]
  
  if (sample_name %in% colnames(count_df_phyto) && taxa %in% rownames(count_df_phyto)) {
    count_df_phyto[taxa, sample_name] <- count_df_phyto[taxa, sample_name] + 1
  }
}

write.table(count_df_phyto,"C:/Users/johan/OneDrive/R/ML/count_table_phytoplankton")
