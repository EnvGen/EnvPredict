## Change working directory (add_argument(p, "-w", help="working directory",...)
## Before running on a new machine
rm(list = ls())
options(warn=-1)




#1. Installing packages

if ( ! "argparser" %in% installed.packages()[,"Package"]) {
  installed.packages("argparser")
}

if ( ! "DiagrammeRsvg" %in% installed.packages()[,"Package"]) {
  devtools::install_github('rich-iannone/DiagrammeRsvg')
}


bioconduc_packages=c("ape")
new.packages <- bioconduc_packages[!(bioconduc_packages %in% installed.packages()[,"Package"])]
if(length(new.packages) >0 ) {
  for (s in new.packages) { suppressPackageStartupMessages(BiocManager::install(s))}
}

#libraries
list_of_packages <- c("tidyverse","anytime","viridis","reshape2","vegan","ggpubr","ape","lubridate",
                      "dplyr","plyr","svglite","gridExtra","geosphere","corrplot","ggthemes",
                      "cowplot", "data.table", "mlr", "tidymodels", "randomForest", "xgboost", "DiagrammeR",
                      "Ckmeans.1d.dp") #,"DiagrammeRsvg", "rsvg")

new.packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new.packages) >0 ) {
  for (s in new.packages) { suppressPackageStartupMessages(install.packages(s))}
}



suppressPackageStartupMessages(library("argparser"))

# arguments
p <- arg_parser("Read datasets")
p <- add_argument(p, "-b", help="16S_ASV sequences", default="../seq_data/combined/16S/asvs_16S.txt")
p <- add_argument(p, "-e", help="18S_ASV sequences", default="../seq_data/combined/18S/asvs_18S.txt")
p <- add_argument(p, "-n", help="16S_count table", default="../seq_data/combined/16S/seqtab_16S.tsv")
p <- add_argument(p, "-y", help="18S_count table", default="../seq_data/combined/18S/seqtab_18S.tsv")
p <- add_argument(p, "-f", help="Barrnap-filered 16S ASV ids", default="../seq_data/combined/barrnap_cleaned_16S_18S/ASV_IDs_16S.txt")
p <- add_argument(p, "-g", help="Barrnap-filered 16S ASV ids", default="../seq_data/combined/barrnap_cleaned_16S_18S/ASV_IDs_18S.txt")
p <- add_argument(p, "-m", help="Metadata", default="../env_data/combined/physical_chemical_processed.tsv")
# p <- add_argument(p, "-a", help="Abiotic parameters", default="contextual_data/physical_chemical_data_2019-2020_2022-03-16.txt")
p <- add_argument(p, "-z", help="taxa_16S table", default="../seq_data/combined/16S/taxa_16S.tsv")
p <- add_argument(p, "-s", help="taxa_16S SILVA table", default="../seq_data/combined/16S/taxa_16S_SILVA.tsv")
p <- add_argument(p, "-v", help="taxa_18S table", default="../seq_data/combined/18S/taxa_18S.tsv")
p <- add_argument(p, "-w", help="working directory", default="//wsl.localhost/Ubuntu/home/krzjur/EnvPredict/code")
# p <- add_argument(p, "-d", help="Minimum depth", default=10)
p <- add_argument(p, "-o", help="output directory", default="../output/")



argv <- parse_args(p)

for (s in list_of_packages) { suppressPackageStartupMessages(library(s, character.only = TRUE))}


setwd(argv$w)
dir.create(argv$o)


## load sequence tables with counts
seqtab_16S = as.matrix(read.delim(argv$n))
seqtab_18S = as.matrix(read.delim(argv$y))
# load asv tables
asv_16S = as.vector(read.csv(argv$b, sep = '\t', header = FALSE)[[2]])
names(asv_16S) = as.vector(read.csv(argv$b, sep = '\t', header = FALSE)[[1]])
asv_16S = toupper(asv_16S)

asv_18S = as.vector(read.csv(argv$e, sep = '\t', header = FALSE)[[2]])
names(asv_18S) = as.vector(read.csv(argv$e, sep = '\t', header = FALSE)[[1]])
asv_18S = toupper(asv_18S)

# load tax tables
taxa_16S = as.matrix(read.table(argv$z, sep = '\t'))
taxa_16S_silva = as.matrix(read.table(argv$s, sep = '\t'))
taxa_18S = as.matrix(read.table(argv$v, sep = '\t'))

## Keep only ASVs identified as 16S/18S by Barrnap

barrnap_16S_ids = as.vector(read.csv(argv$f, sep = '\t', header = FALSE)[[1]])
barrnap_18S_ids = as.vector(read.csv(argv$g, sep = '\t', header = FALSE)[[1]])

taxa_16S = taxa_16S[rownames(taxa_16S) %in% barrnap_16S_ids,]
taxa_16S_silva = taxa_16S_silva[rownames(taxa_16S_silva) %in% barrnap_16S_ids,]
taxa_18S = taxa_18S[rownames(taxa_18S) %in% barrnap_18S_ids,]

seqtab_16S = seqtab_16S[rownames(seqtab_16S) %in% barrnap_16S_ids,]
seqtab_18S = seqtab_18S[rownames(seqtab_18S) %in% barrnap_18S_ids,]

asv_16S = asv_16S[names(asv_16S) %in% barrnap_16S_ids]
asv_18S = asv_18S[names(asv_18S) %in% barrnap_18S_ids]

#load metadata file
metadata = read.delim(paste(argv$m))
#### Other data (phyto,ect)

## Exclude spike reads - not needed with Barrnap!
#these are the spike DNA sequences
# spike_16S = 'CGGGCAGCTCTCGATAACCGGCGGAAGGTGGTAGCCACGGACAGGATCAGAACAATTAGAAGTGCCGCAGGTGGCCAAGTCCCCCGGACACAAGACGAGGCCGGAGGCCTGGTATATACACGTAGCTAAGAAGAGCTCATCCAGACTGGGAACGGTGTGCCAGCAGCCGCGGTAACATCACCACAACGTATTCGGTCACAAATTGATCGGAGGGAGAAATCGTCCGCAGGATCTCAAACTTTAACTAAGGACTAGTACTACATAGGCTCGAGAAGAGCTACCGTTTGCAGGGTCGCCGGGTACCGCTTAACCATAAAAGATCCACTCAGGTAGCCGTCCAGTTTCCTCTGAAATGATGGGGCGAGAAACACGGCTGGGCGTTATACGAGTGCTTTAGAATATGAGGAGAGACAGGGGTATATTCAAGG'
# spike_18S = 'CTTCGTTATCGTCACGGAGAGAACTGCCTTTAGCGATCTGTCTAGAGACCCGCCAAATATAAGCCTTGGGGTTATTAACAAAGACGCCAAACACTAGTGAATATGACAGTGAAGGGGTGTGAGATAGCTTTACTGGGTGAGAAAACACTCGTTAAAAAGAATTAGACCGGATAATCCCCGAGGGGCCGTAGGCATGGACTTGTCGTTGCCACCGAGCATAGCGGTTTCGAAATAGCCGAGATGGGCACTGGCGAATTAACCCACTGGTTTATATGGATCCGATGGGTTCACTTAATAAGCTCGTACCAGGGATGAATAAAGCGTTACGAGAATTATAAACATGGAGTTCCTATTGATTTGAGGTTAATACCGAACGGGAACATTTGTCGATCATGCTTCACATAGAGT'
# 
# # For 16S - identify spike ASVs
# spike_16S_ix = agrep(spike_16S, asv_16S)
# spiketab_16S = seqtab_16S[spike_16S_ix, ]
# dim(spiketab_16S)
# 
# ## Check if spike ASVs are unannotated
# # taxa_16S[spike_16S_ix,]
# 
# # For 18S - identify spike ASVs
# spike_18S_ix = agrep(spike_18S, asv_18S)
# spiketab_18S = seqtab_18S[spike_18S_ix, ]


## 16S exclude mitochondria and chloroplasts

ix_taxa_16S = which(taxa_16S_silva[,4] != 'Chloroplast' | is.na(taxa_16S_silva[,4]))
ix_taxa_16S = intersect(ix_taxa_16S, which(taxa_16S_silva[,5] != 'Mitochondria' | is.na(taxa_16S_silva[,5])))

taxa_16S_silva = taxa_16S_silva[ix_taxa_16S,]

taxa_16S = taxa_16S[ix_taxa_16S,]
seqtab_16S = seqtab_16S[ix_taxa_16S,]

## For 18S - remove Metazoa
ix_taxa_18S = which(taxa_18S[,4] != 'Metazoa' | is.na(taxa_18S[,4]))
taxa_18S = taxa_18S[ix_taxa_18S,]
seqtab_18S = seqtab_18S[ix_taxa_18S,]

### Filter samples ###

## Adjust seqtab colnames


colnames(seqtab_18S) = gsub(colnames(seqtab_18S), pattern = '_18S', replacement = '')
colnames(seqtab_16S) = gsub(colnames(seqtab_16S), pattern = '_16S', replacement = '')

colnames(seqtab_18S) = gsub(colnames(seqtab_18S), pattern = '^X', replacement = '')
colnames(seqtab_16S) = gsub(colnames(seqtab_16S), pattern = '^X', replacement = '')

## Pick only samples which are in the seqtab

ix = intersect(which(metadata$sample_id %in% colnames(seqtab_16S)),
               which(metadata$sample_id %in% colnames(seqtab_18S)))
metadata = metadata[ix,]

## Pick only columns from seqtab which correspond to the monitoring samples

iy = which(colnames(seqtab_16S) %in% metadata$sample_id)
seqtab_16S = seqtab_16S[,iy]

iy = which(colnames(seqtab_18S) %in% metadata$sample_id)
seqtab_18S = seqtab_18S[,iy]

# And set the same order of samples in metadata and seqtab
iy = match(metadata$sample_id, colnames(seqtab_16S))
seqtab_16S = seqtab_16S[,iy]

iy = match(metadata$sample_id, colnames(seqtab_18S))
seqtab_18S = seqtab_18S[,iy]

## Exclude samples which are not high-quality

samples_to_exclude = c('20191015_263953_1') ## suspected contamination
samples_to_exclude = c(samples_to_exclude, '20190820_264450_1') ## undersequenced samples with replicates available

# Choose just one of replicates from the same time and station

dates_stations = c()
meta_ix = setdiff(1:nrow(metadata), which(metadata$sample_id %in% samples_to_exclude))

for(i in meta_ix){
  ds = paste(metadata$date[i], metadata$station_name[i])
  if(! ds %in% dates_stations){
    dates_stations = c(dates_stations, ds)
  }else{
    meta_ix = meta_ix[! meta_ix %in% i]
  }
}

# Pick only selected stations for further analyses
metadata = metadata[meta_ix,]
ix = metadata$sample_id
seqtab_16S = seqtab_16S[,ix]
seqtab_18S = seqtab_18S[,ix]

## Change NB1 / 3 to B7
## Since they are practically the same sampling location

metadata$station_name[metadata$station_name == 'NB1 / B3'] = 'B7'

### Analyze nr of reads per sample ###
### Exclude under/overssequenced samples ###

## 18S

samples_to_exclude = c()

## Different thresholds across the years, since different sequencing technologies were used (wore reads in general in 2015-2017 samples)

ix_2019 = which(metadata$date > as.Date('2019-01-01'))

# png(paste0(plot_folder, 'tot_reads_sorted_', data_type, '_2019.png'),
#     width = 5000, height = 3000, res = 500, units = 'px')
# pdf(paste0(plot_folder, 'tot_reads_sorted_', data_type, '.pdf'),
#     width = 10, height = 6)
par(mar = c(2,2,2,2))
plot(sort(colSums(seqtab_18S)[ix_2019]), xlab = 'Samples', ylab = 'Reads')

# dev.off()

sort(colSums(seqtab_18S)[ix_2019])

ix = which(colSums(seqtab_18S)[ix_2019] > 200000)
samples_to_exclude = c(samples_to_exclude, names(colSums(seqtab_18S)[ix_2019][ix]))

ix_2015 = which(metadata$date < as.Date('2019-01-01'))

# png(paste0(plot_folder, 'tot_reads_sorted_', data_type, '_2015.png'),
#     width = 5000, height = 3000, res = 500, units = 'px')
# pdf(paste0(plot_folder, 'tot_reads_sorted_', data_type, '.pdf'),
#     width = 10, height = 6)
par(mar = c(2,2,2,2))
plot(sort(colSums(seqtab_18S)[ix_2015]), xlab = 'Samples', ylab = 'Reads')

# dev.off()

sort(colSums(seqtab_18S)[ix_2015])

## Pick only the samples with enough reads

ix = which(! colnames(seqtab_18S) %in% samples_to_exclude)
seqtab_18S = seqtab_18S[,ix]

## 16S

samples_to_exclude = c()

## Different thresholds across the years, since different sequencing technologies were used (wore reads in general in 2015-2017 samples)

ix_2019 = which(metadata$date > as.Date('2019-01-01'))

# png(paste0(plot_folder, 'tot_reads_sorted_', data_type, '_2019.png'),
#     width = 5000, height = 3000, res = 500, units = 'px')
# pdf(paste0(plot_folder, 'tot_reads_sorted_', data_type, '.pdf'),
#     width = 10, height = 6)
par(mar = c(2,2,2,2))
plot(sort(colSums(seqtab_16S)[ix_2019]), xlab = 'Samples', ylab = 'Reads')

# dev.off()

sort(colSums(seqtab_16S)[ix_2019])

ix = which(colSums(seqtab_16S)[ix_2019] < 30000)
samples_to_exclude = c(samples_to_exclude, names(colSums(seqtab_16S)[ix_2019][ix]))

ix_2015 = which(metadata$date < as.Date('2019-01-01'))

# png(paste0(plot_folder, 'tot_reads_sorted_', data_type, '_2015.png'),
#     width = 5000, height = 3000, res = 500, units = 'px')
# pdf(paste0(plot_folder, 'tot_reads_sorted_', data_type, '.pdf'),
#     width = 10, height = 6)
par(mar = c(2,2,2,2))
plot(sort(colSums(seqtab_16S)[ix_2015]), xlab = 'Samples', ylab = 'Reads')

# dev.off()

sort(colSums(seqtab_16S)[ix_2015])

ix = which(colSums(seqtab_16S)[ix_2015] > 2000000)
samples_to_exclude = c(samples_to_exclude, names(colSums(seqtab_16S)[ix_2015][ix]))

## Pick only the samples with enough reads

ix = which(! colnames(seqtab_16S) %in% samples_to_exclude)
seqtab_16S = seqtab_16S[,ix]

## Get only samples with quality 16S AND 18S data - discuss

common_samples = intersect(colnames(seqtab_16S),
                           colnames(seqtab_18S))
seqtab_16S = seqtab_16S[,common_samples]
seqtab_18S = seqtab_18S[,common_samples]

ix = which(metadata$sample_id %in% common_samples)
metadata = metadata[ix,]

## Normalize the reads

norm_seqtab_16S = seqtab_16S
for (i in 1:ncol(seqtab_16S)) {
  norm_seqtab_16S[,i] = seqtab_16S[,i]/sum(seqtab_16S[,i])
}

# For 18S
norm_seqtab_18S = seqtab_18S
for (i in 1:ncol(seqtab_18S)) {
  norm_seqtab_18S[,i] = seqtab_18S[,i]/sum(seqtab_18S[,i])
}

## Sum up clade counts at each taxonomic level

## 16S

clade_counts_16S = list()
norm_clade_counts_16S = list()

for (i in 1:ncol(taxa_16S)) {
  matr = norm_matr = NULL
  clade = unique(taxa_16S[,i])
  clade = clade[!is.na(clade)]
  for (j in 1:length(clade)) {
    ix = which(clade[j]==taxa_16S[,i])
    if (length(ix) > 1) {
      matr = rbind(matr, apply(seqtab_16S[ix,], 2, sum, na.rm=TRUE))
      norm_matr = rbind(norm_matr, apply(norm_seqtab_16S[ix,], 2, sum, na.rm=TRUE))
    } else {
      matr = rbind(matr, seqtab_16S[ix,])
      norm_matr = rbind(norm_matr, norm_seqtab_16S[ix,])
    }
  }
  rownames(matr) = rownames(norm_matr) = clade
  colnames(matr) = colnames(norm_matr) = metadata$sample_id
  clade_counts_16S[[i]] = matr
  norm_clade_counts_16S[[i]] = norm_matr
}

## 18S

clade_counts_18S = list()
norm_clade_counts_18S = list()

for (i in 1:ncol(taxa_18S)) {
  matr = norm_matr = NULL
  clade = unique(taxa_18S[,i])
  clade = clade[!is.na(clade)]
  for (j in 1:length(clade)) {
    ix = which(clade[j]==taxa_18S[,i])
    if (length(ix) > 1) {
      matr = rbind(matr, apply(seqtab_18S[ix,], 2, sum, na.rm=TRUE))
      norm_matr = rbind(norm_matr, apply(norm_seqtab_18S[ix,], 2, sum, na.rm=TRUE))
    } else {
      matr = rbind(matr, seqtab_18S[ix,])
      norm_matr = rbind(norm_matr, norm_seqtab_18S[ix,])
    }
  }
  rownames(matr) = rownames(norm_matr) = clade
  colnames(matr) = colnames(norm_matr) = metadata$sample_id
  clade_counts_18S[[i]] = matr
  norm_clade_counts_18S[[i]] = norm_matr
}

write.table(seqtab_16S,
            file = paste('../seq_data/combined/16S/filtered_seqtab_16S.tsv', sep = ''),
            sep = '\t', row.names = TRUE, col.names = TRUE)
write.table(seqtab_18S,
            file = paste('../seq_data/combined/18S/filtered_seqtab_18S.tsv', sep = ''),
            sep = '\t', row.names = TRUE, col.names = TRUE)
write.table(metadata,
            file = paste('../seq_data/combined/filtered_metadata.tsv', sep = ''),
            sep = '\t', row.names = FALSE, col.names = TRUE)
## Save taxa tables

write.table(taxa_16S,
            file = paste('../seq_data/combined/16S/filtered_taxa_16S.tsv', sep = ''),
            sep = '\t', row.names = TRUE, col.names = TRUE)
write.table(taxa_18S,
            file = paste('../seq_data/combined/18S/filtered_taxa_18S.tsv', sep = ''),
            sep = '\t', row.names = TRUE, col.names = TRUE)


## Loops to save each matrix in clade_counts_16S as a separate file

for (i in 1:length(clade_counts_16S)) {
  write.table(clade_counts_16S[[i]],
              file = paste('../seq_data/combined/16S/clade_counts_16S_', i, '.tsv', sep = ''),
              sep = '\t', row.names = TRUE, col.names = TRUE)
  write.table(norm_clade_counts_16S[[i]],
              file = paste('../seq_data/combined/16S/norm_clade_counts_16S_', i, '.tsv', sep = ''),
              sep = '\t', row.names = TRUE, col.names = TRUE)
}

for (i in 1:length(clade_counts_18S)) {
  write.table(clade_counts_18S[[i]],
              file = paste('../seq_data/combined/18S/clade_counts_18S_', i, '.tsv', sep = ''),
              sep = '\t', row.names = TRUE, col.names = TRUE)
  write.table(norm_clade_counts_18S[[i]],
              file = paste('../seq_data/combined/18S/norm_clade_counts_18S_', i, '.tsv', sep = ''),
              sep = '\t', row.names = TRUE, col.names = TRUE)
}

## Save normalized count tables

write.table(norm_seqtab_16S,
            file = paste('../seq_data/combined/16S/norm_seqtab_16S.tsv', sep = ''),
            sep = '\t', row.names = TRUE, col.names = TRUE)
write.table(norm_seqtab_18S,
            file = paste('../seq_data/combined/18S/norm_seqtab_18S.tsv', sep = ''),
            sep = '\t', row.names = TRUE, col.names = TRUE)

save.image('read_data.Rdata')

### FInished reading data ###
