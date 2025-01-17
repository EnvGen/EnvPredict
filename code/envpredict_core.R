
## Set working directory
setwd("~/aquatic/envpredict")

## Load libraries
library(lubridate)
library(ape)
library(phangorn)
library(seqinr)
library(vegan)
library(pheatmap)
library(randomForest)
library(ranger)
library(rsample)
library(gbm)
library(stringi)
library(lubridate)

## Set files
seqtab_file_18S = "../seq_data/combined/18S/filtered_seqtab_18S.tsv" ## count matrix for each ASV & sample
taxa_file_18S = "../seq_data/combined/18S/filtered_taxa_18S.tsv" ## taxonomic annotation for each ASV
seqtab_file_16S = "../seq_data/combined/16S/filtered_seqtab_16S.tsv" ## count matrix for each ASV & sample
taxa_file_16S =  "../seq_data/combined/16S/filtered_taxa_16S.tsv" ## taxonomic annotation for each ASV
phys_chem_file = "../env_data/combined/physical_chemical_processed_translation.tsv"
#bact_plan_file = "env_data/combined/bacterioplankton_processed.tsv"
#pico_plan_file = "env_data/combined/picoplankton_processed.tsv"
phyt_plan_file = "../env_data/combined/phytoplankton_processed.tsv"
zoop_plan_file = "../env_data/combined/zooplankton_processed.tsv"
id_translation_file = "../env_data/combined/physical_chemical_processed_translation.tsv"

# read sequencing data
asv_counts_18S = as.matrix(read.delim(seqtab_file_18S, row.names = 1))
asv_taxa_18S = as.matrix(read.delim(taxa_file_18S, row.names = 1))[,1:9]
asv_counts_16S = as.matrix(read.delim(seqtab_file_16S, row.names = 1))
asv_taxa_16S = as.matrix(read.delim(taxa_file_16S, row.names = 1))[,1:7]

colnames(asv_counts_18S) = gsub("^X", "", colnames(asv_counts_18S))
colnames(asv_counts_16S) = gsub("^X", "", colnames(asv_counts_16S))
rownames(asv_counts_18S) = paste(rownames(asv_counts_18S), "18S", sep = "_")
rownames(asv_taxa_18S) = paste(rownames(asv_taxa_18S), "18S", sep = "_")
rownames(asv_counts_16S) = paste(rownames(asv_counts_16S), "16S", sep = "_")
rownames(asv_taxa_16S) = paste(rownames(asv_taxa_16S), "16S", sep = "_")

identical(rownames(asv_counts_18S), rownames(asv_taxa_18S))
identical(rownames(asv_counts_16S), rownames(asv_taxa_16S))
identical(colnames(asv_counts_16S), colnames(asv_counts_18S))
colnames(asv_counts_16S) = gsub("^X", "", colnames(asv_counts_16S))
samples = colnames(asv_counts_16S)

alt_ids = as.matrix(read.delim(id_translation_file))
alt_ids = alt_ids[,c(1,ncol(alt_ids))]
alt_ids[,2] = gsub("-","", alt_ids[,2])
samples_alt = alt_ids[match(samples, alt_ids[,1]),2]

ix = setdiff(1:length(samples_alt), grep("NA_", samples_alt)) # samples with station id
asv_counts_16S = asv_counts_16S[,ix]
asv_counts_18S = asv_counts_18S[,ix]
samples = samples[ix]
samples_alt = samples_alt[ix]
colnames(asv_counts_16S) = colnames(asv_counts_18S) = samples_alt

norm_asv_counts_16S = t(t(asv_counts_16S)/colSums(asv_counts_16S))
norm_asv_counts_18S = t(t(asv_counts_18S)/colSums(asv_counts_18S))

# read physchem
phys_chem = read.delim(phys_chem_file)
rownames(phys_chem) = phys_chem[,1]
year = year(phys_chem[,4])
yday = yday(phys_chem[,4])
ylength = rep(NA, length(yday))

## TIme since midnight
times = phys_chem$time_h
times = interval(start = strptime("00:00:00", "%H:%M:%S"), end = strptime(times, "%H:%M:%S"))
times = time_length(times, unit = 'minutes')
times_since_midnight = times
time_xcord = cos(2*pi*times/24*60)
time_ycord = sin(2*pi*times/24*60)

ix = which(!is.na(yday))
ylength[ix] = yday(paste(as.character(year[ix]), "-12-31", sep = ""))
yangle = 2*pi*yday/ylength
yday_xcord = cos(yangle)
yday_ycord = sin(yangle)

phys_chem = t(cbind(year, yday, yday_xcord, yday_ycord, times, times_xcord, times_ycord, phys_chem[,c(9:26)]))


ix = match(samples, colnames(phys_chem)) # samples for which we have sequence data
phys_chem = phys_chem[,ix]
colnames(phys_chem) = samples_alt

phys_chem_small = phys_chem[c(3:8,12:18),]
phys_chem_small = phys_chem_small[, which(!is.na(colSums(phys_chem_small)))] # lacks samples with NA in any of the parameters

samples = samples_alt # no need of original sample names any more

#bact_plan = read.delim(bact_plan_file)
#pico_plan = read.delim(pico_plan_file)

# read plankton microscopy data
phyt_plan = read.delim(phyt_plan_file)
colnames(phyt_plan) = gsub("^X", "", colnames(phyt_plan))
colnames(phyt_plan) = gsub("\\.", "-", colnames(phyt_plan))
colnames(phyt_plan) = gsub("-", "", colnames(phyt_plan))
phyt_plan = as.matrix(phyt_plan)
ix = which(!is.na(match(colnames(phyt_plan), samples)))
phyt_plan = phyt_plan[,ix]

zoo_plan = read.delim(zoop_plan_file)
colnames(zoo_plan) = gsub("^X", "", colnames(zoo_plan))
colnames(zoo_plan) = gsub("\\.", "-", colnames(zoo_plan))
colnames(zoo_plan) = gsub("-", "", colnames(zoo_plan))
zoo_plan = as.matrix(zoo_plan)
ix = which(!is.na(match(colnames(zoo_plan), samples)))
zoo_plan = zoo_plan[,ix]

# sum per genus for seq data
sum_seq_per_genus <- function (count_matrix, taxa_matrix) {
  genus = taxa_matrix[,8]
  unique_genus = sort(unique(genus))
  genus_count_matrix = matrix(ncol = ncol(count_matrix), nrow = length(unique_genus))
  rownames(genus_count_matrix) = unique_genus
  colnames(genus_count_matrix) = colnames(count_matrix)
  for (i in 1:length(unique_genus)) {
    ix = which(genus == unique_genus[i])
    if (length(ix) > 1) {
      genus_count_matrix[i,] = colSums(count_matrix[ix,])
    } else {
      genus_count_matrix[i,] = count_matrix[ix,]
    }
  }
  return(genus_count_matrix)
}
norm_asv_counts_18S_genus = sum_seq_per_genus(norm_asv_counts_18S, asv_taxa_18S)

# sum per genus for microscopy data
sum_micr_per_genus <- function (count_matrix) {
  a<-strsplit(rownames(count_matrix), " ")
  res<-as.data.frame(t(stri_list2matrix(a)))
  genus = res[,1]
  unique_genus = sort(unique(genus))
  genus_count_matrix = matrix(ncol = ncol(count_matrix), nrow = length(unique_genus))
  rownames(genus_count_matrix) = unique_genus
  colnames(genus_count_matrix) = colnames(count_matrix)
  for (i in 1:length(unique_genus)) {
    ix = which(genus == unique_genus[i])
    if (length(ix) > 1) {
      genus_count_matrix[i,] = colSums(count_matrix[ix,])
    } else {
      genus_count_matrix[i,] = count_matrix[ix,]
    }
  }
  return(genus_count_matrix)
}
zoo_plan_genus = sum_micr_per_genus(zoo_plan)
phyt_plan_genus = sum_micr_per_genus(phyt_plan)


## Functions for preparing feature and response matrices
extract_shared_samples <- function(features_matrix_full, responses_matrix_full) { # selecting only samples shared by features and responses
  return_matrices = list()
  shared_samples = intersect(colnames(features_matrix_full), colnames(responses_matrix_full))
  if (length(grep("NA_", shared_samples)) > 0) {
    shared_samples = shared_samples[-grep("NA_", shared_samples)]
  }
  ix = match(shared_samples, colnames(features_matrix_full))
  features_matrix = features_matrix_full[,ix]
  ix = match(shared_samples, colnames(responses_matrix_full))
  responses_matrix = responses_matrix_full[,ix]
  responses_matrix = responses_matrix[which(rowSums(responses_matrix, na.rm = T) > 0),]
  return_matrices[[1]] = features_matrix
  return_matrices[[2]] = responses_matrix
  names(return_matrices) = c("features_matrix", "responses_matrix")
  return(return_matrices)
}

do_feature_selection <- function(features_matrix, occupancy) { # occupancy = proportion of samples the feature should be present in
  binary_features_matrix = features_matrix
  binary_features_matrix[which(features_matrix > 0)] = 1
  ix = which(rowSums(binary_features_matrix) > occupancy * ncol(features_matrix))
  features_matrix = features_matrix[ix,]
}

## Functions for running randomforest predictions
run_randomforest_out_of_bag <- function(features_matrix, responses_matrix) {
  predicted_responses_matrix = matrix(nrow = nrow(responses_matrix), ncol = ncol(responses_matrix))
  for (i in 1:nrow(responses_matrix)) {
    #rf = randomForest(t(features_matrix), y=responses_matrix[i,], importance=T, ntree = 1000)
    #predicted_responses_matrix[i,] = rf$predicted
    ix = which(!is.na(responses_matrix[i,]))
    df = as.data.frame(cbind(responses_matrix[i,ix], t(features_matrix[,ix])))
    rf = ranger(V1 ~ ., data = df, num.trees = 2000, importance = 'impurity')
    predicted_responses_matrix[i,ix] = rf$predictions
  }
  return(predicted_responses_matrix)
}

run_randomforest <- function(features_matrix, responses_matrix, numfolds, min_samples) {
  min_samples_in_fold = 1
  predicted_responses_matrix = matrix(nrow = nrow(responses_matrix), ncol = ncol(responses_matrix))
  for (i in 1:nrow(responses_matrix)) {
    response = responses_matrix[i,]
    if (length(which(!is.na(response))) < min_samples) { next } # skip this parameter if too few non-NA samples
    df = as.data.frame(cbind(response, t(features_matrix)))
    folds <- vfold_cv(df, v = numfolds, strata = response)
    #folds <- vfold_cv(df, v = numfolds, strata = NULL)
    for (j in 1:nrow(foldxs)) {
      trainingIndex <- folds$splits[[j]]$in_id
      testingIndex <- setdiff(1:length(response), trainingIndex)
      trainingIndex = trainingIndex[which(!is.na(response[trainingIndex]))] # remove NA response samples
      if (length(trainingIndex) >= min_samples_in_fold) {
        rf = ranger(response ~ ., data = df[trainingIndex,], num.trees = 2000, importance = 'impurity')
        predicted_responses_matrix[i,testingIndex] = predict(rf, df[testingIndex,])$predictions  
      }
    }
  }
  return(predicted_responses_matrix)
}

run_xgboost <- function(features_matrix, responses_matrix, numfolds, min_samples) {
  min_samples_in_fold = 1
  predicted_responses_matrix = matrix(nrow = nrow(responses_matrix), ncol = ncol(responses_matrix))
  for (i in 1:nrow(responses_matrix)) {
    response = responses_matrix[i,]
    if (length(which(!is.na(response))) < min_samples) { next } # skip this parameter if too few non-NA samples
    df = as.data.frame(cbind(response, t(features_matrix)))
    folds <- vfold_cv(df, v = numfolds, strata = response)
    for (j in 1:nrow(folds)) {
      trainingIndex <- folds$splits[[j]]$in_id 
      testingIndex <- setdiff(1:length(response), trainingIndex)
      trainingIndex = trainingIndex[which(!is.na(response[trainingIndex]))] # remove NA response samples
      if (length(trainingIndex) >= min_samples_in_fold) {
        gbm <- gbm(response ~ ., data = df[trainingIndex,], distribution = "gaussian", n.trees = 10000, shrinkage = 0.001, interaction.depth = 2, n.minobsinnode = 1)
        predicted_responses_matrix[i,testingIndex] = predict(gbm, df[testingIndex,], n.trees = 10000)
      }
    }
  }
  return(predicted_responses_matrix)
}

## Run ML predictions

## use one of the below as features_matrix_full
#features_matrix_full = asv_counts_16S
#features_matrix_full = norm_asv_counts_16S
#features_matrix_full = asv_counts_18S
#features_matrix_full = norm_asv_counts_18S
#features_matrix_full = rbind(norm_asv_counts_16S, norm_asv_counts_18S)
#features_matrix_full = phyt_plan_genus
#features_matrix_full = zoo_plan_genus
#features_matrix_full = phys_chem_small
  
## use one of the below as responseses_matrix_full
#responses_matrix_full = phys_chem
#responses_matrix_full = zoo_plan
#responses_matrix_full = zoo_plan_genus
#responses_matrix_full = phyt_plan
#responses_matrix_full = phyt_plan_genus
#responses_matrix_full = norm_asv_counts_18S_genus
  

## for using 16S or 18S as features and physchem as responses
features_matrix_full = norm_asv_counts_16S
#features_matrix_full = norm_asv_counts_18S
responses_matrix_full = phys_chem
features_matrix = extract_shared_samples(features_matrix_full, responses_matrix_full)$features_matrix
responses_matrix = extract_shared_samples(features_matrix_full, responses_matrix_full)$responses_matrix
features_matrix = do_feature_selection(features_matrix, 0.1)
predicted_responses_matrix_rf_ob = run_randomforest_out_of_bag(features_matrix, responses_matrix)

## for using plankton as features and physchem as responses
features_matrix_full = zoo_plan_genus
#features_matrix_full = phyt_plan_genus
responses_matrix_full = phys_chem
features_matrix = extract_shared_samples(features_matrix_full, responses_matrix_full)$features_matrix
responses_matrix = extract_shared_samples(features_matrix_full, responses_matrix_full)$responses_matrix
features_matrix = do_feature_selection(features_matrix, 0.1)
predicted_responses_matrix_rf_ob = run_randomforest_out_of_bag(features_matrix, responses_matrix)

## for using 16S or 18S as features and plankton as responses
features_matrix_full = norm_asv_counts_16S
#features_matrix_full = norm_asv_counts_18S
responses_matrix_full = phyt_plan_genus
features_matrix = extract_shared_samples(features_matrix_full, responses_matrix_full)$features_matrix
responses_matrix = extract_shared_samples(features_matrix_full, responses_matrix_full)$responses_matrix
features_matrix = do_feature_selection(features_matrix, 0.1)
predicted_responses_matrix_rf_ob = run_randomforest_out_of_bag(features_matrix, responses_matrix)

## for using physchem as features and plankton genera as responses
features_matrix_full = phys_chem
responses_matrix_full = phyt_plan_genus
features_matrix = extract_shared_samples(features_matrix_full, responses_matrix_full)$features_matrix
responses_matrix = extract_shared_samples(features_matrix_full, responses_matrix_full)$responses_matrix
#features_matrix = do_feature_selection(features_matrix, 0.1)
predicted_responses_matrix_rf_ob = run_randomforest_out_of_bag(features_matrix, responses_matrix)

## for using physchem as features and 18S genera as responses
features_matrix_full = phys_chem
responses_matrix_full = norm_asv_counts_18S_genus
features_matrix = extract_shared_samples(features_matrix_full, responses_matrix_full)$features_matrix
responses_matrix = extract_shared_samples(features_matrix_full, responses_matrix_full)$responses_matrix
#features_matrix = do_feature_selection(features_matrix, 0.1)
predicted_responses_matrix_rf_ob = run_randomforest_out_of_bag(features_matrix, responses_matrix)


## examples on how to run the actual predictions
#predicted_responses_matrix_rf_5f = run_randomforest(features_matrix, responses_matrix, 5, 10)
predicted_responses_matrix_rf_10f = run_randomforest(features_matrix, responses_matrix, 10, 10)
#predicted_responses_matrix_rf_ob = run_randomforest_out_of_bag(features_matrix, responses_matrix)
#predicted_responses_matrix_xgb = run_xgboost(features_matrix, responses_matrix, 5, 10)

####

cor_matr = matrix(ncol = 6, nrow = nrow(responses_matrix))
rownames(cor_matr) = rownames(responses_matrix)
colnames(cor_matr) = c("#samp","#pos_samp","RF_5f", "RF_10f", "RF_OB","XGB_5f")
for (i in 1:nrow(responses_matrix)) {
  cor_matr[i,1] = length(which(!is.na(responses_matrix[i,])))
  cor_matr[i,2] = length(which(responses_matrix[i,] > 0))
  #cor_matr[i,3] = cor.test(responses_matrix[i,], predicted_responses_matrix_rf_5f[i,])$est
  #cor_matr[i,4] = cor.test(responses_matrix[i,], predicted_responses_matrix_rf_10f[i,])$est
  cor_matr[i,5] = cor.test(responses_matrix[i,], predicted_responses_matrix_rf_ob[i,])$est
  #if (!is.na(predicted_responses_matrix_xgb[i,1])) {
  #  cor_matr[i,6] = cor.test(responses_matrix[i,], predicted_responses_matrix_xgb[i,])$est
  #}
}
summary(cor_matr)
#cbind(cor_matr[,1], cor_matr[,2], cor_matr[,3], cor_matr[,4], cor_matr[,5])
cbind(cor_matr[,1], cor_matr[,2], cor_matr[,5])
plot(cor_matr[,2], cor_matr[,5])

par(mfrow = c(2,2))
par(mfrow = c(1,1))
plot(rowSums(responses_matrix), cor_matr[,2], col = "darkgreen", xlab = "micr. tot. count", ylab = "R-value") # RF_5f
points(rowSums(responses_matrix), cor_matr[,3], col = "blue") # RF_10f
points(rowSums(responses_matrix), cor_matr[,4], col = "black") # RF_ob
points(rowSums(responses_matrix), cor_matr[,5], col = "red") # XGB

plot(log10(rowSums(responses_matrix)), cor_matr[,2], col = "darkgreen", xlab = "micr. tot. count", ylab = "R-value") # RF_5f
points(log10(rowSums(responses_matrix)), cor_matr[,3], col = "blue") # RF_10f
points(log10(rowSums(responses_matrix)), cor_matr[,4], col = "black") # RF_ob
points(log10(rowSums(responses_matrix)), cor_matr[,5], col = "red") # XGB

plot(cor_matr[,3], cor_matr[,4], xlab = "RF_5f", ylab = "RF_10f", xlim = c(0,1), ylim = c(0,1))
lines(c(0,1),c(0,1))

plot(cor_matr[,4], cor_matr[,5], xlab = "RF_10f", ylab = "RF_out-of-bag", xlim = c(0,1), ylim = c(0,1))
lines(c(0,1),c(0,1))

plot(cor_matr[,4], cor_matr[,6], xlab = "RF_10f", ylab = "XGB_5f", xlim = c(0,1), ylim = c(0,1))
lines(c(0,1),c(0,1))

### for plankton
par(mfcol=c(1,1), mar=c(4,4,1,1))
plot(cor_matr[,2], cor_matr[,4], xlab = "Presence in #samples", ylab = "Prediction R-value", ylim = c(0,1))

# example cor plots
par(mfcol = c(3,3), mar=c(2,3,1,1))
taxon = "Tripos"
ix = which(rownames(responses_matrix) == taxon)
cor.test(responses_matrix[ix,], predicted_responses_matrix_rf_10f[ix,], xlab = "", ylab = "")
plot(responses_matrix[ix,], predicted_responses_matrix_rf_10f[ix,], xlab = "", ylab = "")
taxon = "Dinophysis"
ix = which(rownames(responses_matrix) == taxon)
cor.test(responses_matrix[ix,], predicted_responses_matrix_rf_10f[ix,], xlab = "", ylab = "")
plot(responses_matrix[ix,], predicted_responses_matrix_rf_10f[ix,], xlab = "", ylab = "")
taxon = "Teleaulax"
ix = which(rownames(responses_matrix) == taxon)
cor.test(responses_matrix[ix,], predicted_responses_matrix_rf_10f[ix,], xlab = "", ylab = "")
plot(responses_matrix[ix,], predicted_responses_matrix_rf_10f[ix,], xlab = "", ylab = "")

## for physiochem
par(mfrow = c(5,4), mar=c(2,2,1,1), xpd = TRUE) # good size of figure: 750 x 530
parameters_to_plot = c(
  "Longitude",     
  "Latitude",      
  "Salinity",      
  "Temperature",   
  "pH",            
  "Alkalinity",
  "Secchi_depth",  
  "SiO3",          
  "N_tot",         
  "DIN",           
  "NH4",           
  "NO3_NO2",       
  "P_tot",         
  "Phosphate",     
  "DOC",            
  "Humus",         
  "Chl")
ix = match(parameters_to_plot, rownames(responses_matrix))
for (i in 1:length(ix)) {
  plot(responses_matrix[ix[i],], predicted_responses_matrix_rf_ob[ix[i],], xlab = "Observed", ylab = "Predicted", main = paste(parameters_to_plot[i], " r=", round(cor_matr[ix[i],5], 2), sep = ""), pch = 1, col = "#3182bd", cex = 0.8)
  #plot(responses_matrix[i,], predicted_responses_matrix_rf_10f[i,], xlab = "Observed", ylab = "Predicted", main = paste(rownames(cor_matr)[i], " r=", round(cor_matr[i,4], 2), sep = ""), pch = 1, col = "#3182bd", cex = 0.8)
  #plot(responses_matrix[i,], predicted_responses_matrix_xgb[i,], xlab = "Observed", ylab = "Predicted", main = paste(rownames(cor_matr)[i], " r=", round(cor_matr[i,6], 2), sep = ""), pch = 1, col = "#3182bd", cex = 0.8)
}

par(mfrow = c(1,1), mar=c(4,4,1,1), xpd = TRUE) # good size of figure: 750 x 530
plot(cor_matr_physchem_16S[,4], cor_matr_physchem_18S[,4], xlim = c(0,1), ylim = c(0,1), xlab = "ML based on 16S", ylab = "ML based on 18S")
lines(c(0,1),c(0,1))
wilcox.test(cor_matr_physchem_16S[,4], cor_matr_physchem_18S[,4], paired = T)


### run 'predictions' based on sequence matching

features_matrix_full = norm_asv_counts_18S_genus
#responses_matrix_full = zoo_plan_genus
responses_matrix_full = phyt_plan_genus

use_alt_samplenames = 1 # for plankton
if (use_alt_samplenames == 1) {
  shared_samples = intersect(samples_alt, colnames(responses_matrix_full))
  ix = match(shared_samples, samples_alt)
} else {
  shared_samples = intersect(samples, colnames(responses_matrix_full))
  ix = match(shared_samples, samples)
}
features_matrix = features_matrix_full[,ix]
features_matrix = features_matrix[which(rowSums(features_matrix) > 0),]
ix = match(shared_samples, colnames(responses_matrix_full))
responses_matrix = responses_matrix_full[,ix]
responses_matrix = responses_matrix[which(rowSums(responses_matrix, na.rm = T) > 0),]

c(nrow(features_matrix), nrow(responses_matrix), length(intersect(rownames(features_matrix), rownames(responses_matrix))))

shared_genus = sort(intersect(rownames(features_matrix), rownames(responses_matrix)))
cor_matr_seq_match = matrix(ncol = 3, nrow = length(shared_genus))
colnames(cor_matr_seq_match) = c("pos_samp_micr","pos_samp_barc","cor-R")
rownames(cor_matr_seq_match) = shared_genus
for (i in 1:length(shared_genus)) {
  ixf = match(shared_genus[i], rownames(features_matrix))
  ixr = match(shared_genus[i], rownames(responses_matrix))
  cor_matr_seq_match[i,1] = length(which(responses_matrix[ixr,] > 0))
  cor_matr_seq_match[i,2] = length(which(features_matrix[ixf,] > 0))
  cor_matr_seq_match[i,3] = cor.test(responses_matrix[ixr,], features_matrix[ixf,])$est
}

plot(cor_matr_seq_match[,2], cor_matr_seq_match[,3], xlab = "Presence in #MB samples", ylab = "Prediction R-value", ylim = c(0,1))
plot(cor_matr_seq_match[,1], cor_matr_seq_match[,3], xlab = "Presence in #Mic samples", ylab = "Prediction R-value", ylim = c(0,1))
plot(cor_matr_seq_match[,1]*cor_matr_seq_match[,2], cor_matr_seq_match[,3], xlab = "Presence in #Mic samples", ylab = "Prediction R-value", ylim = c(0,1))

ix1 = intersect(which(cor_matr_seq_match[,1] >= 10), which(cor_matr_seq_match[,2] >= 10))
ix2 = match(shared_genus[ix1], rownames(cor_matr))
plot(cor_matr_seq_match[ix1,3], cor_matr[ix2,4], xlim = c(0,1), ylim = c(0,1), xlab = c("Directly from sequence annotation"), ylab = c("From machine-learning"), main = "Correlation coefficient per genus")
lines(c(0,1),c(0,1))

# example cor plots
par(mfcol = c(3,3), mar=c(2,3,1,1))
taxon = "Tripos"
ix1 = which(rownames(features_matrix) == taxon)
ix2 = which(rownames(responses_matrix) == taxon)
cor.test(responses_matrix[ix2,], features_matrix[ix1,])
plot(responses_matrix[ix2,], features_matrix[ix1,])
taxon = "Dinophysis"
ix1 = which(rownames(features_matrix) == taxon)
ix2 = which(rownames(responses_matrix) == taxon)
cor.test(responses_matrix[ix2,], features_matrix[ix1,])
plot(responses_matrix[ix2,], features_matrix[ix1,])
taxon = "Teleaulax"
ix1 = which(rownames(features_matrix) == taxon)
ix2 = which(rownames(responses_matrix) == taxon)
cor.test(responses_matrix[ix2,], features_matrix[ix1,])
plot(responses_matrix[ix2,], features_matrix[ix1,])


### plot gaussion distributions for presentation

par(mfcol = c(3,3))
lw = 3
curve(dnorm(x, 0, 1), from=-4, to=4, axes = F, xlab = "", ylab = "", col = "lightblue", lw = lw)
box()
curve(dnorm(x, 0, 1), from=-5, to=3, axes = F, xlab = "", ylab = "", col = "lightblue", lw = lw)
box()
curve(dnorm(x, 0, 1), from=-3, to=5, axes = F, xlab = "", ylab = "", col = "lightblue", lw = lw)
box()
curve(dnorm(x, 0, 1), from=-4, to=4, axes = F, xlab = "", ylab = "", col = "lightgreen", lw = lw)
box()
curve(dnorm(x, 0, 1), from=-5, to=3, axes = F, xlab = "", ylab = "", col = "lightgreen", lw = lw)
box()
curve(dnorm(x, 0, 1), from=-4, to=4, axes = F, xlab = "", ylab = "", col = "lightgreen", lw = lw)
box()
curve(dnorm(x, 0, 1), from=-5, to=3, axes = F, xlab = "", ylab = "", col = "orange", lw = lw)
box()
curve(dnorm(x, 0, 1), from=-4, to=4, axes = F, xlab = "", ylab = "", col = "orange", lw = lw)
box()
curve(dnorm(x, 0, 1), from=-3, to=5, axes = F, xlab = "", ylab = "", col = "orange", lw = lw)
box()



### shouldn't be needed anymore ###
# adjusting sample names
if (use_alt_samplenames_features == 1) {
  features_samples_full = alt_ids[match(colnames(features_matrix_full), alt_ids[,1]),2]
} else {
  features_samples_full = colnames(features_matrix_full)
}
if (use_alt_samplenames_responses == 0) {
  responses_samples_full = colnames(responses_matrix_full)
} else {
  responses_samples_full = alt_ids[match(colnames(responses_matrix_full), alt_ids[,1]),2]
}
colnames(features_matrix_full) = features_samples_full
colnames(responses_matrix_full) = responses_samples_full








