
## Set working directory
setwd("~/aquatic/envpredict/code")

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
library(caret)

## Set files
seqtab_file_18S = "../seq_data/combined/18S/filtered_seqtab_18S.tsv" ## count matrix for each ASV & sample
taxa_file_18S = "../seq_data/combined/18S/filtered_taxa_18S.tsv" ## taxonomic annotation for each ASV
## Set files
seqtab_file_18S_with_metazoa = "../seq_data/combined/18S/seqtab_18S.tsv" ## count matrix for each ASV & sample
taxa_file_18S_with_metazoa = "../seq_data/combined/18S/taxa_18S.tsv" ## taxonomic annotation for each ASV
seqtab_file_16S = "../seq_data/combined/16S/filtered_seqtab_16S.tsv" ## count matrix for each ASV & sample
taxa_file_16S =  "../seq_data/combined/16S/filtered_taxa_16S.tsv" ## taxonomic annotation for each ASV
phys_chem_file = "../env_data/combined/physical_chemical_processed_translation.tsv"
#bact_plan_file = "env_data/combined/bacterioplankton_processed.tsv"
#pico_plan_file = "env_data/combined/picoplankton_processed.tsv"
phyt_plan_file = "../env_data/combined/phytoplankton_processed.tsv"
zoop_plan_file = "../env_data/combined/zooplankton_processed.tsv"
# id_translation_file = "../env_data/combined/physical_chemical_processed_translation.tsv"

# read sequencing data
asv_counts_18S = as.matrix(read.delim(seqtab_file_18S, row.names = 1))
asv_taxa_18S = as.matrix(read.delim(taxa_file_18S, row.names = 1))[,1:9]
#asv_counts_18S_with_metazoa = as.matrix(read.delim(seqtab_file_18S_with_metazoa, row.names = 1))
#asv_taxa_18S_with_metazoa = as.matrix(read.delim(taxa_file_18S_with_metazoa, row.names = 1))[,1:9]

asv_counts_18S_with_metazoa = as.matrix(read.delim(seqtab_file_18S, row.names = 1)) # wrong !
asv_taxa_18S_with_metazoa = as.matrix(read.delim(taxa_file_18S, row.names = 1))[,1:9] # wrong !

asv_counts_16S = as.matrix(read.delim(seqtab_file_16S, row.names = 1))
asv_taxa_16S = as.matrix(read.delim(taxa_file_16S, row.names = 1))[,1:7]

colnames(asv_counts_18S) = gsub("^X", "", colnames(asv_counts_18S))
colnames(asv_counts_18S_with_metazoa) = gsub("^X", "", colnames(asv_counts_18S_with_metazoa))
colnames(asv_counts_16S) = gsub("^X", "", colnames(asv_counts_16S))
rownames(asv_counts_18S) = paste(rownames(asv_counts_18S), "18S", sep = "_")
rownames(asv_taxa_18S) = paste(rownames(asv_taxa_18S), "18S", sep = "_")
rownames(asv_counts_18S_with_metazoa) = paste(rownames(asv_counts_18S_with_metazoa), "18S", sep = "_")
rownames(asv_taxa_18S_with_metazoa) = paste(rownames(asv_taxa_18S_with_metazoa), "18S", sep = "_")
rownames(asv_counts_16S) = paste(rownames(asv_counts_16S), "16S", sep = "_")
rownames(asv_taxa_16S) = paste(rownames(asv_taxa_16S), "16S", sep = "_")

identical(rownames(asv_counts_18S), rownames(asv_taxa_18S))
identical(rownames(asv_counts_16S), rownames(asv_taxa_16S))
identical(colnames(asv_counts_18S), colnames(asv_counts_18S_with_metazoa))
identical(colnames(asv_counts_16S), colnames(asv_counts_18S))
identical(colnames(asv_counts_16S), colnames(asv_counts_18S_with_metazoa))
colnames(asv_counts_16S) = gsub("^X", "", colnames(asv_counts_16S))
samples = colnames(asv_counts_16S)


# read physchem
phys_chem = read.delim(phys_chem_file)

## Match the original and alternative sample ideas
alt_ids = cbind(phys_chem$sample_id, phys_chem$station_id_date)
samples_alt = alt_ids[match(samples, alt_ids[,1]),2]

## Remove without station ID
ix = setdiff(1:length(samples_alt), grep("NA_", samples_alt)) # samples with station id
asv_counts_16S = asv_counts_16S[,ix]
asv_counts_18S = asv_counts_18S[,ix]
asv_counts_18S_with_metazoa = asv_counts_18S_with_metazoa[,ix]
samples = samples[ix]
samples_alt = samples_alt[ix]
colnames(asv_counts_16S) = colnames(asv_counts_18S) = colnames(asv_counts_18S_with_metazoa) = samples_alt


norm_asv_counts_16S = t(t(asv_counts_16S)/colSums(asv_counts_16S))
norm_asv_counts_18S = t(t(asv_counts_18S)/colSums(asv_counts_18S))
norm_asv_counts_18S_with_metazoa = t(t(asv_counts_18S_with_metazoa)/colSums(asv_counts_18S_with_metazoa))

## Prepare physchem data
rownames(phys_chem) = phys_chem$sample_id
year = year(phys_chem$date)
yday = yday(phys_chem$date)
ylength = rep(NA, length(yday))

time = phys_chem$time_h
time = interval(start = strptime("00:00:00", "%H:%M:%S"), end = strptime(time, "%H:%M:%S"))
time = time_length(time, unit = 'minutes')
min_since_midnight = time
time_xcord = cos(2*pi*time/(24*60))
time_ycord = sin(2*pi*time/(24*60))

ix = which(!is.na(yday))
ylength[ix] = yday(paste(as.character(year[ix]), "-12-31", sep = ""))
yangle = 2*pi*yday/ylength
yday_xcord = cos(yangle)
yday_ycord = sin(yangle)

phys_chem = t(cbind(
  year, 
  yday, 
  yday_xcord, 
  yday_ycord, 
  min_since_midnight,
  time_xcord, 
  time_ycord,
  phys_chem[,c('Longitude','Latitude','Salinity','Temperature','pH','Alkalinity','Secchi_depth','SiO3','N_tot','DIN','NH4','NO2','NO3','NO3_NO2','P_tot','Phosphate','DOC','Humus','Chl')]
))

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
use_alternative_sample_names <- function(matrix) {
  ix = match(colnames(matrix), alt_ids[,1])
  colnames(matrix) = alt_ids[ix,2]
  return(matrix)
}

extract_shared_samples <- function(features_matrix_full, responses_matrix_full) { # selecting only samples shared by features and responses
  return_matrices = list()
  shared_samples = intersect(colnames(responses_matrix_full), colnames(features_matrix_full))
  if (length(grep("NA_", shared_samples)) > 0) {
    shared_samples = shared_samples[-grep("NA_", shared_samples)]
  }
  ix = match(shared_samples, colnames(features_matrix_full))
  features_matrix = features_matrix_full[,ix]
  ix = match(shared_samples, colnames(responses_matrix_full))
  responses_matrix = responses_matrix_full[,ix]
  responses_matrix_bin = responses_matrix
  responses_matrix_bin[which(responses_matrix > 0)] = 1
  responses_matrix_bin[which(responses_matrix < 0)] = 1
  responses_matrix = responses_matrix[which(rowSums(responses_matrix_bin, na.rm = T) > 0),] # only include rows with at least one non-NA and non-0 values
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

## Functions for running ML predictions
run_randomforest_out_of_bag <- function(features_matrix, responses_matrix) {
  predicted_responses_matrix = responses_matrix
  predicted_responses_matrix[,] = NA
  rownames(features_matrix) = paste("f", 1:nrow(features_matrix), sep= "")
  for (i in 1:nrow(responses_matrix)) {
    ix = which(!is.na(responses_matrix[i,]))
    df = as.data.frame(cbind(responses_matrix[i,ix], t(features_matrix[,ix])))
    #rf = ranger(V1 ~ ., data = df, num.trees = 2000, importance = 'impurity')
    rf = ranger(V1 ~ ., data = df, num.trees = 2000, importance = 'none')
    predicted_responses_matrix[i,ix] = rf$predictions
  }
  return(predicted_responses_matrix)
}

run_randomforest <- function(features_matrix, responses_matrix, numfolds, min_samples) {
  min_samples_in_fold = 1
  predicted_responses_matrix = responses_matrix
  predicted_responses_matrix[,] = NA
  rownames(features_matrix) = paste("f", 1:nrow(features_matrix), sep= "")
  for (i in 1:nrow(responses_matrix)) {
    response = responses_matrix[i,]
    if (length(which(!is.na(response))) < min_samples) { next } # skip this parameter if too few non-NA samples
    df = as.data.frame(cbind(response, t(features_matrix)))
    folds <- vfold_cv(df, v = numfolds, strata = response)
    #folds <- vfold_cv(df, v = numfolds, strata = NULL)
    #for (j in 1:nrow(foldxs)) {
    for (j in 1:nrow(folds)) {
      trainingIndex <- folds$splits[[j]]$in_id
      testingIndex <- setdiff(1:length(response), trainingIndex)
      trainingIndex = trainingIndex[which(!is.na(response[trainingIndex]))] # remove NA response samples
      if (length(trainingIndex) >= min_samples_in_fold) {
        #rf = ranger(response ~ ., data = df[trainingIndex,], num.trees = 2000, importance = 'impurity')
        rf = ranger(response ~ ., data = df[trainingIndex,], num.trees = 2000, importance = 'none')
        predicted_responses_matrix[i,testingIndex] = predict(rf, df[testingIndex,])$predictions  
      }
    }
  }
  return(predicted_responses_matrix)
}

## come here 1
predict_randomforest <- function(features_matrix_train, responses_matrix_train, features_matrix_target, min_samples = 10) {
  if (!identical(rownames(features_matrix_train), rownames(features_matrix_target))) {
    stop("Rownames of the feature matrices (i.e., ASVs/taxa) do not match!")
  }
  predicted_responses_matrix = matrix(ncol = ncol(features_matrix_target), nrow = nrow(responses_matrix_train))
  colnames(predicted_responses_matrix) = colnames(features_matrix_target)
  rownames(predicted_responses_matrix) = rownames(responses_matrix_train)
  for (i in 1:nrow(responses_matrix_train)) {
    response = responses_matrix_train[i,]
    if (length(which(!is.na(response))) < min_samples) { next } # skip this parameter if too few non-NA samples
    df = as.data.frame(cbind(response, t(features_matrix_train)))
    df = df[which(!is.na(response)),]
    #rf = ranger(response ~ ., data = df[trainingIndex,], num.trees = 2000, importance = 'impurity')
    rf = ranger(response ~ ., data = df, num.trees = 2000, importance = 'none')
    predicted_responses_matrix[i,] = predict(rf, t(features_matrix_target))$predictions  
  }
  return(predicted_responses_matrix)
}


run_gbm <- function(features_matrix, responses_matrix, numfolds, min_samples) {
  min_samples_in_fold = 1
  predicted_responses_matrix = responses_matrix
  predicted_responses_matrix[,] = NA
  rownames(features_matrix) = paste("f", 1:nrow(features_matrix), sep= "")
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

run_xgboost_caret <- function(features_matrix, responses_matrix, numfolds_outer, numfolds_inner, min_samples) {
  output_files_path = output_files_path = "../output/physchem_based_on_seqdata_xgboost"
  outfile_actual = "norm_seqtab_16S_xgboost_Actual.tsv"
  outfile_predicted = "norm_seqtab_16S_xgboost_Predictions.tsv"
  
  train_control = trainControl(method = "cv", number = numfolds_inner)
  min_samples_in_fold = 1
  predicted_responses_matrix = responses_matrix
  predicted_responses_matrix[,] = NA
  rownames(features_matrix) = paste("f", 1:nrow(features_matrix), sep= "")
  for (i in 1:nrow(responses_matrix)) {
  #for (i in 1:3) {  
    response = responses_matrix[i,]
    if (length(which(!is.na(response))) < min_samples) { next } # skip this parameter if too few non-NA samples
    df = as.data.frame(cbind(response, t(features_matrix)))
    folds <- vfold_cv(df, v = numfolds_outer, strata = response)
    df = as.data.frame(t(features_matrix))
    for (j in 1:nrow(folds)) { # outer fold
      trainingIndex <- folds$splits[[j]]$in_id 
      testingIndex <- setdiff(1:length(response), trainingIndex)
      trainingIndex = trainingIndex[which(!is.na(response[trainingIndex]))] # remove NA response samples
      if (length(trainingIndex) >= min_samples_in_fold) {
        xgb <- caret::train(
          x = df[trainingIndex,],
          y = response[trainingIndex],
          trControl = train_control,
          #tuneGrid = grid_default,
          #method = "xgbTree",
          method = "ranger",
          modelType = "regression",
          verbose = TRUE
        )
        predicted_responses_matrix[i,testingIndex] = predict(xgb, df[testingIndex,])
      }
      #write.table(responses_matrix, paste(output_files_path, outfile_actual, sep = "/"), sep="\t")
      #write.table(predicted_responses_matrix, paste(output_files_path, outfile_predicted, sep = "/"), sep="\t")
    }
  }
  return(predicted_responses_matrix)
}

## Function for getting correlations between predicted and meaured variables
get_correlations_predictions <- function(responses_matrix, predicted_responses_matrix) {
  cor_matr = matrix(ncol = 4, nrow = nrow(responses_matrix))
  rownames(cor_matr) = rownames(responses_matrix)
  colnames(cor_matr) = c("#non-NA_samples","#non-0_samples", "pearson", "spearman")
  for (i in 1:nrow(responses_matrix)) {
    cor_matr[i,1] = length(which(!is.na(responses_matrix[i,])))
    cor_matr[i,2] = length(which(responses_matrix[i,] > 0))
    cor_matr[i,3] = round(cor.test(responses_matrix[i,], predicted_responses_matrix[i,], method = "pearson")$est , 3)
    cor_matr[i,4] = round(cor.test(responses_matrix[i,], predicted_responses_matrix[i,], method = "spearman")$est , 3)
  }
  return(cor_matr)
}

make_scatterplots_actual_vs_predicted <- function(responses_matrix, predicted_responses_matrix) {
  parameters_to_plot = c(
    "Longitude", "Latitude", 
    "yday_xcord", "yday_ycord",
    "time_xcord", "time_ycord",
    "Salinity", "Temperature",   
    "pH", "Alkalinity",
    "Secchi_depth", "SiO3",          
    "N_tot", "DIN",           
    "NH4", "NO3_NO2",       
    "P_tot", "Phosphate",     
    "DOC","Humus",         
    "Chl"
  )
  len = ceiling(sqrt(length(parameters_to_plot)))
  par(mfrow = c(len,len), mar=c(2,2,1,1), xpd = TRUE) # good size of figure: 750 x 530
  ix = match(parameters_to_plot, rownames(responses_matrix))
  for (i in 1:length(ix)) {
    n = length(which(!is.na(responses_matrix[ix[i],])))
    c = round(cor.test(responses_matrix[ix[i],], predicted_responses_matrix_rf_ob[ix[i],])$est, 2)
    plot(responses_matrix[ix[i],], predicted_responses_matrix_rf_ob[ix[i],], xlab = "Observed", ylab = "Predicted", main = paste(parameters_to_plot[i], "r=", c, "n=", n), pch = 1, col = "#3182bd", cex = 0.8)
  }
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
#features_matrix_full = rbind(norm_asv_counts_16S, norm_asv_counts_18S)
responses_matrix_full = phys_chem
#responses_matrix_full = phys_chem[c(10,12,26),]
features_matrix = extract_shared_samples(features_matrix_full, responses_matrix_full)$features_matrix
responses_matrix = extract_shared_samples(features_matrix_full, responses_matrix_full)$responses_matrix
features_matrix = do_feature_selection(features_matrix, 0.1)

predicted_responses_matrix_rf_ob = run_randomforest_out_of_bag(features_matrix, responses_matrix)
predicted_responses_matrix_rf_10f = run_randomforest(features_matrix, responses_matrix, 10, 10)
predicted_responses_matrix_xgb = run_xgboost_caret(features_matrix, responses_matrix, numfolds_outer = 10, numfolds_inner = 5, min_samples = 10)

cor_ob = get_correlations_predictions(responses_matrix, predicted_responses_matrix_rf_ob)
cor_10f = get_correlations_predictions(responses_matrix, predicted_responses_matrix_rf_10f)
cor_xgb = get_correlations_predictions(responses_matrix, predicted_responses_matrix_xgb)

## for using plankton as features and physchem as responses
#features_matrix_full = zoo_plan_genus
features_matrix_full = phyt_plan_genus
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



### Running predictions to be used in paper ###

## Running physiochem predictions on deep feature representation files
output_files_path = "../output/RepresentationsFromDeepMicro"
if (!dir.exists(output_files_path)) { dir.create(output_files_path) }
features_files_path = "../seq_data/combined/RepresentationsFromDeepMicro"
list.files(features_files_path)
infiles = sort(list.files(features_files_path))
for (i in 1:length(infiles)) {
  outfile_actual = paste(gsub(".csv$","",infiles[i]),"RF10fold_Actual.tsv",sep ="_")
  outfile_predicted = paste(gsub(".csv$","",infiles[i]),"RF10fold_Predictions.tsv",sep ="_")
  responses_matrix_full = phys_chem
  features_matrix_full = as.matrix(read.delim(paste(features_files_path, infiles[i], sep="/"), row.names = 1))
  features_matrix_full = t(features_matrix_full)
  features_matrix_full = use_alternative_sample_names(features_matrix_full)
  features_matrix = extract_shared_samples(features_matrix_full, responses_matrix_full)$features_matrix
  responses_matrix = extract_shared_samples(features_matrix_full, responses_matrix_full)$responses_matrix
  #features_matrix = do_feature_selection(features_matrix, 0.1)
  predicted_responses_matrix = run_randomforest(features_matrix, responses_matrix, 10, 10)
  write.table(responses_matrix, paste(output_files_path, outfile_actual, sep = "/"), sep="\t")
  write.table(predicted_responses_matrix, paste(output_files_path, outfile_predicted, sep = "/"), sep="\t")
}

## Running physiochem predictions on seq files for different taxonomic levels
# 18S levels: Domain	Supergroup	Division	Subdivision	Class	Order	Family	Genus	Species
# 16S levels: Domain	Phylum	Class	Order	Family Genus	Species
output_files_path = "../output/DifferentTaxonomicLevels"
if (!dir.exists(output_files_path)) { dir.create(output_files_path) }
features_files_path_16S = "../seq_data/combined/16S"
features_files_path_18S = "../seq_data/combined/18S"
infiles = sort(c(list.files(features_files_path_16S, pattern="norm_.+tsv$", full.names = TRUE), list.files(features_files_path_18S, pattern="norm_.+tsv$", full.names = TRUE))) # only include files starting with norm_
infiles = infiles[-grep("_1\\.tsv$", infiles)]
for (i in 1:length(infiles)) {
  outfile_actual = gsub(".tsv$","_RF10fold_Actual.tsv",basename(infiles[i])) ## Do we really neeed a separate actual table for each option?
  outfile_predicted = gsub(".tsv$","_RF10fold_Predictions.tsv",basename(infiles[i]))
  responses_matrix_full = phys_chem
  features_matrix_full = as.matrix(read.delim(infiles[i], row.names = 1))
  colnames(features_matrix_full) = gsub("^X", "", colnames(features_matrix_full))
  #features_matrix_full = t(features_matrix_full)
  features_matrix_full = use_alternative_sample_names(features_matrix_full)
  features_matrix = extract_shared_samples(features_matrix_full, responses_matrix_full)$features_matrix
  responses_matrix = extract_shared_samples(features_matrix_full, responses_matrix_full)$responses_matrix
  features_matrix = do_feature_selection(features_matrix, 0.1)
  predicted_responses_matrix = run_randomforest(features_matrix, responses_matrix, 10, 10)
  write.table(responses_matrix, paste(output_files_path, outfile_actual, sep = "/"), sep="\t")
  write.table(predicted_responses_matrix, paste(output_files_path, outfile_predicted, sep = "/"), sep="\t")
}

## Running physiochem predictions on phytoplankton (based on all and based on genera only) and on seq data, using the same set of samples
output_files_path = "../output/physchem_based_on_seqdata_or_phytoplan"
responses_matrix_full = phys_chem
# based on phyt_plan
features_matrix_full_1 = phyt_plan
features_matrix_full_2 = norm_asv_counts_16S
these_sample = intersect(colnames(features_matrix_full_1), colnames(features_matrix_full_2))
ix = match(these_sample, colnames(features_matrix_full_1))
features_matrix_full_1 = features_matrix_full_1[, ix]
features_matrix = extract_shared_samples(features_matrix_full_1, responses_matrix_full)$features_matrix
features_matrix = do_feature_selection(features_matrix, 0.1)
responses_matrix = extract_shared_samples(features_matrix_full_1, responses_matrix_full)$responses_matrix
features_matrix = do_feature_selection(features_matrix, 0.1)
predicted_responses_matrix = run_randomforest(features_matrix, responses_matrix, 10, 10)
outfile_actual = "phytoplan-based_physchem_Actual.tsv"
write.table(responses_matrix, paste(output_files_path, outfile_actual, sep = "/"), sep="\t")
outfile_predicted = "phytoplan-based_physchem_Predictions.tsv"
write.table(predicted_responses_matrix, paste(output_files_path, outfile_predicted, sep = "/"), sep="\t")
# based on phyt_plan_genus
features_matrix_full_1 = phyt_plan_genus
features_matrix_full_2 = norm_asv_counts_16S
these_sample = intersect(colnames(features_matrix_full_1), colnames(features_matrix_full_2))
ix = match(these_sample, colnames(features_matrix_full_1))
features_matrix_full_1 = features_matrix_full_1[, ix]
features_matrix = extract_shared_samples(features_matrix_full_1, responses_matrix_full)$features_matrix
features_matrix = do_feature_selection(features_matrix, 0.1)
responses_matrix = extract_shared_samples(features_matrix_full_1, responses_matrix_full)$responses_matrix
features_matrix = do_feature_selection(features_matrix, 0.1)
predicted_responses_matrix = run_randomforest(features_matrix, responses_matrix, 10, 10)
outfile_actual = "phytoplan_genus-based_physchem_Actual.tsv"
write.table(responses_matrix, paste(output_files_path, outfile_actual, sep = "/"), sep="\t")
outfile_predicted = "phytoplan_genus-based_physchem_Predictions.tsv"
write.table(predicted_responses_matrix, paste(output_files_path, outfile_predicted, sep = "/"), sep="\t")
# based on 16S ASVs
ix = match(these_sample, colnames(features_matrix_full_2))
features_matrix_full_2 = features_matrix_full_2[, ix]
features_matrix = extract_shared_samples(features_matrix_full_2, responses_matrix_full)$features_matrix
features_matrix = do_feature_selection(features_matrix, 0.1)
responses_matrix = extract_shared_samples(features_matrix_full_2, responses_matrix_full)$responses_matrix
features_matrix = do_feature_selection(features_matrix, 0.1)
predicted_responses_matrix = run_randomforest(features_matrix, responses_matrix, 10, 10)
outfile_actual = "16S-based_physchem_Actual.tsv"
write.table(responses_matrix, paste(output_files_path, outfile_actual, sep = "/"), sep="\t")
outfile_predicted = "16S-based_physchem_Predictions.tsv"
write.table(predicted_responses_matrix, paste(output_files_path, outfile_predicted, sep = "/"), sep="\t")
ncol(responses_matrix) # 328

## Running physiochem predictions on zooplankton (based on all and based on genera only) and on seq data, using the same set of samples
output_files_path = "../output/physchem_based_on_seqdata_or_zooplan"
responses_matrix_full = phys_chem
# based on zoo_plan
features_matrix_full_1 = zoo_plan
features_matrix_full_2 = norm_asv_counts_16S
these_sample = intersect(colnames(features_matrix_full_1), colnames(features_matrix_full_2))
ix = match(these_sample, colnames(features_matrix_full_1))
features_matrix_full_1 = features_matrix_full_1[, ix]
features_matrix = extract_shared_samples(features_matrix_full_1, responses_matrix_full)$features_matrix
features_matrix = do_feature_selection(features_matrix, 0.1)
responses_matrix = extract_shared_samples(features_matrix_full_1, responses_matrix_full)$responses_matrix
features_matrix = do_feature_selection(features_matrix, 0.1)
predicted_responses_matrix = run_randomforest(features_matrix, responses_matrix, 10, 10)
outfile_actual = "zooplan-based_physchem_Actual.tsv"
write.table(responses_matrix, paste(output_files_path, outfile_actual, sep = "/"), sep="\t")
outfile_predicted = "zooplan-based_physchem_Predictions.tsv"
write.table(predicted_responses_matrix, paste(output_files_path, outfile_predicted, sep = "/"), sep="\t")
# based on zoo_plan_genus
features_matrix_full_1 = zoo_plan_genus
features_matrix_full_2 = norm_asv_counts_16S
these_sample = intersect(colnames(features_matrix_full_1), colnames(features_matrix_full_2))
ix = match(these_sample, colnames(features_matrix_full_1))
features_matrix_full_1 = features_matrix_full_1[, ix]
features_matrix = extract_shared_samples(features_matrix_full_1, responses_matrix_full)$features_matrix
features_matrix = do_feature_selection(features_matrix, 0.1)
responses_matrix = extract_shared_samples(features_matrix_full_1, responses_matrix_full)$responses_matrix
features_matrix = do_feature_selection(features_matrix, 0.1)
predicted_responses_matrix = run_randomforest(features_matrix, responses_matrix, 10, 10)
outfile_actual = "zooplan_genus-based_physchem_Actual.tsv"
write.table(responses_matrix, paste(output_files_path, outfile_actual, sep = "/"), sep="\t")
outfile_predicted = "zooplan_genus-based_physchem_Predictions.tsv"
write.table(predicted_responses_matrix, paste(output_files_path, outfile_predicted, sep = "/"), sep="\t")
# based on 16S ASVs
ix = match(these_sample, colnames(features_matrix_full_2))
features_matrix_full_2 = features_matrix_full_2[, ix]
features_matrix = extract_shared_samples(features_matrix_full_2, responses_matrix_full)$features_matrix
features_matrix = do_feature_selection(features_matrix, 0.1)
responses_matrix = extract_shared_samples(features_matrix_full_2, responses_matrix_full)$responses_matrix
features_matrix = do_feature_selection(features_matrix, 0.1)
predicted_responses_matrix = run_randomforest(features_matrix, responses_matrix, 10, 10)
outfile_actual = "16S-based_physchem_Actual.tsv"
write.table(responses_matrix, paste(output_files_path, outfile_actual, sep = "/"), sep="\t")
outfile_predicted = "16S-based_physchem_Predictions.tsv"
write.table(predicted_responses_matrix, paste(output_files_path, outfile_predicted, sep = "/"), sep="\t")
ncol(responses_matrix) # 241


## Running predictions of phytoplankton genera in microscopy and from 18S data, respectively, using physiochem data
output_files_path = "../output/phytoplankton_genus_micr_and_18S_based_on_physchem"
features_matrix_full = phys_chem
responses_matrix_full_1 = norm_asv_counts_18S_genus
responses_matrix_full_2 = phyt_plan_genus
shared_samples = intersect(colnames(responses_matrix_full_1), colnames(responses_matrix_full_2))
responses_matrix_full_1 = responses_matrix_full_1[,shared_samples]
responses_matrix_full_2 = responses_matrix_full_2[,shared_samples]
features_matrix = features_matrix_full[,shared_samples]
features_to_include = c("yday_xcord",
  "yday_ycord",
  "time_xcord",
  "time_ycord",
  "Longitude",
  "Latitude",
  "Salinity",
  "Temperature",
  "SiO3", 
  "N_tot", 
  "DIN", 
  "NH4",
  "NO2",
  "NO3",
  "P_tot",
  "Phosphate",
  "Chl"
)
features_matrix = features_matrix[features_to_include,]
features_matrix = features_matrix[, which(complete.cases(t(features_matrix[features_to_include,])))]
pheatmap(features_matrix)
shared_samples = colnames(features_matrix)
responses_matrix_1 = responses_matrix_full_1[,shared_samples]
responses_matrix_2 = responses_matrix_full_2[,shared_samples]
length(shared_samples) # 

predicted_responses_matrix = run_randomforest(features_matrix, responses_matrix_1, 10, 10)
outfile_actual = "physchem-based_18S-genus_Actual.tsv"
write.table(responses_matrix_1, paste(output_files_path, outfile_actual, sep = "/"), sep="\t")
outfile_predicted = "physchem-based_18S-genus_Predictions.tsv"
write.table(predicted_responses_matrix, paste(output_files_path, outfile_predicted, sep = "/"), sep="\t")

predicted_responses_matrix = run_randomforest(features_matrix, responses_matrix_2, 10, 10)
outfile_actual = "physchem-based_phytoplan-genus_Actual.tsv"
write.table(responses_matrix_2, paste(output_files_path, outfile_actual, sep = "/"), sep="\t")
outfile_predicted = "physchem-based_phytoplan-genus_Predictions.tsv"
write.table(predicted_responses_matrix, paste(output_files_path, outfile_predicted, sep = "/"), sep="\t")



## Running phytoplankton predictions (genera only?) based on seq data, seq data sequence matching, and physchem data, using the same set of samples

output_files_path = "../output/phytoplankton_predicted/"
if (!dir.exists(output_files_path)) { dir.create(output_files_path) }


feature_tables = list(norm_asv_counts_18S, norm_asv_counts_16S, phys_chem)
names(feature_tables) = c("norm_asv_counts_18S", "norm_asv_counts_16S", "phys_chem")
response_tables = list(phyt_plan, phyt_plan_genus)
names(response_tables) = c("phyt_plan", "phyt_plan_genus")

## Get common cols across all the tables

feature_col_names = lapply(feature_tables, colnames)
response_col_names = lapply(response_tables, colnames)
all_col_names = c(feature_col_names, response_col_names)
cols_to_keep = Reduce(intersect, all_col_names)
for(i in 1:length(feature_tables)) {
  feature_tables[[i]] = feature_tables[[i]][,cols_to_keep]
}
for(i in 1:length(response_tables)) {
  response_tables[[i]] = response_tables[[i]][,cols_to_keep]
}

## Choose only the physicochemical paramters available for more than 70% of the samples
## Threshold based on the hustogram plotted using the first code line below
phys_chem_complete = feature_tables[['phys_chem']]

hist(rowSums(! is.na(phys_chem_complete))/ncol(phys_chem_complete), breaks = 20)
sort(rowSums(! is.na(phys_chem_complete))/ncol(phys_chem_complete))
ix_phys_chem = which(rowSums(! is.na(phys_chem_complete))/ncol(phys_chem_complete) > 0.80)
iy_phys_chem = which(complete.cases(t(phys_chem_complete[ix_phys_chem,])))

phys_chem_complete = phys_chem_complete[ix_phys_chem,iy_phys_chem]

ncol(phys_chem_complete)/ncol(phys_chem)

feature_tables[['phys_chem']] = phys_chem_complete

cols_to_keep = colnames(phys_chem_complete)
for(i in 1:length(feature_tables)) {
  feature_tables[[i]] = feature_tables[[i]][,cols_to_keep]
}
for(i in 1:length(response_tables)) {
  response_tables[[i]] = response_tables[[i]][,cols_to_keep]
}


for (i in 1:length(response_tables)) {
  response_matrix = response_tables[[i]]
  response_name = names(response_tables)[i]
  for (j in 1:length(feature_tables)) {
    feature_matrix = feature_tables[[j]]
    feature_name = names(feature_tables)[j]
    feature_matrix = do_feature_selection(feature_matrix, 0.1)
    predicted_responses_matrix_rf_ob = run_randomforest_out_of_bag(feature_matrix, response_matrix)
    write.table(response_matrix, paste(output_files_path, paste(feature_name, "_", response_name,  "_RF10fold_Actual.tsv", sep = ""), sep = "/"), sep="\t")
    write.table(predicted_responses_matrix_rf_ob, paste(output_files_path, paste(feature_name, "_", response_name, "_RF10fold_Predictions.tsv", sep = ""), sep = "/"), sep="\t")
  }
}

## Relative abundance based on matching - only matching ASVs left

features_matrix = norm_asv_counts_18S_genus[,cols_to_keep]

responses_matrix = phyt_plan_genus[,cols_to_keep]

c(nrow(features_matrix), nrow(responses_matrix), length(intersect(rownames(features_matrix), rownames(responses_matrix))))

shared_genera = sort(intersect(rownames(features_matrix), rownames(responses_matrix)))

predicted_matrix = features_matrix[shared_genera,]
responses_matrix = responses_matrix[shared_genera,]

write.table(responses_matrix, paste(output_files_path, 'direct_matching_norm_clade_counts_18S_8_phyt_plan_genus_Actual.tsv', sep = "/"), sep="\t")
write.table(predicted_matrix, paste(output_files_path, 'direct_matching_norm_clade_counts_18S_8_phyt_plan_genus_Predictions.tsv', sep = "/"), sep="\t")

features_matrix = features_matrix[shared_genera,]
renorm_matching_abundance = features_matrix

for(k in 1:ncol(renorm_matching_abundance)){
  renorm_matching_abundance[,k] = renorm_matching_abundance[,k]/sum(renorm_matching_abundance[,k])
}

predicted_matrix = renorm_matching_abundance

write.table(responses_matrix, paste(output_files_path, 'renomralized_direct_matching_norm_clade_counts_18S_8_phyt_plan_genus_Actual.tsv', sep = "/"), sep="\t")
write.table(predicted_matrix, paste(output_files_path, 'renomralized_direct_matching_norm_clade_counts_18S_8_phyt_plan_genus_Predictions.tsv', sep = "/"), sep="\t")

## Running zooplankton predictions (genera only?) based on seq data and physchem data, using the same set of samples

output_files_path = "../output/zooplankton_predicted/"
if (!dir.exists(output_files_path)) { dir.create(output_files_path) }


feature_tables = list(norm_asv_counts_18S, norm_asv_counts_16S, phys_chem)
names(feature_tables) = c("norm_asv_counts_18S", "norm_asv_counts_16S", "phys_chem")
response_tables = list(zoo_plan, zoo_plan_genus)
names(response_tables) = c("zoo_plan", "zoo_plan_genus")

## Get common cols across all the tables

feature_col_names = lapply(feature_tables, colnames)
response_col_names = lapply(response_tables, colnames)
all_col_names = c(feature_col_names, response_col_names)
cols_to_keep = Reduce(intersect, all_col_names)
for(i in 1:length(feature_tables)) {
  feature_tables[[i]] = feature_tables[[i]][,cols_to_keep]
}
for(i in 1:length(response_tables)) {
  response_tables[[i]] = response_tables[[i]][,cols_to_keep]
}

## Choose only the physicochemical paramters available for more than 70% of the samples
## Threshold based on the hustogram plotted using the first code line below
phys_chem_complete = feature_tables[['phys_chem']]

hist(rowSums(! is.na(phys_chem_complete))/ncol(phys_chem_complete), breaks = 20)
sort(rowSums(! is.na(phys_chem_complete))/ncol(phys_chem_complete))
ix_phys_chem = which(rowSums(! is.na(phys_chem_complete))/ncol(phys_chem_complete) > 0.80)
iy_phys_chem = which(complete.cases(t(phys_chem_complete[ix_phys_chem,])))

phys_chem_complete = phys_chem_complete[ix_phys_chem,iy_phys_chem]

ncol(phys_chem_complete)/ncol(phys_chem)

feature_tables[['phys_chem']] = phys_chem_complete

cols_to_keep = colnames(phys_chem_complete)
for(i in 1:length(feature_tables)) {
  feature_tables[[i]] = feature_tables[[i]][,cols_to_keep]
}
for(i in 1:length(response_tables)) {
  response_tables[[i]] = response_tables[[i]][,cols_to_keep]
}


for (i in 1:length(response_tables)) {
  response_matrix = response_tables[[i]]
  response_name = names(response_tables)[i]
  for (j in 1:length(feature_tables)) {
    feature_matrix = feature_tables[[j]]
    feature_name = names(feature_tables)[j]
    feature_matrix = do_feature_selection(feature_matrix, 0.1)
    predicted_responses_matrix_rf_ob = run_randomforest_out_of_bag(feature_matrix, response_matrix)
    write.table(response_matrix, paste(output_files_path, paste(feature_name, "_", response_name,  "_RF10fold_Actual.tsv", sep = ""), sep = "/"), sep="\t")
    write.table(predicted_responses_matrix_rf_ob, paste(output_files_path, paste(feature_name, "_", response_name, "_RF10fold_Predictions.tsv", sep = ""), sep = "/"), sep="\t")
  }
}

## Read the asv_counts with metazoa

matching_norm_asv_counts_18S_with_metazoa = norm_asv_counts_18S_with_metazoa[,cols_to_keep]
identical(colnames(matching_norm_asv_counts_18S_with_metazoa), colnames(response_tables[['zoo_plan_genus']]))

## Get normalized genus counts

with_metazoa_genus = NULL
i = (ncol(asv_taxa_18S_with_metazoa)-1)
clade = unique(asv_taxa_18S_with_metazoa[,i])
clade = clade[!is.na(clade)]
for (j in 1:length(clade)) {
  ix = which(clade[j]==asv_taxa_18S_with_metazoa[,i])
  if (length(ix) > 1) {
    with_metazoa_genus = rbind(with_metazoa_genus, apply(matching_norm_asv_counts_18S_with_metazoa[ix,], 2, sum, na.rm=TRUE))
  } else {
    with_metazoa_genus = rbind(with_metazoa_genus, matching_norm_asv_counts_18S_with_metazoa[ix,])
  }
}
rownames(with_metazoa_genus) = clade
colnames(with_metazoa_genus) = colnames(matching_norm_asv_counts_18S_with_metazoa)

## Relative abundance based on matching - only matching ASVs left, only Eukaryotes?

features_matrix = with_metazoa_genus

responses_matrix = zoo_plan_genus[,cols_to_keep]

c(nrow(features_matrix), nrow(responses_matrix), length(intersect(rownames(features_matrix), rownames(responses_matrix))))

shared_genera = sort(intersect(rownames(features_matrix), rownames(responses_matrix)))

predicted_matrix = features_matrix[shared_genera,]
responses_matrix = responses_matrix[shared_genera,]

write.table(responses_matrix, paste(output_files_path, 'direct_matching_norm_clade_counts_18S_8_zoo_plan_genus_Actual.tsv', sep = "/"), sep="\t")
write.table(predicted_matrix, paste(output_files_path, 'direct_matching_norm_clade_counts_18S_8_zoo_plan_genus_Predictions.tsv', sep = "/"), sep="\t")

features_matrix = features_matrix[shared_genera,]
renorm_matching_abundance = features_matrix

for(k in 1:ncol(renorm_matching_abundance)){
  renorm_matching_abundance[,k] = renorm_matching_abundance[,k]/sum(renorm_matching_abundance[,k])
}

predicted_matrix = renorm_matching_abundance

write.table(responses_matrix, paste(output_files_path, 'renomralized_direct_matching_norm_clade_counts_18S_8_zoo_plan_genus_18S_Actual.tsv', sep = "/"), sep="\t")
write.table(predicted_matrix, paste(output_files_path, 'renomralized_direct_matching_norm_clade_counts_18S_8_zoo_plan_genus_18S_Predictions.tsv', sep = "/"), sep="\t")

## Build predictors based on 2019-2020 dataset and predict 2015-2017

output_files_path = "../output/predict_2015_2017/"
if (!dir.exists(output_files_path)) { dir.create(output_files_path) }

ix_2019_2020 = which(phys_chem['year',] %in% c(2019, 2020))
ix_2015_2017 = which(phys_chem['year',] %in% c(2015, 2016, 2017))

samples_2019_2020 = colnames(phys_chem)[ix_2019_2020]
samples_2015_2017 = colnames(phys_chem)[ix_2015_2017]

## phys_chem_prediction based on 16S

features_matrix_train = norm_asv_counts_16S[,samples_2019_2020]
responses_matrix_train = phys_chem[,samples_2019_2020]
features_matrix_target = norm_asv_counts_16S[,samples_2015_2017]
responses_matrix_target = phys_chem[,samples_2015_2017]

features_matrix_train = do_feature_selection(features_matrix_train, 0.1)
features_matrix_target = features_matrix_target[rownames(features_matrix_train),]

identical(colnames(features_matrix_train),colnames(responses_matrix_train))

predicted_responses = predict_randomforest(features_matrix_train, responses_matrix_train, features_matrix_target)

write.table(responses_matrix_target, paste(output_files_path, 'norm_asv_counts_16S_phys_chem_2015_2017_Actual.tsv', sep = "/"), sep="\t")
write.table(predicted_responses, paste(output_files_path, 'norm_asv_counts_16S_phys_chem_2015_2017_Predictions.tsv', sep = "/"), sep="\t")

## phyt_plan_genus_prediction based on 18S

output_files_path = "../output/predict_2015_2017/"
if (!dir.exists(output_files_path)) { dir.create(output_files_path) }

features_matrix_train = norm_asv_counts_18S[,samples_2019_2020]
responses_matrix_train = phyt_plan_genus[,which(colnames(phyt_plan_genus) %in% samples_2019_2020)]
features_matrix_target = norm_asv_counts_18S[,samples_2015_2017]
responses_matrix_target = phyt_plan_genus[,which(colnames(phyt_plan_genus) %in% samples_2015_2017)]

features_matrix_train = extract_shared_samples(features_matrix_train, responses_matrix_train)$features_matrix
responses_matrix_train = extract_shared_samples(features_matrix_train, responses_matrix_train)$responses_matrix

identical(colnames(features_matrix_train),colnames(responses_matrix_train))

features_matrix_target = extract_shared_samples(features_matrix_target, responses_matrix_target)$features_matrix
responses_matrix_target = extract_shared_samples(features_matrix_target, responses_matrix_target)$responses_matrix

identical(colnames(features_matrix_target),colnames(responses_matrix_target))

features_matrix_train = do_feature_selection(features_matrix_train, 0.1)
features_matrix_target = features_matrix_target[rownames(features_matrix_train),]

predicted_responses = predict_randomforest(features_matrix_train, responses_matrix_train, features_matrix_target)

write.table(responses_matrix_target, paste(output_files_path, 'norm_asv_counts_18S_phyt_plan_genus_2015_2017_Actual.tsv', sep = "/"), sep="\t")
write.table(predicted_responses, paste(output_files_path, 'norm_asv_counts_18S_phyt_plan_genus_2015_2017_Predictions.tsv', sep = "/"), sep="\t")

## zoo_plan_genus_prediction based on 16S

output_files_path = "../output/predict_2015_2017/"
if (!dir.exists(output_files_path)) { dir.create(output_files_path) }

features_matrix_train = norm_asv_counts_16S[,samples_2019_2020]
responses_matrix_train = zoo_plan_genus[,which(colnames(zoo_plan_genus) %in% samples_2019_2020)]
features_matrix_target = norm_asv_counts_16S[,samples_2015_2017]
responses_matrix_target = zoo_plan_genus[,which(colnames(zoo_plan_genus) %in% samples_2015_2017)]

features_matrix_train = extract_shared_samples(features_matrix_train, responses_matrix_train)$features_matrix
responses_matrix_train = extract_shared_samples(features_matrix_train, responses_matrix_train)$responses_matrix

identical(colnames(features_matrix_train),colnames(responses_matrix_train))

features_matrix_target = extract_shared_samples(features_matrix_target, responses_matrix_target)$features_matrix
responses_matrix_target = extract_shared_samples(features_matrix_target, responses_matrix_target)$responses_matrix

identical(colnames(features_matrix_target),colnames(responses_matrix_target))

features_matrix_train = do_feature_selection(features_matrix_train, 0.1)
features_matrix_target = features_matrix_target[rownames(features_matrix_train),]

predicted_responses = predict_randomforest(features_matrix_train, responses_matrix_train, features_matrix_target)

write.table(responses_matrix_target, paste(output_files_path, 'norm_asv_counts_16S_zoo_plan_genus_2015_2017_Actual.tsv', sep = "/"), sep="\t")
write.table(predicted_responses, paste(output_files_path, 'norm_asv_counts_16S_zoo_plan_genus_2015_2017_Predictions.tsv', sep = "/"), sep="\t")


## Cross-validation based predicitons on 2019/2020 dataset

ix_2019_2020 = which(phys_chem['year',] %in% c(2019, 2020))

samples_2019_2020 = colnames(phys_chem)[ix_2019_2020]

## phys_chem_prediction based on 16S

output_files_path = "../output/only_2019_2020/"
if (!dir.exists(output_files_path)) { dir.create(output_files_path) }

outfile_actual = "norm_asv_counts_16S_phys_chem_2019_2020_Actual.tsv"
outfile_predicted = "norm_asv_counts_16S_phys_chem_2019_2020_Predictions.tsv"
responses_matrix_full = phys_chem[,samples_2019_2020]
features_matrix_full = norm_asv_counts_16S[,samples_2019_2020]

features_matrix = extract_shared_samples(features_matrix_full, responses_matrix_full)$features_matrix
responses_matrix = extract_shared_samples(features_matrix_full, responses_matrix_full)$responses_matrix
features_matrix = do_feature_selection(features_matrix, 0.1)
predicted_responses_matrix = run_randomforest(features_matrix, responses_matrix, 10, 10)
write.table(responses_matrix, paste(output_files_path, outfile_actual, sep = "/"), sep="\t")
write.table(predicted_responses_matrix, paste(output_files_path, outfile_predicted, sep = "/"), sep="\t")

## phyt_plan_genus_prediction based on 18S

output_files_path = "../output/only_2019_2020/"
if (!dir.exists(output_files_path)) { dir.create(output_files_path) }

outfile_actual = "norm_asv_counts_18S_phyt_plan_genus_2019_2020_Actual.tsv"
outfile_predicted = "norm_asv_counts_18S_phyt_plan_genus_2019_2020_Predictions.tsv"
responses_matrix_full = phyt_plan_genus[,which(colnames(phyt_plan_genus) %in% samples_2019_2020)]
features_matrix_full = norm_asv_counts_18S[,samples_2019_2020]

features_matrix = extract_shared_samples(features_matrix_full, responses_matrix_full)$features_matrix
responses_matrix = extract_shared_samples(features_matrix_full, responses_matrix_full)$responses_matrix
features_matrix = do_feature_selection(features_matrix, 0.1)
predicted_responses_matrix = run_randomforest(features_matrix, responses_matrix, 10, 10)
write.table(responses_matrix, paste(output_files_path, outfile_actual, sep = "/"), sep="\t")
write.table(predicted_responses_matrix, paste(output_files_path, outfile_predicted, sep = "/"), sep="\t")

## zoo_plan_genus_prediction based on 16S

output_files_path = "../output/only_2019_2020/"
if (!dir.exists(output_files_path)) { dir.create(output_files_path) }

outfile_actual = "norm_asv_counts_16S_zoo_plan_genus_2019_2020_Actual.tsv"
outfile_predicted = "norm_asv_counts_16S_zoo_plan_genus_2019_2020_Predictions.tsv"
responses_matrix_full = zoo_plan_genus[,which(colnames(zoo_plan_genus) %in% samples_2019_2020)]
features_matrix_full = norm_asv_counts_16S[,samples_2019_2020]

features_matrix = extract_shared_samples(features_matrix_full, responses_matrix_full)$features_matrix
responses_matrix = extract_shared_samples(features_matrix_full, responses_matrix_full)$responses_matrix
features_matrix = do_feature_selection(features_matrix, 0.1)
predicted_responses_matrix = run_randomforest(features_matrix, responses_matrix, 10, 10)
write.table(responses_matrix, paste(output_files_path, outfile_actual, sep = "/"), sep="\t")
write.table(predicted_responses_matrix, paste(output_files_path, outfile_predicted, sep = "/"), sep="\t")

#####################
#### other stuff ####

summary(cor_matr)
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


par(mfrow = c(1,1), mar=c(4,4,1,1), xpd = TRUE) # good size of figure: 750 x 530
plot(cor_matr_physchem_16S[,4], cor_matr_physchem_18S[,4], xlim = c(0,1), ylim = c(0,1), xlab = "ML based on 16S", ylab = "ML based on 18S")
lines(c(0,1),c(0,1))
wilcox.test(cor_matr_physchem_16S[,4], cor_matr_physchem_18S[,4], paired = T)


### Do 'predictions' based on sequence matching

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

shared_genera = sort(intersect(rownames(features_matrix), rownames(responses_matrix)))
cor_matr_seq_match = matrix(ncol = 3, nrow = length(shared_genera))
colnames(cor_matr_seq_match) = c("pos_samp_micr","pos_samp_barc","cor-R")
rownames(cor_matr_seq_match) = shared_genera
for (i in 1:length(shared_genera)) {
  ixf = match(shared_genera[i], rownames(features_matrix))
  ixr = match(shared_genera[i], rownames(responses_matrix))
  cor_matr_seq_match[i,1] = length(which(responses_matrix[ixr,] > 0))
  cor_matr_seq_match[i,2] = length(which(features_matrix[ixf,] > 0))
  cor_matr_seq_match[i,3] = cor.test(responses_matrix[ixr,], features_matrix[ixf,])$est
}

plot(cor_matr_seq_match[,2], cor_matr_seq_match[,3], xlab = "Presence in #MB samples", ylab = "Prediction R-value", ylim = c(0,1))
plot(cor_matr_seq_match[,1], cor_matr_seq_match[,3], xlab = "Presence in #Mic samples", ylab = "Prediction R-value", ylim = c(0,1))
plot(cor_matr_seq_match[,1]*cor_matr_seq_match[,2], cor_matr_seq_match[,3], xlab = "Presence in #Mic samples", ylab = "Prediction R-value", ylim = c(0,1))

ix1 = intersect(which(cor_matr_seq_match[,1] >= 10), which(cor_matr_seq_match[,2] >= 10))
ix2 = match(shared_genera[ix1], rownames(cor_matr))
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


##







