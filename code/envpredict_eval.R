
library(pheatmap)
par(mgp = c(2.2, 1, 0))

##############################################################################
### Physche#m predictions from different taxonomic levels  - Random Forest ###
results_folder = "../output/DifferentTaxonomicLevels"
prediction_files = list.files(results_folder, pattern = "Predictions.tsv")
prediction_files = prediction_files[c(grep("_16S_", prediction_files), grep("_18S_", prediction_files))]
observation_files = gsub("Predictions", "Actual", prediction_files)
pred_matr = read.delim(paste(results_folder, prediction_files[1], sep = "/"), row.names = 1)
num_samples = c()
cor_matr = matrix(nrow = nrow(pred_matr), ncol = length(prediction_files))
colnames(cor_matr) = gsub("_Predictions.tsv","", prediction_files)
colnames(cor_matr) = gsub("_RF10fold","", colnames(cor_matr))
colnames(cor_matr) = gsub("norm_clade_counts_","", colnames(cor_matr))
colnames(cor_matr) = gsub("norm_seqtab_","", colnames(cor_matr))
colnames(cor_matr) = gsub(".tsv","", colnames(cor_matr))
rownames(cor_matr) = rownames(pred_matr)
for (i in 1:length(prediction_files)) {
  pred_matr = as.matrix(read.delim(paste(results_folder, prediction_files[i], sep = "/"), row.names = 1))
  obs_matr = as.matrix(read.delim(paste(results_folder, observation_files[i], sep = "/"), row.names = 1))
  for (j in 1:nrow(pred_matr)) {
    cor_matr[j,i] = round(cor.test(pred_matr[j,], obs_matr[j,])$est , 3)
    num_samples[j] = length(which(!is.na(obs_matr[j,])))
  }
}
sample_labels = paste(rownames(cor_matr), " [", num_samples, "]", sep = "") 
par(mar = c(2,10,2,2))
barplot(t(cor_matr), beside = T, horiz = T, las = 1, legend = T, names.arg = sample_labels, cex.names = 0.8, args.legend = list(x = "topright"))
par(mar = c(5,5,1,1))
barplot(colMeans(cor_matr), las = 2, ylab = "mean cor")
pheatmap(cor_matr, cluster_cols = F, cluster_rows = F, labels_row = sample_labels)
cor_matr_rf = cor_matr

#####################################################################
### ML Physchem predictions from different taxonomic levels - TabPNF ###
results_folder = "../output/TabPNF/DifferentTaxonomicLevels"
prediction_files = list.files(results_folder, pattern = "Predictions.tsv")
prediction_files = prediction_files[c(grep("_16S_", prediction_files), grep("_18S_", prediction_files))]
observation_files = gsub("Predictions", "Actual", prediction_files)
pred_matr = read.delim(paste(results_folder, prediction_files[1], sep = "/"), row.names = 1)
num_samples = c()
cor_matr = matrix(nrow = nrow(pred_matr), ncol = length(prediction_files))
colnames(cor_matr) = gsub("_Predictions.tsv","", prediction_files)
colnames(cor_matr) = gsub("_RF10fold","", colnames(cor_matr))
colnames(cor_matr) = gsub("norm_clade_counts_","", colnames(cor_matr))
colnames(cor_matr) = gsub("norm_seqtab_","", colnames(cor_matr))
colnames(cor_matr) = gsub(".tsv","", colnames(cor_matr))
rownames(cor_matr) = rownames(pred_matr)
for (i in 1:length(prediction_files)) {
  pred_matr = as.matrix(read.delim(paste(results_folder, prediction_files[i], sep = "/"), row.names = 1))
  #obs_matr = as.matrix(read.delim(paste(results_folder, observation_files[i], sep = "/"), row.names = 1))
  obs_matr = as.matrix(read.delim("../output/DifferentTaxonomicLevels/norm_clade_counts_16S_2_RF10fold_Actual.tsv", row.names = 1))
  for (j in 1:nrow(pred_matr)) {
    cor_matr[j,i] = round(cor.test(pred_matr[j,], obs_matr[j,])$est , 3)
    num_samples[j] = length(which(!is.na(obs_matr[j,])))
  }
}
sample_labels = paste(rownames(cor_matr), " [", num_samples, "]", sep = "") 
par(mar = c(2,10,2,2))
barplot(t(cor_matr), beside = T, horiz = T, las = 1, legend = T, names.arg = sample_labels, cex.names = 0.8, args.legend = list(x = "topright"))
par(mar = c(5,5,1,1))
barplot(colMeans(cor_matr), las = 2, ylab = "mean cor")
pheatmap(cor_matr, cluster_cols = F, cluster_rows = F, labels_row = sample_labels)
cor_matr_tpnf = cor_matr

# plotting results from RF vs TablePNF
# par(mfrow=c(2,2), mar = c(4,4,1,1))
# plot(cor_matr_rf, cor_matr_tpnf, xlim = c(-0.3, 1), ylim = c(-0.3, 1))
# lines(c(-0.2,1), c(-0.2,1))
# wilcox.test(cor_matr_rf, cor_matr_tpnf, paired = T)
# 
# plot(cor_matr_rf, cor_matr_tpnf, xlim = c(0.9, 1), ylim = c(0.9, 1))
# lines(c(0.8,1), c(0.8,1))


##############################################################
### ML Physchem predictions from RepresentationsFromDeepMicro ###
results_folder = "../output/RepresentationsFromDeepMicro"
#results_folder = "~/aquatic/envpredict/output/RepresentationsFromDeepMicro_prefilt_features"
prediction_files = list.files(results_folder, pattern = "Predictions.tsv")
observation_files = gsub("Predictions", "Actual", prediction_files)
pred_matr = read.delim(paste(results_folder, prediction_files[1], sep = "/"), row.names = 1)
num_samples = c()
cor_matr = matrix(nrow = nrow(pred_matr), ncol = length(prediction_files))
colnames(cor_matr) = gsub("_Predictions.tsv","", prediction_files)
colnames(cor_matr) = gsub("_RF10fold","", colnames(cor_matr))
colnames(cor_matr) = gsub("norm_clade_counts_","", colnames(cor_matr))
colnames(cor_matr) = gsub("norm_seqtab_","", colnames(cor_matr))
colnames(cor_matr) = gsub("_filtered","", colnames(cor_matr))
colnames(cor_matr) = gsub("_filt","", colnames(cor_matr))
colnames(cor_matr) = gsub("_pivoted_march_2025_rep","", colnames(cor_matr))
colnames(cor_matr) = gsub(".tsv","", colnames(cor_matr))
rownames(cor_matr) = rownames(pred_matr)
for (i in 1:length(prediction_files)) {
  pred_matr = as.matrix(read.delim(paste(results_folder, prediction_files[i], sep = "/"), row.names = 1))
  obs_matr = as.matrix(read.delim(paste(results_folder, observation_files[i], sep = "/"), row.names = 1))
  for (j in 1:nrow(pred_matr)) {
    cor_matr[j,i] = round(cor.test(pred_matr[j,], obs_matr[j,])$est , 3)
    num_samples[j] = length(which(!is.na(obs_matr[j,])))
  }
}
sample_labels = paste(rownames(cor_matr), " [", num_samples, "]", sep = "") 
par(mar = c(2,10,2,2))
barplot(t(cor_matr), beside = T, horiz = T, las = 1, legend = T, names.arg = sample_labels, cex.names = 0.8, args.legend = list(x = "topright"))
par(mar = c(5,5,1,1))
barplot(colMeans(cor_matr), las = 2, ylab = "mean cor")
pheatmap(cor_matr, cluster_cols = F, cluster_rows = F, labels_row = sample_labels)
cor_matr_de = cor_matr
#cor_matr_depf = cor_matr

########################################################################
### ML Physchem predictions from either 16S-ASVs or plankton microscopy ###
plot_folder = '../output/plots_ML Physchem_preds/'
if (!dir.exists(plot_folder)){
  dir.create(plot_folder)
}

group_colors <- c('#018571', '#a6611a', '#ffffb2')
names(group_colors) <- c('Phytoplankton', 'Zooplankton', '16S')

pdf(paste0(plot_folder, "16S_vs_zooplankton_vs_phytoplankton.pdf"),
    width = 8*0.8, height = 10*0.8)

layout(matrix(c(1,1,1,1,4,4,4,4,2,3,5,6), 4, 3, byrow = F))
results_folders = c() 
results_folders[1] = "../output/ML Physchem_based_on_seqdata_or_phytoplan"
results_folders[2] = "../output/ML Physchem_based_on_seqdata_or_zooplan"
for (res_fol in 1:2) {
  results_folder = results_folders[res_fol]
  prediction_files = list.files(results_folder, pattern = "Predictions.tsv")
  observation_files = gsub("Predictions", "Actual", prediction_files)
  pred_matr = read.delim(paste(results_folder, prediction_files[1], sep = "/"), row.names = 1)
  cor_matr = matrix(nrow = nrow(pred_matr), ncol = length(prediction_files))
  colnames(cor_matr) = gsub("_Predictions.tsv","", prediction_files)
  colnames(cor_matr) = gsub("-based_ML Physchem","", colnames(cor_matr))
  rownames(cor_matr) = rownames(pred_matr)
  num_samples = c()
  for (i in 1:length(prediction_files)) {
    pred_matr = as.matrix(read.delim(paste(results_folder, prediction_files[i], sep = "/"), row.names = 1))
    obs_matr = as.matrix(read.delim(paste(results_folder, observation_files[i], sep = "/"), row.names = 1))
    for (j in 1:nrow(pred_matr)) {
      cor_matr[j,i] = cor.test(pred_matr[j,], obs_matr[j,])$est
      #cor_matr[j,i] = (cor_matr[j,i])^2 # r-squared
      cor_matr[j,i] = round(cor_matr[j,i], 3)
      num_samples[j] = length(which(!is.na(obs_matr[j,])))
    }
  }
  if(res_fol == 1){
    cor_matr = cor_matr[,c("16S","phytoplan_genus","phytoplan")]
    colnames(cor_matr) = c("16S", "Phytoplankton_genus", "Phytoplankton")
  }else{
    cor_matr = cor_matr[,c("16S","zooplan_genus","zooplan")]
    colnames(cor_matr) = c("16S", "Zooplankton_genus", "Zooplankton")
  }
  
  cor_matr = cor_matr[-c(1:7),]
  num_samples = num_samples[-c(1:7)]
  sample_labels = paste(rownames(cor_matr), " [", num_samples, "]", sep = "") 
  #pheatmap(cor_matr, cluster_cols = F, cluster_rows = F, labels_row = sample_labels)
  #layout(matrix(c(1,2,1,3), 2, 2, byrow = T))
  par(mar = c(3,7,4,2), xpd = T)
  #barplot(t(cor_matr[,]), beside = T, horiz = T, las = 1, legend = T, names.arg = sample_labels, cex.names = 0.8, args.legend = list(x = "top"))
  barplot(t(cor_matr[,c(1,3)]), beside = T, horiz = T, las = 1, names.arg = sample_labels, cex.names = 0.8, xlab = "Pearson cor. (R)", col = rev(group_colors[c(res_fol,3)]), xlim = c(0,1))
  
  legend("top", inset = c(0, -0.05),
         legend = colnames(cor_matr)[c(1,3)],
         fill = rev(group_colors[c(res_fol,3)]),
         title = "Model trained on:")
  
  par(mar = c(2,4,2,2))
  
  #boxplot(cor_matr[,1], cor_matr[,2], cor_matr[,3],  names = colnames(cor_matr), las = 2, ylim = c(0.2, 1))
  boxplot(cor_matr[,1], cor_matr[,3],  names = colnames(cor_matr)[c(1,3)], ylim = c(0.2, 1), col = rev(group_colors[c(res_fol,3)]),
          angle = 45)
  par(mar = c(4,4,2,2))
  plot(cor_matr[,1], cor_matr[,3], xlim = c(0,1), ylim = c(0,1), xlab = colnames(cor_matr)[1], ylab = colnames(cor_matr)[3], col  = group_colors[res_fol], bg = group_colors[3], alpha = 0.5, pch = 21)
  lines(c(0,1), c(0,1))
  print(colnames(cor_matr)[3])
  print(wilcox.test(cor_matr[,1], cor_matr[,2], paired = T))
  print(wilcox.test(cor_matr[,1], cor_matr[,3], paired = T))
}

dev.off()

#############################################################################################################
### Phytoplankton predictions from either 16S-ASVs (ML) 18S-ASVs (ML), 18S-seq-matching, or ML Physchem data ###
results_folder = "../output/phytoplankton_predicted"
prediction_files = list.files(results_folder, pattern = "Predictions.tsv")
prediction_files = prediction_files[c(grep("phyt_plan_genus", prediction_files))] # grep("direct_matching_", prediction_files))]
observation_files = gsub("Predictions", "Actual", prediction_files)
pred_matr = read.delim(paste(results_folder, prediction_files[1], sep = "/"), row.names = 1)
shared_genus = rownames(pred_matr)
for (i in 1:length(prediction_files)) {
  pred_matr = as.matrix(read.delim(paste(results_folder, prediction_files[i], sep = "/"), row.names = 1))
  shared_genus = intersect(shared_genus, rownames(pred_matr))
}
cor_matr = matrix(nrow = length(shared_genus), ncol = length(prediction_files))
colnames(cor_matr) = gsub("_Predictions.tsv","", prediction_files)
colnames(cor_matr) = gsub("phyt_plan_genus_","", colnames(cor_matr))
colnames(cor_matr) = gsub("direct_matching_norm_clade_counts_18S_8_phyt_plan_genus","match_18S", colnames(cor_matr))
colnames(cor_matr) = gsub("norm_asv_counts_","ML_", colnames(cor_matr))
colnames(cor_matr) = gsub("_RF10fold","", colnames(cor_matr))
colnames(cor_matr) = gsub("renomralized_matching_18","renorm_matching_18", colnames(cor_matr))
rownames(cor_matr) = shared_genus
length(shared_genus)
num_samples = c()
for (i in 1:length(prediction_files)) {
  pred_matr = as.matrix(read.delim(paste(results_folder, prediction_files[i], sep = "/"), row.names = 1))[shared_genus,]
  obs_matr = as.matrix(read.delim(paste(results_folder, observation_files[i], sep = "/"), row.names = 1))[shared_genus,]
  for (j in 1:nrow(pred_matr)) {
    num_samples[j] = NA 
    if (length(unique(obs_matr[j,])) > 0) { # if not all observations are the same
      cor_matr[j,i] = round(cor.test(pred_matr[j,], obs_matr[j,])$est , 3)
      #num_samples[j] = length(which(!is.na(obs_matr[j,])))
      num_samples[j] = length(which(obs_matr[j,] != 0))
    }
  }
}
ix = which(rowSums(is.na(cor_matr)) == 0)
cor_matr = cor_matr[ix,]
num_samples = num_samples[ix]
cor_matr_phytoplan = cor_matr
num_samples_phytoplan = num_samples
sample_labels = paste(rownames(cor_matr), " [", num_samples, "]", sep = "") 
par(mar = c(2,10,2,2))
barplot(t(cor_matr), beside = T, horiz = T, las = 1, legend = T, names.arg = sample_labels, cex.names = 0.8, args.legend = list(x = "topright"))
pheatmap(cor_matr, cluster_cols = F, cluster_rows = F, labels_row = sample_labels)
# plotting cor vs number of micr samples where genus was observed

group_colors <- c('#99d8c9','#66c2a4','#2ca25f','#006d2c')
names(group_colors) <- c('ML Physchem', 'ML 16S', 'ML 18S', 'Match 18S')
colnames(cor_matr) = c('Match 18S', 'ML 16S', 'ML 18S', 'ML Physchem', "Renormalized match 18s")

pdf(paste0(plot_folder, "predicting_phytoplankton.pdf"),
    width = 9*0.8, height = 6*0.8)

layout(matrix(c(1,2,5,3,4,6), 2, 3, byrow = T))
ix2 = c(4,2,3,1)
par(mar = c(4,4,1,1))
for (i in 1:length(ix2)) {
  plot(num_samples, cor_matr[,ix2[i]], main = colnames(cor_matr)[ix2[i]], ylim = c(0,1), xlab = "Detected in nr samples", ylab = "Pearson correlation coefficient (r)", bg = group_colors[i], col = group_colors[i], pch = 21)
}
ix = which(num_samples > 50)
par(mar = c(4,4,1,1))
boxplot(
  cor_matr[ix,ix2[1]], cor_matr[ix,ix2[2]], cor_matr[ix,ix2[3]], cor_matr[ix,ix2[4]], 
  names = colnames(cor_matr)[ix2], las = 3, ylab = "Pearson correlation coefficient (r)", col = group_colors 
)
par(mar = c(4,4,1,1))
plot(cor_matr[ix,ix2[4]], cor_matr[ix,ix2[3]], xlim = c(0,1), ylim = c(0,1), xlab = colnames(cor_matr)[ix2[4]], ylab = colnames(cor_matr)[ix2[3]], col = group_colors[4], bg = group_colors[3], pch = 21)
lines(c(0,1), c(0,1))

dev.off()

ix = which(num_samples > 50)
wilcox.test(cor_matr[ix,ix2[1]], cor_matr[ix,ix2[2]], paired = T)
wilcox.test(cor_matr[ix,ix2[1]], cor_matr[ix,ix2[3]], paired = T)

wilcox.test(cor_matr[ix,'ML Physchem'], cor_matr[ix,'ML 16S'], paired = T)
wilcox.test(cor_matr[ix,'ML Physchem'], cor_matr[ix,'ML 18S'], paired = T)
wilcox.test(cor_matr[ix,'Match 18S'], cor_matr[ix,'ML 18S'], paired = T)
wilcox.test(cor_matr[ix,'Match 18S'], cor_matr[ix,'ML 16S'], paired = T)
wilcox.test(cor_matr[ix,'Match 18S'], cor_matr[ix,'ML Physchem'], paired = T)
wilcox.test(cor_matr[ix,'ML 18S'], cor_matr[ix,'ML 16S'], paired = T)

## Zooplankton predictions from either 16S-ASVs (ML) 18S-ASVs (ML), 18S-seq-matching, or physiochem data
results_folder = "../output/zooplankton_predicted"
prediction_files = list.files(results_folder, pattern = "Predictions.tsv")
prediction_files = prediction_files[c(grep("zoo_plan_genus", prediction_files))] # grep("direct_matching_", prediction_files))]
observation_files = gsub("Predictions", "Actual", prediction_files)
pred_matr = read.delim(paste(results_folder, prediction_files[1], sep = "/"), row.names = 1)
shared_genus = rownames(pred_matr)
for (i in 1:length(prediction_files)) {
  pred_matr = as.matrix(read.delim(paste(results_folder, prediction_files[i], sep = "/"), row.names = 1))
  shared_genus = intersect(shared_genus, rownames(pred_matr))
}
cor_matr = matrix(nrow = length(shared_genus), ncol = length(prediction_files))
colnames(cor_matr) = gsub("_Predictions.tsv","", prediction_files)
colnames(cor_matr) = gsub("zoo_plan_genus_","", colnames(cor_matr))
colnames(cor_matr) = gsub("direct_matching_norm_clade_counts_18S_8_zoo_plan_genus","match_18S", colnames(cor_matr))
colnames(cor_matr) = gsub("norm_asv_counts_","ML_", colnames(cor_matr))
colnames(cor_matr) = gsub("_RF10fold","", colnames(cor_matr))
colnames(cor_matr) = gsub("renomralized_matching_18","renorm_matching_18", colnames(cor_matr))
rownames(cor_matr) = shared_genus
length(shared_genus)
num_samples = c()
for (i in 1:length(prediction_files)) {
  pred_matr = as.matrix(read.delim(paste(results_folder, prediction_files[i], sep = "/"), row.names = 1))[shared_genus,]
  obs_matr = as.matrix(read.delim(paste(results_folder, observation_files[i], sep = "/"), row.names = 1))[shared_genus,]
  for (j in 1:nrow(pred_matr)) {
    num_samples[j] = NA 
    if (length(unique(obs_matr[j,])) > 0) { # if not all observations are the same
      cor_matr[j,i] = round(cor.test(pred_matr[j,], obs_matr[j,])$est , 3)
      #num_samples[j] = length(which(!is.na(obs_matr[j,])))
      num_samples[j] = length(which(obs_matr[j,] != 0))
    }
  }
}
ix = which(rowSums(is.na(cor_matr)) == 0)
cor_matr = cor_matr[ix,]
num_samples = num_samples[ix]
cor_matr_zooplan = cor_matr
num_samples_zooplan = num_samples
sample_labels = paste(rownames(cor_matr), " [", num_samples, "]", sep = "") 
par(mfrow = c(1,1), mar = c(2,10,2,2))
barplot(t(cor_matr), beside = T, horiz = T, las = 1, legend = T, names.arg = sample_labels, cex.names = 0.8, args.legend = list(x = "topright"))
pheatmap(cor_matr, cluster_cols = F, cluster_rows = F, labels_row = sample_labels)
layout(matrix(c(1,2,5,3,4,6), 2, 3, byrow = T))
ix2 = c(4,2,3,1)
par(mar = c(4,4,1,1))
for (i in 1:length(ix2)) {
  plot(num_samples, cor_matr[,ix2[i]], main = colnames(cor_matr)[ix2[i]], ylim = c(0,1), xlab = "Detected in num samples", ylab = "Pearson cor (R)")
}
ix = which(num_samples > 50)
par(mar = c(6,4,1,1))
boxplot(
  cor_matr[ix,ix2[1]], cor_matr[ix,ix2[2]], cor_matr[ix,ix2[3]], cor_matr[ix,ix2[4]], 
  names = colnames(cor_matr)[ix2], las = 3, ylab = "Pearson cor (R)"
)
par(mar = c(4,4,1,1))
plot(cor_matr[ix,ix2[4]], cor_matr[ix,ix2[3]], xlim = c(0,1), ylim = c(0,1), xlab = colnames(cor_matr)[ix2[4]], ylab = colnames(cor_matr)[ix2[3]])
lines(c(0,1), c(0,1))
ix = which(num_samples > 50)
wilcox.test(cor_matr[ix,ix2[1]], cor_matr[ix,ix2[2]], paired = T)
wilcox.test(cor_matr[ix,ix2[1]], cor_matr[ix,ix2[3]], paired = T)


group_colors <- c('#fec44f','#fe9929','#d95f0e','#993404')
names(group_colors) <- c('ML Physchem', 'ML 16S', 'ML 18S', 'Match 18S')
colnames(cor_matr) = c('Match 18S', 'ML 16S', 'ML 18S', 'ML Physchem', "Renormalized match 18s")

pdf(paste0(plot_folder, "predicting_zooplankton.pdf"),
    width = 9*0.8, height = 6*0.8)

layout(matrix(c(1,2,5,3,4,6), 2, 3, byrow = T))
ix2 = c(4,2,3,1)
par(mar = c(4,4,1,1))
for (i in 1:length(ix2)) {
  plot(num_samples, cor_matr[,ix2[i]], main = colnames(cor_matr)[ix2[i]], ylim = c(0,1), xlab = "Detected in nr samples", ylab = "Pearson correlation coefficient (r)", bg = group_colors[i], col = group_colors[i], pch = 21)
}
ix = which(num_samples > 50)
par(mar = c(4,4,1,1))
boxplot(
  cor_matr[ix,ix2[1]], cor_matr[ix,ix2[2]], cor_matr[ix,ix2[3]], cor_matr[ix,ix2[4]], 
  names = colnames(cor_matr)[ix2], las = 3, ylab = "Pearson correlation coefficient (r)", col = group_colors 
)
par(mar = c(4,4,1,1))
plot(cor_matr[ix,ix2[4]], cor_matr[ix,ix2[3]], xlim = c(0,1), ylim = c(0,1), xlab = colnames(cor_matr)[ix2[4]], ylab = colnames(cor_matr)[ix2[3]], col = group_colors[4], bg = group_colors[3], pch = 21)
lines(c(0,1), c(0,1))

dev.off()

ix = which(num_samples > 50)
wilcox.test(cor_matr[ix,ix2[1]], cor_matr[ix,ix2[2]], paired = T)
wilcox.test(cor_matr[ix,ix2[1]], cor_matr[ix,ix2[3]], paired = T)

wilcox.test(cor_matr[ix,'ML Physchem'], cor_matr[ix,'ML 16S'], paired = T)
wilcox.test(cor_matr[ix,'ML Physchem'], cor_matr[ix,'ML 18S'], paired = T)
wilcox.test(cor_matr[ix,'Match 18S'], cor_matr[ix,'ML 18S'], paired = T)
wilcox.test(cor_matr[ix,'Match 18S'], cor_matr[ix,'ML 16S'], paired = T)
wilcox.test(cor_matr[ix,'Match 18S'], cor_matr[ix,'ML Physchem'], paired = T)
wilcox.test(cor_matr[ix,'ML 18S'], cor_matr[ix,'ML 16S'], paired = T)
# plotting the phyto and zooplankton predictions in one figure (requires the above 2 sections to be run first)
layout(matrix(c(1,2,5,3,4,6,7,8,11,9,10,12), 4, 3, byrow = T))
ix2 = c(4,2,3,1)
num_samples = num_samples_phytoplan
cor_matr = cor_matr_phytoplan
par(mar = c(5,4,1,1), mgp = c(2,1,0))
for (i in 1:length(ix2)) {
  plot(num_samples, cor_matr[,ix2[i]], main = colnames(cor_matr)[ix2[i]], ylim = c(0,1), xlab = "Detected in num samples", ylab = "Pearson cor (R)")
}
ix = which(num_samples > 50)
boxplot(
  cor_matr[ix,ix2[1]], cor_matr[ix,ix2[2]], cor_matr[ix,ix2[3]], cor_matr[ix,ix2[4]], 
  horizontal = F, names = colnames(cor_matr)[ix2], las = 3, ylab = "Pearson cor (R)"
)
par(mar = c(4,4,1,1))
plot(cor_matr[ix,ix2[4]], cor_matr[ix,ix2[3]], xlim = c(0,1), ylim = c(0,1), xlab = colnames(cor_matr)[ix2[4]], ylab = colnames(cor_matr)[ix2[3]])
lines(c(0,1), c(0,1))
wilcox.test(cor_matr[ix,'phys_chem'], cor_matr[ix,'ML_16S'], paired = T)
wilcox.test(cor_matr[ix,'phys_chem'], cor_matr[ix,'ML_18S'], paired = T)
wilcox.test(cor_matr[ix,'match_18S'], cor_matr[ix,'ML_18S'], paired = T)

num_samples = num_samples_zooplan
cor_matr = cor_matr_zooplan
par(mar = c(5,4,1,1), mgp = c(2,1,0))
for (i in 1:length(ix2)) {
  plot(num_samples, cor_matr[,ix2[i]], main = colnames(cor_matr)[ix2[i]], ylim = c(0,1), xlab = "Detected in num samples", ylab = "Pearson cor (R)")
}
ix = which(num_samples > 50)
boxplot(
  cor_matr[ix,ix2[1]], cor_matr[ix,ix2[2]], cor_matr[ix,ix2[3]], cor_matr[ix,ix2[4]], 
  horizontal = F, names = colnames(cor_matr)[ix2], las = 3, ylab = "Pearson cor (R)"
)
par(mar = c(4,4,1,1))
plot(cor_matr[ix,ix2[4]], cor_matr[ix,ix2[3]], xlim = c(0,1), ylim = c(0,1), xlab = colnames(cor_matr)[ix2[4]], ylab = colnames(cor_matr)[ix2[3]])
lines(c(0,1), c(0,1))

############################################################################
### Phytoplankton-genus and 18S-genus predictions based on ML Physchem data ###
results_folder = "~/aquatic/envpredict/output/phytoplankton_genus_micr_and_18S_based_on_ML Physchem"
prediction_files = list.files(results_folder, pattern = "Predictions.tsv")
observation_files = gsub("Predictions", "Actual", prediction_files)
pred_matr = read.delim(paste(results_folder, prediction_files[1], sep = "/"), row.names = 1)
shared_genus = rownames(pred_matr)
for (i in 1:length(prediction_files)) {
  pred_matr = as.matrix(read.delim(paste(results_folder, prediction_files[i], sep = "/"), row.names = 1))
  shared_genus = intersect(shared_genus, rownames(pred_matr))
}
length(shared_genus) # 153
cor_matr = matrix(nrow = length(shared_genus), ncol = 2*length(prediction_files))
colnames(cor_matr) = c(paste(prediction_files[1], "cor"), paste(prediction_files[1], "non-zero"), paste(prediction_files[2], "cor"), paste(prediction_files[2], "non-zero"))
colnames(cor_matr) = gsub("_Predictions.tsv","", colnames(cor_matr))
colnames(cor_matr) = gsub("ML Physchem-based_","", colnames(cor_matr))
rownames(cor_matr) = shared_genus

for (i in 1:length(prediction_files)) {
  pred_matr = as.matrix(read.delim(paste(results_folder, prediction_files[i], sep = "/"), row.names = 1))[shared_genus,]
  obs_matr = as.matrix(read.delim(paste(results_folder, observation_files[i], sep = "/"), row.names = 1))[shared_genus,]
  for (j in 1:nrow(pred_matr)) {
    if (length(unique(obs_matr[j,])) > 0) { # if not all observations are the same
      #cor_matr[j,(2*i) - 1] = round(cor.test(pred_matr[j,], obs_matr[j,])$est , 3)
      cor_matr[j,(2*i) - 1] = round(cor.test(pred_matr[j,], obs_matr[j,], method = "spearman")$est , 3)
      #cor_matr[j,2*i] = length(which(!is.na(obs_matr[j,])))
      cor_matr[j,2*i] = length(which(obs_matr[j,] != 0))
    }
  }
}
ix = which(rowSums(is.na(cor_matr)) == 0)
cor_matr = cor_matr[ix,]

par(mfcol = c(1,1), mar = c(4,4,1,1))
plot(cor_matr[,2], cor_matr[,1], xlab = "Detected in number samples", ylab = "Cor")
points(cor_matr[,4], cor_matr[,3], col = "blue")

ix = intersect(which(cor_matr[,2] > 10), which(cor_matr[,4] > 10))
plot(cor_matr[ix,1], cor_matr[ix,3], xlab = "18S", ylab = "Micr", xlim = c(0,1), ylim = c(0,1))
lines(c(0,1), c(0,1))
ix = intersect(which(cor_matr[,2] > 50), which(cor_matr[,4] > 50))
plot(cor_matr[ix,1], cor_matr[ix,3], xlab = "18S", ylab = "Micr", xlim = c(0,1), ylim = c(0,1))
lines(c(0,1), c(0,1))


pred_matr_18S = as.matrix(read.delim(paste(results_folder, prediction_files[1], sep = "/"), row.names = 1))[shared_genus,]
obs_matr_18S = as.matrix(read.delim(paste(results_folder, observation_files[1], sep = "/"), row.names = 1))[shared_genus,]
pred_matr_mic = as.matrix(read.delim(paste(results_folder, prediction_files[2], sep = "/"), row.names = 1))[shared_genus,]
obs_matr_mic = as.matrix(read.delim(paste(results_folder, observation_files[2], sep = "/"), row.names = 1))[shared_genus,]


plot(pred_matr_18S["Alexandrium",], obs_matr_18S["Alexandrium",])
plot(pred_matr_mic["Alexandrium",], obs_matr_mic["Alexandrium",])

cor.test(pred_matr_18S["Alexandrium",], obs_matr_18S["Alexandrium",])
cor.test(pred_matr_mic["Alexandrium",], obs_matr_mic["Alexandrium",])









