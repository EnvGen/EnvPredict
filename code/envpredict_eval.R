
library(pheatmap)

# physchem predictions from different taxonomic levels  - Random Forest
results_folder = "~/aquatic/envpredict/output/DifferentTaxonomicLevels"
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
pheatmap(cor_matr, cluster_cols = F, cluster_rows = F, labels_row = sample_labels)
cor_matr_rf = cor_matr


# physchem predictions from different taxonomic levels  - TabPNF
results_folder = "~/aquatic/envpredict/output/DifferentTaxonomicLevels_TabPNF"
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
  obs_matr = as.matrix(read.delim("~/aquatic/envpredict/output/DifferentTaxonomicLevels/norm_clade_counts_16S_2_RF10fold_Actual.tsv", row.names = 1))
  for (j in 1:nrow(pred_matr)) {
    cor_matr[j,i] = round(cor.test(pred_matr[j,], obs_matr[j,])$est , 3)
    num_samples[j] = length(which(!is.na(obs_matr[j,])))
  }
}
sample_labels = paste(rownames(cor_matr), " [", num_samples, "]", sep = "") 
par(mar = c(2,10,2,2))
barplot(t(cor_matr), beside = T, horiz = T, las = 1, legend = T, names.arg = sample_labels, cex.names = 0.8, args.legend = list(x = "topright"))
pheatmap(cor_matr, cluster_cols = F, cluster_rows = F, labels_row = sample_labels)
cor_matr_tpnf = cor_matr


# plotting results from RF vs TablePNF
par(mfrow=c(2,2), mar = c(4,4,1,1))
plot(cor_matr_rf, cor_matr_tpnf, xlim = c(-0.3, 1), ylim = c(-0.3, 1))
lines(c(-0.2,1), c(-0.2,1))
wilcox.test(cor_matr_rf, cor_matr_tpnf, paired = T)

plot(cor_matr_rf, cor_matr_tpnf, xlim = c(0.9, 1), ylim = c(0.9, 1))
lines(c(0.8,1), c(0.8,1))


# physchem predictions from RepresentationsFromDeepMicro
results_folder = "~/aquatic/envpredict/output/RepresentationsFromDeepMicro"
prediction_files = list.files(results_folder, pattern = "Predictions.tsv")
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
pheatmap(cor_matr, cluster_cols = F, cluster_rows = F, labels_row = sample_labels)


# physchem predictions from either 16S-ASVs or plankton microscopy
results_folder = "~/aquatic/envpredict/output/physchem_based_on_seqdata_or_phytoplan"
prediction_files = list.files(results_folder, pattern = "Predictions.tsv")
observation_files = gsub("Predictions", "Actual", prediction_files)
pred_matr = read.delim(paste(results_folder, prediction_files[1], sep = "/"), row.names = 1)
cor_matr = matrix(nrow = nrow(pred_matr), ncol = length(prediction_files))
colnames(cor_matr) = gsub("_Predictions.tsv","", prediction_files)
colnames(cor_matr) = gsub("-based_physchem","", colnames(cor_matr))
rownames(cor_matr) = rownames(pred_matr)
num_samples = c()
for (i in 1:length(prediction_files)) {
  pred_matr = as.matrix(read.delim(paste(results_folder, prediction_files[i], sep = "/"), row.names = 1))
  obs_matr = as.matrix(read.delim(paste(results_folder, observation_files[i], sep = "/"), row.names = 1))
  for (j in 1:nrow(pred_matr)) {
    cor_matr[j,i] = round(cor.test(pred_matr[j,], obs_matr[j,])$est , 3)
    num_samples[j] = length(which(!is.na(obs_matr[j,])))
  }
}
cor_matr
sample_labels = paste(rownames(cor_matr), " [", num_samples, "]", sep = "") 
par(mfcol = c(1,1), mar = c(2,10,2,20), xpd = T)
barplot(t(cor_matr), beside = T, horiz = T, las = 1, legend = T, names.arg = sample_labels, cex.names = 0.8, args.legend = list(x = "bottom", inset = c(1,-0.25)))
pheatmap(cor_matr, cluster_cols = F, cluster_rows = F, labels_row = sample_labels)
boxplot(cor_matr[,1], cor_matr[,2], cor_matr[,3],  names = colnames(cor_matr))
par(mfcol = c(1,1), mar = c(4,4,2,2))
plot(cor_matr[,1], cor_matr[,3], xlim = c(0,1), ylim = c(0,1), xlab = colnames(cor_matr)[1], ylab = colnames(cor_matr)[3])
lines(c(0,1), c(0,1))



# plankton predictions from either 16S-ASVs (ML) 18S-ASVs (ML), 18S-seq-matching, or physiochem data
results_folder = "~/aquatic/envpredict/output/phytoplankton_predicted"
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
colnames(cor_matr) = gsub("renomralized_direct_matching_norm_asv_counts_18S","rn_direct_matching_18S", colnames(cor_matr))
colnames(cor_matr) = gsub("direct_matching_norm_asv_counts_18S","direct_matching_18S", colnames(cor_matr))
colnames(cor_matr) = gsub("norm_asv_counts_","ML_", colnames(cor_matr))
colnames(cor_matr) = gsub("_complete","", colnames(cor_matr))
rownames(cor_matr) = shared_genus
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
sample_labels = paste(rownames(cor_matr), " [", num_samples, "]", sep = "") 

par(mfrow = c(1,1), mar = c(2,10,2,2))
barplot(t(cor_matr), beside = T, horiz = T, las = 1, legend = T, names.arg = sample_labels, cex.names = 0.8, args.legend = list(x = "topright"))
pheatmap(cor_matr, cluster_cols = F, cluster_rows = F, labels_row = sample_labels)
# plotting cor vs number of micr samples where genus was observed
par(mfrow = c(3,2), mar = c(4,4,1,1))
for (i in 1:length(prediction_files)) {
  plot(num_samples, cor_matr[,i], main = colnames(cor_matr)[i])
}
# boxplots of cor per method, including genera observed in at least of 50 samples
ix = which(num_samples > 50)
ix2 = c(4,2,3,1,5)
par(mfrow = c(1,1), mar = c(10,4,1,1))
boxplot(
  cor_matr[ix,ix2[1]], cor_matr[ix,ix2[2]], cor_matr[ix,ix2[3]], cor_matr[ix,ix2[4]], cor_matr[ix,ix2[5]], 
  names = colnames(cor_matr)[ix2], las = 3, ylab = "Pearson cor (R)"
)
wilcox.test(cor_matr[,"phys_chem"], cor_matr[,"ML_16S"], paired = T)
wilcox.test(cor_matr[,"phys_chem"], cor_matr[,"ML_18S"], paired = T)
plot(cor_matr[ix,ix2[1]], cor_matr[ix,ix2[3]])
lines(c(0.3,1), c(0.3,1))







prediction_files = list.files(folder, pattern = "Predictions.tsv")
observation_files = gsub("Predictions", "Actual", prediction_files)
pred_matr = read.delim(paste(folder, prediction_files[1], sep = "/"), row.names = 1)
cor_matr = matrix(nrow = nrow(pred_matr), ncol = length(prediction_files))
colnames(cor_matr) = gsub("_Predictions.tsv","", prediction_files)
colnames(cor_matr) = gsub("_RF10fold","", colnames(cor_matr))
colnames(cor_matr) = gsub("norm_clade_counts_","", colnames(cor_matr))
colnames(cor_matr) = gsub("norm_seqtab_","", colnames(cor_matr))
colnames(cor_matr) = gsub(".tsv","", colnames(cor_matr))
rownames(cor_matr) = rownames(pred_matr)
for (i in 1:length(prediction_files)) {
  pred_matr = as.matrix(read.delim(paste(folder, prediction_files[i], sep = "/"), row.names = 1))
  obs_matr = as.matrix(read.delim(paste(folder, observation_files[i], sep = "/"), row.names = 1))
  for (j in 1:nrow(pred_matr)) {
    cor_matr[j,i] = round(cor.test(pred_matr[j,], obs_matr[j,])$est , 3)
  }
}
cor_matr
cor_matr_taxa = cor_matr


barplot(colMeans(cor_matr_deep), ylim = c(0,1), las = 3)
sort(colMeans(cor_matr_deep), decreasing = T)[1:10]

barplot(colMeans(cor_matr_taxa), ylim = c(0,1), las = 3)
sort(colMeans(cor_matr_taxa), decreasing = T)[1:10]


ix_16S = grep("16S", colnames(cor_matr))
ix_18S = grep("18S", colnames(cor_matr))

barplot(colMeans(cor_matr)[ix_16S], ylim = c(0,1))
barplot(colMeans(cor_matr)[ix_18S], ylim = c(0,1))

barplot(colMeans(cor_matr)[ix_16S], ylim = c(0,1))









