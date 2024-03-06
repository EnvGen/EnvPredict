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
p <- arg_parser("clean_dataset.R")
p <- add_argument(p, "-b", help="16S_ASV sequences", default="../seq_data/combined/16S/asvs_16S.txt")
p <- add_argument(p, "-e", help="18S_ASV sequences", default="../seq_data/combined/18S/asvs_18S.txt")
p <- add_argument(p, "-n", help="16S_count table", default="../seq_data/combined/16S/seqtab_16S.tsv")
p <- add_argument(p, "-y", help="18S_count table", default="../seq_data/combined/18S/seqtab_18S.tsv")
p <- add_argument(p, "-m", help="Metadata", default="../env_data/combined/physical_chemical_processed.tsv")
# p <- add_argument(p, "-a", help="Abiotic parameters", default="contextual_data/physical_chemical_data_2019-2020_2022-03-16.txt")
p <- add_argument(p, "-z", help="taxa_16S table", default="../seq_data/combined/16S/taxa_16S.tsv")
p <- add_argument(p, "-v", help="taxa_18S table", default="../seq_data/combined/18S/taxa_18S.tsv")
p <- add_argument(p, "-w", help="working directory", default="//wsl.localhost/Ubuntu/home/krzjur/EnvPredict/code")
# p <- add_argument(p, "-d", help="Minimum depth", default=10)
p <- add_argument(p, "-o", help="output directory", default="../output/")



argv <- parse_args(p)

for (s in list_of_packages) { suppressPackageStartupMessages(library(s, character.only = TRUE))}


setwd(argv$w)
dir.create(argv$o)

##Functions

export_table<-function(dfcount,dftax,file_name,tax_level ) {
  #dfcount=norm_seqtab_16S;dftax=taxa_16S;file_name="norm_seqtab_16S.tsv" 
  dataset_VAE <- data.frame(dfcount)
  taxa_VAE<-data.frame(dftax)
  
  taxa_VAE_tmp = taxa_VAE %>% filter(!is.na(.[[tax_level]]))
  
  dataset_VAE_temp <- dataset_VAE[row.names(dataset_VAE) %in% row.names(taxa_VAE_tmp),]
  
  new_row_names=rep("", nrow(dataset_VAE_temp))
  tokeep=rep(0, nrow(dataset_VAE_temp))
  for ( i in 1:nrow(taxa_VAE_tmp)) {
    tempotext=row.names(taxa_VAE_tmp)[i]
    if (tempotext %in% row.names(dataset_VAE_temp)) {
      for (c in 1:ncol(taxa_VAE_tmp)) {
        if (!is.na(taxa_VAE_tmp[i,c])) tempotext=paste(tempotext,taxa_VAE_tmp[i,c],sep="_")
      }
      new_row_names[i]  <-tempotext
      tokeep[i]=i
    } 
    
    
  }
  
  
  new_row_names_clean=new_row_names[tokeep[tokeep>0]]
  
  names(dataset_VAE_temp)<-gsub("X","",names(dataset_VAE_temp))
  row.names(dataset_VAE_temp) <-new_row_names_clean
  
  
  dataset_VAE_final <- data.frame(Samples=names(dataset_VAE_temp),    # Append new column to front of data
                                  data.frame(t(dataset_VAE_temp)))
  
  fwrite(dataset_VAE_final, file=paste(argv$o, paste(tax_level,file_name,sep="_"), sep="/"), quote=FALSE, sep='\t', row.names = F)
  
  
  fwrite(dataset_VAE_temp, file=paste(argv$o, paste("T",tax_level,file_name,sep="_"), sep="/"), quote=FALSE, sep='\t', row.names = F)
  
  
  return(dataset_VAE_final)
}

RF_analysis<-function(MDF, target, myNtree, dirout, with_CV, with_opt_mtry) {
  prefi=paste0("RF_",target)
  #MDF=DF; target="Salinity"
  dir.create(dirout, recursive = T)
  set.seed(123)
  
  patern=paste("ASV",target,  sep="|")
  X0=MDF[,grep(patern,names(MDF),ignore.case = T)]
  X0=na.omit(X0)
  
  full_name=names(MDF)[grep(target,names(MDF),ignore.case = T)]
  names(X0)[grep(full_name,names(X0) )]<-"Response"
  
  df_split <- initial_split(X0, strata = Response)
  X <- training(df_split)
  df_test <- testing(df_split)
  
  
  if (with_opt_mtry) {
    try(bestmtry <- tuneRF(X[, !names(X) %in% "Response"], X$Response,ntreeTry = myNtree, stepFactor = 1.2, improve = 0.01, trace=T, plot= T), silent = TRUE)
    if (exists("bestmtry")) {
      cat("INFO: Parameter mtry has been optimised\n")
      Mtry=bestmtry[match(min(bestmtry[,2]), bestmtry[,2]),1]
      
      rf2 <- randomForest(X[,!names(X) %in% "Response"],X$Response, mtry=Mtry,ntree = myNtree, keep.forest=TRUE, importance=TRUE )
    } else {
      cat("INFO: Parameter mtry couldn't be optimised\n")
      rf2 <- randomForest(X[,!names(X) %in% "Response"],X$Response,ntree = myNtree, keep.forest=TRUE, importance=TRUE ) 
      }
  } else {
    cat("INFO: Using default mtry value\n")
    rf2 <- randomForest(X[,!names(X) %in% "Response"],X$Response,ntree = myNtree, keep.forest=TRUE, importance=TRUE )   
  }
  
  result<-rfcv(X[!names(X) %in% "Response"], X$Response)
  NV<-as.integer(names(which.min(result$error.cv)))
  total_features=length(names(X))-1
  if (with_CV) { 
 # if (NV >= nrow(X) & NV < total_features && with_CV) {  
    cat("INFO: Total features will be reduced from ", total_features, "to ", NV, "\n")
    imp<-as.data.frame(importance(rf2))
    imp<-imp[order(imp[[1]], decreasing = TRUE), ]
    imp<-imp[1:NV,]
    
    selected_seqs=c("Response",row.names(imp))
    
    Xs<-X[, (names(X) %in% selected_seqs)]
    if (with_opt_mtry) {
      try(bestmtry2 <- tuneRF(Xs[, !names(Xs) %in% "Response"], Xs$Response,ntreeTry = myNtree,stepFactor = 1.2, improve = 0.01, trace=T, plot= T), silent = TRUE)
      if (exists("bestmtry2")) {
        cat("INFO: Parameter mtry has been optimised\n")
        Mtry2=bestmtry2[match(min(bestmtry2[,2]), bestmtry2[,2]),1]
        rf <- randomForest(Xs[,!names(Xs) %in% "Response"],Xs$Response, mtry=Mtry2,ntree = myNtree, keep.forest=TRUE, importance=TRUE )
      } else {  
        cat("INFO: Parameter mtry couldn't be optimised\n")
        rf <- randomForest(Xs[,!names(Xs) %in% "Response"],Xs$Response,ntree = myNtree, keep.forest=TRUE, importance=TRUE )
      }
    } else {
      cat("INFO: Using default mtry value\n")
      rf <- randomForest(Xs[,!names(Xs) %in% "Response"],Xs$Response,ntree = myNtree, keep.forest=TRUE, importance=TRUE )
    }
    
  } else {
    cat("INFO: Total number of features used ", total_features,". Cross-validation suggests ",NV," features could be enough\n")
    rf<-rf2
    
  }
  
  seq_import<-as.data.frame(varImpPlot(rf))
  
  seq_import$Type="ASV"
  #seq_import$Type[grep("ASV.*",row.names(seq_import))]<-"ASV"
  seq_import$Type[grep("Archaea",row.names(seq_import))]<-"Archaea"
  seq_import$Type[grep("Bacteria",row.names(seq_import))]<-"Bacteria"
  
  newnames=c()
  vec=row.names(seq_import)
  for (i in 1:nrow(seq_import)) { 
    if (length(grep("ASV_",vec)) != 0) {
      m=strsplit(vec[i], split = "\\.")[[1]]
      n=paste(m[1],m[length(m)-1], sep=" ")
    } else {
      n=vec[i]
    }
    newnames<-c(newnames,n)
  }
  
  row.names(seq_import)<-newnames
  
  seq_import=seq_import[order(seq_import$`%IncMSE`, decreasing = TRUE),]
  #SEL=row.names(seq_import)
  
  #re=list(SELEC=SEL, data=Xs, RF=rf)
  #Exporting resutls
  sink(paste(dirout,paste(prefi,"txt",sep="."), sep = "/"))
  print(rf)
  print("Importance")
  print(seq_import)
  sink()
  ###
  
  #Creating figures
  seq_import=seq_import[(seq_import$`%IncMSE` > 4 | seq_import$`%IncMSE` < -4),]
  gC <- ggplot(data = seq_import, aes(x = reorder(row.names(seq_import),`%IncMSE`), y = `%IncMSE`)) +
    geom_bar(stat = "identity",aes(fill=Type)) + ggtitle("Random Forest importance - %IncMSE > |4|") +
    theme_bw()+       
    guides(fill=guide_legend(title=""))+ theme(legend.title = element_text(size=12)) +
    scale_color_manual(values = c("#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666"))+
    theme(axis.title = element_blank(), axis.ticks.y = element_blank(), panel.border = element_blank(),axis.text.x = element_text(size=6)) +
    coord_flip()
  
  ggsave(paste(dirout,paste(prefi,"pdf",sep="."),sep = "/"), gC, width = 82, height = 42, units = "cm")
  
  predicted_tr<-predict(rf, X, type="response")
  ggDF=data.frame(Actual=X$Response,Predicted=predicted_tr)
  wTr<-ggplot(data=ggDF, aes(x=Actual,y=Predicted))+geom_point()+geom_smooth(method =lm )+stat_cor(label.y = max(X$Response))+theme_bw()
  
  ggsave(paste(dirout,paste("train",prefi,"pdf",sep="."),sep = "/"), wTr, width = 42, height = 42, units = "cm")
  
  
  predicted_te<-predict(rf, df_test[,grep("ASV",names(df_test))], type="response")
  ggDFtest=data.frame(Actual=df_test$Response,Predicted=predicted_te)
  wTe<-ggplot(data=ggDFtest, aes(x=Actual,y=Predicted))+geom_point()+geom_smooth(method =lm )+stat_cor(label.y = max(df_test$Response))+theme_bw()
  
  ggsave(paste(dirout,paste("test",prefi,"pdf",sep="."),sep = "/"), wTe, width = 42, height = 42, units = "cm")
  
  
  #return(re)
}

xgb_analysis<-function(MDF, target, dirout) {
  #MDF=DF; target="salinity"; dirout=Dout2
  dir.create(dirout, recursive = T)
  patern=paste("ASV",target,  sep="|")
  
  
  X0=MDF[,grep(patern,names(MDF),ignore.case = T)]
  X0=na.omit(X0)
  
  full_name=names(MDF)[grep(target,names(MDF),ignore.case = T)]
  names(X0)[grep(full_name,names(X0) )]<-"Response"
  
  set.seed(123)
  
  df_split <- initial_split(X0, strata = Response)
  X <- training(df_split)
  df_test <- testing(df_split)
  
  
  prefi=paste0("xgb_",target)
  dtrain <- xgb.DMatrix(data = as.matrix(X[,!names(X) %in% "Response"]), label=X$Response)
  dtest <- xgb.DMatrix(data = as.matrix(df_test[,!names(df_test) %in% "Response"]), label=df_test$Response)
  
  df_split2 <- initial_split(X, strata = Response)
  X2 <- training(df_split2)
  df_valid <- testing(df_split2)
  dtrain2 <- xgb.DMatrix(data = as.matrix(X2[,!names(X2) %in% "Response"]), label=X2$Response)
  dvalid <- xgb.DMatrix(data = as.matrix(df_valid[,!names(df_valid) %in% "Response"]), label=df_valid$Response)
  
  watchlist <- list(train=dtrain2, val=dvalid)
  #watchlist <- list(train=dtrain)
  ###
  params <- list(objective = "reg:squarederror", gamma=5) #booster = "gbtree", 
  xgbcv <- xgb.cv( params = params, data = dtrain, nrounds = 100, nfold = 5, showsd = T, verbose=F, early.stop.round = 20, maximize = F)
  #bst <- xgb.train(params = params, data = dtrain, nrounds = xgbcv$best_iteration, watchlist =watchlist , print.every.n = 10, early.stop.round = 10, maximize = F , eval_metric = "error")
  #list(val=dtest,train=dtrain)
  ###
  
  bst <- xgb.train(data=dtrain, nrounds=xgbcv$best_iteration, watchlist=watchlist, eval.metric = "error", #eval.metric = "logloss", 
                   objective = "reg:squarederror", verbose=0)
  
  cat("INFO: Optimised nrounds: ", xgbcv$best_iteration,"\n")
  
  predtr <- predict(bst, dtrain)
  ggDF=data.frame(Actual=X$Response,Predicted=predtr)
  wTr<-ggplot(data=ggDF, aes(x=Actual,y=Predicted))+geom_point()+geom_smooth(method =lm )+stat_cor(label.y = max(X$Response))+theme_bw()
  
  ggsave(paste(dirout,paste("train",prefi,"pdf",sep="."),sep = "/"), wTr, width = 42, height = 42, units = "cm")
  
  
  
  pred <- predict(bst, dtest)
  
  ggDFtest=data.frame(Actual=df_test$Response,Predicted= pred)
  wTe<-ggplot(data=ggDFtest, aes(x=Actual,y=Predicted))+geom_point()+geom_smooth(method =lm )+stat_cor(label.y = max(df_test$Response))+theme_bw()
  
  ggsave(paste(dirout,paste("test",prefi,"pdf",sep="."),sep = "/"), wTe, width = 42, height = 42, units = "cm")
  
  
  importance_matrix <- xgb.importance(model = bst)
  
  importance_matrix2 <- importance_matrix[order(importance_matrix$Gain),]
  
  sink(paste(dirout,paste(prefi,"txt",sep="."), sep = "/"))
  print(importance_matrix2)
  sink()
  
  
  xgb_ggplot <- xgb.ggplot.importance(importance_matrix = importance_matrix[1:60]) +theme_bw() #,rel_to_first = TRUE, xlab = "Relative importance")
  
  ggsave(paste(dirout,paste(prefi,"pdf",sep="."),sep = "/"), xgb_ggplot, width = 82, height = 42, units = "cm")
  gr<-xgb.plot.tree(model = bst,render=FALSE) #, fname=paste(dirout,paste(prefi,"tree2","pdf",sep="."),sep = "/"))
  
  export_graph(gr, paste(dirout,paste(prefi,"tree","pdf",sep="."),sep = "/"), width=1500, height=1900)
}

#### Anders code
# load sequence tables with counts
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
taxa_18S = as.matrix(read.table(argv$v, sep = '\t'))
#### remove Barnp taxonomy affected 
#load metadata file
metadata = read.delim(paste(argv$m))
#### Other data (phyto,ect)
#these are the spike DNA sequences
spike_16S = 'CGGGCAGCTCTCGATAACCGGCGGAAGGTGGTAGCCACGGACAGGATCAGAACAATTAGAAGTGCCGCAGGTGGCCAAGTCCCCCGGACACAAGACGAGGCCGGAGGCCTGGTATATACACGTAGCTAAGAAGAGCTCATCCAGACTGGGAACGGTGTGCCAGCAGCCGCGGTAACATCACCACAACGTATTCGGTCACAAATTGATCGGAGGGAGAAATCGTCCGCAGGATCTCAAACTTTAACTAAGGACTAGTACTACATAGGCTCGAGAAGAGCTACCGTTTGCAGGGTCGCCGGGTACCGCTTAACCATAAAAGATCCACTCAGGTAGCCGTCCAGTTTCCTCTGAAATGATGGGGCGAGAAACACGGCTGGGCGTTATACGAGTGCTTTAGAATATGAGGAGAGACAGGGGTATATTCAAGG'
spike_18S = 'CTTCGTTATCGTCACGGAGAGAACTGCCTTTAGCGATCTGTCTAGAGACCCGCCAAATATAAGCCTTGGGGTTATTAACAAAGACGCCAAACACTAGTGAATATGACAGTGAAGGGGTGTGAGATAGCTTTACTGGGTGAGAAAACACTCGTTAAAAAGAATTAGACCGGATAATCCCCGAGGGGCCGTAGGCATGGACTTGTCGTTGCCACCGAGCATAGCGGTTTCGAAATAGCCGAGATGGGCACTGGCGAATTAACCCACTGGTTTATATGGATCCGATGGGTTCACTTAATAAGCTCGTACCAGGGATGAATAAAGCGTTACGAGAATTATAAACATGGAGTTCCTATTGATTTGAGGTTAATACCGAACGGGAACATTTGTCGATCATGCTTCACATAGAGT'

# For 16S - identify spike ASVs
spike_16S_ix = agrep(spike_16S, asv_16S)
spiketab_16S = seqtab_16S[spike_16S_ix, ]
dim(spiketab_16S)

## Check if spike ASVs are unannotated
# taxa_16S[spike_16S_ix,]

# For 18S - identify spike ASVs
spike_18S_ix = agrep(spike_18S, asv_18S)
spiketab_18S = seqtab_18S[spike_18S_ix, ]

## Exclude spike reads and Metazoa from taxa 

ix_taxa_16S =  setdiff(1:nrow(taxa_16S), spike_16S_ix)
taxa_16S = taxa_16S[ix_taxa_16S,]
seqtab_16S = seqtab_16S[ix_taxa_16S,]

ix_taxa_18S =  setdiff(1:nrow(taxa_18S), spike_18S_ix)
ix_taxa_18S = intersect(ix_taxa_18S, which(taxa_18S[,4] != 'Metazoa' | is.na(taxa_18S[,4])))
ix_taxa_18S = union(ix_taxa_18S, which(is.na(taxa_18S[,3])))
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
metadata = metadata[ix,]
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

## Discuss - one more sample to remove?

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
metadata = metadata[ix,]
seqtab_16S = seqtab_16S[,ix]

## Get only samples with quality 16S AND 18S data - discuss

common_samples = intersect(colnames(seqtab_16S),
                           colnames(seqtab_18S))
seqtab_16S = seqtab_16S[,common_samples]
seqtab_18S = seqtab_16S[,common_samples]

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

# Rarefying reads per sample to the same counts - not needed (discuss)
# For 16S
m=min(colSums(seqtab_16S))  # should we use alfat/beta diversity to select m?
r_seqtab_16S = t(rrarefy(x = t(seqtab_16S), sample = m))
colSums(r_seqtab_16S)
r_seqtab_16S = t(r_seqtab_16S)
# which(colSums(seqtab) == min(colSums(seqtab)))

r_norm_seqtab_16S = r_seqtab_16S
for (i in 1:ncol(r_seqtab_16S)) {
  r_norm_seqtab_16S[,i] = r_seqtab_16S[,i]/sum(r_seqtab_16S[,i])
}

# For 18S
m=min(colSums(seqtab_18S))
r_seqtab_18S = rrarefy(x = t(seqtab_18S), sample = m)
colSums(r_seqtab_18S)
r_seqtab_18S = t(r_seqtab_18S)
# which(colSums(seqtab) == min(colSums(seqtab)))

r_norm_seqtab_18S = r_seqtab_18S
for (i in 1:ncol(r_seqtab_18S)) {
  r_norm_seqtab_18S[,i] = r_seqtab_18S[,i]/sum(r_seqtab_18S[,i])
}

### FInished reading data ###

### Anders code ends here
#Average of Temp,Salinity, Lat, Lon per station
#abiotics=metadata 

## I do not understand this part (KTJ) - discuss

colnames(metadata)
abiotics=metadata[,c("Salinity", "Temperature", "pH", "Alkalinity",
                     "SiO3", "N_tot", "DIN", "NH4",
                     "NO3_NO2", "P_tot", "Phosphate",
                     "DOC", "Humus", "Chl", "Chl_hose", "Chl_combined")]

rest_metadata=metadata[,c("Longitude", "Latitude","Secchi_depth")]
  
rest_metadata=rest_metadata %>% mutate_if(is.character, function(x) as.numeric(x))

ab_factors=c("Salinity", "Temperature", "pH", "Alkalinity",
             "SiO3", "N_tot", "DIN", "NH4",
             "NO3_NO2", "P_tot", "Phosphate",
             "DOC", "Humus", "Chl", "Chl_hose", "Chl_combined")

for (patern in ab_factors) {
  item=paste0(patern,"_average")
  abiotics[[item]]=apply(rest_metadata[,grep(patern,names(rest_metadata),ignore.case = T)],1, function(x) mean(x,
    na.rm=T))
}

# abiotics$salinity_average=as.numeric(metadata$salinity_average)

#abiotics=abiotics[,c(5,16,17,28,40,45,50,53,56,59,62,65,68,71,74,77,80,83,86)]

names(abiotics)[1] <- "Samples"

abiotics$Samples<-factor(abiotics$Samples, levels = unique(abiotics$Samples))
#str(abiotics)
abiotics=abiotics %>% mutate_if(is.character, function(x) as.numeric(x))

abiotics_means=aggregate(abiotics, list(abiotics$Samples), function(x) mean(x,na.rm = T))
abiotics_means$Samples=abiotics_means$Group.1
abiotics_means=abiotics_means[,!names(abiotics_means) %in% c("Group.1")]
fwrite(abiotics_means, file=paste(argv$o, "abiotic_avg_factors.tsv", sep="/"), quote=FALSE, sep='\t', row.names = F)

####



#"Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"   "Species"
Tax_Level="Genus"

DF_norm_16S=export_table(norm_seqtab_16S,taxa_16S,"norm_seqtab_16S.tsv",Tax_Level)

DF_r_norm_16S=export_table(r_norm_seqtab_16S,taxa_16S,"r_norm_seqtab_16S.tsv",Tax_Level)

missing_data_stations=abiotics_means$Samples[! abiotics_means$Samples %in% DF_r_norm_16S$Samples]

DF=merge(DF_r_norm_16S,abiotics_means, by="Samples")
DF=DF[,!names(DF) %in% c("Samples")]

extra="with_default"
Dout=paste(argv$o,paste("ML_average",Tax_Level,extra,sep="_"),"RF",sep="/")

ab_factors=c("Salinity", "Temperature",  "oxygen",  "ammonium",  "phosphate", "phosphoru",  "nitrogen",  "silicate",  "humic_substance", "chlorophyll")

for (fc in ab_factors) { #fc="Salinity"
  cat("RF with ",fc, " \n")
  RF_analysis(DF, fc, 1000, Dout, with_CV=F, with_opt_mtry=F) 
}


Dout2=paste(argv$o,paste("ML_avg",Tax_Level,sep="_"),"xgb",sep="/")

ab_factors=c("Salinity", "Temperature",  "oxygen",  "ammonium",  "phosphate", "phosphoru",  "nitrogen",  "silicate",  "humic_substance", "chlorophyll")
#fc="oxygen"
for (fc in ab_factors) {
  cat("xgb with ",fc, " \n")
#xgb_analysis(DF, fc, 7, 7,Dout2)
#  xgb_analysis(DF, fc, Dout2)
  }

# save model to binary local file
#xgb.save(bst, "xgboost.model")

# load binary model to R
#bst2 <- xgb.load("xgboost.model")
#pred2 <- predict(bst2, test$data)


