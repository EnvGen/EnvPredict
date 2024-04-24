## Change working directory (add_argument(p, "-w", help="working directory",...)
## Before running on a new machine

#Other RF (e.g., mtry) and XGBOOST (e.g., nround) hyperparameters can be also tuned 

rm(list = ls())
options(warn=-1)

#1. Installing packages

if ( ! "argparser" %in% installed.packages()[,"Package"]) {
  installed.packages("argparser")
}

#if ( ! "DiagrammeRsvg" %in% installed.packages()[,"Package"]) {
#  devtools::install_github('rich-iannone/DiagrammeRsvg')
#}


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


# 2. Setting up argparser and paths

suppressPackageStartupMessages(library("argparser"))

# arguments
p <- arg_parser("RFandXGB.R [options]")
p <- add_argument(p, "-p", help="Directory with 16S_norm_clade_counts_taxlevel.stv", default="../seq_data/combined/16S")
p <- add_argument(p, "-e", help="Directory with 18S_norm_clade_counts_taxlevel.stv", default="../seq_data/combined/18S")
p <- add_argument(p, "-m", help="Metadata", default="/Users/luisdelgado/Documents/EnvPredict/seq_data/combined/filtered_metadata.tsv")
p <- add_argument(p, "-a", help="Directory with VAE latent features files", default="/Users/luisdelgado/Documents/EnvPredict/seq_data/combined/RepresentationsFromDeepMicro")
p <- add_argument(p, "-b", help="plankton factors", default="/Users/luisdelgado/Documents/EnvPredict/env_data/combined/zooplankton_processed.tsv")
p <- add_argument(p, "-t", help="Biotic factors", default="/Users/luisdelgado/Documents/EnvPredict/env_data/combined/physical_chemical_processed_translation.tsv")
p <- add_argument(p, "-w", help="working directory", default="/Users/luisdelgado/Documents/EnvPredict/code")
p <- add_argument(p, "-o", help="output directory", default="../output")


argv <- parse_args(p)

for (s in list_of_packages) { suppressPackageStartupMessages(library(s, character.only = TRUE))}

setwd(argv$w)
dir.create(argv$o)


#Extra settings
only_RF_part=F
RF_extra="with_CV"
rf_with_CV=T
rf_with_opt_mtry=F
##Functions


RF_analysis<-function(counts,Abiotic, target, myNtree, dirout, with_CV, with_opt_mtry) {
#counts=DF_norm;Abiotic=abiotics;target=fc; myNtree=1000; dirout=Dout; with_CV=rf_with_CV; with_opt_mtry=rf_with_opt_mtry
  
  Fact=Abiotic[names(Abiotic) %in% c("sample_id", target) ]
  names(Fact)<- c("sample_id", "Response")
  MDF=merge(counts,Fact, by="sample_id")
  row.names(MDF)<-MDF$sample_id
  MDF=MDF[,!names(MDF) %in% c("sample_id")]
  
  prefi=paste0("RF_",target)
  dir.create(dirout, recursive = T)
  set.seed(123)
  
  X0=na.omit(MDF)
  cat("INFO: Total number of samples ", nrow(X0), "\n")


  train_and_test<-function(fold){ #fold=folds$splits[[1]]
    extratext=as.character(fold$id)
    cat("INFO: Working on", extratext, "\n")
    
    X <- fold %>% analysis() 
    df_test <- fold %>% assessment()
    cat("INFO: Number of samples (training)", nrow(X), "- Number of samples (testing)", nrow(df_test), "\n")
  
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
  
  
  if (argv$b != "") {seq_import$Type=""} else if (argv$a != "") {seq_import$Type="Latent feature"} else { seq_import$Type="Amplicon"}
  seq_import$Type[grep("Archaea",row.names(seq_import))]<-"Archaea"
  seq_import$Type[grep("Bacteria",row.names(seq_import))]<-"Bacteria"
  seq_import$Type[grep("Eukaryote",row.names(seq_import))]<-"Eukaryote"
  
  
  seq_import=seq_import[order(seq_import$`%IncMSE`, decreasing = TRUE),]
  
  #Exporting resutls
  sink(paste(dirout,paste(paste(prefi,extratext, sep="_"),"txt",sep="."), sep = "/"))
  print(rf)
  print("Importance")
  print(seq_import)
  sink()
  ###
  
  #Creating figures
  seq_import=seq_import[(seq_import$`%IncMSE` > 4),]
  gC <- ggplot(data = seq_import, aes(x = reorder(row.names(seq_import),`%IncMSE`), y = `%IncMSE`)) +
    geom_bar(stat = "identity",aes(fill=Type)) + ggtitle("Random Forest importance - %IncMSE > |4|") +
    theme_bw()+       
    guides(fill=guide_legend(title=""))+ theme(legend.title = element_text(size=12)) +
    scale_color_manual(values = c("#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666"))+
    theme(axis.title = element_blank(), axis.ticks.y = element_blank(), panel.border = element_blank(),axis.text.x = element_text(size=6)) +
    coord_flip()
  
  ggsave(paste(dirout,paste(paste(prefi,extratext, sep="_"),"pdf",sep="."),sep = "/"), gC, width = 82, height = 42, units = "cm")
  
  predicted_tr<-predict(rf, X[,!names(X) %in% "Response"], type="response")
  ggDF=data.frame(Actual=X$Response,Predicted=predicted_tr)
  wTr<-ggplot(data=ggDF, aes(x=Actual,y=Predicted))+geom_point()+geom_smooth(method =lm, formula = 'y ~ x' )+stat_cor(label.y = max(X$Response))+theme_bw()
  
  ggsave(paste(dirout,paste("train",paste(prefi,extratext, sep="_"),"pdf",sep="."),sep = "/"), wTr, width = 42, height = 42, units = "cm")
  
  predicted_te<-predict(rf, df_test[,!names(df_test) %in% "Response"], type="response")
  
  #---- save prediction
  table_pred<-data.frame(sample_id=row.names(df_test), target=predicted_te) #, Actual=df_test$Response)
  names(table_pred)[2]<-paste0("Predicted_",target)
  #---
  
  ggDFtest=data.frame(Actual=df_test$Response,Predicted=predicted_te)
  wTe<-ggplot(data=ggDFtest, aes(x=Actual,y=Predicted))+geom_point()+geom_smooth(method =lm, formula = 'y ~ x' )+stat_cor(label.y = max(df_test$Response))+theme_bw()
  
  ggsave(paste(dirout,paste("test",paste(prefi,extratext, sep="_"),"pdf",sep="."),sep = "/"), wTe, width = 42, height = 42, units = "cm")
  
  
  return(table_pred)
  }
  
  folds <- vfold_cv(X0, v = 5, strata = Response)
  predictions_list <- map(folds$splits, train_and_test)
  predictions <- bind_rows(predictions_list)
  
  return(predictions)
}

xgb_analysis<-function(counts,Abiotic,target, dirout) {
  
  #----
  Fact=Abiotic[names(Abiotic) %in% c("sample_id", target) ]
  names(Fact)<- c("sample_id", "Response")
  MDF=merge(counts,Fact, by="sample_id")
  row.names(MDF)<-MDF$sample_id
  MDF=MDF[,!names(MDF) %in% c("sample_id")]
  X0=na.omit(MDF)
  cat("INFO: Total number of samples ", nrow(X0), "\n")
  #---
  
  dir.create(dirout, recursive = T)
  set.seed(123)
  
  
  ###
  train_and_test_XGB<-function(fold){ #fold=folds$splits[[1]]
    extratext=as.character(fold$id)
    cat("INFO: Working on", extratext, "\n")
    
    X <- fold %>% analysis() 
    df_test <- fold %>% assessment()
    cat("INFO: Number of samples (training)", nrow(X), "- Number of samples (testing)", nrow(df_test), "\n")
  
  
  ###
  prefi=paste0("xgb_",target)
  dtrain <- xgb.DMatrix(data = as.matrix(X[,!names(X) %in% "Response"]), label=X$Response)
  dtest <- xgb.DMatrix(data = as.matrix(df_test[,!names(df_test) %in% "Response"]), label=df_test$Response)
  
  watchlist <- list(train=dtrain)
  ###
  params <- list(objective = "reg:squarederror", gamma=5)  
  xgbcv <- xgb.cv( params = params, data = dtrain, nrounds = 100, nfold = 5, showsd = T, verbose=F, early.stop.round = 20, maximize = F)
  ###
  
  bst <- xgb.train(data=dtrain, nrounds=xgbcv$best_iteration, watchlist=watchlist, eval.metric = "error", #eval.metric = "logloss", 
                   objective = "reg:squarederror", verbose=0)
  
  cat("INFO: Optimised nrounds: ", xgbcv$best_iteration,"\n")
  
  predtr <- predict(bst, dtrain)
  ggDF=data.frame(Actual=X$Response,Predicted=predtr)
  wTr<-ggplot(data=ggDF, aes(x=Actual,y=Predicted))+geom_point()+geom_smooth(method =lm, formula = 'y ~ x' )+stat_cor(label.y = max(X$Response))+theme_bw()
  
  ggsave(paste(dirout,paste("train",paste(prefi,extratext, sep="_"),"pdf",sep="."),sep = "/"), wTr, width = 42, height = 42, units = "cm")
  
  pred <- predict(bst, dtest)
  
  #---- save prediction
  table_pred<-data.frame(sample_id=row.names(df_test), target=pred) #, Actual=df_test$Response)
  names(table_pred)[2]<-paste0("Predicted_",target)
  #---
  
  ggDFtest=data.frame(Actual=df_test$Response,Predicted= pred)
  wTe<-ggplot(data=ggDFtest, aes(x=Actual,y=Predicted))+geom_point()+geom_smooth(method =lm, formula = 'y ~ x' )+stat_cor(label.y = max(df_test$Response))+theme_bw()
  
  ggsave(paste(dirout,paste("test",paste(prefi,extratext, sep="_"),"pdf",sep="."),sep = "/"), wTe, width = 42, height = 42, units = "cm")
  
  importance_matrix <- xgb.importance(model = bst)
  
  importance_matrix2 <- importance_matrix[order(importance_matrix$Gain),]
  
  sink(paste(dirout,paste(paste(prefi,extratext, sep="_"),"txt",sep="."), sep = "/"))
  print(importance_matrix2)
  sink()
  
  
  xgb_ggplot <- xgb.ggplot.importance(importance_matrix = importance_matrix[1:60]) +theme_bw() 
  
  ggsave(paste(dirout,paste(paste(prefi,extratext, sep="_"),"pdf",sep="."),sep = "/"), xgb_ggplot, width = 82, height = 42, units = "cm")
  gr<-xgb.plot.tree(model = bst,render=FALSE) 
  
  export_graph(gr, paste(dirout,paste(paste(prefi,extratext, sep="_"),"tree","pdf",sep="."),sep = "/"), width=1500, height=1900)
  
  return(table_pred)
  }
  
  ###
  folds <- vfold_cv(X0, v = 5, strata = Response)
  predictions_list <- map(folds$splits, train_and_test_XGB)
  predictions <- bind_rows(predictions_list)
  
  return(predictions)
  ###
}

## Reading input data

norm_clade_counts_16S_files=list.files(argv$p, pattern = "norm_clade_counts_16S_\\d.tsv"  )
norm_clade_counts_18S_files=list.files(argv$e, pattern = "norm_clade_counts_18S_\\d.tsv"  )

if (argv$b != "") { # processing abiotic data
  abiot=as.data.frame(t(read.delim(argv$b, header = TRUE)))
  rownames(abiot)<-gsub("^X","", rownames(abiot))  
  samples_names=as.data.frame(read.delim(argv$t, header = TRUE))
  samples_names$sample<-gsub("-",".",samples_names$sample)
  snames<-samples_names %>% select(sample, sample_id)
  abiot$sample=rownames(abiot)
  abiotics=merge(abiot, snames)
  abiotics<-abiotics[,!names(abiotics) %in% "sample"]
  abiotics <- abiotics %>% select("sample_id", everything())
  
  ab_factors=names(abiotics)[!names(abiotics) %in% "sample_id"][1:2]
  add_sufix="Biotic"
  
} else {    
  abiotics=read_tsv(argv$m, show_col_types = FALSE) 
  ab_factors=names(abiotics)[10:21]
  add_sufix="Abiotic"
}


for (rRNA in c("16S","18S")) { #rRNA="16S"
  if (rRNA == "16S")  tax_level_numbers<-list("Specie"=7, "Genus"=6,"Family"=5, "Order"=4, "Class"=3)
  if (rRNA == "18S")  tax_level_numbers<-list("Specie"=9, "Genus"=8,"Family"=7, "Order"=6, "Class"=5)
  
  for (Tax_Level in names(tax_level_numbers)){ #Tax_Level="Specie"
    num=tax_level_numbers[[Tax_Level]]
    
    if (rRNA == "16S") DF_norm=read.delim(paste(argv$p,norm_clade_counts_16S_files[num], sep="/"), header = TRUE)
    if (rRNA == "18S") DF_norm=read.delim(paste(argv$e,norm_clade_counts_18S_files[num], sep="/"), header = TRUE)
    
    DF_norm <- as.data.frame(t(DF_norm))
    DF_norm$sample_id<-row.names(DF_norm)
    
    DF_norm$sample_id <- sub("^X", "", DF_norm$sample_id)
    
    Dout=paste(argv$o,paste("RF",RF_extra,sep="_"),paste(Tax_Level,rRNA, sep="_"),sep="/")
    
    table_pred<-vector(mode = "list", length = length(ab_factors))
    i=1
    for (fc in ab_factors) { #fc="Amphipoda"                       
      cat("*** RF with ",Tax_Level," - ",rRNA," - ",fc, " \n")
      
      table_pred[[i]]<-RF_analysis(DF_norm,abiotics, fc, 1000, Dout, with_CV=rf_with_CV, with_opt_mtry=rf_with_opt_mtry)
      i=i+1
    }
    
    final_pred_table<-table_pred %>% reduce(full_join, by='sample_id') %>% arrange(sample_id)
    
    prefix_files=paste("RF",RF_extra,Tax_Level,rRNA, add_sufix, sep="_")
    doutpred=paste(argv$o,"Summary", sep="/")
    dir.create(doutpred, recursive = T)
    write.table(final_pred_table, file=paste(doutpred,paste(prefix_files,"Predictions.tsv",sep="_"),sep = "/"), sep = '\t', row.names = F, col.names = T )
    
    final_actual_table <- abiotics %>% filter(sample_id %in% final_pred_table$sample_id) %>% select(sub("Predicted_","",names(final_pred_table))) %>% arrange(sample_id)
    write.table(final_actual_table, file=paste(doutpred,paste(prefix_files,"Actual_values.tsv",sep="_"),sep = "/"), sep = '\t', row.names = F, col.names = T )
    
    
    if (only_RF_part == FALSE) {
      Dout2=paste(argv$o,"XGB",paste(Tax_Level,rRNA,sep="_"),sep="/")
      table_pred_xgb<-vector(mode = "list", length = length(ab_factors))
      j=1
      for (fc in ab_factors) {
        cat("*** xgb with ",Tax_Level," - ",rRNA, " - ",fc, " \n")
        table_pred_xgb[[j]]<-xgb_analysis(DF_norm,abiotics, fc, Dout2)
        j=j+1
      }
      
      final_pred_table_xgb<- table_pred_xgb %>% reduce(full_join, by='sample_id') %>% arrange(sample_id)
      
      prefix_files=paste("XGB",Tax_Level,rRNA, add_sufix,sep="_")
      write.table(final_pred_table_xgb, file=paste(doutpred,paste(prefix_files,"Predictions.tsv",sep="_"),sep = "/"), sep = '\t', row.names = F, col.names = T )
      
      final_actual_table_xgb <- abiotics %>% filter(sample_id %in% final_pred_table_xgb$sample_id) %>% select(sub("Predicted_","",names(final_pred_table_xgb))) %>% arrange(sample_id)
      write.table(final_actual_table_xgb, file=paste(doutpred,paste(prefix_files,"Actual_values.tsv",sep="_"),sep = "/"), sep = '\t', row.names = F, col.names = T )
      
    }
  }
  
  if (argv$a != "") {
    if (rRNA == "16S") {
      VAE_features_files=list.files(argv$a, pattern = "16S_filt.csv")
      archt_numbers=as.vector(sapply(VAE_features_files,function(x) sub("_norm_seqtab_16S_filt.csv","", x)))
    }
    if (rRNA == "18S") {
      VAE_features_files=list.files(argv$a, pattern = "18S_filt.csv")
      archt_numbers=as.vector(sapply(VAE_features_files,function(x) sub("_norm_seqtab_18S_filt.csv","", x)))
    }
    
    
    for (archt in archt_numbers ) {
      num=grep(archt,VAE_features_files)
      DF_norm=read.delim(paste(argv$a,VAE_features_files[num], sep="/"), header = TRUE)
      names(DF_norm)[1]<-"sample_id"
      DF_norm$sample_id <- sub("^X", "", DF_norm$sample_id)
      
      Dout3=paste(argv$o,"VAE", "RF",paste(archt,rRNA,sep="_"),sep="/")
      
      table_pred<-vector(mode = "list", length = length(ab_factors))
      i=1
      for (fc in ab_factors) { 
        cat("*** RF with VAE ",archt, " - ",rRNA," - ", fc, " \n")
        table_pred[[i]]<-RF_analysis(DF_norm,abiotics, fc, 1000, Dout3, with_CV=rf_with_CV, with_opt_mtry=rf_with_opt_mtry)
        i=i+1
      }
      
      final_pred_table<-table_pred %>% reduce(full_join, by='sample_id') %>% arrange(sample_id)
      
      prefix_files=paste("VAE",archt,"RF",rRNA,add_sufix,sep="_")
      write.table(final_pred_table, file=paste(doutpred,paste(prefix_files,"Predictions.tsv",sep="_"),sep = "/"), sep = '\t', row.names = F, col.names = T )
      
      final_actual_table <- abiotics %>% filter(sample_id %in% final_pred_table$sample_id) %>% select(sub("Predicted_","",names(final_pred_table))) %>% arrange(sample_id)
      write.table(final_actual_table, file=paste(doutpred,paste(prefix_files,"Actual_values.tsv",sep="_"),sep = "/"), sep = '\t', row.names = F, col.names = T )
      
      
      if (only_RF_part == FALSE) {
        Dout2A=paste(argv$o,"VAE","XGB",paste(archt,rRNA,sep="_"),sep="/")
        table_pred_xgb<-vector(mode = "list", length = length(ab_factors))
        j=1
        for (fc in ab_factors) {
          cat("*** xgb with VAE ",archt," - ",rRNA, " - ",fc, " \n")
          table_pred_xgb[[j]]<-xgb_analysis(DF_norm,abiotics, fc, Dout2A)
          j=j+1
        }
        
        final_pred_table_xgb<- table_pred_xgb %>% reduce(full_join, by='sample_id') %>% arrange(sample_id)
        
        prefix_files=paste("VAE",archt,"XGB",rRNA,add_sufix,sep="_")
        write.table(final_pred_table_xgb, file=paste(doutpred,paste(prefix_files,"Predictions.tsv",sep="_"),sep = "/"), sep = '\t', row.names = F, col.names = T )
        
        final_actual_table_xgb <- abiotics %>% filter(sample_id %in% final_pred_table_xgb$sample_id) %>% select(sub("Predicted_","",names(final_pred_table_xgb))) %>% arrange(sample_id)
        write.table(final_actual_table_xgb, file=paste(doutpred,paste(prefix_files,"Actual_values.tsv",sep="_"),sep = "/"), sep = '\t', row.names = F, col.names = T )
        
      }  
    }
  }
  
}    

