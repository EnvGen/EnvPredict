rm(list = ls())
options(warn=-1)

#libraries
list_of_packages <- c("tidyverse","plyr","ggpubr","DiagrammeR","tidymodels", "randomForest", "xgboost","future") 

#1. Setting up argparser and paths

suppressPackageStartupMessages(library("argparser"))

# arguments
p <- arg_parser("RF_XGB_ab.R [options]")
p <- add_argument(p, "-p", help="Directory with 16S_norm_clade_counts_taxlevel.stv", default="../seq_data/combined/16S")
p <- add_argument(p, "-e", help="Directory with 18S_norm_clade_counts_taxlevel.stv", default="../seq_data/combined/18S")
p <- add_argument(p, "-a", help="Directory with VAE latent features files", default="../RepresentationsFromDeepMicro")
p <- add_argument(p, "-m", help="Physicochemical parameters to be predicted", default="../env_data/combined/physical_chemical_processed.tsv")
p <- add_argument(p, "-b", help="plankton factors to be predicted", default="")
p <- add_argument(p, "-t", help="Translation samples IDs, required when importing data from plankton factors", default="../env_data/combined/sample_id_translation.tsv")
p <- add_argument(p, "-w", help="working directory", default="")
p <- add_argument(p, "-o", help="output directory", default="../output")

argv <- parse_args(p)

for (s in list_of_packages) { suppressPackageStartupMessages(library(s, character.only = TRUE))}

setwd(argv$w)
dir.create(argv$o)
plan(multisession)
#Extra settings
only_RF_part=F
RF_extra="withCV_mtryopt"
rf_with_CV=T
rf_with_opt_mtry=T
minSamp=10 #Features present in less than 100 samples are not included in the analysis
##Functions


RF_analysis<-function(counts,Abiotic, target, myNtree, dirout, with_CV, with_opt_mtry,type_variable) {
  
  train_and_test<-function(fold){ #fold=folds$splits[[1]]
    extratext=as.character(fold$id)
    cat("INFO: Working on", extratext, "\n")
    
    X <- fold %>% analysis() 
    df_test <- fold %>% assessment()
    cat("   Number of samples (training)", nrow(X), "- Number of samples (testing)", nrow(df_test), "\n")
    
    if (with_opt_mtry && with_CV == F) {
      tryCatch(bestmtry <- tuneRF(X[, !names(X) %in% "Response"], X$Response,ntreeTry = myNtree, stepFactor = 1.2, improve = 0.01, trace=F, plot= F), error = function(e){cat("   Parameter mtry couldn't be optimised\n")})
      if (exists("bestmtry")) {
        cat("   Parameter mtry has been optimised\n")
        Mtry=bestmtry[match(min(bestmtry[,2]), bestmtry[,2]),1]
        
        rf2 %<-% {randomForest(X[,!names(X) %in% "Response"],X$Response, mtry=Mtry,ntree = myNtree, keep.forest=TRUE, importance=TRUE )}
      }else {
        rf2 %<-% {randomForest(X[,!names(X) %in% "Response"],X$Response,ntree = myNtree, keep.forest=TRUE, importance=TRUE ) }
      }
    } else {
      cat("   Using default mtry value\n")
      rf2 %<-% {randomForest(X[,!names(X) %in% "Response"],X$Response,ntree = myNtree, keep.forest=TRUE, importance=TRUE ) }  
    }
    
    total_features=length(names(X))-1
    
    if (with_CV) { 
      result %<-% {rfcv(X[!names(X) %in% "Response"], X$Response, cv.fold=5)}
      NV<-as.integer(names(which.min(result$error.cv)))
      cat("   Total features will be reduced from ", total_features, "to ", NV, "\n")
      imp<-as.data.frame(importance(rf2))
      imp<-imp[order(imp[[1]], decreasing = TRUE), ]
      imp<-imp[1:NV,]

      selected_seqs=c("Response",row.names(imp))
      Xs<-X[, (names(X) %in% selected_seqs)]
      Xs_no_response <- as.data.frame(Xs[, !names(Xs) %in% "Response"]) ## For it not to crash if NV = 1
      if(NV == 1){
        colnames(Xs_no_response) <- selected_seqs[2]
      }
      if (with_opt_mtry) {
        tryCatch(bestmtry2 <- tuneRF(Xs_no_response, Xs$Response,ntreeTry = myNtree,stepFactor = 1.2, improve = 0.01, trace=F, plot= F), error = function(e){cat("   Parameter mtry couldn't be optimised\n")})
        if (exists("bestmtry2")) {
          cat("   Parameter mtry has been optimised\n")
          Mtry2=bestmtry2[match(min(bestmtry2[,2]), bestmtry2[,2]),1]
          rf %<-% {randomForest(Xs_no_response,Xs$Response, mtry=Mtry2,ntree = myNtree, keep.forest=TRUE, importance=TRUE )}
        } else {  
          rf %<-% {randomForest(Xs_no_response,Xs$Response,ntree = myNtree, keep.forest=TRUE, importance=TRUE )}
        }
      } else {
        cat("   Using default mtry value\n")
        rf %<-% {randomForest(Xs_no_response,Xs$Response,ntree = myNtree, keep.forest=TRUE, importance=TRUE )}
      }
      
    } else {
      # cat("INFO: Total number of features used ", total_features,". Cross-validation suggests ",NV," features could be enough\n")
      rf<-rf2 ## WHY?
      cat("   INFO: Total number of features used ", total_features,"\n")
    }
    
    seq_import<-as.data.frame(varImpPlot(rf))
    seq_import$Type=type_variable
    
    seq_import$Type[grep("Archaea",row.names(seq_import))]<-"Archaea"
    seq_import$Type[grep("Bacteria",row.names(seq_import))]<-"Bacteria"
    seq_import$Type[grep("Eukaryote",row.names(seq_import))]<-"Eukaryote"
    
    
    seq_import=seq_import[order(seq_import$`%IncMSE`, decreasing = TRUE),]
    
    #Exporting results
    sink(paste(dirout,paste(paste(prefi,extratext, sep="_"),"txt",sep="."), sep = "/"))
    print(rf)
    print("Importance")
    print(seq_import)
    sink()
    
    
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
    
    
    ggDFtest=data.frame(Actual=df_test$Response,Predicted=predicted_te)
    wTe<-ggplot(data=ggDFtest, aes(x=Actual,y=Predicted))+geom_point()+geom_smooth(method =lm, formula = 'y ~ x' )+stat_cor(label.y = max(df_test$Response))+theme_bw()
    
    ggsave(paste(dirout,paste("test",paste(prefi,extratext, sep="_"),"pdf",sep="."),sep = "/"), wTe, width = 42, height = 42, units = "cm")
    
    
    return(table_pred)
  }
  
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
  
  folds <- vfold_cv(X0, v = 5, strata = Response)
  predictions_list <- map(folds$splits, train_and_test)
  predictions <- bind_rows(predictions_list)
  return(predictions)
}

xgb_analysis<-function(counts,Abiotic,target, dirout) {
  
  train_and_test_XGB<-function(fold){ #fold=folds$splits[[1]]
    extratext=as.character(fold$id)
    cat("INFO: Working on", extratext, "\n")
    
    X <- fold %>% analysis() 
    df_test <- fold %>% assessment()
    cat("   Number of samples (training)", nrow(X), "- Number of samples (testing)", nrow(df_test), "\n")
    
    prefi=paste0("xgb_",target)
    dtrain <- xgb.DMatrix(data = as.matrix(X[,!names(X) %in% "Response"]), label=X$Response)
    dtest <- xgb.DMatrix(data = as.matrix(df_test[,!names(df_test) %in% "Response"]), label=df_test$Response)
    
    watchlist <- list(train=dtrain)
    
    params <- list(objective = "reg:squarederror", gamma=5)  
    if(length(unique(X[,'Response'])) >= 5){
      xgbcv <- xgb.cv( params = params, data = dtrain, nrounds = 100, nfold = 5, showsd = T, verbose=F, early.stop.round = 20, maximize = F)}else{
      xgbcv <-  xgb.cv( params = params, data = dtrain, nrounds = 100, nfold = 5, showsd = T, verbose=F, early.stop.round = 20, maximize = F, stratified = FALSE)
      }
    
    
    bst <- xgb.train(data=dtrain, nrounds=xgbcv$best_iteration, watchlist=watchlist, eval.metric = "error", #eval.metric = "logloss", 
                     objective = "reg:squarederror", verbose=0)
    
    cat("   Optimised nrounds: ", xgbcv$best_iteration,"\n")
    
    predtr <- predict(bst, dtrain)
    ggDF=data.frame(Actual=X$Response,Predicted=predtr)
    wTr<-ggplot(data=ggDF, aes(x=Actual,y=Predicted))+geom_point()+geom_smooth(method =lm, formula = 'y ~ x' )+stat_cor(label.y = max(X$Response))+theme_bw()
    
    ggsave(paste(dirout,paste("train",paste(prefi,extratext, sep="_"),"pdf",sep="."),sep = "/"), wTr, width = 42, height = 42, units = "cm")
    
    pred <- predict(bst, dtest)
    
    #---- save prediction
    table_pred<-data.frame(sample_id=row.names(df_test), target=pred) #, Actual=df_test$Response)
    names(table_pred)[2]<-paste0("Predicted_",target)
    
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
    tryCatch(gr<-xgb.plot.tree(model = bst,render=FALSE), error = function(e) {cat("   Non-tree model (no tree plot)")})
    if(exists("gr")){
      export_graph(gr, paste(dirout,paste(paste(prefi,extratext, sep="_"),"tree","pdf",sep="."),sep = "/"), width=1500, height=1900)
    }
    
    return(table_pred)
  }
  
  
  Fact=Abiotic[names(Abiotic) %in% c("sample_id", target) ]
  names(Fact)<- c("sample_id", "Response")
  MDF=merge(counts,Fact, by="sample_id")
  row.names(MDF)<-MDF$sample_id
  MDF=MDF[,!names(MDF) %in% c("sample_id")]
  X0=na.omit(MDF)
  cat("INFO: Total number of samples ", nrow(X0), "\n")
  
  dir.create(dirout, recursive = T)
  set.seed(123)
  print(unique(X0$Response))
  if(length(unique(X0$Response)) >= 5){
    cat("   Stratified cross-validation")
    folds <- vfold_cv(X0, v = 5, strata = Response)
  }else{
    cat("   Non-stratified cross-validation (too few unique response values)")
    folds <- vfold_cv(X0, v = 5)
  }
  predictions_list <- map(folds$splits, train_and_test_XGB)
  predictions <- bind_rows(predictions_list)
  
  return(predictions)
  
}

## Reading input data

norm_clade_counts_16S_files=list.files(argv$p, pattern = "norm_clade_counts_16S_\\d.tsv|norm_seqtab_16S.tsv")
norm_clade_counts_18S_files=list.files(argv$e, pattern = "norm_clade_counts_18S_\\d.tsv|norm_seqtab_18S.tsv")

if (argv$b != "") { 
  abiot=as.data.frame(t(read.delim(argv$b, header = TRUE)))
  rownames(abiot)<-gsub("^X","", rownames(abiot))  
  samples_names=as.data.frame(read.delim(argv$t, header = TRUE))
  samples_names$station_id_date <-gsub("-",".",samples_names$station_id_date)
  abiot$station_id_date=rownames(abiot)
  abiotics=merge(abiot, samples_names, by="station_id_date")
  abiotics<-abiotics[,!names(abiotics) %in% "station_id_date"]
  abiotics <- abiotics %>% select("sample_id", everything())
  ab_factors=names(abiotics)[!names(abiotics) %in% "sample_id"]
  add_sufix="Biotic"
} else {    
  abiotics=read_tsv(argv$m, show_col_types = FALSE) 
  ab_factors=names(abiotics)[8:28] 
  add_sufix="Abiotic"
}

## Pick paramters for troubleshooting
# rRNA = '16S'
# Tax_Level = 'Class'
# fc = 'Thalassiosira rotula'

for (rRNA in c("18S","16S")) { #rRNA="16S"
  if (rRNA == "16S")  tax_level_numbers<-list("Class"=3, "Order"=4,"Family"=5, "Genus"=6, "Species"=7,"ASV"=8 )
  if (rRNA == "18S")  tax_level_numbers<-list("Class"=5, "Order"=6,"Family"=7, "Genus"=8, "Species"=9,"ASV"=10 )
  
  for (Tax_Level in names(tax_level_numbers)){ #Tax_Level="Species"
    num=tax_level_numbers[[Tax_Level]]
    
    if (rRNA == "16S") DF_norm=read.delim(paste(argv$p,norm_clade_counts_16S_files[num], sep="/"), header = TRUE)
    if (rRNA == "18S") DF_norm=read.delim(paste(argv$e,norm_clade_counts_18S_files[num], sep="/"), header = TRUE)
    
    DF_norm <- as.data.frame(t(DF_norm))
    DF_norm$sample_id<-row.names(DF_norm)
    
    DF_norm$sample_id <- sub("^X", "", DF_norm$sample_id)
    
    Dout=paste(argv$o,paste("RF",RF_extra,sep="_"),paste(Tax_Level,rRNA, sep="_"),sep="/")
    
    table_pred<-vector(mode = "list", length = length(ab_factors))
    i=0
    for (fc in ab_factors) { #fc="Amphipoda"   fc="DOC"                    
      check<-na.omit(abiotics[[fc]])
      if (length(check) > minSamp ) {
        cat("*** RF with ",Tax_Level," - ",rRNA," - ",fc, " \n")
        i=i+1
        table_pred[[i]]<-RF_analysis(DF_norm,abiotics, fc, 1000, Dout, with_CV=rf_with_CV, with_opt_mtry=rf_with_opt_mtry,  type_variable="Amplicon")
        
      } else {cat("*** WARNING:", fc, "has less than ",minSamp," samples and won't be include in the analysis\n")}
    }
    
    table_pred = table_pred[lapply(table_pred,length)>0] #removing empty lists (factors that were not predicted)
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
      j=0
      for (fc in ab_factors) {
        check<-na.omit(abiotics[[fc]])
        if (length(check) > minSamp) {
          cat("*** xgb with ",Tax_Level," - ",rRNA, " - ",fc, " \n")
          j=j+1
          table_pred_xgb[[j]]<-xgb_analysis(DF_norm,abiotics, fc, Dout2)
          
        } else {cat("*** WARNING:", fc, "has less than ",minSamp," samples and won't be include in the analysis\n")}
      }
      
      table_pred_xgb= table_pred_xgb[lapply(table_pred_xgb,length)>0] #removing empty lists (factors that were not predicted)
      final_pred_table_xgb<-table_pred_xgb %>% reduce(full_join, by='sample_id') %>% arrange(sample_id)
      
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
      num=grep(paste0("^",archt, '_'),VAE_features_files)
      DF_norm=read.delim(paste(argv$a,VAE_features_files[num], sep="/"), header = TRUE)
      names(DF_norm)[1]<-"sample_id"
      DF_norm$sample_id <- sub("^X", "", DF_norm$sample_id)
      
      Dout3=paste(argv$o,"VAE", "RF",paste(archt,rRNA,sep="_"),sep="/")
      
      table_pred<-vector(mode = "list", length = length(ab_factors))
      i=0
      for (fc in ab_factors) { 
        check<-na.omit(abiotics[[fc]])
        if (length(check) > minSamp) {
          cat("*** RF with VAE ",archt, " - ",rRNA," - ", fc, " \n")
          i=i+1
          table_pred[[i]]<-RF_analysis(DF_norm,abiotics, fc, 1000, Dout3, with_CV=rf_with_CV, with_opt_mtry=rf_with_opt_mtry,  type_variable="Latent feature")
        } else {cat("*** WARNING:", fc, "has less than ",minSamp," samples and won't be include in the analysis\n")}
      }
      table_pred= table_pred[lapply(table_pred,length)>0]
      final_pred_table<-table_pred %>% reduce(full_join, by='sample_id') %>% arrange(sample_id)
      
      prefix_files=paste("VAE",archt,"RF",rRNA,add_sufix,sep="_")
      write.table(final_pred_table, file=paste(doutpred,paste(prefix_files,"Predictions.tsv",sep="_"),sep = "/"), sep = '\t', row.names = F, col.names = T )
      
      final_actual_table <- abiotics %>% filter(sample_id %in% final_pred_table$sample_id) %>% select(sub("Predicted_","",names(final_pred_table))) %>% arrange(sample_id)
      write.table(final_actual_table, file=paste(doutpred,paste(prefix_files,"Actual_values.tsv",sep="_"),sep = "/"), sep = '\t', row.names = F, col.names = T )
      
      
      if (only_RF_part == FALSE) {
        Dout2A=paste(argv$o,"VAE","XGB",paste(archt,rRNA,sep="_"),sep="/")
        table_pred_xgb<-vector(mode = "list", length = length(ab_factors))
        j=0
        for (fc in ab_factors) {
          check<-na.omit(abiotics[[fc]])
          if (length(check) > minSamp) {
            cat("*** xgb with VAE ",archt," - ",rRNA, " - ",fc, " \n")
            j=j+1
            table_pred_xgb[[j]]<-xgb_analysis(DF_norm,abiotics, fc, Dout2A)
          } else {cat("*** WARNING:", fc, "has less than ",minSamp," samples and won't be include in the analysis\n")}
          
        }
        
        table_pred_xgb= table_pred_xgb[lapply(table_pred_xgb,length)>0]
        final_pred_table_xgb<- table_pred_xgb %>% reduce(full_join, by='sample_id') %>% arrange(sample_id)
        
        prefix_files=paste("VAE",archt,"XGB",rRNA,add_sufix,sep="_")
        write.table(final_pred_table_xgb, file=paste(doutpred,paste(prefix_files,"Predictions.tsv",sep="_"),sep = "/"), sep = '\t', row.names = F, col.names = T )
        
        final_actual_table_xgb <- abiotics %>% filter(sample_id %in% final_pred_table_xgb$sample_id) %>% select(sub("Predicted_","",names(final_pred_table_xgb))) %>% arrange(sample_id)
        write.table(final_actual_table_xgb, file=paste(doutpred,paste(prefix_files,"Actual_values.tsv",sep="_"),sep = "/"), sep = '\t', row.names = F, col.names = T )
        
      }  
    }
  }
  
}    

## Get a name for the RData file
if(argv$b != ''){n
  ame = strsplit(argv$b, '/')
}else{
  name = strsplit(argv$m, '/')
}
name = name[[1]][length(name[[1]])]
name = strsplit(name, '.', fixed = TRUE)[[1]][1]

## Save the session
save.image(paste('prediction_', name ,'.RData', sep = ''))
