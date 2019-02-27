
#########
#INITIALIZATION


my_packages <- c("data.table", "tidyverse", "dplyr", "devtools","car","ggpubr","glmnet","summarytools","knitr","htmltools",
                 "corrplot","factoextra","Metrics","readr","gplots","stringr","readxl","e1071","ggplot2","reshape2",
                 "ROCR","caret","gridExtra","doMC","parallel","multtest")




library(data.table)
library(tidyverse)
library(dplyr)

library(devtools)

library(car)
library(ggpubr)
library(glmnet)


library(knitr)
library(htmltools)
library(factoextra)
library(readr)
library(gplots)
library(stringr)
library(readxl)
library(e1071)

library(ggplot2)
library(reshape2)



library(ROCR)
library(gridExtra)
library(SciencesPo)
library(caret)







#multtest not available for this new version - apparently can go from one version to the other. 


#packages(my_packages) #first time
#library(my_packages) #second time

# source("https://bioconductor.org/biocLite.R")
# biocLite()
# biocLite("sva")
# BiocInstaller::biocLite(c("gmodels","caret","caretEnsemble ","doParallel","plyr","DESeq2","STRINGdb","igraph","hash","rain","pheatmap","genefilter",
#                           "RColorBrewer","ggplot2","pheatmap","RColorBrewer","ggplot2","gridExtra","cowplot","grid","gridExtra","dplyr","missForest","Amelia",
#                           "AnnotationDbi","AnnotationHub","ArrayExpress","pROC","AUC","BaySeq","biclust","bibtex","BiocGenerics","BiocParallel","biomaRt","Biostrings",
#                           "BMA","BSgenome","car","caTools","chrom","clusterProfiler","corrplot","COSINE","curl","cvTools","data.table","DiffCorr","DO.db","doMC","DOSE","e1071",
#                           "edgeR","Epi","foreach","foreign","GeneCycle","geneplotter","GenomicFeatures","GenomicRanges","ggthemes","git2r","GO.db","gplots","gProfileR","graph","Gviz","haven","h20","hda",
#                           "Hmisc","httr","igraph","igraphdata","IRanges","jsonlite","KEGG.db","KEGGREST","KernSmooth","kernlab","knitr","lattice","limma","lme4","LogicReg","loo","magrittr","mboost","mclust","mcmc",
#                           "MCMCpack","metafor","metap","mi","mice","mlogit","NbClust","NHANES","openxlsx","org.Dm.eg.db","org.hs.eb.db","org.Mm.eg.db","pamr","pcaMethods","penalized","perm","permute","pheatmap","plotly",
#                           "plot3D","pwr","psych","pvclust","qvalue","rain","randomForest","rARPACK","RCircos","RColorBrewer","Rccp","RccpArmadillo","RccpEigen","RCurl","RCytoscape","readstata13","reshape2","Rgraphviz","RJSONIO",
#                           "rmeta","robust","ROCR","rpart","rpart.plot","rPython","Rsamtools","rtracklayer","Rtsne","samplesize","shiny","shynyBS","snow","STRINGdb","stringr","TDARACNE","tidyr","tseries","tsne","VenDiagram","WGCNA","xtable","zoo"))
#  

########################

settings <- list(
  
  number_of_cores = 4,
  
  #to save model
  save_models = F,
  
  #save the datasets
  save_datasets = F,
  
  #to load from the saved models
  load_and_run_from_models = F,
  
  seed=123
)


settingsPipe <- list(
  
 loops = 10, 
  traintestdiv = 0.75,
 nlambda = 100
  
)

settingsAUC <- list(
  
  N =1000
)

if (Part1 ==1 && Part2 ==0){
ErrorsFinEN<-vector(mode="double", length=settingsPipe$loops)
BetasFinEN<-vector(mode="character", length=settingsPipe$loops)
LambdaFinEN<-vector(mode="double", length=settingsPipe$loops)
BNumFinEN<-vector(mode="double", length=settingsPipe$loops)
see2EN<-data.frame(All="All")
LauCoef1<-data.frame(Coeff="See",stringsAsFactors=FALSE)
BetasTodo<-data.frame(Features="Name",Coefficients=1)



ErrorsFin<-vector(mode="double", length=settingsPipe$loops)
BetasFin<-vector(mode="character", length=settingsPipe$loops)
LambdaFin<-vector(mode="double", length=settingsPipe$loops)
BNumFin<-vector(mode="double", length=settingsPipe$loops)
see2<-data.frame(All="All")
LauCoef1L<-data.frame(Coeff="See",stringsAsFactors=FALSE)
BetasTodoL<-data.frame(Features="Name",Coefficients=1)
}

if (Part1 ==0 && Part2 ==1){
  
  multipleAUC<-matrix(rnorm(2),1,N) 
  multipleAUCR<-matrix(rnorm(2),1,N) 
  multipleAUCNB<-matrix(rnorm(2),1,N) 
  multipleAUCNBR<-matrix(rnorm(2),1,N) 
  
  multipleROC<-matrix(as.list(rnorm(2)),1,N)  
  multipleROCR<-matrix(as.list(rnorm(2)),1,N)  
  multipleNBROC<-matrix(as.list(rnorm(2)),1,N) 
  multipleNBROCR<-matrix(as.list(rnorm(2)),1,N) 
  
  singleROC<-list()
  doubleROC<-list()
  singleROCR<-list()
  doubleROCR<-list()
  
  doublePlus<-list()
  singlePlus<-list()
  
  singleAUC<-matrix(rnorm(2),NumVar-1,N)    
  doubleAUC<-matrix(rnorm(2),(NumVar-1),N) 
  singleAUCR<-matrix(rnorm(2),NumVar-1,N)    
  doubleAUCR<-matrix(rnorm(2),(NumVar-1),N)
  
  doubleAUCSVMR<-matrix(rnorm(2),(NumVar-1),N) 
  doubleAUCSVM<-matrix(rnorm(2),(NumVar-1),N) 
  doubleAUCSVMCrossR<-matrix(rnorm(2),(NumVar-1),N) 
  doubleAUCSVMCross<-matrix(rnorm(2),(NumVar-1),N) 
  doubleAUCRFCross<-matrix(rnorm(2),(NumVar-1),N) 
  doubleAUCRFCrossR<-matrix(rnorm(2),(NumVar-1),N) 
  doubleAUCNBR<-matrix(rnorm(2),(NumVar-1),N) 
  doubleAUCNB<-matrix(rnorm(2),(NumVar-1),N) 
  
  MatsingleROC<-matrix(as.list(rnorm(2)),NumVar-1,N)  
  MatsingleROCR<-matrix(as.list(rnorm(2)),NumVar-1,N)  
  MatdoubleROC<-matrix(as.list(rnorm(2)),(NumVar-1),N) 
  MatdoubleROCR<-matrix(as.list(rnorm(2)),(NumVar-1),N) 
  
  MatsinglePlus<-matrix(as.list(rnorm(2)),NumVar-1,N) 
  MatdoublePlus<-matrix(as.list(rnorm(2)),(NumVar-1),N)
  MatsinglePlusR<-matrix(as.list(rnorm(2)),NumVar-1,N) 
  MatdoublePlusR<-matrix(as.list(rnorm(2)),(NumVar-1),N)
  multiplePlus<-matrix(as.list(rnorm(2)),1,N)
  multiplePlusR<-matrix(as.list(rnorm(2)),1,N)
  
  
  
  
  
}

