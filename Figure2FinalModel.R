#File needed to create Figure 2 of paper

Part1 <- 1
Part2 <- 0


WD <- paste(path.expand("~"), "/Documents/PhD/R/Pipeline/", sep="")
setwd(WD)


#load settings and set working directory
source("SettingsPipe.R")

CF_ISS_1 <- read.csv("CF_ISS1.csv")
CF_ISS_2 <- read.csv("CF_ISS2.csv")
CF_ISS_3 <- read.csv("CF_ISS3.csv")

names(CF_ISS_1) <- c("","ID","Label","White Blood Cells [109/L]","Neutrophil [109/L]","Lymphocytes [109/L]","Monocytes [109/L]","Immature Granulocytes [109/L]","Neutrophils %","Lymphocytes %","Monocytes %","Immature Granulocytes %","ROS TO PMA (%)","ROS TO PMA (MFI)","CD62L (MedFI)","CD88 (MedFI)","CD16 (MedFI)","CD11B (MedFI)","CXCR1 (MedFI)","CXCR2 (MedFI)","CD11B after FMLF (MedFI)","Fold increase in CD11b after fMLf","CD63 (MedFI)","CD62L after FMLF (MedFI)","Fold decrease in CD62L (MedFI)","CD16BRIGHT CD62LDIM % ","CD16BRIGHT CD62LDIM [106/L]","CD16+ %","CD16 (MedFI).1","CD14+ 16- %","CD14+ 16+ % ","CD14+ 16- [x106/L]","CD14+ 16+ [x106/L]","HLA-DR+ %","HLA-DR (MedFI)","CD14+ HLA-DRLow/Neg %","CD14+ HLA-DRLow/Neg [x106/L]","TLR4+ %","TLR4 (MedFI)","TLR2+ %","TLR2 (MedFI)","CD86+ %","CD86 (MedFI)","B Cells [x106/L]","Natural Killler Cells [x106/L]","CD56DIM Natural Killer Cells [x106/L]","CD56BRIGHT Natural Killer cells [x106/L]","Natural Killer T Cells [x106/L]","CD3+ [106/L]","CD3+ 4+ [106/L]","Terminal Helper T Cells","Naive Helper T Cells","Central Memory Helper T Cells","Effector Memory Helper T Cells","CD3+ 8+ [106/L]","Terminal Cytotoxic T Cells","Naive Cytotoxic T Cells","Central Memory Cytotoxic T Cells","Effector Memory Cytotoxic T Cells","IL-6 1LPS","IL-8 1LPS","IL-10 1LPS","TNF-A 1LPS","MCP-1 1LPS","IL-6 10LPS","IL-8 10LPS","IL-10 10LPS","TNF-A 10LPS","MCP-1 10LPS","Serum IL1.Ra","Serum IL.6","Serum IL.8","Serum IL.10","Serum G.CSF","Serum MCP.1","Serum TNF-A","Cortisol [ng/ml]","Cell Free DNA","HMGB1","ISS","NISS")
names(CF_ISS_2) <- c("","ID","Label","White Blood Cells [109/L]","Neutrophil [109/L]","Lymphocytes [109/L]","Monocytes [109/L]","Immature Granulocytes [109/L]","Neutrophils %","Lymphocytes %","Monocytes %","Immature Granulocytes %","ROS TO PMA (%)","ROS TO PMA (MFI)","CD62L (MedFI)","CD88 (MedFI)","CD16 (MedFI)","CD11B (MedFI)","CXCR1 (MedFI)","CXCR2 (MedFI)","CD11B after FMLF (MedFI)","Fold increase in CD11b after fMLf","CD63 (MedFI)","CD62L after FMLF (MedFI)","Fold decrease in CD62L (MedFI)","CD16BRIGHT CD62LDIM % ","CD16BRIGHT CD62LDIM [106/L]","CD16+ %","CD16 (MedFI).1","CD14+ 16- %","CD14+ 16+ % ","CD14+ 16- [x106/L]","CD14+ 16+ [x106/L]","HLA-DR+ %","HLA-DR (MedFI)","CD14+ HLA-DRLow/Neg %","CD14+ HLA-DRLow/Neg [x106/L]","TLR4+ %","TLR4 (MedFI)","TLR2+ %","TLR2 (MedFI)","CD86+ %","CD86 (MedFI)","B Cells [x106/L]","Natural Killler Cells [x106/L]","CD56DIM Natural Killer Cells [x106/L]","CD56BRIGHT Natural Killer cells [x106/L]","Natural Killer T Cells [x106/L]","CD3+ [106/L]","CD3+ 4+ [106/L]","Terminal Helper T Cells","Naive Helper T Cells","Central Memory Helper T Cells","Effector Memory Helper T Cells","CD3+ 8+ [106/L]","Terminal Cytotoxic T Cells","Naive Cytotoxic T Cells","Central Memory Cytotoxic T Cells","Effector Memory Cytotoxic T Cells","IL-6 1LPS","IL-8 1LPS","IL-10 1LPS","TNF-A 1LPS","MCP-1 1LPS","IL-6 10LPS","IL-8 10LPS","IL-10 10LPS","TNF-A 10LPS","MCP-1 10LPS","Nuclear DNA [ng/ml]","Mitochondrial DNA [ng/ml]","Serum IL1.Ra","Serum IL.6","Serum IL.8","Serum IL.10","Serum G.CSF","Serum MCP.1","Serum TNF-A","Cortisol [ng/ml]","Cell Free DNA","HMGB1","ISS","NISS")
names(CF_ISS_3)<- c("","ID","Label","White Blood Cells [109/L]","Neutrophil [109/L]","Lymphocytes [109/L]","Monocytes [109/L]","Immature Granulocytes [109/L]","Neutrophils %","Lymphocytes %","Monocytes %","Immature Granulocytes %","ROS TO PMA (%)","ROS TO PMA (MFI)","CD62L (MedFI)","CD88 (MedFI)","CD16 (MedFI)","CD11B (MedFI)","CXCR1 (MedFI)","CXCR2 (MedFI)","CD11B after FMLF (MedFI)","Fold increase in CD11b after fMLf","CD63 (MedFI)","CD62L after FMLF (MedFI)","Fold decrease in CD62L (MedFI)","CD16BRIGHT CD62LDIM % ","CD16BRIGHT CD62LDIM [106/L]","CD16+ %","CD16 (MedFI).1","CD14+ 16- %","CD14+ 16+ % ","CD14+ 16- [x106/L]","CD14+ 16+ [x106/L]","HLA-DR+ %","HLA-DR (MedFI)","CD14+ HLA-DRLow/Neg %","CD14+ HLA-DRLow/Neg [x106/L]","TLR4+ %","TLR4 (MedFI)","TLR2+ %","TLR2 (MedFI)","CD86+ %","CD86 (MedFI)","B Cells [x106/L]","Natural Killler Cells [x106/L]","CD56DIM Natural Killer Cells [x106/L]","CD56BRIGHT Natural Killer cells [x106/L]","Natural Killer T Cells [x106/L]","CD3+ [106/L]","CD3+ 4+ [106/L]","Terminal Helper T Cells","Naive Helper T Cells","Central Memory Helper T Cells","Effector Memory Helper T Cells","CD3+ 8+ [106/L]","Terminal Cytotoxic T Cells","Naive Cytotoxic T Cells","Central Memory Cytotoxic T Cells","Effector Memory Cytotoxic T Cells","IL-6 1LPS","IL-8 1LPS","IL-10 1LPS","TNF-A 1LPS","MCP-1 1LPS","IL-6 10LPS","IL-8 10LPS","IL-10 10LPS","TNF-A 10LPS","MCP-1 10LPS","Nuclear DNA [ng/ml]","Mitochondrial DNA [ng/ml]","Serum IL1.Ra","Serum IL.6","Serum IL.8","Serum IL.10","Serum G.CSF","Serum MCP.1","Serum TNF-A","Cortisol [ng/ml]","Cell Free DNA","HMGB1","ISS","NISS")




settings$number_of_cores <- detectCores()
test_name <- "First" #Make sure this is changed between time points: Second, Third
settings$seed <- 132


set.seed(settings$seed)


#########FeatureSelection

#1# Load Data Set 

Betas_select <- data.frame(CF_ISS_1$`Monocytes %`,CF_ISS_2$`CD11B (MedFI)`,CF_ISS_3$`CD62L after FMLF (MedFI)`,CF_ISS_1$NISS)

Betas_select["Label"]<- CF_ISS_1$Label

names(Betas_select) <- c("First Time Monocytes %","Second Time CD11B (MedFI)","Third Time CD62L after FMLF (MedFI)","NISS","Label")

Betas_select2 <- Betas_select

Part1 <- 0
Part2 <- 1

NumVar <- length(Betas_select)
names <-colnames(Betas_select)
names1 <- as.character(names)
names1 <- strsplit(names1,", ")
names <- as.data.frame(names1)

N <- settingsAUC$N

source("SettingsPipe.R")

for (j in 1:N){ #N different measurements of AUC values, mean done at the end. 
  
  smp_size <- floor(0.65 * nrow(Betas_select2))
  #set.seed(907)
  train_ind <- sample(seq_len(nrow(Betas_select2)), size = smp_size)
  
  # Training set
  train <- Betas_select2[train_ind, ]
  
  # Test set
  test <- Betas_select2[-train_ind, ]
  
  xtrain <- model.matrix(Label~. -1, data = train)
  
  ytrain <- train$Label
  xtest <- model.matrix(Label~. -1, data = test)
  ytest <- test$Label
  
  xtest <- data.frame(xtest)
  xtrain <- data.frame(xtrain)
  names(xtrain) <- names(train)[1:dim(train)[2]-1]
  names(xtest) <- names(test)[1:dim(train)[2]-1]
  
  source('FunctionsAUC.R')
  
  
  multipleAUCNB[1,j] <- multipleAUCfunNB(xtrain, ytrain,xtest,ytest)
  multipleAUC[1,j] <- multipleAUCfun(xtrain, ytrain,xtest,ytest)
  
  multipleNBROC[[j]] <- multipleROCfunNB(xtrain, ytrain,xtest,ytest)
  
  multipleROC[[j]] <- multipleROCfun(xtrain, ytrain,xtest,ytest)
  multiplePlus[1,j] <- multiplePlusfun(xtrain, ytrain,xtest,ytest)
  
  #separate
  
  for (s in (1:(NumVar-1))){
    s <- as.numeric(s)
    singleAUC[s,j] <- singleAUCfun(xtrain, ytrain,xtest,ytest,s) 
    singleROC[[s]] <- singleROCfun(xtrain, ytrain,xtest,ytest,s) 
    MatsinglePlus[s,j] <- singlePlusfun(xtrain, ytrain,xtest,ytest,s) 
  }
  MatsingleROC[,j] <- matrix(singleROC)
  
  
  #Null Hypothesis #### 
  
  
  # Training set
  
  train$Label <- sample(train$Label)
  test$Label <- sample(test$Label)
  #Permuted data, will make sure that are models are really valid as randomizing the label should yield around 0.5 AUC values. The same testing and training arrangements for the real per model are used. 
  
  # Test set
  
  
  xtrain <- model.matrix(Label~. -1, data = train)
  
  ytrain <- train$Label
  xtest <- model.matrix(Label~. -1, data = test)
  ytest <- test$Label
  
  xtest <- data.frame(xtest)
  xtrain <- data.frame(xtrain)
  names(xtrain) <- names(train)[1:dim(train)[2]-1]
  names(xtest) <- names(test)[1:dim(train)[2]-1]
  source('FunctionsAUC.R')
  
  multipleAUCNBR[1,j] <- multipleAUCfunNB(xtrain, ytrain,xtest,ytest)
  multipleAUCR[1,j] <- multipleAUCfun(xtrain, ytrain,xtest,ytest)
  multiplePlusR[1,j] <- multiplePlusfun(xtrain, ytrain,xtest,ytest)
  
  multipleNBROCR[[j]] <- multipleROCfunNB(xtrain, ytrain,xtest,ytest)
  multipleROCR[[j]] <- multipleROCfun(xtrain, ytrain,xtest,ytest)
  
  
  for (s in (1:(NumVar-1))){
    s <- as.numeric(s)
    singleAUCR[s,j] <- singleAUCfun(xtrain, ytrain,xtest,ytest,s) 
    singleROCR[[s]] <- singleROCfun(xtrain, ytrain,xtest,ytest,s) 
    MatsinglePlusR[s,j] <- singlePlusfun(xtrain, ytrain,xtest,ytest,s) 
  }
  MatsingleROCR[,j] <- matrix(singleROCR)
  
}


trial <- NULL
source("FunctionsAUC.R")

#Means 

singleAUC <- as.data.frame(singleAUC)
singleAUC <- mutate(singleAUC, Means=rowMeans(singleAUC))
row.names(singleAUC) <- names(xtrain)

singleAUCR <- as.data.frame(singleAUCR)
singleAUCR <- mutate(singleAUCR, Means=rowMeans(singleAUCR))
row.names(singleAUCR) <- names(xtrain)

multipleAUCNBR <- as.data.frame(multipleAUCNBR)
multipleAUCNBR["Means"] <- rowMeans(multipleAUCNBR)
multipleAUCR <- as.data.frame(multipleAUCR)
multipleAUCR["Means"] <- rowMeans(as.data.frame(multipleAUCR))

multipleAUCNB <- as.data.frame(multipleAUCNB)
multipleAUCNB["Means"] <- rowMeans(as.data.frame(multipleAUCNB))
multipleAUC <- as.data.frame(multipleAUC)
multipleAUC["Means"] <- rowMeans(as.data.frame(multipleAUC))

Final <- data.frame(MultiNB=t(multipleAUCNB),MultiNBRand=t(multipleAUCNBR),Multi=t(multipleAUC),MultiRand=t(multipleAUCR) )
FinalMeans <- data.frame(MultiNB=multipleAUCNB$Means,MultiNBRand=multipleAUCNBR$Means,Multi=multipleAUC$Means,MultiRand=multipleAUCR$Means )

write.csv(FinalMeans, paste("Results/",test_name, "_FinalMeansAUC.csv"))

MonoSingle <- plotAUCSingle(singleAUC, singleAUCR,NumVar,test_name)

n <- length(plot)

# ###Accuracy, Precision etc
#
plot3 <-list()
plot4 <-list()
Data <-c("Accuracy","Sensitivity","Specificity","Precision")

MonoParameters <- plotParamSingle(MatsinglePlus, MatsinglePlusR,NumVar,Data,N,test_name)


names <-names(Betas_select2)
names1 <-as.character(names)
names1 <-strsplit(names1,", ")
names1 <-as.data.frame(names1)


MonoMultipleNB <- plotAUCMultiple(multipleAUCNB, multipleAUCNBR,1,test_name)
MonoMultiple <- plotAUCMultiple(multipleAUC, multipleAUCR,2,test_name)

Final <-data.frame(MultiNB=t(multipleAUCNB),MultiNBRand=t(multipleAUCNBR),Multi=t(multipleAUC),MultiRand=t(multipleAUCR) )
FinalMeans <-data.frame(MultiNB=multipleAUCNB$Means,MultiNBRand=multipleAUCNBR$Means,Multi=multipleAUC$Means,MultiRand=multipleAUCR$Means )

Total <- data.frame(Single=singleAUC$Means,SingleRandom=singleAUCR$Means)
m <-length(names)
rownames(Total) <- colnames(Betas_select)[1:(length(names)-1)]


####multiplePlusfun final accuracy for multiple, all together
MonoParamMult <- multipleParam(multiplePlus,multiplePlusR,N,test_name)

write.csv(FinalMeans,paste("Results/",test_name, "_FinalMeansMultiple.csv"))
write.csv(Total,paste("Results/",test_name, "FinalMeanssingle.csv"))
write.csv(Betas_select2,paste("Results/",test_name, "Betas_select.csv"))


###########

#DUMBELL


MeansMultiple <- FinalMeans
MeansSingle <- Total


colnames(MeansSingle) <- c("Real", "Permuted_Data")
Total <- rbind(MeansSingle, Multivariate=data.frame(Real=MeansMultiple$Multi,Permuted_Data=MeansMultiple$MultiRand))
Total["Names"]<- rownames(Total)

plot1 <- ggplot(Total, aes(y=Names, x=Permuted_Data,xend=Real))+geom_dumbbell(color="cornflowerblue", size=1.5, point.colour.l="goldenrod1")+theme(axis.text.x = element_text(angle = 90), axis.text=element_text(size=12),axis.title=element_text(size=18),strip.text.x = element_text(size = 8, colour = "black"),axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))+labs(x="AUC values",y="Selected Features")+xlim(0.48,0.9)
saveRDS(plot1, paste("Results/",test_name, "_Finalplot.rds"))

plot2 <- plotAUCMultipleFig(multipleAUC, multipleAUCR,2)


save5 <- grid.arrange(arrangeGrob(plot1,plot2, ncol=2), ncol=1,name=test_name)

saveRDS(save5, paste("Results/","TotalFinalplot.rds"))



#######


# get the table with p-values and standard deviation

dd<- data.frame(Random=as.numeric(multipleAUCR[1,]),Normal=as.numeric(multipleAUC[1,]))
ddsd <- format(round(sd(dd$Random), 3), nsmall = 3)
ddsd2 <- format(round(sd(dd$Normal), 3), nsmall = 3)

count2 <- 0
count <- list()
gg <- 0 
ss <- 0

for (i in 1:length(as.numeric(multipleAUCR[1,]))){
  ss <- ss+1 
  for (j in 1:length(as.numeric(multipleAUCR[1,]))){
  gg <- gg+1
   if (dd$Random[i] > dd$Normal[j]){
    count2 <- count2 + 1
    
   }
  }
  count[i] <- count2
  count2 <- 0
}


s <- data.frame(table(as.numeric(count)))
pval <- length(s$Var1)/length(as.numeric(multipleAUCR[1,]))
pval <- format(round(pval, 3), nsmall = 3)

MeansMultiple$Multi <- format(round(MeansMultiple$Multi, 3), nsmall = 3)
MeansMultiple$MultiRand <- format(round(MeansMultiple$MultiRand, 3), nsmall = 3)

FinalVal <- paste(MeansMultiple$Multi,"+/-",ddsd2)
FinalValRand <- paste(MeansMultiple$MultiRand,"+/-",ddsd)
FinalTable<- data.frame(FinalVal,FinalValRand,pval)
rownames(FinalTable) <- c("Multivariate")
names(FinalTable) <- c("  Model AUC Value  ","  Permuted AUC Value  ","  p value  ")


mytheme <- gridExtra::ttheme_default(
  core = list(fg_params=list(cex = 1.5)),
  colhead = list(fg_params=list(cex = 1.5)),
  rowhead = list(fg_params=list(cex = 1.5)))

tbl <- tableGrob(FinalTable, theme=mytheme)
pdf(paste("Figures/", "Total_Figure2.pdf"),width=10, height=5)
grid.arrange(save5, tbl,
             nrow=2,
             as.table=TRUE,
             heights=c(2,1))
graphics.off()


###Violin Plot


#open saved files so this has to follow from Main Pipe directly. 

library("tibble")

All1 <- read.csv("CF_ISS_normal1.csv")
names(All1) <- c("","MODS","ID","Age","Gender","GCS","ISS","NISS","Base Excess","Lactate","Mech_Injury","Injury_State","Deceased","Hospital Free Days","ICU Free Days","White Blood Cells [109/L]","Neutrophil [109/L]","Lymphocytes [109/L]","Monocytes [109/L]","Immature Granulocytes [109/L]","Neutrophils %","Lymphocytes %","Monocytes %","Immature Granulocytes %","ROS TO PMA (%)","ROS TO PMA (MFI)","CD62L (MedFI)","CD88 (MedFI)","CD16 (MedFI)","CD11B (MedFI)","CXCR1 (MedFI)","CXCR2 (MedFI)","CD11B after FMLF (MedFI)","Fold increase in CD11b after fMLf","CD63 (MedFI)","CD62L after FMLF (MedFI)","Fold decrease in CD62L (MedFI)","CD16BRIGHT CD62LDIM % ","CD16BRIGHT CD62LDIM [106/L]","CD16+ %","CD16 (MedFI).1","CD14+ 16- %","CD14+ 16+ % ","CD14+ 16- [x106/L]","CD14+ 16+ [x106/L]","HLA-DR+ %","HLA-DR (MedFI)","CD14+ HLA-DRLow/Neg %","CD14+ HLA-DRLow/Neg [x106/L]","TLR4+ %","TLR4 (MedFI)","TLR2+ %","TLR2 (MedFI)","CD86+ %","CD86 (MedFI)","B Cells [x106/L]","Natural Killler Cells [x106/L]","CD56DIM Natural Killer Cells [x106/L]","CD56BRIGHT Natural Killer cells [x106/L]","Natural Killer T Cells [x106/L]","CD3+ [106/L]","CD3+ 4+ [106/L]","Terminal Helper T Cells","Naive Helper T Cells","Central Memory Helper T Cells","Effector Memory Helper T Cells","CD3+ 8+ [106/L]","Terminal Cytotoxic T Cells","Naive Cytotoxic T Cells","Central Memory Cytotoxic T Cells","Effector Memory Cytotoxic T Cells","IL-6 1LPS","IL-8 1LPS","IL-10 1LPS","TNF-A 1LPS","MCP-1 1LPS","IL-6 10LPS","IL-8 10LPS","IL-10 10LPS","TNF-A 10LPS","MCP-1 10LPS","Nuclear DNA [ng/ml]","Mitochondrial DNA [ng/ml]","Serum IL1.Ra","Serum IL.6","Serum IL.8","Serum IL.10","Serum G.CSF","Serum MCP.1","Serum TNF-α","Cortisol [ng/ml]","Cell Free DNA","NISS","HMGB1")

All2 <- read.csv("CF_ISS_normal2.csv")                              
names(All2) <- c("","ID","White Blood Cells [109/L]","Neutrophil [109/L]","Lymphocytes [109/L]","Monocytes [109/L]","Immature Granulocytes [109/L]","Neutrophils %","Lymphocytes %","Monocytes %","Immature Granulocytes %","ROS TO PMA (%)","ROS TO PMA (MFI)","CD62L (MedFI)","CD88 (MedFI)","CD16 (MedFI)","CD11B (MedFI)","CXCR1 (MedFI)","CXCR2 (MedFI)","CD11B after FMLF (MedFI)","Fold increase in CD11b after fMLf","CD63 (MedFI)","CD62L after FMLF (MedFI)","Fold decrease in CD62L (MedFI)","CD16BRIGHT CD62LDIM % ","CD16BRIGHT CD62LDIM [106/L]","CD16+ %","CD16 (MedFI).1","CD14+ 16- %","CD14+ 16+ % ","CD14+ 16- [x106/L]","CD14+ 16+ [x106/L]","HLA-DR+ %","HLA-DR (MedFI)","CD14+ HLA-DRLow/Neg %","CD14+ HLA-DRLow/Neg [x106/L]","TLR4+ %","TLR4 (MedFI)","TLR2+ %","TLR2 (MedFI)","CD86+ %","CD86 (MedFI)","B Cells [x106/L]","Natural Killler Cells [x106/L]","CD56DIM Natural Killer Cells [x106/L]","CD56BRIGHT Natural Killer cells [x106/L]","Natural Killer T Cells [x106/L]","CD3+ [106/L]","CD3+ 4+ [106/L]","Terminal Helper T Cells","Naive Helper T Cells","Central Memory Helper T Cells","Effector Memory Helper T Cells","CD3+ 8+ [106/L]","Terminal Cytotoxic T Cells","Naive Cytotoxic T Cells","Central Memory Cytotoxic T Cells","Effector Memory Cytotoxic T Cells","IL-6 1LPS","IL-8 1LPS","IL-10 1LPS","TNF-A 1LPS","MCP-1 1LPS","IL-6 10LPS","IL-8 10LPS","IL-10 10LPS","TNF-A 10LPS","MCP-1 10LPS","Nuclear DNA [ng/ml]","Mitochondrial DNA [ng/ml]","Serum IL1.Ra","Serum IL.6","Serum IL.8","Serum IL.10","Serum G.CSF","Serum MCP.1","Serum TNF-α","Cortisol [ng/ml]","Cell Free DNA","MODS","ISS","NISS","HMGB1") 

All3 <- read.csv("CF_ISS_normal3.csv")  
names(All3) <- c("","ID","White Blood Cells [109/L]","Neutrophil [109/L]","Lymphocytes [109/L]","Monocytes [109/L]","Immature Granulocytes [109/L]","Neutrophils %","Lymphocytes %","Monocytes %","Immature Granulocytes %","ROS TO PMA (%)","ROS TO PMA (MFI)","CD62L (MedFI)","CD88 (MedFI)","CD16 (MedFI)","CD11B (MedFI)","CXCR1 (MedFI)","CXCR2 (MedFI)","CD11B after FMLF (MedFI)","Fold increase in CD11b after fMLf","CD63 (MedFI)","CD62L after FMLF (MedFI)","Fold decrease in CD62L (MedFI)","CD16BRIGHT CD62LDIM % ","CD16BRIGHT CD62LDIM [106/L]","CD16+ %","CD16 (MedFI).1","CD14+ 16- %","CD14+ 16+ % ","CD14+ 16- [x106/L]","CD14+ 16+ [x106/L]","HLA-DR+ %","HLA-DR (MedFI)","CD14+ HLA-DRLow/Neg %","CD14+ HLA-DRLow/Neg [x106/L]","TLR4+ %","TLR4 (MedFI)","TLR2+ %","TLR2 (MedFI)","CD86+ %","CD86 (MedFI)","B Cells [x106/L]","Natural Killler Cells [x106/L]","CD56DIM Natural Killer Cells [x106/L]","CD56BRIGHT Natural Killer cells [x106/L]","Natural Killer T Cells [x106/L]","CD3+ [106/L]","CD3+ 4+ [106/L]","Terminal Helper T Cells","Naive Helper T Cells","Central Memory Helper T Cells","Effector Memory Helper T Cells","CD3+ 8+ [106/L]","Terminal Cytotoxic T Cells","Naive Cytotoxic T Cells","Central Memory Cytotoxic T Cells","Effector Memory Cytotoxic T Cells","IL-6 1LPS","IL-8 1LPS","IL-10 1LPS","TNF-A 1LPS","MCP-1 1LPS","IL-6 10LPS","IL-8 10LPS","IL-10 10LPS","TNF-A 10LPS","MCP-1 10LPS","Nuclear DNA [ng/ml]","Mitochondrial DNA [ng/ml]","Serum IL1.Ra","Serum IL.6","Serum IL.8","Serum IL.10","Serum G.CSF","Serum MCP.1","Serum TNF-α","Cortisol [ng/ml]","Cell Free DNA","MODS","ISS","NISS","HMGB1")


#These are the selected features, and so the figures we are generating
Names_to_plot2 <- c("Monocytes %","CD11B (MedFI)","CD62L after FMLF (MedFI)","MODS","ID")


Betas_plot_1h <- All1[,colnames(All1[,intersect(gsub("`", "", Names_to_plot2), colnames(All1))])]
Betas_plot_2h <- All2[,colnames(All2[,intersect(gsub("`", "",  Names_to_plot2), colnames(All2))])]
Betas_plot_3h <- All3[,colnames(All3[,intersect(gsub("`", "",  Names_to_plot2), colnames(All3))])]

#Before merging add a column with time points
Betas_plot_1h["Time"] <- as.data.frame(as.factor(c(rep(1,dim(Betas_plot_1h)[1]))))
Betas_plot_2h["Time"] <- as.data.frame(as.factor(c(rep(2,dim(Betas_plot_2h)[1]))))
Betas_plot_3h["Time"] <- as.data.frame(as.factor(c(rep(3,dim(Betas_plot_3h)[1]))))

All <- list(Betas_plot_1h,Betas_plot_2h,Betas_plot_3h)
Betas_plot <- do.call("rbind",All)


plot <- list()
plot2 <- list()
plot3 <- list()


names(Betas_plot) <-  c("MonocytesPercentage","CD11BMedFI","CD62LafterFMLFMedFI","MODS","ID","Time")
nm <- names(Betas_plot)
levels(Betas_plot$MODS) <- c("No","Yes")
levels(Betas_plot$Time) <- c("T<1h","4<T<12h","48<T<72h")
names <- c("MODS vs non MODS for Monocytes %","MODS vs non MODS for CD11B (MedFI)","MODS vs non MODS for CD62L after FMLF (MedFI)")
axis <- c("Monocytes %","CD11B (MedFI)","CD62L after FMLF (MedFI)")

counts <- 0
for (i in nm){
  counts <- counts +1
  dodge <- position_dodge(width = 0.9)
  print(i)
  p1 <- ggplot(Betas_plot,aes_string(x="MODS",y=i,fill="Time")) + geom_violin(position=dodge)+geom_boxplot(position=dodge, width=0.2)+labs(y=paste(axis[counts]),x="MODS")+geom_hline(yintercept = 0, linetype = 2)+labs(x="MODS")
  #print(Betas_plot,aes_string(x="MODS",y=i[i],fill="TimePoint")) + geom_violin(position=dodge)+geom_boxplot(position=dodge, width=0.2)+labs(title= paste("Count of", R[i]),y=paste("Real",R[i],"count"),x="Development of MODS")+geom_hline(yintercept = 0, linetype = 2)
  p2 <- ggplot(Betas_plot,aes_string(x="Time", i,group="ID",label="ID"))+geom_point()+geom_line()+facet_grid(. ~ MODS)+geom_text(check_overlap = TRUE,vjust = 0, nudge_x = 0.2)+scale_x_discrete(labels=c("1" = "T<=1h", "2" = "4<T<12","3" = "48<T<72"))+labs(y=paste(axis[counts]))+
    labs(x="TimePoints")
  p <- ggboxplot(Betas_plot,x="Time", y=i,color = "MODS", palette = "jco")+labs(title= paste(names[counts]),y=paste(axis[counts]),x=("Time"))
  p3 <- p+geom_hline(yintercept = 0, linetype = 2)+stat_compare_means(aes(group = MODS), label = "p.format",method = "t.test") 
  
  plot[[i]] <- p1
  #print(plot[[i]])
  plot2[[i]] <- p2
  #print(plot2[[i]])
  plot3[[i]] <- p3
  #print(plot3[[i]])
}

pdf(paste("Figures/","ViolinFigure3.pdf"),width=13, height=8)
grid.arrange(plot3[[1]],plot[[1]],plot3[[2]],plot[[2]],plot3[[3]],plot[[3]], nrow=3, ncol=2)
graphics.off()







