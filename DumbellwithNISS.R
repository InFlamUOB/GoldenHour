Part1 <- 1
Part2 <- 0


WD <- paste(path.expand("~"), "/Documents/PhD/R/Pipeline/", sep="")
setwd(WD)



#load settings and set working directory
source("SettingsPipe.R")
CF_ISS_1 <- read.csv("CF_ISS1.csv")
CF_ISS_2 <- read.csv("CF_ISS2.csv")
CF_ISS_3 <- read.csv("CF_ISS3.csv")

names(CF_ISS_1) <- c("","ID","Label","White Blood Cells [109/L]","Neutrophil [109/L]","Lymphocytes [109/L]","Monocytes [109/L]","Immature Granulocytes [109/L]","Neutrophils %","Lymphocytes %","Monocytes %","Immature Granulocytes %","ROS TO PMA (%)","ROS TO PMA (MFI)","CD62L (MedFI)","CD88 (MedFI)","CD16 (MedFI)","CD11B (MedFI)","CXCR1 (MedFI)","CXCR2 (MedFI)","CD11B after FMLF (MedFI)","Fold increase in CD11b after fMLf","CD63 (MedFI)","CD62L after FMLF (MedFI)","Fold decrease in CD62L (MedFI)","CD16BRIGHT CD62LDIM % ","CD16BRIGHT CD62LDIM [106/L]","CD16+ %","CD16 (MedFI)__1","CD14+ 16- %","CD14+ 16+ % ","CD14+ 16- [x106/L]","CD14+ 16+ [x106/L]","HLA-DR+ %","HLA-DR (MedFI)","CD14+ HLA-DRLow/Neg %","CD14+ HLA-DRLow/Neg [x106/L]","TLR4+ %","TLR4 (MedFI)","TLR2+ %","TLR2 (MedFI)","CD86+ %","CD86 (MedFI)","B Cells [x106/L]","Natural Killler Cells [x106/L]","CD56DIM Natural Killer Cells [x106/L]","CD56BRIGHT Natural Killer cells [x106/L]","Natural Killer T Cells [x106/L]","CD3+ [106/L]","CD3+ 4+ [106/L]","Terminal Helper","Naive Helper","Central Memory Helper","Effector Memory Helper","CD3+ 8+ [106/L]","Terminal Cytotoxic","Naive Cytotoxic","Central Memory Cytotoxic","Effector Memory Cytotoxic","IL-6 1LPS","IL-8 1LPS","IL-10 1LPS","TNF-A 1LPS","MCP-1 1LPS","IL-6 10LPS","IL-8 10LPS","IL-10 10LPS","TNF-A 10LPS","MCP-1 10LPS","Serum IL1.Ra","Serum IL.6","Serum IL.8","Serum IL.10","Serum G.CSF","Serum MCP.1","Serum TNF-α","Cortisol [ng/ml]","Cell Free DNA","HMGB1","ISS","NISS")
names(CF_ISS_2) <- c("","ID","Label","White Blood Cells [109/L]","Neutrophil [109/L]","Lymphocytes [109/L]","Monocytes [109/L]","Immature Granulocytes [109/L]","Neutrophils %","Lymphocytes %","Monocytes %","Immature Granulocytes %","ROS TO PMA (%)","ROS TO PMA (MFI)","CD62L (MedFI)","CD88 (MedFI)","CD16 (MedFI)","CD11B (MedFI)","CXCR1 (MedFI)","CXCR2 (MedFI)","CD11B after FMLF (MedFI)","Fold increase in CD11b after fMLf","CD63 (MedFI)","CD62L after FMLF (MedFI)","Fold decrease in CD62L (MedFI)","CD16BRIGHT CD62LDIM % ","CD16BRIGHT CD62LDIM [106/L]","CD16+ %","CD16 (MedFI)__1","CD14+ 16- %","CD14+ 16+ % ","CD14+ 16- [x106/L]","CD14+ 16+ [x106/L]","HLA-DR+ %","HLA-DR (MedFI)","CD14+ HLA-DRLow/Neg %","CD14+ HLA-DRLow/Neg [x106/L]","TLR4+ %","TLR4 (MedFI)","TLR2+ %","TLR2 (MedFI)","CD86+ %","CD86 (MedFI)","B Cells [x106/L]","Natural Killler Cells [x106/L]","CD56DIM Natural Killer Cells [x106/L]","CD56BRIGHT Natural Killer cells [x106/L]","Natural Killer T Cells [x106/L]","CD3+ [106/L]","CD3+ 4+ [106/L]","Terminal Helper","Naive Helper","Central Memory Helper","Effector Memory Helper","CD3+ 8+ [106/L]","Terminal Cytotoxic","Naive Cytotoxic","Central Memory Cytotoxic","Effector Memory Cytotoxic","IL-6 1LPS","IL-8 1LPS","IL-10 1LPS","TNF-A 1LPS","MCP-1 1LPS","IL-6 10LPS","IL-8 10LPS","IL-10 10LPS","TNF-A 10LPS","MCP-1 10LPS","Nuclear DNA [ng/ml]","Mitochondrial DNA [ng/ml]","Serum IL1.Ra","Serum IL.6","Serum IL.8","Serum IL.10","Serum G.CSF","Serum MCP.1","Serum TNF-α","Cortisol [ng/ml]","Cell Free DNA","HMGB1","ISS","NISS")
names(CF_ISS_3)<- c("","ID","Label","White Blood Cells [109/L]","Neutrophil [109/L]","Lymphocytes [109/L]","Monocytes [109/L]","Immature Granulocytes [109/L]","Neutrophils %","Lymphocytes %","Monocytes %","Immature Granulocytes %","ROS TO PMA (%)","ROS TO PMA (MFI)","CD62L (MedFI)","CD88 (MedFI)","CD16 (MedFI)","CD11B (MedFI)","CXCR1 (MedFI)","CXCR2 (MedFI)","CD11B after FMLF (MedFI)","Fold increase in CD11b after fMLf","CD63 (MedFI)","CD62L after FMLF (MedFI)","Fold decrease in CD62L (MedFI)","CD16BRIGHT CD62LDIM % ","CD16BRIGHT CD62LDIM [106/L]","CD16+ %","CD16 (MedFI)__1","CD14+ 16- %","CD14+ 16+ % ","CD14+ 16- [x106/L]","CD14+ 16+ [x106/L]","HLA-DR+ %","HLA-DR (MedFI)","CD14+ HLA-DRLow/Neg %","CD14+ HLA-DRLow/Neg [x106/L]","TLR4+ %","TLR4 (MedFI)","TLR2+ %","TLR2 (MedFI)","CD86+ %","CD86 (MedFI)","B Cells [x106/L]","Natural Killler Cells [x106/L]","CD56DIM Natural Killer Cells [x106/L]","CD56BRIGHT Natural Killer cells [x106/L]","Natural Killer T Cells [x106/L]","CD3+ [106/L]","CD3+ 4+ [106/L]","Terminal Helper","Naive Helper","Central Memory Helper","Effector Memory Helper","CD3+ 8+ [106/L]","Terminal Cytotoxic","Naive Cytotoxic","Central Memory Cytotoxic","Effector Memory Cytotoxic","IL-6 1LPS","IL-8 1LPS","IL-10 1LPS","TNF-A 1LPS","MCP-1 1LPS","IL-6 10LPS","IL-8 10LPS","IL-10 10LPS","TNF-A 10LPS","MCP-1 10LPS","Nuclear DNA [ng/ml]","Mitochondrial DNA [ng/ml]","Serum IL1.Ra","Serum IL.6","Serum IL.8","Serum IL.10","Serum G.CSF","Serum MCP.1","Serum TNF-α","Cortisol [ng/ml]","Cell Free DNA","HMGB1","ISS","NISS")



settings$number_of_cores <- detectCores()
test_name <- "First" #Make sure this is changed between time points: Second, Third
settings$seed <- 132


set.seed(settings$seed)


#########FeatureSelection

#1# Load Data Set 


s <- CF_ISS_1 #obtained from source file with correct column names
#clean data
s <- s[,3:dim(s)[2]]



#Make sure the column to classify with has Label as column name
#s$Label<-s$Group
#s$Group <- NULL

s$Label <- as.factor(s$Label)
levels(s$Label) <- c("N","Y")# has to be N and Y 




load_set <- function() {
  ###
  # load the dataset into the variable s
  
  return(list(data=s))
}


#Load functions to use for feature selection. Plus set variables
source("/Users/laura/Documents/PhD/R/Pipeline/RFANN.R") #Randm Forest and ANN functions load

settingsPipe$loops <- 100 #Number of loops. Per loop a model with features selected. Then filter accordingly.
settingsPipe$traintestdiv <- 0.75
settingsPipe$nlambda <- 100 #cross validation kambda for regularization lambda optimization
settingsPipe$filter <- 55
settingsAUC$N<-1000 #number of models produced  (both permuted (random) and real) - for all univariate (Each feature), bivariate and multivariate. 




All2 <- load_set()$data;


for (i in 1:settingsPipe$loops){
  
  smp_size = floor(settingsPipe$traintestdiv * nrow(All2))
  
  train_ind = sample(seq_len(nrow(All2)), size = smp_size)
  
  #Training set
  train = All2[train_ind, ]
  
  #Test set
  test = All2[-train_ind, ]
  
  #Creates matrices for independent and dependent variables.
  xtrain = model.matrix(Label~. -1, data = train)
  ytrain = train$Label
  xtest = model.matrix(Label~. -1, data = test)
  ytest = test$Label
  
  
  #Choose lambda value that minimize missclassification error. 
  #0.5 as elastic nets, all variables with EN are based on ElasticNets analysis. 100 lambdas sampled with 10 cross validation for each, already internalized in method
  
  alpha <- 0.5
  nlambda <- settingsPipe$nlambda
  
  
  
  Betas2EN <- RegMeGeneral(xtrain, ytrain, alpha, nlambda,xtest,ytest)
  ValuesEN <- RegMeFeatures(Betas2EN,BetasFinEN,BNumFinEN)
  BetasFinEN <- as.character(ValuesEN$BetasFinEN) #this in no characater gives problems, as factor not able to append?
  BNumFinEN <- ValuesEN$BNumFinEN
  see2EN <- RegMeFrequencies(Betas2EN,see2EN)
  
  
  # LASSO : alpha of 1 
  
  alpha <- 1
  nlambda <- settingsPipe$nlambda
  
  
  Betas2 <- RegMeGeneral(xtrain, ytrain, alpha, nlambda,xtest,ytest)
  Values <- RegMeFeatures(Betas2,BetasFin,BNumFin)
  BetasFin <- as.character(Values$BetasFin) #this in no characater gives problems, as factor not able to append?
  BNumFin <- ValuesEN$BNumFin
  see2 <- RegMeFrequencies(Betas2,see2)
  
  #Beta coefficients
  BetasTodo <- RegMeCoefficients(Betas2EN,BetasTodo)
  BetasTodoL <- RegMeCoefficients(Betas2,BetasTodoL)
  
}


#Select features LASSO and EN

#obtain in a data frame all error, betas names, number and lamda for the N models for each lasso and EN
All_info <- data.frame(BetasNames=BetasFin, BetasNum=BNumFin)
All_infoEN <- data.frame(BetasNames=BetasFinEN, BetasNum=BNumFinEN)

saveRDS(All_info, paste("Results/",test_name, "_AllinfoLASSO.rds")) 
saveRDS(All_infoEN, paste("Results/",test_name, "_AllinfoEN.rds")) 

#same as before for LASSO and EN

Final_LASSO1 <- processRegularization(see2,settingsPipe$filter)
Final_EN1 <- processRegularization(see2EN,settingsPipe$filter)

m<-count(see2, All)
m<-m[-1,]
Final_LASSO<-m[order(-m$n),]
Freqs<-Final_LASSO
Freqs$All <- factor(Freqs$All, levels = Freqs$All[order(-Freqs$n)]) #plot in a bar graph the frequencies of ocurrance of the betas and order from highest to smalles

Freqs$All<-gsub("`","",Freqs$All)

plot2<-ggplot(Freqs, aes(reorder(All,-n),n))+geom_bar(stat = "identity")+theme(axis.text=element_text(angle = 90,size=12),axis.title=element_text(size=14),strip.text.x = element_text(size = 14, colour = "black"), axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))+labs(x="Features",y="Frequency")+geom_hline(yintercept=55, linetype="dashed", color = "red")
                                                                                                                                                                                                                                                                                                       
                                                                                                                                                              


#Boxplot with Betas and its coefficients

Boxplot1 <- BetasTodo[BetasTodo$Features %in% Final_EN1$All,] #see which features appear in the filtered features (80%) and obtain the coefficients associated. 
Boxplot1["Method"]<-as.factor("EN")

Boxplot2 <- BetasTodoL[BetasTodoL$Features %in% Final_LASSO1$All,]
Boxplot2["Method"] <- as.factor("LASSO")

Fin_Boxplot <- rbind(Boxplot1,Boxplot2) #Unite both boxplots LASSO and EN 
Fin_Boxplot$Features<-gsub("`","",Fin_Boxplot$Features)

pdf(paste("Figures/", test_name, "_BoxplotAll.pdf"))
plot1<-ggplot(Fin_Boxplot,aes(Fin_Boxplot$Features,Fin_Boxplot$Coefficients))+geom_boxplot(aes(color=Method))+geom_jitter()+theme(axis.text.x = element_text(size=12, angle=90))+ggtitle("Beta coefficients EN and LASSO")+ labs(x = "Features",y="Beta Coefficients")
graphics.off()

pdf(paste("Figures/", test_name, "_BoxplotAll.pdf"),width=11,height = 10)
grid.arrange(plot2,plot1,nrow=2)
graphics.off()

grid.arrange(plot1,plot2,nrow=2)
#################################################################################################

All_Feat <- rbind(Final_LASSO1,Final_EN1)
All_Feat2 <- unique(All_Feat$All) #select the filtered betas in both. 
Betas_select <- All2[,intersect(gsub("`", "", All_Feat2), colnames(All2))]
Betas_select["Label"] <- All2$Label
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
ISS <- which(rownames(MeansSingle) %in% "ISS")
MeansSingle <- MeansSingle[-ISS,]


colnames(MeansSingle) <- c("Real", "Permuted_Data")
Total <- rbind(MeansSingle, Multivariate=data.frame(Real=MeansMultiple$Multi,Permuted_Data=MeansMultiple$MultiRand))
Total["Names"]<- rownames(Total)


plot1 <- ggplot(Total, aes(y=Names, x=Permuted_Data,xend=Real))+geom_dumbbell(color="cornflowerblue", size=1.5, point.colour.l="goldenrod1")+theme(axis.text.x = element_text(angle = 90), axis.text=element_text(size=12),axis.title=element_text(size=18),strip.text.x = element_text(size = 8, colour = "black"),axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))+labs(x="AUC values",y="Selected Features")+xlim(0.48,0.9)
saveRDS(plot1, paste("Results/",test_name, "_Finalplot.rds"))


plot2 <- plotAUCMultipleFig(multipleAUC, multipleAUCR,2)


pdf(paste("FigAUC/", test_name,"_FinalPlot.pdf"),width=10, height=4)
save5 <- grid.arrange(arrangeGrob(plot1,plot2, ncol=2), ncol=1,name=test_name)
graphics.off()

saveRDS(save5, paste("Results/",test_name, "_Finalplot.rds"))
plot(save5)
