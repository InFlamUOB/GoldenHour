
DumbellwithoutNISSFun <- function(test_nameNISS,s,i){
  

Part1 <- 1
Part2 <- 0
source("SettingsPipe.R")


s <- s[,3:dim(s)[2]]
d <- which(names(s) %in% c("ISS","NISS"))
s <- s[,-d]
#Here take away ISS and NISS columns


#Make sure the column to classify with has Label as column name
#s$Label<-s$Group
#s$Group <- NULL

s$Label <- as.factor(s$Label)
levels(s$Label) <- c("N","Y")# has to be N and Y 
test_name <- test_nameNISS




#Load functions to use for feature selection. Plus set variables
source("RFANN.R") #Randm Forest and ANN functions load

settingsPipe$loops <- 100 #Number of loops. Per loop a model with features selected. Then filter accordingly.
settingsPipe$traintestdiv <- 0.75
settingsPipe$nlambda <- 100 #cross validation kambda for regularization lambda optimization
filter <- c(75,75,75)
settingsPipe$filter <- filter[i] #with NISS set it higher
settingsAUC$N<-1000 #number of models produced  (both permuted (random) and real) - for all univariate (Each feature), bivariate and multivariate. 




All2 <- s


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

saveRDS(see2, paste("Results/",test_name, "_AllinfoLASSO.rds")) 
saveRDS(see2EN, paste("Results/",test_name, "_AllinfoEN.rds")) 

#same as before for LASSO and EN

Final_LASSO1 <- processRegularization(see2,settingsPipe$filter)
Final_EN1 <- processRegularization(see2EN,settingsPipe$filter)

#figure plot with histogram (This information is included inside the function processRegularization)
m<-count(see2, All)
m<-m[-1,]
Final_LASSO<-m[order(-m$n),]
Freqs<-Final_LASSO
Freqs$All <- factor(Freqs$All, levels = Freqs$All[order(-Freqs$n)]) #plot in a bar graph the frequencies of ocurrance of the betas and order from highest to smalles

Freqs$All<-gsub("`","",Freqs$All)
plot2<-ggplot(Freqs, aes(reorder(All,-n),n))+geom_bar(stat = "identity")+theme(axis.text=element_text(angle = 90,size=12),axis.title=element_text(size=14),strip.text.x = element_text(size = 14, colour = "black"), axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))+labs(x="Features",y="Frequency")+geom_hline(yintercept=settingsPipe$filter, linetype="dashed", color = "red")+ggtitle("Frequency of Features in LASSO")









#Boxplot with Betas and its coefficients

Boxplot1 <- BetasTodo[BetasTodo$Features %in% Final_EN1$All,] #see which features appear in the filtered features (80%) and obtain the coefficients associated. 
Boxplot1["Method"]<-as.factor("EN")

Boxplot2 <- BetasTodoL[BetasTodoL$Features %in% Final_LASSO1$All,]
Boxplot2["Method"] <- as.factor("LASSO")



Fin_Boxplot <- rbind(Boxplot1,Boxplot2) #Unite both boxplots LASSO and EN 

Fin_Boxplot$Features<-gsub("`","",Fin_Boxplot$Features)
plot1<- ggplot(Fin_Boxplot,aes(Fin_Boxplot$Features,Fin_Boxplot$Coefficients))+geom_boxplot(aes(color=Method))+geom_jitter()+theme(axis.text.x = element_text(size=12, angle=90))+ggtitle("Beta coefficients Elastic Net and LASSO")+labs(x = "Features",y="Beta Coefficients")

pdf(paste("Figures/",test_name,"_BoxplotAll.pdf"),width=9,height=6)
grid.arrange(plot1)
graphics.off()

pdf(paste("Figures/",test_name,"_Figure4.pdf"),width=14,height=11)
grid.arrange(plot2,plot1,nrow=2)
graphics.off()


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

source("SettingsPipe.R",local=TRUE)





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


source("FunctionsAUC.R",local=TRUE)


MonoSingle <- plotAUCSingle(singleAUC, singleAUCR,NumVar,test_name)

n <- length(plot)

# ###Accuracy, Precision etc
#
plot3 <-list()
plot4 <-list()
Data <-c("Accuracy","Sensitivity","Specificity","Precision")

MonoParameters <- plotParamSingle(MatsinglePlus, MatsinglePlusR,NumVar,Data,N,test_name)

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
write.csv(Total,paste("Results/",test_name, "_FinalMeanssingle.csv"))
write.csv(Betas_select2,paste("Results/",test_name, "Betas_select.csv"))


###########

#DUMBELL



MeansMultiple <- FinalMeans
MeansSingle <- Total


colnames(MeansSingle) <- c("Real", "Permuted_Data")
Total <- rbind(MeansSingle, Multivariate=data.frame(Real=MeansMultiple$Multi,Permuted_Data=MeansMultiple$MultiRand))
Total["Names"]<- rownames(Total)


plot1 <- ggplot(Total, aes(y=Names, x=Permuted_Data,xend=Real))+geom_dumbbell(color="cornflowerblue", size=1.5, point.colour.l="goldenrod1")+theme(axis.text.x = element_text(angle = 90), axis.text=element_text(size=14),axis.title=element_text(size=18),strip.text.x = element_text(size = 8, colour = "black"),axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))+labs(x="AUC values",y="Selected Features")+xlim(0.45,0.9)
saveRDS(plot1, paste("Results/",test_name, "_Finalplot.rds"))


plot2 <- plotAUCMultipleFig(multipleAUC, multipleAUCR,2)


pdf(paste("FigAUC/", test_name,"_FinalPlot.pdf"),width=10, height=4)
save5 <- grid.arrange(arrangeGrob(plot1,plot2, ncol=2), ncol=1,name=test_name)
graphics.off()

saveRDS(save5, paste("Results/",test_name,"_Finalplot.rds"))
plot(save5)


return(save5)

}





