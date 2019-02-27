


TP <- c("First", "Second", "Third")
#TP <- c("First","First","First")

plota <- readRDS(paste("Results/",TP[1], "_Finalplot.rds"))
plotb <- readRDS(paste("Results/",TP[2], "_Finalplot.rds"))
plotc <- readRDS(paste("Results/",TP[3], "_Finalplot.rds"))

FinalFinal <- grid.arrange(plota,plotb,plotc , nrow=3)

#Figure 1 plot 

TPNISS <- c("FirstNISS", "SecondNISS", "ThirdNISS")
#TPNISS <- c("FirstNISS", "FirstNISS", "FirstNISS")



FinalMeansMultipleNISS <- lapply(1:3, function(x){paste("Results/",TPNISS[x], "_FinalMeansMultiple.csv")})
FinalMeansMultiple <- lapply(1:3, function(x){paste("Results/",TP[x], "_FinalMeansMultiple.csv")})
FMM <- lapply(FinalMeansMultiple, read.csv)
FMMNISS <- lapply(FinalMeansMultipleNISS, read.csv)

FinalMeansSingle <- lapply(1:3, function(x){paste("Results/",TP[x], "_FinalMeanssingle.csv")})
FinalMeansSingleNISS <- lapply(1:3, function(x){paste("Results/",TPNISS[x], "_FinalMeanssingle.csv")})
FMS <- lapply(FinalMeansSingle, read.csv)
FMSNISS <- lapply(FinalMeansSingleNISS, read.csv)



#######



AllFilesPerTime <- lapply(1:3, function(x) list(FMM[[x]],FMMNISS[[x]],FMS[[x]],FMSNISS[[x]])) #one per time point with ISS and no ISS

#s <- do.call("rbind", AllFilesPerTime)

#Here NISS is the one without NISS, if NISS here 1. 

#CHOOSE A TIEM POINT!!!

PlotTP<- 1

FirstNISS <- data.frame(Names=unique(AllFilesPerTime[[PlotTP]][[4]]$X),Single=AllFilesPerTime[[PlotTP]][[4]]$Single,NISS=rep(0),Rand=rep(0))
FirstNISSRand <- data.frame(Names=unique(AllFilesPerTime[[PlotTP]][[4]]$X),Single=AllFilesPerTime[[PlotTP]][[4]]$SingleRandom,NISS=rep(0),Rand=rep(1))
First <- data.frame(Names=unique(AllFilesPerTime[[PlotTP]][[3]]$X),Single=AllFilesPerTime[[PlotTP]][[3]]$Single,NISS=rep(1),Rand=rep(0))
FirstRand <- data.frame(Names=unique(AllFilesPerTime[[PlotTP]][[3]]$X),Single=AllFilesPerTime[[PlotTP]][[3]]$SingleRandom,NISS=rep(1),Rand=rep(1))
FirstMultipleNISS <- data.frame(Names="Multiple",Single=AllFilesPerTime[[PlotTP]][[2]]$Multi,NISS=rep(0),Rand=rep(0))
FirstMultipleNISSRand <- data.frame(Names="Multiple",Single=AllFilesPerTime[[PlotTP]][[2]]$MultiRand,NISS=rep(0),Rand=rep(1))
FirstMultiple <- data.frame(Names="Multiple",Single=AllFilesPerTime[[PlotTP]][[1]]$Multi,NISS=rep(1),Rand=rep(0))
FirstMultipleRand <- data.frame(Names="Multiple",Single=AllFilesPerTime[[PlotTP]][[1]]$MultiRand,NISS=rep(1),Rand=rep(1))


FirstTotal <- rbind(First,FirstRand, FirstNISS,FirstNISSRand, FirstMultiple,FirstMultipleRand, FirstMultipleNISS,FirstMultipleNISSRand)

FirstTotal$NISS <- as.factor(FirstTotal$NISS)
FirstTotal$Rand <- as.factor(FirstTotal$Rand)

FirstTotal["Levels"] <- rep(5)


FirstTotal[FirstTotal$NISS == 0 & FirstTotal$Rand ==0, "Levels"] <- 0
FirstTotal[FirstTotal$NISS == 0 & FirstTotal$Rand ==1, "Levels"] <- 1
FirstTotal[FirstTotal$NISS == 1 & FirstTotal$Rand ==0, "Levels"] <- 2
FirstTotal[FirstTotal$NISS == 1 & FirstTotal$Rand ==1, "Levels"] <- 3

FirstTotal$Levels <- as.factor( FirstTotal$Levels)

End <- FirstTotal %>% filter(FirstTotal$Rand==0) 
Start <- FirstTotal %>% filter(FirstTotal$Rand==1) 
End["SingleRand"] <- Start$Single 

levels(End$NISS)[levels(End$NISS)==0] <- "NISS included"
levels(End$NISS)[levels(End$NISS)==1] <- "NISS not included"


ggplot(End, aes(y=Names, x=SingleRand,xend=Single))+facet_grid(. ~ NISS)+geom_dumbbell(color="cornflowerblue", size=1.5, point.colour.l="goldenrod1")+theme(axis.text.x = element_text(angle = 90), axis.text=element_text(size=12),axis.title=element_text(size=18),strip.text.x = element_text(size = 8, colour = "black"),axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))+labs(x="AUC values",y="Selected Features")+xlim(0.4,0.9)

#store in different and then do grid arrange!

