WD <- paste(path.expand("~"), "/Documents/PhD/R/Pipeline/", sep="")
setwd(WD)
source("DumbellwithoutNISSFun.R")
source("DumbellwithNISSFun.R")
Part1 <- 1
Part2 <- 0
source("SettingsPipe.R")

settings$seed <- 132


set.seed(settings$seed)




CF_ISS_1 <- read.csv("CF_ISS1.csv")
CF_ISS_2 <- read.csv("CF_ISS2.csv")
CF_ISS_3 <- read.csv("CF_ISS3.csv")

names(CF_ISS_1) <- c("","ID","Label","White Blood Cells [109/L]","Neutrophil [109/L]","Lymphocytes [109/L]","Monocytes [109/L]","Immature Granulocytes [109/L]","Neutrophils %","Lymphocytes %","Monocytes %","Immature Granulocytes %","ROS TO PMA (%)","ROS TO PMA (MFI)","CD62L (MedFI)","CD88 (MedFI)","CD16 (MedFI)","CD11B (MedFI)","CXCR1 (MedFI)","CXCR2 (MedFI)","CD11B after FMLF (MedFI)","Fold increase in CD11b after fMLf","CD63 (MedFI)","CD62L after FMLF (MedFI)","Fold decrease in CD62L (MedFI)","CD16BRIGHT CD62LDIM % ","CD16BRIGHT CD62LDIM [106/L]","CD16+ %","CD16 (MedFI).1","CD14+ 16- %","CD14+ 16+ % ","CD14+ 16- [x106/L]","CD14+ 16+ [x106/L]","HLA-DR+ %","HLA-DR (MedFI)","CD14+ HLA-DRLow/Neg %","CD14+ HLA-DRLow/Neg [x106/L]","TLR4+ %","TLR4 (MedFI)","TLR2+ %","TLR2 (MedFI)","CD86+ %","CD86 (MedFI)","B Cells [x106/L]","Natural Killler Cells [x106/L]","CD56DIM Natural Killer Cells [x106/L]","CD56BRIGHT Natural Killer cells [x106/L]","Natural Killer T Cells [x106/L]","CD3+ [106/L]","CD3+ 4+ [106/L]","Terminal Helper T Cells","Naive Helper T Cells","Central Memory Helper T Cells","Effector Memory Helper T Cells","CD3+ 8+ [106/L]","Terminal Cytotoxic T Cells","Naive Cytotoxic T Cells","Central Memory Cytotoxic T Cells","Effector Memory Cytotoxic T Cells","IL-6 1LPS","IL-8 1LPS","IL-10 1LPS","TNF-A 1LPS","MCP-1 1LPS","IL-6 10LPS","IL-8 10LPS","IL-10 10LPS","TNF-A 10LPS","MCP-1 10LPS","Serum IL1.Ra","Serum IL.6","Serum IL.8","Serum IL.10","Serum G.CSF","Serum MCP.1","Serum TNF-A","Cortisol [ng/ml]","Cell Free DNA","HMGB1","ISS","NISS")
names(CF_ISS_2) <- c("","ID","Label","White Blood Cells [109/L]","Neutrophil [109/L]","Lymphocytes [109/L]","Monocytes [109/L]","Immature Granulocytes [109/L]","Neutrophils %","Lymphocytes %","Monocytes %","Immature Granulocytes %","ROS TO PMA (%)","ROS TO PMA (MFI)","CD62L (MedFI)","CD88 (MedFI)","CD16 (MedFI)","CD11B (MedFI)","CXCR1 (MedFI)","CXCR2 (MedFI)","CD11B after FMLF (MedFI)","Fold increase in CD11b after fMLf","CD63 (MedFI)","CD62L after FMLF (MedFI)","Fold decrease in CD62L (MedFI)","CD16BRIGHT CD62LDIM % ","CD16BRIGHT CD62LDIM [106/L]","CD16+ %","CD16 (MedFI).1","CD14+ 16- %","CD14+ 16+ % ","CD14+ 16- [x106/L]","CD14+ 16+ [x106/L]","HLA-DR+ %","HLA-DR (MedFI)","CD14+ HLA-DRLow/Neg %","CD14+ HLA-DRLow/Neg [x106/L]","TLR4+ %","TLR4 (MedFI)","TLR2+ %","TLR2 (MedFI)","CD86+ %","CD86 (MedFI)","B Cells [x106/L]","Natural Killler Cells [x106/L]","CD56DIM Natural Killer Cells [x106/L]","CD56BRIGHT Natural Killer cells [x106/L]","Natural Killer T Cells [x106/L]","CD3+ [106/L]","CD3+ 4+ [106/L]","Terminal Helper T Cells","Naive Helper T Cells","Central Memory Helper T Cells","Effector Memory Helper T Cells","CD3+ 8+ [106/L]","Terminal Cytotoxic T Cells","Naive Cytotoxic T Cells","Central Memory Cytotoxic T Cells","Effector Memory Cytotoxic T Cells","IL-6 1LPS","IL-8 1LPS","IL-10 1LPS","TNF-A 1LPS","MCP-1 1LPS","IL-6 10LPS","IL-8 10LPS","IL-10 10LPS","TNF-A 10LPS","MCP-1 10LPS","Nuclear DNA [ng/ml]","Mitochondrial DNA [ng/ml]","Serum IL1.Ra","Serum IL.6","Serum IL.8","Serum IL.10","Serum G.CSF","Serum MCP.1","Serum TNF-A","Cortisol [ng/ml]","Cell Free DNA","HMGB1","ISS","NISS")
names(CF_ISS_3)<- c("","ID","Label","White Blood Cells [109/L]","Neutrophil [109/L]","Lymphocytes [109/L]","Monocytes [109/L]","Immature Granulocytes [109/L]","Neutrophils %","Lymphocytes %","Monocytes %","Immature Granulocytes %","ROS TO PMA (%)","ROS TO PMA (MFI)","CD62L (MedFI)","CD88 (MedFI)","CD16 (MedFI)","CD11B (MedFI)","CXCR1 (MedFI)","CXCR2 (MedFI)","CD11B after FMLF (MedFI)","Fold increase in CD11b after fMLf","CD63 (MedFI)","CD62L after FMLF (MedFI)","Fold decrease in CD62L (MedFI)","CD16BRIGHT CD62LDIM % ","CD16BRIGHT CD62LDIM [106/L]","CD16+ %","CD16 (MedFI).1","CD14+ 16- %","CD14+ 16+ % ","CD14+ 16- [x106/L]","CD14+ 16+ [x106/L]","HLA-DR+ %","HLA-DR (MedFI)","CD14+ HLA-DRLow/Neg %","CD14+ HLA-DRLow/Neg [x106/L]","TLR4+ %","TLR4 (MedFI)","TLR2+ %","TLR2 (MedFI)","CD86+ %","CD86 (MedFI)","B Cells [x106/L]","Natural Killler Cells [x106/L]","CD56DIM Natural Killer Cells [x106/L]","CD56BRIGHT Natural Killer cells [x106/L]","Natural Killer T Cells [x106/L]","CD3+ [106/L]","CD3+ 4+ [106/L]","Terminal Helper T Cells","Naive Helper T Cells","Central Memory Helper T Cells","Effector Memory Helper T Cells","CD3+ 8+ [106/L]","Terminal Cytotoxic T Cells","Naive Cytotoxic T Cells","Central Memory Cytotoxic T Cells","Effector Memory Cytotoxic T Cells","IL-6 1LPS","IL-8 1LPS","IL-10 1LPS","TNF-A 1LPS","MCP-1 1LPS","IL-6 10LPS","IL-8 10LPS","IL-10 10LPS","TNF-A 10LPS","MCP-1 10LPS","Nuclear DNA [ng/ml]","Mitochondrial DNA [ng/ml]","Serum IL1.Ra","Serum IL.6","Serum IL.8","Serum IL.10","Serum G.CSF","Serum MCP.1","Serum TNF-A","Cortisol [ng/ml]","Cell Free DNA","HMGB1","ISS","NISS")


#### Function 

test_nameNISS <- c("FirstNISS","SecondNISS","ThirdNISS")
test_name <- c("First","Second","Third")
s <-list(CF_ISS_1,CF_ISS_2,CF_ISS_3)


listNISS <- list()
list <- list()

for (i in 1:3){
  listNISS[[i]] <- DumbellwithoutNISSFun(test_nameNISS[i], s[[i]],i)
  list[[i]] <- DumbellwithNISSFun(test_name[i], s[[i]],i)
}




################################

#CHANGED ALL SCRIPTS NAMES - SEE FROM WHAT COMES OUT.
# 
# 


# For first plot - maybe just as simple as this, plus send it to the figure file. 

pdf(paste("Figures/","WithoutNISS_FinalPlot.pdf"),width=15, height=15)
FinalFin<- grid.arrange(listNISS[[1]],listNISS[[2]],listNISS[[3]],nrow=3)
graphics.off()


pdf(paste("Figures/","WithNISS_Figure1.pdf"),width=14, height=11)
FinalFin<- grid.arrange(list[[1]],list[[2]],list[[3]],nrow=3)
graphics.off()

 
TPNISS <- test_nameNISS
TP <- test_name 

FinalMeansMultipleNISS <- lapply(1:3, function(x){paste("Results/",TPNISS[x],"_FinalMeansMultiple.csv")})
FinalMeansMultiple <- lapply(1:3, function(x){paste("Results/",TP[x],"_FinalMeansMultiple.csv")})
FMM <- lapply(FinalMeansMultiple, read.csv)
FMMNISS <- lapply(FinalMeansMultipleNISS, read.csv)

FinalMeansSingle <- lapply(1:3, function(x){paste("Results/",TP[x], "FinalMeanssingle.csv")})
FinalMeansSingleNISS <- lapply(1:3, function(x){paste("Results/",TPNISS[x], "_FinalMeanssingle.csv")})
FMS <- lapply(FinalMeansSingle, read.csv)
FMSNISS <- lapply(FinalMeansSingleNISS, read.csv)



#######



AllFilesPerTime <- lapply(1:3, function(x) list(FMM[[x]],FMMNISS[[x]],FMS[[x]],FMSNISS[[x]])) #one per time point with ISS and no ISS

#s <- do.call("rbind", AllFilesPerTime)

#Here NISS is the one without NISS, if NISS here 1.

#CHOOSE A TIME POINT!!!

Final <- list()
for (PlotTP in 1:3){

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

levels(End$NISS)[levels(End$NISS)==0] <- "NISS not included"
levels(End$NISS)[levels(End$NISS)==1] <- "NISS included"


Final[[PlotTP]]<- ggplot(End, aes(y=Names, x=SingleRand,xend=Single))+facet_grid(. ~ NISS)+geom_dumbbell(color="cornflowerblue", size=1.5, point.colour.l="goldenrod1")+theme(axis.text.x = element_text(angle = 90), axis.text=element_text(size=12),axis.title=element_text(size=18),strip.text.x = element_text(size = 8, colour = "black"),axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))+labs(x="AUC values",y="Selected Features")+xlim(0.47,0.9)

#store in different and then do grid arrange!

}

pdf(paste("Figures/","WithoutNISS_SupplementaryFig1.pdf"),width=9, height=12)
FinalFin<- grid.arrange(Final[[1]],Final[[2]],Final[[3]],nrow=3)
graphics.off()

