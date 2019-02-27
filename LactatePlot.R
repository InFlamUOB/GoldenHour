

LactateDataA <- read_xlsx("Lactate.xlsx")
LactateData <- LactateDataA[1:1585,] #take out blank space from data base

LactateData$`Date First Screened` <- as.Date(LactateData$`Date First Screened`)

LactateData["TimePoint"] <-data.frame(rep(0,length(LactateData$DateTime)))

LactateData$DateTime <- as.POSIXct(LactateData$DateTime)   


LactateData3 <- LactateData[ LactateData$HospitalNo %like% "GH",]
LactateData3 <- LactateData3[!(LactateData3$Lactate %like% "NA"),]

LengthMat<-length(LactateData3$HospitalNo)
LactateData3["TimeDiff"] <- data.frame(rep(0,length(LactateData3$HospitalNo)))


#loop that establishes the time differnce between time of accident and blood sample taken in minutes. For each individual new time. 
set<- 1
count <-1 
for (i in 2:LengthMat){
  count=count+1
  
  if (LactateData3$HospitalNo[i-1] != LactateData3$HospitalNo[i]){
    LactateData3$TimeDiff[count] <- 0
    set <- i
  }
  
  else{
    LactateData3$TimeDiff[count] <- difftime(LactateData3$DateTime[i],LactateData3$DateTime[set],units='hours')
  }
}

#0, 4-12,48-72

s <- data.frame(table(LactateData3$HospitalNo)) #can see that each patient has difefrent readings.
MV <- data.frame(IDs = s$Var1)

MV3 <- LactateData3

counts <- 0
MV4 <- data.frame(NULL)
Lactate9 <- data.frame(matrix(0,3,3))

for (i in 1:dim(MV)[1]){
  
  counts <- counts + 1;
  # LactateData5 <- MV3[ MV3$HospitalNo %like% MV2$IDs[i],]
  
  
  LactateData5 <- MV3[ MV3$HospitalNo %like% MV$IDs[i],]
  
  firstPoint <- as.numeric(LactateData5[1,4])
  Lactate7 <- filter(LactateData5 , TimeDiff<=12 & TimeDiff>=4)
  secondPoint <- mean(as.numeric(Lactate7$Lactate))
  Lactate8 <- filter(LactateData5 , TimeDiff<=72 & TimeDiff>=48)
  thirdPoint <- mean(as.numeric(Lactate8$Lactate))
  
  Lactate9["ID"] <- data.frame(rep(MV$IDs[i],3))
  Lactate9["Lactate"] <- as.data.frame(rbind(firstPoint,secondPoint,thirdPoint))
  Lactate9["Time"] <- as.data.frame(c(1,2,3))
  
  
  MV4 <- rbind(MV4,Lactate9)
  
  
}

MV4 <- MV4[,-(1:3)]

MV5 <- MV4[complete.cases(MV4),]


#Here not filtered with those that have not loads of measurements per time point so more NANs, key decide what to do with these. 
#Plot all lactates together, with those filtered and lines unitingthem to get to see the dynamics. Exploit data.

#Now have to write MODS and nonMODS

CF_ISS_1 <- read.csv("CF_ISS1.csv")
names(CF_ISS_1) <- c("","ID","Label","White Blood Cells [109/L]","Neutrophil [109/L]","Lymphocytes [109/L]","Monocytes [109/L]","Immature Granulocytes [109/L]","Neutrophils %","Lymphocytes %","Monocytes %","Immature Granulocytes %","ROS TO PMA (%)","ROS TO PMA (MFI)","CD62L (MedFI)","CD88 (MedFI)","CD16 (MedFI)","CD11B (MedFI)","CXCR1 (MedFI)","CXCR2 (MedFI)","CD11B after FMLF (MedFI)","Fold increase in CD11b after fMLf","CD63 (MedFI)","CD62L after FMLF (MedFI)","Fold decrease in CD62L (MedFI)","CD16BRIGHT CD62LDIM % ","CD16BRIGHT CD62LDIM [106/L]","CD16+ %","CD16 (MedFI).1","CD14+ 16- %","CD14+ 16+ % ","CD14+ 16- [x106/L]","CD14+ 16+ [x106/L]","HLA-DR+ %","HLA-DR (MedFI)","CD14+ HLA-DRLow/Neg %","CD14+ HLA-DRLow/Neg [x106/L]","TLR4+ %","TLR4 (MedFI)","TLR2+ %","TLR2 (MedFI)","CD86+ %","CD86 (MedFI)","B Cells [x106/L]","Natural Killler Cells [x106/L]","CD56DIM Natural Killer Cells [x106/L]","CD56BRIGHT Natural Killer cells [x106/L]","Natural Killer T Cells [x106/L]","CD3+ [106/L]","CD3+ 4+ [106/L]","Terminal Helper T Cells","Naive Helper T Cells","Central Memory Helper T Cells","Effector Memory Helper T Cells","CD3+ 8+ [106/L]","Terminal Cytotoxic T Cells","Naive Cytotoxic T Cells","Central Memory Cytotoxic T Cells","Effector Memory Cytotoxic T Cells","IL-6 1LPS","IL-8 1LPS","IL-10 1LPS","TNF-A 1LPS","MCP-1 1LPS","IL-6 10LPS","IL-8 10LPS","IL-10 10LPS","TNF-A 10LPS","MCP-1 10LPS","Serum IL1.Ra","Serum IL.6","Serum IL.8","Serum IL.10","Serum G.CSF","Serum MCP.1","Serum TNF-A","Cortisol [ng/ml]","Cell Free DNA","HMGB1","ISS","NISS")


Hazeldine_1h <- CF_ISS_1
MV6 <- MV5[MV5$Time==1,]

e <- data.frame(hello=1)
e["hello"] <- data.frame(NULL)

see1 <- strsplit(as.character(MV6$ID),"H")
for (i in 1:length(see1)){
  
  e[i] <- as.numeric(see1[[i]][2])
}


tt2 <- match(Hazeldine_1h$ID,e)

Hazeldine_1h["Lactate"] <- MV6$Lactate[tt2]
Hazeldine_1h <- Hazeldine_1h[complete.cases(Hazeldine_1h),]



#names(Hazeldine_1h)<-sapply(1:dim(Hazeldine_1h)[2], function(i){paste0("First ",names(Hazeldine_1h)[i])})

Hazeldine_1h["Time"] <- as.data.frame(rep(1,dim(Hazeldine_1h)[1]))


###Script were merge with Hazeldine 1h and see if running it with 41 patients makes the feature come up as important or not- it doesnt. 
#Now plot. 



MV5["ID"] <- as.numeric(gsub("GH","",MV5$ID))
MODS <- data.frame(ID=Hazeldine_1h$ID,Label=Hazeldine_1h$Label)
FinalPlot <- inner_join(MV5,MODS)

FinalPlot$Time <- as.factor(FinalPlot$Time)

ggplot(FinalPlot, aes(Time, Lactate, fill=Time, palette = "jco"))+geom_boxplot(aes(fill=Label))+scale_x_discrete(labels=c("1" = "T<=1h", "2" = "4<T<12","3" = "48<T<72"))+labs(title= paste("MODS vs NonMODS of Lactate at Three Different Times"))
d <- list()
for ( i in 1:3){
  d[i] <- dim(FinalPlot[FinalPlot$Time==i & FinalPlot$Label=="No"  ,])[1]
  count[i] <- dim(FinalPlot[FinalPlot$Time==i & FinalPlot$Label=="Yes"  ,])[1]
}

plot1<-ggplot(FinalPlot, aes(Time, Lactate, fill=Time, palette = "jco"))+geom_boxplot(aes(fill=Label))+scale_x_discrete(labels=c("1" = paste("T<1h [No:",d[1],",Yes:",count[1],"]"),"2" = paste("4<T<12h [No:",d[2],",Yes:",count[2],"]"),"3" = paste("48<T<72h [No:",d[3],",Yes:",count[3],"]")))+
  labs(title= paste("MODS vs NonMODS of Lactate at Three Different Times"))+stat_compare_means(aes(group = Label), label = "p.format", paired=FALSE,method="t.test")+theme(axis.text.x = element_text(size=12))+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=16),strip.text.x = element_text(size = 16, colour = "black"))+
  theme(legend.text=element_text(size=16),legend.title = element_text(size=16))


pdf(paste("Figures/Lactate_SummaryFigure4.pdf"),width=9,height=4)
grid.arrange(plot1)
graphics.off()


