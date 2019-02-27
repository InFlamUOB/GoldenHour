
#open saved files so this has to follow from Main Pipe directly. 

library("tibble")
setwd("/Users/laura/Documents/PhD/R/GoldenHourGit")

source("/Users/laura/Documents/PhD/R/Pipeline/Names_better.R")


#These are the selected features, and so the figures we are generating
Names_to_plot2 <- c("CD62L after FMLF (MedFI)","CD11B (MedFI)","Monocytes %","MODS","ID")

 

Betas_plot_1h <- All1[,colnames(All1[,intersect(gsub("`", "",  Names_to_plot2), colnames(All1))])]
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

names(Betas_plot) <- c("CD62LafterFMLFMedFI","CD11BMedFI","MonocytesPercentage","MODS","ID","Time")

nm <- names(Betas_plot)
levels(Betas_plot$MODS) <- c("No","Yes")
levels(Betas_plot$Time) <- c("T<1h","4<T<12h","48<T<72h")



for (i in nm){
 
  dodge <- position_dodge(width = 0.9)
    print(i)
    p1 <- ggplot(Betas_plot,aes_string(x="MODS",y=i,fill="Time")) + geom_violin(position=dodge)+geom_boxplot(position=dodge, width=0.2)+labs(title= paste("MODS vs Non MODS of ", i),y=paste(i,"count"),x="MODS")+geom_hline(yintercept = 0, linetype = 2)+labs(x="MODS")
    #print(Betas_plot,aes_string(x="MODS",y=i[i],fill="TimePoint")) + geom_violin(position=dodge)+geom_boxplot(position=dodge, width=0.2)+labs(title= paste("Count of", R[i]),y=paste("Real",R[i],"count"),x="Development of MODS")+geom_hline(yintercept = 0, linetype = 2)
    p2 <- ggplot(Betas_plot,aes_string(x="Time", i,group="ID",label="ID"))+geom_point()+geom_line()+facet_grid(. ~ MODS)+geom_text(check_overlap = TRUE,vjust = 0, nudge_x = 0.2)+scale_x_discrete(labels=c("1" = "T<=1h", "2" = "4<T<12","3" = "48<T<72"))+labs(title= paste("MODS vs NonMODS of", i),y=paste("Real",i,"count"))+
      labs(x="TimePoints")
    p <- ggboxplot(Betas_plot,x="Time", y=i,color = "MODS", palette = "jco")+labs(title= paste(i),y=paste(i,"count"),x=("Time"))
    p3 <- p+geom_hline(yintercept = 0, linetype = 2)+stat_compare_means(aes(group = MODS), label = "p.format",method = "t.test") 
    
    plot[[i]] <- p1
    #print(plot[[i]])
    plot2[[i]] <- p2
    #print(plot2[[i]])
    plot3[[i]] <- p3
    #print(plot3[[i]])
}


grid.arrange(plot3[[2]],plot[[2]],plot3[[3]],plot[[3]],plot3[[1]],plot[[1]], nrow=3, ncol=2)

#Check p.values are the correct ones
 
s<- Betas_plot %>% filter(Time==3) 
p <- with(s, t.test(CD62LafterFMLFMedFI[MODS == "Y"], CD62LafterFMLFMedFI[MODS == "N"],paired= FALSE))$p.value
print(p)


#Normally Lactate springs up as significant - include Lactate here: Betas_plot_1h: obtain p-value of:0.6420177
#p <- with(Betas_plot_1h, t.test(Lactate[MODS == "Y"], Lactate[MODS == "N"],paired= FALSE))$p.value
#print(p)

#p value works - stick with it

#cannot compute wilcoxon as cannot compute exact p-value with ties
#mm <- with(s, wilcox.test(CD62L_FMLF_MedFI[MODS == "Y"], CD62L_FMLF_MedFI[MODS == "N"],alternative = "two.sided"),paired=FALSE)


