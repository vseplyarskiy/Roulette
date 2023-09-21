require(ggplot2)

setwd("C:/Vova/Vova/study/R_scripts_data/R_data/NGM//")
data<-read.table(file="mutation_regional_ATA_C");
data1<-read.table(file="mutation_regional_ATA_C_expected");

data<-cbind(data,data1[,3])
ggplot(data,aes(y=V2*2,x=data$`data1[, 3]`,color=V5))+
  geom_point() +
ylab("Observed mutation rate")+
    xlab("Expected mutation rate bin")+
  theme_bw()