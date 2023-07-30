
library(ggplot2)

setwd("C:/Vova/Vova/study/R_scripts_data/R_data/NGM//")
data<-read.table(file="ADRB2_Primat_AI",header = FALSE,  sep=" ")
cor.test(data$V4,(data$V5), method="pearson")
ggplot(data, aes(x=V4, y=V5)) +
  #geom_bar(stat="identity",  width = 1,position = "dodge2")+
  geom_point()+
  theme_bw()