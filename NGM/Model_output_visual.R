require(ggplot2)

setwd("C:/Vova/Vova/study/R_scripts_data/R_data/NGM//")
data<-read.table(file="visual");

ggplot(data,aes(x=V2, y=V4,color=V3))+geom_point(size=.2,alpha=0.08)+#geom_point(size=0.8,alpha=1/3)+#geom_smooth()+
  geom_smooth(span = 0.03,method = 'loess',n=200,size=.8,alpha=0.7, se=FALSE)+
  # ylab("mutation rate") + 
  ylab("Expected mutation rate") + 
  xlab("Coordinate") + theme_bw()
  
 