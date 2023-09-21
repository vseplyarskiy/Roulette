require(ggplot2)
require(reshape2)
require(PropCIs)
library(scales)

setwd("C:/Vova/study/R_scripts_data/R_data/NGM//")
data<-read.table(file="rate_Histones_agrig", sep=" ");

#Cis<-riskscoreci(x1=as.numeric(data$V4), n1=as.numeric(data$V6), x2=as.numeric(data$V5), n2=as.numeric(data$V7), conf.level=0.95)



#data$V1 <- factor(data$V1, levels = data$V1)
p<-ggplot(data,aes(x=V2, fill=V1,y=V3))+ geom_bar(position="dodge", stat="identity")+ 
  #geom_violin(trim=FALSE, fill="gray", adjust=2)+
  
  scale_y_continuous(trans = log2_trans())+#,,limits=c(0.5, 3))+
  
  #facet_wrap(~V8)+
 # geom_errorbar(aes(ymin=lower_timing, ymax=upper_timing),stat="identity", position=position_dodge())+
  #scale_fill_manual(values=c("red", adjustcolor("red",0.4),"blue"),name="")+
  #geom_point(alpha=0.1, size=0.4)+
  
  #scale_color_manual(name = "Mutation type", values =colorConverter(data$V2))+
  ylab("observed/expected") + 
  xlab("") + theme_bw()+#xlim(c(genelb,generb))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme( axis.text.y=element_text(size=11),axis.text.x=element_text(size=10,angle = 90,vjust = .25),axis.title=element_text(size=9),
         legend.text=element_text(size=11))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())#+ #theme(legend.position = "none")
#labs(fill = "Origin")
p