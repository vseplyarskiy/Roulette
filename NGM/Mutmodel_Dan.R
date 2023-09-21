require(ggplot2)
require(reshape2)
require(PropCIs)
library(scales)


setwd("C:/Vova/study/R_scripts_data/R_data/NGM")
data<-read.table(file="Grant_figure_new_sv3_pol_dan", header = F, sep=" ")

x=1000000
x=1
data<-cbind(data,-log(1-data$V2))
names(data)[6]<-'V6'
data.lm<-lm(data=data,V6~V5+0)
p<-ggplot(data,aes(x = V5,y=V6))+#scale_x_discrete(limits = positions)+#geom_bar()
  # geom_errorbar(aes(ymin=lower_timing, ymax=upper_timing), size=0.5, width=.25,position=position_dodge(.9))+
  geom_point(size = 1.2)+
  theme_bw()+
  geom_abline(slope = coef(data.lm)[[1]], intercept = 0)+
#ylab("Observed rate (de novo)") +
  ylab("Fracrion of segregated sites, SNP") +
  ylab("-log(Monomorphic)") +
  
  #ggtitle("CpG>ApG  (Process 9)")+
#  scale_color_discrete(name="Mutation\nType") + 
  xlab("Predicted mutation rate") +
  #ylim(c(0,.3))+xlim(c(0,500))+
  #ylim(c(0,.029))+
  # xlim(0,10)+ylim(0,10)+
  #scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2")+
  theme(legend.text = element_text(size = 15),legend.title = element_text( size=15, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(  axis.text=element_text(size = 15),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          #axis.text.y=element_blank(),
          #axis.ticks.x=element_blank(),
          #axis.ticks.y=element_blank(),
          #axis.line.x = ,
          #axis.line.y = black,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          #panel.border = element_blank(),
          panel.background = element_blank())+theme(legend.position = "none")
#theme(legend.justification=c(1,1), legend.position=c(.92,.95))

p
