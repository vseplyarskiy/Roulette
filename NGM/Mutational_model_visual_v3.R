require(ggplot2)
require(epitools)
library(epitools)
setwd("C:/Vova/study/R_scripts_data/R_data/NGM//")

data<-read.table(file="Grant_figure_new_sv3_dan_coding_tri",sep=" ")
data<-read.table(file="Grant_figure_new_sv3_dan_coding",sep=" ")
#<-data[data$V1 != 'ACG_T' & data$V1 != 'TCG_T' &  data$V1 != 'CCG_T' & data$V1 != 'GCG_T',]
data<-data[!((data$V1 == 'ACG_T' | data$V1 == 'TCG_T' |  data$V1 == 'CCG_T' | data$V1 == 'GCG_T') & (data$V5 <= 500000)),]

nrow(data)
data$V5<-data$V5*10^-6#(data$V3/max(data$V3))*200
data$V2<-data$V2*1#(data$V6/max(data$V6))*200
data1<-data[data$V1 == 'TTA_A' | data$V1=='ATA_C' |data$V1=='ATT_C' |data$V1=='ATT_A',]
data1<-data[data$V1 == 'CTC_A' ,]
data1<-data[data$V1 == 'TTT_C' ,]


CI<-pois.approx(data$V3, pt = 1, conf.level = 0.95)


p<-ggplot(data,aes(x = V5,y=V2))+#scale_x_discrete(limits = positions)+#geom_bar()
  # geom_errorbar(aes(ymin=CI[,4]/data$V4, ymax=CI[,5]/data$V4), size=0.5, width=.25,position=position_dodge(.9))+
  geom_pointrange(aes(ymin=CI[,4]/data$V4, ymax=CI[,5]/data$V4))+

 # geom_point(size = 1.2)+
  theme_bw()+
  ylab("Observed rate (de novo)") +
  #ylab("Fracrion of segregated sites") +
  
  #ggtitle("CpG>ApG  (Process 9)")+
 # scale_color_discrete(name="Mutation\nType") + 
  xlab("Predicted mutation rate") +
  #ylim(c(0.9,11))+xlim(c(0.9,11))+
 #ylim(c(0,.029))+
 # xlim(0,10)+ylim(0,10)+
  scale_x_continuous(trans = "log2")+scale_y_continuous(trans = "log2")+
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
          panel.background = element_blank())+theme(legend.position = "none")+
  geom_abline(intercept = 0, slope = 1, color="grey", 
                linetype="dashed", size=1.5)
#theme(legend.justification=c(1,1), legend.position=c(.92,.95))

p



p<-ggplot(data,aes(x = V6*10**4,y=V6))+#scale_x_discrete(limits = positions)+#geom_bar()
  # geom_errorbar(aes(ymin=lower_timing, ymax=upper_timing), size=0.5, width=.25,position=position_dodge(.9))+
  geom_point(size = 1.2)+
  theme_bw()+
  #ylab("Observed rate (de novo)") +
  ylab("Fracrion of segregated sites") +
  
  #ggtitle("CpG>ApG  (Process 9)")+
  scale_color_discrete(name="Mutation\nType") + 
  xlab("Predicted mutation rate") +
  #ylim(c(0.9,11))+xlim(c(0.9,11))+
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
data[data$V6>100 & data$V3<.03,]


data1<-data[data$V1 == 'TTC_C',]
data[data$V3<15 & data$V6>25,]
data[data$V6>130,]