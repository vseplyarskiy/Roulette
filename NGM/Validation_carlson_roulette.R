require(ggplot2)
require(epitools)


data<-read.table(file="rate_denovo_Carlson", sep=" ");
data[1,1]<-0.05
c=1.29286405764761e-05
mil=1000000;
data<-cbind(data,data[,4])
names(data)[5]<-"V5"
data[,1]*c->data[,1]
data[,1]*mil->data[,1]
data[,2]*mil->data[,2]
data[,4]/mil->data[,4]



CI<-pois.approx(data$V3, pt = 1, conf.level = 0.95)


p<-ggplot(data,aes(x = V1,y=V2,size=V5))+#scale_x_discrete(limits = positions)+#geom_bar()
  # geom_errorbar(aes(ymin=CI[,4]/data$V4, ymax=CI[,5]/data$V4), size=0.5, width=.25,position=position_dodge(.9))+
  geom_pointrange(aes(ymin=CI[,4]/data$V4, ymax=CI[,5]/data$V4))+
   scale_size(range = c(0.1, 1.0))+
  # geom_point(size = 1.2)+
  theme_bw()+
  ylab("Observed rate (de novo)") +
  ggtitle("Carlson et al")+
  # scale_color_discrete(name="Mutation\nType") + 
  #geom_smooth(method='lm')+
  xlab("Predicted mutation rate") +
  #ylim(c(0.9,11))+xlim(c(0.9,11))+
  #ylim(c(0,.029))+
 #  xlim(0.05,100)+
#+ylim(0,10)+
  scale_x_continuous(trans = "log2",limits = c(1,1000),breaks = c(1,10,100,1000))+
  scale_y_continuous(trans = "log2",limits = c(1,1000),breaks = c(1,10,100,1000))+
  theme(legend.text = element_text(size = 15),legend.title = element_text( size=15, face="bold"))+
  theme(plot.title = element_text(hjust = 0.45))+
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
              linetype="dashed", size=1)
#theme(legend.justification=c(1,1), legend.position=c(.92,.95))

p



data<-read.table(file="rate_denovo_our_model", sep=" ");
#data[1,1]<-0.05
c=0.000653909515334231
mil=1000000;
data<-cbind(data,data[,4])
names(data)[5]<-"V5"
data[,1]*c->data[,1]
data[,1]*mil->data[,1]
data[,2]*mil->data[,2]
data[,4]/mil->data[,4]

CI<-pois.approx(data$V3, pt = 1, conf.level = 0.95)


p<-ggplot(data,aes(x = V1,y=V2,size=V5))+#scale_x_discrete(limits = positions)+#geom_bar()
  # geom_errorbar(aes(ymin=CI[,4]/data$V4, ymax=CI[,5]/data$V4), size=0.5, width=.25,position=position_dodge(.9))+
  geom_pointrange(aes(ymin=CI[,4]/data$V4, ymax=CI[,5]/data$V4))+
  # geom_point(size = 1.2)+
  scale_size(range = c(0.1, 1.0))+
  
  theme_bw()+
  ylab("Observed rate (de novo)") +
  ggtitle("Roulette")+
  # scale_color_discrete(name="Mutation\nType") + 
  #geom_smooth(method='lm')+
  xlab("Predicted mutation rate") +
  #ylim(c(0.9,11))+xlim(c(0.9,11))+
  #ylim(c(0,.029))+
  #  xlim(0.05,100)+
  #+ylim(0,10)+
  scale_x_continuous(trans = "log2",limits = c(1,1000),breaks = c(1,10,100,1000))+
  scale_y_continuous(trans = "log2",limits = c(1,1000),breaks = c(1,10,100,1000))+
  theme(legend.text = element_text(size = 15),legend.title = element_text( size=15, face="bold"))+
  theme(plot.title = element_text(hjust = 0.45))+
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
              linetype="dashed", size=1)
#theme(legend.justification=c(1,1), legend.position=c(.92,.95))

p
