require(ggplot2)
require(Hmisc)
setwd("C:/Vova/Vova/study/R_scripts_data/R_data/NGM//")

data<-read.table(file="Grant_figure",sep=" ")


#positions <- c("CpG outside CGI","CpG outside CGI","CpG outside CGI","CpN inside CGI","CpN inside CGI","CpN inside CGI","CpN outside CGI","CpN outside CGI","CpN outside CGI");
conf_timing<-binom.confint(data$V4, data$V5, conf.level = 0.95,methods = "asymptotic")
lower_timing<-conf_timing[,4]
upper_timing<-conf_timing[,5]


p<-ggplot(data,aes(x = V2,y=V3))+#scale_x_discrete(limits = positions)+#geom_bar()
  # geom_errorbar(aes(ymin=lower_timing, ymax=upper_timing), size=0.5, width=.25,position=position_dodge(.9))+
  geom_bar(colour="black",position="dodge",stat="identity")+
  scale_fill_discrete(name="") + 
  ylab("Mutation density in validation cohort (de novo)") +
  #ggtitle("CpG>ApG  (Process 9)")+
  facet_wrap(~V1,scale="free",ncol = 4, nrow = 24)+
  
  
  xlab("Decile of predicted mutation rate") +
  
  theme(legend.text = element_text(size = 30),legend.title = element_text( size=30, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(  axis.text.x=element_blank(),
          axis.title.x = element_text(size = 50),
          axis.title.y = element_text(size = 50),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+theme(legend.position = "none")
#theme(legend.justification=c(1,1), legend.position=c(.92,.95))

#p
ggsave(filename="C:/Vova/Vova/study/Signature germline/Shamil_grant/Simple_decile.svg", plot = p,device = "svg", path = NULL,
       scale = 1, width = 600, height = 1000, units = c("mm"),
       dpi = 300, limitsize = TRUE)



data<-read.table(file="Grant_figure_small",sep=" ")


#positions <- c("CpG outside CGI","CpG outside CGI","CpG outside CGI","CpN inside CGI","CpN inside CGI","CpN inside CGI","CpN outside CGI","CpN outside CGI","CpN outside CGI");
conf_timing<-binom.confint(data$V4, data$V5, conf.level = 0.95,methods = "asymptotic")
lower_timing<-conf_timing[,4]
upper_timing<-conf_timing[,5]


p<-ggplot(data,aes(x = V2,y=V3))+#scale_x_discrete(limits = positions)+#geom_bar()
  # geom_errorbar(aes(ymin=lower_timing, ymax=upper_timing), size=0.5, width=.25,position=position_dodge(.9))+
  geom_bar(colour="black",position="dodge",stat="identity")+
  scale_fill_discrete(name="") + 
  ylab("Mutation density in validation cohort (de novo)") +
  #ggtitle("CpG>ApG  (Process 9)")+
  facet_wrap(~V1,scale="free",ncol = 3, nrow = 4)+
  
  
  xlab("Decile of predicted mutation rate") +
  
  theme(legend.text = element_text(size = 30),legend.title = element_text( size=30, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(  axis.text.x=element_blank(),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          axis.text.y=element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.line.x = element_blank(),
          axis.line.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())+theme(legend.position = "none")
#theme(legend.justification=c(1,1), legend.position=c(.92,.95))

p
ggsave(filename="C:/Vova/Vova/study/Signature germline/Shamil_grant/Simple_decile_small.svg", plot = p,device = "svg", path = NULL,
       scale = 1, width = 360, height = 200, units = c("mm"),
       dpi = 300, limitsize = TRUE)




data<-read.table(file="Grant_figure_new_short",sep=" ")
data<-data[data$V1 != 'TTT_C',]
data$V3<-(data$V3/min(data$V3))#*11
data$V6<-(data$V6/min(data$V6))#*11

p<-ggplot(data,aes(x = V6,y=V3,color=V1))+#scale_x_discrete(limits = positions)+#geom_bar()
  # geom_errorbar(aes(ymin=lower_timing, ymax=upper_timing), size=0.5, width=.25,position=position_dodge(.9))+
  geom_point(size = 2.5)+
  theme_bw()+
  ylab("Observed rate (de novo)") +
  #ggtitle("CpG>ApG  (Process 9)")+
   scale_color_discrete(name="Mutation\nType") + 
    xlab("Predicted mutation rate") +
  #ylim(c(0.9,11))+xlim(c(0.9,11))+
  ylim(c(1,7))+xlim(1,7)+
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
          panel.background = element_blank())#+theme(legend.position = "none")
#theme(legend.justification=c(1,1), legend.position=c(.92,.95))

p
ggsave(filename="C:/Vova/Vova/study/Signature germline/Shamil_grant/Simple_nice_small.png", plot = p,device = "png", path = NULL,
       scale = 1, width = 100, height = 70, units = c("mm"),
       dpi = 300, limitsize = TRUE)

r_mut=c('TTT_A','TTT_G','TTT_C','CCC_A','CCC_G','CCC_T','ACT_T','ACT_G','ACT_A','GTA_A','GTA_C','GTA_G');



