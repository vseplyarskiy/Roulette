require(ggplot2)

setwd("C:/Vova/study/R_scripts_data/R_data/NGM/")

data<-read.table(file="rate_evan_penta_outliers_v3",sep=" ",header=F)
sum (data$V7)
data$V2 <- factor(data$V2,levels = unique(data$V2))

p<-ggplot(data, aes(x=V2, y=V4,fill=V3)) +
  geom_bar(stat="identity",  width = 1,position = "dodge2")+
  theme_bw()+#ggtitle("rep")+#ylim(c(0,0.031))+
  #scale_fill_discrete(name="Mutation\ntype")+
  #scale_fill_manual("legend",values = cols)+ 
  #scale_fill_manual(values=cols2, name="mutation\ntype")+
  
  #theme(title=element_text(size=25))+
  ylab("SFS distortion (AF 0.05-0.2)")+xlab("mutation rate")+
  #theme( axis.text.y=element_text(size=25),axis.text.x=element_text(size=4),axis.title.y =element_text(size=25))+theme(axis.text.x = element_text(angle = 90, vjust = 0.35,size=8))+

 # theme(legend.key.width = unit(0.5, "cm"),legend.key.height = unit(0.5, "cm"))+
 # ggtitle(compon)+theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())#+
  #theme(legend.position = "none")#+
#guides(fill=guide_legend(ncol=3))+ theme(legend.title.align=0.5,legend.title =  element_text(size=25))+
#scale_fill_manual(name="", labels=c("b"="remaning genome", "r"="high intensity"))+

#theme(legend.position = "none")+
#theme(axis.title.x=element_blank(),
#     axis.text.x=element_blank(),
#    axis.ticks.x=element_blank(),axis.line.x=element_blank())
p