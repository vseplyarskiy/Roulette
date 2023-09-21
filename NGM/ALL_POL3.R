library(ggplot2)

setwd("C:/Vova/study/R_scripts_data/R_data/NGM/")
#data<-read.table(file="distance_low_FDR", header = F, sep=" ")
#data<-read.table(file="rate_Rayn", header = F, sep=" ")
data<-read.table(file="ALLP3", header = F, sep=" ")

data1<-data#data1<-data[data[,2]<5000,]
#data1<-cbind(data1,data1[,4]/data1[,3])
p<-ggplot(data1, aes(x=log(data1[,3]), y=data1[,4]/data1[,3],label=data1[,6], color=as.factor(data1[,7]))) + geom_point() + 
  geom_text(label=data1[,6],nudge_x = .5,size=4)+

  ylab("observed/expected")+xlab("expected")+
 # geom_abline(intercept = 0, slope = 1, color="gray", linetype="dashed", size=1)+
  theme_classic()+
  #scale_x_continuous(trans = 'log2')+
  scale_color_manual(name = "",values=c("red","black"),
                     labels=c("active gene", "pseudogene"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme( axis.text.y=element_text(size=15),axis.text.x=element_text(size=15),
         
         axis.title=element_text(size=17),
         legend.text=element_text(size=14))+
  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border.y = element_blank()
  )+
  ggtitle("Other Pol III transcripts") +
  theme(legend.position=c(.75,.9))


p


ggsave(filename="C:/Vova/study/NGM/Figures/Figure2/ALLPolIII.svg", plot = p,device = "svg", path = NULL,
       scale = 1, width = 150, height = 90, units = c("mm"),
       dpi = 300, limitsize = TRUE)