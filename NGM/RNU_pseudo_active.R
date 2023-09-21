library(ggplot2)

setwd("C:/Vova/study/R_scripts_data/R_data/NGM/")
#data<-read.table(file="distance_low_FDR", header = F, sep=" ")
data<-read.table(file="rate_RNU_ind", header = F, sep=" ")
#data<-read.table(file="SINE_summ", header = F, sep=" ")


p<-ggplot(data, aes(x=data[,4], y=data[,5], color=V2)) + geom_point() + 
  
  ylab("observed")+xlab("expected")+
  geom_abline(intercept = 0, slope = 1, color="gray", 
              linetype="dashed", size=1)+
  theme_classic()+
  scale_color_manual(name = "",values=c("red","black"),
  labels=c("active gene", "pseudogene"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme( axis.text.y=element_text(size=15),axis.text.x=element_text(size=15),axis.title=element_text(size=17),
         legend.text=element_text(size=14))+
  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border.y = element_blank()
  )+
  ggtitle("Effect of active transcription") +
 theme(legend.position=c(.75,.9))
p

ggsave(filename="C:/Vova/study/NGM/Figures/Figure2/RNU_pseudo.svg", plot = p,device = "svg", path = NULL,
       scale = 1, width = 150, height = 90, units = c("mm"),
       dpi = 300, limitsize = TRUE)
