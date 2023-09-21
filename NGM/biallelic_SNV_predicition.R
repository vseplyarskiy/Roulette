library(ggplot2)

setwd("C:/Vova/study/R_scripts_data/R_data/NGM/")
#data<-read.table(file="distance_low_FDR", header = F, sep=" ")
data<-read.table(file="rate_tri", header = F, sep=" ")
#data<-read.table(file="SINE_summ", header = F, sep=" ")

data<-data[data$V1*data$V3>300,]
p<-ggplot(data, aes(x=data[,1], y=data[,2])) + geom_point() + 
  
  ylab("rate of triallelic SNVs")+xlab("expectation under independence")+
  geom_abline(intercept = 0, slope = 1, color="gray", 
              linetype="dashed", size=1)+
  theme_classic()+
  scale_color_manual(name = "",#values=c("red","black"),
                     #labels=c("active gene", "pseudogene")
                     )+
  theme(legend.text = element_text(size = 12),axis.title =element_text(size = 16),
        axis.text  =element_text(size = 12),
        legend.title = element_text( size=12, face="bold"))+
  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border.y = element_blank()
  )+
  ggtitle("No cryptic variation in mutation rate") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position=c(.75,.9))
p

ggsave(filename="C:/Vova/study/NGM/Figures/Figure2/No_cryptic.svg", plot = p,device = "svg", path = NULL,
       scale = 1, width = 150, height = 90, units = c("mm"),
       dpi = 300, limitsize = TRUE)
