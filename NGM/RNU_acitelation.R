library(ggplot2)

setwd("C:/Vova/study/R_scripts_data/R_data/NGM/")
#data<-read.table(file="distance_low_FDR", header = F, sep=" ")
data<-read.table(file="RNU_Acetelation_rate", header = F, sep=" ")
#data<-read.table(file="SINE_summ", header = F, sep=" ")

data[,6]+1->data[,6]
data1<-data[data[,2] == 'pseudogene' ,]
p<-ggplot(data, aes(y=data[,6], x=data[,5]/data[,4], color=V2)) + 
  geom_hex() + scale_fill_viridis_c(trans = "log10")+
  ylab("H3k27ac")+xlab("observed/expected")+
  #geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed", size=1)+
  theme_classic()+
  scale_color_manual(values=c("red", "green"),name="",labels  = c("acitve gene","pseudogene"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme( axis.text.y=element_text(size=15),axis.text.x=element_text(size=15),axis.title=element_text(size=17),
         legend.text=element_text(size=14))+
  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border.y = element_blank()
  )+
  ggtitle("RNU acetylation")
p

ggsave(filename="C:/Vova/study/NGM/Figures/Supplementary_Fig/RNU_H3K27ac.svg", plot = p,device = "svg", path = NULL,
       scale = 1, width = 150, height = 90, units = c("mm"),
       dpi = 300, limitsize = TRUE)
