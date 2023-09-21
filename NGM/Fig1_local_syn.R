require(ggplot2)
library(DescTools)
require(PropCIs)
setwd("C:/Vova/study/R_scripts_data/R_data/NGM/")
data<-read.table(file="rate_local", header = F, sep=" ")
data1<-data[data[,1] == 'C>G' | data[,1] == 'G>C' & data[,3]>=10,]
data1<-data1[ data1[,3]>=10,]




p<-ggplot(data1, aes( x=V2, y=V7,color=V8)) + 
  geom_point(size=.5) +
  scale_color_manual(name = "",values=(c("cyan","orange")))+
geom_smooth(method = 'loess',n=10,size=1, se=FALSE)+
  # scale_alpha_manual(name="strand",values = c(1,.3))+
  ylab("mutation rate")+
  xlab("coordinate")+
  theme_classic() +

  #theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1))+
  ggtitle("C>G rate along Chr. 8")+theme(plot.title = element_text(hjust = 0.5))+
  #theme(legend.position = "none") + 
  theme(plot.title = element_text(hjust = 0.5))+
  theme( axis.text.y=element_text(size=15),axis.text.x=element_text(size=15),axis.title=element_text(size=17),
         legend.text=element_text(size=14))+
  #theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + theme(legend.position=c(.75,.9))#theme(legend.justification=c(1,1), legend.position=c(1,1))

p
ggsave(filename="C:/Vova/study/NGM/Figures/Figure1/Local_variation.svg", plot = p,device = "svg", path = NULL,
       scale = 1, width = 160, height = 90, units = c("mm"),
       dpi = 300, limitsize = TRUE) 