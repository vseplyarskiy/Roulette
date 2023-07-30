require(ggplot2)

setwd("C:/Vova/study/R_scripts_data/R_data/NGM//")


require(ggplot2)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(hrbrthemes)
library(viridis)
data<-read.table(file="rate_mut_v2", sep=" ")
data[data$V2=='other',2]<-'no overlap'

#data<-data[data$V7>=50 & data$V3=='C>A',]
data<-data[data$V7>=50,]



p1<-ggplot(data, aes(as.factor(V4*20+50), y=V5, fill=V2)) +
  geom_boxplot(width=0.8) +
  
  #geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  scale_fill_viridis(discrete = TRUE, name="Tissue specificity") +
  theme_classic() +
  theme(
    #legend.position="none",
    plot.title = element_text(size=17)
  ) +
  
  ggtitle("mutations at active TFBS") +theme(plot.title = element_text(hjust = .5))+
  
  scale_y_continuous(trans='log2', breaks =c(0.5,1,2),limits = c(0.5,2.5))+

  facet_wrap(~V3,strip.position="top",scales = "free")+
  xlab("Distance from TFBS center")+ ylab("Observed/Expected") +
  theme(strip.text.x = element_text(size=17,hjust = .5),
        axis.text.y=element_text(size=15),axis.text.x=element_text(size=13,angle = 90,vjust = .25),axis.title=element_text(size=17),
        legend.text=element_text(size=17))

p1


ggsave(filename="C:/Vova/study/NGM/Figures/Figure2/TFBS_facet.svg", plot = p1,device = "svg", path = NULL,
       scale = 1, width = 480, height = 280, units = c("mm"),
       dpi = 300, limitsize = TRUE)


data<-read.table(file="rate_mut_v2", sep=" ")
data[data$V2=='other',2]<-'no overlap'

data<-data[data$V7>=50 & data$V3=='T>G',]
#data<-data[data$V7>=50,]
median(data[data[,4]==0 & data[,2]=='Testis',5])


p1<-ggplot(data, aes(as.factor(V4*20+50), y=V5, fill=V2)) +
  geom_boxplot(width=0.8,outlier.alpha = 0.1) +
  
  #geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  scale_fill_viridis(discrete = TRUE, name="Tissue specificity") +
  theme_classic() +

  
  ggtitle("T>G rate at TFBS") +theme(plot.title = element_text(hjust = .5,size=24))+
  
  scale_y_continuous(trans='log2', breaks =c(0.5,1,2),limits = c(0.5,2.5))+
  
 # facet_wrap(~V3,strip.position="top",scales = "free")+
  xlab("Distance from TFBS center")+ ylab("Observed/Expected") +
  theme(#strip.text.x = element_text(size=20,hjust = .5),
        axis.text.y=element_text(size=20),axis.text.x=element_text(size=20,angle = 90,vjust = .25),axis.title=element_text(size=20),
        legend.text=element_text(size=13),
        legend.title = element_text(size=15)        #axis.title.y = element_text(vjust = .5),
)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border.y = element_blank()
  )

p1


ggsave(filename="C:/Vova/study/NGM/Figures/Figure2/TFBS_T_G.svg", plot = p1,device = "svg", path = NULL,
       scale = 1, width = 200, height = 150, units = c("mm"),
       dpi = 300, limitsize = TRUE)