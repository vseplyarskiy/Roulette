setwd("C:/Vova/study/R_scripts_data/R_data/NGM//")

library(DescTools)
require(ggplot2)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(hrbrthemes)
library(viridis)
setwd("C:/Vova/study/R_scripts_data/R_data/NGM//")

data<-read.table(file="rate_TFBS_R_T_efff_prom_mel", sep=" ")
#data[data$V2=='other',2]<-'no overlap'

#data<-data[data$V7>=50 & data$V3=='C>A',]
#data<-data[data$V7>=50,]
PoissonCI(data[,5], conf.level = 0.95, sides = c("two.sided"),
          method = c("exact"))[,c(2,3)]->CI


p1<-ggplot(data, aes(x=V2, y = V4, ymin = CI[,1]/V6, ymax = CI[,2]/V6,color=V3)) +
  geom_pointrange(size=0.3) +geom_line(size=1.2)+
  
  #geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  #scale_color_manual(values = c("#FFFF66","#CCCC00"), name="") +
  scale_color_manual(values = c("black","grey"), name="") +
  theme_classic() +
  ggtitle("Effect of a promoter on mutation rate") +theme(plot.title = element_text(hjust = .5,size=24))+
  
  scale_y_continuous(trans='log2', breaks =c(1,1.5,2),limits = c(0.9,2))+
  
  
  
    theme(
    #legend.position="none",
    plot.title = element_text(size=17)
  ) +
  
 
  facet_wrap(~V1,strip.position="top",scales = "free")+
  xlab("Distance from TFBS center")+ ylab("Observed/Expected") +
  theme(strip.text.x = element_text(size=17,hjust = .5),
        axis.text.y=element_text(size=15),axis.text.x=element_text(size=13,angle = 90,vjust = .25),axis.title=element_text(size=17),
        legend.text=element_text(size=17))

p1


ggsave(filename="C:/Vova/study/NGM/Figures/Figure5/TFBS_promoterSM_facet.svg", plot = p1,device = "svg", path = NULL,
       scale = 1, width = 200, height = 150, units = c("mm"),
       dpi = 300, limitsize = TRUE)


data<-read.table(file="rate_TFBS_R_T_efff_prom_mel", sep=" ")
#data[data$V2=='other',2]<-'no overlap'

#data<-data[data$V1=='T>G',]
#data<-data[data$V7>=50,]
#median(data[data[,4]==0 & data[,2]=='Testis',5])

PoissonCI(data[,5], conf.level = 0.95, sides = c("two.sided"),
          method = c("exact"))[,c(2,3)]->CI


p1<-ggplot(data, aes(x=V2, y = V4, ymin = CI[,1]/V6, ymax = CI[,2]/V6,color=V3)) +
  geom_pointrange(size=0.3) +geom_line(size=1.2)+
  
  #geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  scale_color_manual(values = c("grey","black"), name="") +
  theme_classic() +
  
  
  ggtitle("TCC>T mutations in melanoma") +theme(plot.title = element_text(hjust = .5,size=24))+
  
  scale_y_continuous(trans='log2',breaks =c(.008,.016),limits = c(0.006,.024))+
  
  # facet_wrap(~V3,strip.position="top",scales = "free")+
  xlab("Distance from TFBS center")+ ylab("mutation rate") +
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


ggsave(filename="C:/Vova/study/NGM/Figures/Figure5/TFBS_promoter_melanoma.svg", plot = p1,device = "svg", path = NULL,
       scale = 1, width = 250, height = 180, units = c("mm"),
       dpi = 300, limitsize = TRUE)