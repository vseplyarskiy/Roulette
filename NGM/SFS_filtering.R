require(ggplot2)

setwd("C:/Vova/study/R_scripts_data/R_data/NGM//")


require(ggplot2)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(hrbrthemes)
library(viridis)
data<-read.table(file="SFS_problem_visual_v2", sep=" ")
#data[data$V2=='other',2]<-'no overlap'

#data<-data[data$V7>=50 & data$V3=='C>A',]
#data<-data[data$V7>=50,]



p1<-ggplot(data, aes(x=V4,  color=V11))+#,y = stat(density))) +
  #geom_line(stat = "density", size=1.5) +
  stat_ecdf(size=1.8,alpha=0.6)+
  #geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  scale_color_viridis(discrete = TRUE, name="Filtering") +
  theme_classic() +
  theme(
    #legend.position="none",
    plot.title = element_text(size=17)
  ) +
  ggtitle("SFS distortion") +theme(plot.title = element_text(hjust = .5))+
  
 # scale_y_continuous(trans='log2', breaks =c(0.5,1,2),limits = c(0.5,2.5))+
  
  facet_wrap(~V2,strip.position="top",scales = "free")+
  scale_x_continuous(limits=c(0,.3))+
  xlab("Fraction of high frequency variants")+ 
  ylab("CDF") +
  theme(strip.text.x = element_text(size=17,hjust = .5),
        axis.text.x=element_text(size=13,angle = 90,vjust = .25),axis.title=element_text(size=17),
        legend.text=element_text(size=17),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        #axis.title.y = element_text(vjust = .5),
        #panel.background = element_rect(fill = "grey95",color = NA)
        )

p1


ggsave(filename="C:/Vova/study/NGM/Figures/Supplementary_Fig/SFS_facet.svg", plot = p1,device = "svg", path = NULL,
       scale = 1, width = 480, height = 280, units = c("mm"),
       dpi = 300, limitsize = TRUE)