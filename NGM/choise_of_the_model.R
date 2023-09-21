library(ggplot2)
library(jtools) # Load jtools
library(viridis)

setwd("C:/Vova/study/R_scripts_data/R_data/NGM");
data<-read.table(file="models_distribution", sep="\t")

name='AS_pab_max'
data$V1 <- factor(data$V1, levels = unique(data$V1))
  data<-cbind(data,c("Intercept","log-link", "logistic no interaction", "logistic pairwise"))
names(data)[3]<-'V3'
data$V3 <- factor(data$V3, levels = unique(data$V3))

p1<-ggplot(data, aes(x=V3, y=V2)) +
  geom_bar(position="dodge", stat="identity") +
  #geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  #scale_color_discrete( name="") +
  theme_classic() +
  theme(
    #legend.position="none",
    plot.title = element_text(size=17)
  ) +
  
  #ggtitle("Hypermutable regions, AS_pab_max") +  theme(plot.title = element_text(hjust = .5))+
  xlab("Model")+ylab("Fraction of genome")+

theme(strip.text.x = element_text(size=17,hjust = .5),
      axis.text.y=element_text(size=17),axis.text.x=element_text(size=10),#,angle = 90,vjust = .25),
      axis.title=element_text(size=17),
      legend.text=element_text(size=17),
      legend.title = element_text(size=19))
#axis.title.y = element_text(vjust = .5),
#    panel.background = element_rect(fill = "grey95",color = NA)

p1
fi<-paste("C:/Vova/study/NGM/Figures/Supplementary_Fig/","Model_selection",".svg",sep ="")

ggsave(filename=fi, plot = p1,device = "svg", path = NULL,
       scale = 1, width = 215, height = 85, units = c("mm"),
       dpi = 300, limitsize = TRUE)