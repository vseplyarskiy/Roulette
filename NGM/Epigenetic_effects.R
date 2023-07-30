library(ggplot2)
library(lattice)
library(hrbrthemes)
library(viridis)
setwd("C:/Vova/study/R_scripts_data/R_data/NGM/")

data<-read.table(file="Epi_Figure", sep=" ")




p1<-ggplot(data, aes(x=V1, y=V3, fill=V4)) +
  geom_bar(position="dodge", stat="identity") +
  
  #geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  scale_fill_viridis(discrete = TRUE, name="Tissue specificity") +
  theme_classic() +
  theme(
    #legend.position="none",
    plot.title = element_text(size=17)
  ) +
  
  ggtitle("Tracing epigenetic features") +theme(plot.title = element_text(hjust = .5))+
  
  scale_y_continuous(trans='log2', breaks =c(0.5,1,2),limits = c(0.42,2.1))+
  
  facet_wrap(~V2,strip.position="top",scales = "free")+
  xlab("Epigenetic bin")+ ylab("Normalized mutation rate") +
  theme(strip.text.x = element_text(size=17,hjust = .5),
        axis.text.y=element_text(size=15),axis.text.x=element_text(size=13),axis.title=element_text(size=17),
        legend.text=element_text(size=17))
        #axis.title.y = element_text(vjust = .5),
    #    panel.background = element_rect(fill = "grey95",color = NA)

p1


ggsave(filename="C:/Vova/study/NGM/Figures/Supplementary_Fig/Epigenetic_effects.svg", plot = p1,device = "svg", path = NULL,
       scale = 1, width = 480, height = 280, units = c("mm"),
       dpi = 300, limitsize = TRUE)