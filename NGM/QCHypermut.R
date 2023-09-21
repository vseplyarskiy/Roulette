library(ggplot2)
library(jtools) # Load jtools
library(viridis)

setwd("C:/Vova/study/R_scripts_data/R_data/NGM");
data<-read.table(file="rate_qc_control",  col.names = paste0("V",seq_len(8)), fill = TRUE)

#name='AS_pab_max'
name='AS_VQLSOD'

p1<-ggplot(data, aes(x=V2, color=V1)) +
  geom_density(size=1) +
  
  #geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  scale_color_discrete( name="") +
  theme_classic() +
  theme(
    #legend.position="none",
    plot.title = element_text(size=17)
  ) +
  
  #ggtitle("Hypermutable regions, AS_pab_max") + 
 # ggtitle("VQLSOD") + 
  
  theme(plot.title = element_text(hjust = .5))+
  xlab(name)+
xlim (c(-11,11))+
  theme(strip.text.x = element_text(size=17,hjust = .5),
        axis.text.y=element_text(size=15),axis.text.x=element_text(size=13,angle = 90,vjust = .25),
        axis.title=element_text(size=17),
        legend.text=element_text(size=17),
        legend.title = element_text(size=19))
#axis.title.y = element_text(vjust = .5),
#    panel.background = element_rect(fill = "grey95",color = NA)

p1
fi<-paste("C:/Vova/study/NGM/Figures/Supplementary_Fig/","QC",name,".svg",sep ="")

ggsave(filename=fi, plot = p1,device = "svg", path = NULL,
       scale = 1, width = 215, height = 85, units = c("mm"),
       dpi = 300, limitsize = TRUE)