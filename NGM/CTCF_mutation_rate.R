setwd("C:/Vova/study/R_scripts_data/R_data/NGM/")


data<-read.table(file="CTCF_rate", sep=" ")

p1<-ggplot(data, aes(x=V1, y=V3, color=V2)) +
  geom_point(size=0.3) +
  geom_smooth(span=0.1,se=F) +
  #geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  scale_color_discrete( name="Mutation\ntype") +
  theme_classic() +
  
  
  ggtitle("CTCF") +theme(plot.title = element_text(hjust = .5,size=24))+
  
  #scale_y_continuous(trans='log2', breaks =c(0.5,1,2),limits = c(0.5,2.5))+
  
  # facet_wrap(~V3,strip.position="top",scales = "free")+
  xlab("Distance from CTCF center")+ ylab("Observed/Expected") +
  theme(#strip.text.x = element_text(size=20,hjust = .5),
    axis.text.y=element_text(size=20),axis.text.x=element_text(size=20,angle = 90,vjust = .25),axis.title=element_text(size=20),
    legend.text=element_text(size=13),
    legend.title = element_text(size=15), 
    
    #axis.title.y = element_text(vjust = .5),
    )+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border.y = element_blank()
  )

p1


ggsave(filename="C:/Vova/study/NGM/Figures/Figure2/CTCF_v2.svg", plot = p1,device = "svg", path = NULL,
       scale = 1, width = 150, height = 90, units = c("mm"),
       dpi = 300, limitsize = TRUE)


data<-read.table(file="rate_ench", sep=" ")

p1<-ggplot(data, aes(x=V1, y=V2, fill=V1)) +
  geom_bar(stat="identity") +
#  geom_smooth(span=0.1,se=F) +
  #geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  scale_fill_viridis(discrete = T, name="") +
  theme_classic() +
  
  
  ggtitle("ENCODE features") +theme(plot.title = element_text(hjust = .5,size=24))+
  
  scale_y_continuous( breaks =c(0,0.5,1))+
  
  # facet_wrap(~V3,strip.position="top",scales = "free")+
  xlab("ENCODE annotation")+ ylab("Observed/Expected") +
  theme(#strip.text.x = element_text(size=20,hjust = .5),
    axis.text.y=element_text(size=20),
    axis.text.x=element_blank(),
    axis.title=element_text(size=20),
    legend.text=element_text(size=13),
    legend.title = element_text(size=15), 
    
    #axis.title.x = element_blank()#element_text(vjust = .5),
  )

p1
ggsave(filename="C:/Vova/study/NGM/Figures/Supplementary_Fig//Encode.svg", plot = p1,device = "svg", path = NULL,
       scale = 1, width = 150, height = 90, units = c("mm"),
       dpi = 300, limitsize = TRUE)