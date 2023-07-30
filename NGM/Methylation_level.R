library(ggplot2)
setwd("C:/Vova/Vova/study/R_scripts_data/R_data/NGM/");
#table<-read.table(file="5_visual_8_7_anacor_3mut_anacor_3mutr_3mut_3mut",sep=" ");
table<-read.table(file="methylation_bins",sep=" ");
#table<-read.table(file="DNM_visual",sep=" ");

ggplot(table, aes(x = V2, fill = V2,y=V3)) +
  geom_bar(position="dodge", stat="identity")+  xlab("Mutation type") +
  ylab("Fraction") +  
  scale_fill_discrete(name="Orientation")+
  facet_wrap(~V1)+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = .5, 
                                   size = 10, hjust =0))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+theme(plot.title = element_text(hjust = 0.5))
