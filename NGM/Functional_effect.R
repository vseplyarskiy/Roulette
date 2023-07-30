require(ggplot2)
library(scales)


setwd("C:/Vova/Vova/study/R_scripts_data/R_data/NGM/")
data<-read.table(file="Fig_values",  sep=" ")
data1<-read.table(file="Fig_values_line",  sep=" ")

cor.test(data1$V3,data1$V6, method="spearman")
ggplot(data,aes(x=V1, color=V2,y=V3))+ geom_point()+ 
  #geom_violin(trim=FALSE, fill="gray", adjust=2)+
  scale_y_continuous(trans = log2_trans())+
 # scale_y_continuous(trans = log2_trans(),limits=c(0.5, 2))+
  
  #facet_wrap(~V8)+
  #geom_errorbar(aes(ymin=lower_timing, ymax=upper_timing),stat="identity", position=position_dodge())+
  
  #geom_point(alpha=0.1, size=0.4)+
  
  #scale_color_manual(name = "Mutation type", values =colorConverter(data$V2))+
  ylab("Activity") + 
  xlab("") + theme_bw()+#xlim(c(genelb,generb))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme( axis.text.y=element_text(size=11),axis.text.x=element_text(size=10,angle = 90,vjust = .25),axis.title=element_text(size=11),
         legend.text=element_text(size=11))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+ #theme(legend.position = "none")
  labs(fill = "Origin")