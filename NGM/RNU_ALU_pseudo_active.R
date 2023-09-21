library(ggplot2)
library(svglite)
setwd("C:/Vova/study/R_scripts_data/R_data/NGM/")
#data<-read.table(file="distance_low_FDR", header = F, sep=" ")
data<-read.table(file="rate_RNU_ind", header = F, sep=" ")
#data<-read.table(file="SINE_summ", header = F, sep=" ")

data2<-read.table(file="RNU_Acetelation_rate", header = F, sep=" ")
merge(data, data2, by = "V1", all.x=TRUE, all.y=TRUE)->datam
quantile(datam[datam[,2]=='pseudogene',10],0.99)->q
datam[datam[,2]=='pseudogene' & datam[,10]>q,2]<-'pseudogene + H3K27ac'
p<-ggplot(datam, aes(x=V4.x, y=V5.x, color=V2.x)) + geom_point() + 
  
  ylab("observed")+xlab("expected")+
  geom_abline(intercept = 0, slope = 1, color="gray", 
              linetype="dashed", size=1)+
  theme_classic()+
  scale_color_manual(name = "",values=c("red","black","blue"),
  labels=c("active gene", "pseudogene",'pseudogene + H3K27ac'))+
 
  theme(plot.title = element_text(hjust = 0.5))+
  theme( axis.text.y=element_text(size=15),axis.text.x=element_text(size=15),axis.title=element_text(size=17),
         legend.text=element_text(size=14))+
  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border.y = element_blank()
  )+
  ggtitle("Effect of active transcription") +
 theme(legend.position=c(.75,.9))
p

ggsave(filename="C:/Vova/study/NGM/Figures/Supplementary_Fig/RNU_pseudo_v2.svg", plot = p,device = "svg", path = NULL,
       scale = 1, width = 150, height = 90, units = c("mm"),
       dpi = 300, limitsize = TRUE)

setwd("C:/Vova/study/R_scripts_data/R_data/NGM/")
data<-read.table(file="RNU_pseudo_genes_barplot", header = F, sep=" ")


p<-ggplot(data,aes(x=V6, y=data$V5/data$V4,color=V2,fill=V2,size=3))+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=35,binwidth=.005, 
               position=position_dodge(0.8))+
  scale_fill_manual(values=c("red","black"),name='',
                     labels=c("active gene", "pseudogene"))+
  scale_color_manual(values=c("red","black"),name='',labels=c("", "") )+
  #scale_color_manual(name = "Mutation type", values =colorConverter(data$V2))+
  ylab('Obs/Exp') + guides(colour=FALSE)+
  xlab("Gene type") + theme_classic()+#xlim(c(genelb,generb))+
  #theme(plot.title = element_text(hjust = 0.5))+
  theme( axis.text.x=element_text(size=11),axis.text.y=element_text(size=25),axis.title=element_text(size=25),
         legend.text=element_text(size=25))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
  p



  ggsave(filename="C:/Vova/study/NGM/Figures/Supplementary_Fig/RNU_pseudo_barplot.svg", plot = p,device = "svg", path = NULL,
         scale = 1, width = 400, height = 150, units = c("mm"),
         dpi = 300, limitsize = TRUE)



library(ggplot2)

setwd("C:/Vova/study/R_scripts_data/R_data/NGM/")

data<-read.table(file="SINE_summ", header = F, sep=" ")


p<-ggplot(data, aes(x=V4, y=V3)) + geom_point() + 
  
  ylab("observed")+xlab("expected")+
  geom_abline(intercept = 0, slope = 1, color="gray", 
              linetype="dashed", size=1)+
  theme_classic()+
 # scale_color_manual(name = "",values=c("red","black"),
 #                    labels=c("active gene", "pseudogene"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme( axis.text.y=element_text(size=15),axis.text.x=element_text(size=15),axis.title=element_text(size=17),
         legend.text=element_text(size=14))+
  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border.y = element_blank()
  )+
  ggtitle("ALU repeats transcribed by Polymerase III") +
  theme(legend.position=c(.75,.9))
p

ggsave(filename="C:/Vova/study/NGM/Figures/Supplementary_Fig/ALU.svg", plot = p,device = "svg", path = NULL,
       scale = 1, width = 150, height = 90, units = c("mm"),
       dpi = 300, limitsize = TRUE)