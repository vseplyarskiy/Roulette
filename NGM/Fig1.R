require(ggplot2)
require(reshape2)
require(TTR)
library(flexmix)
setwd("C:/Vova/study/R_scripts_data/R_data/NGM/")
#data<-read.table(file="distance_low_FDR", header = F, sep=" ")
data<-read.table(file="decile_denovo", header = F, sep=" ")
p<-ggplot(data, aes(x=V2)) + 
  geom_histogram(color="black", fill="white",binwidth=.07)+
  scale_x_continuous(trans='log2')+
  xlab("mutation rate in mutable vs non-mutable sites")+
  ggtitle("Mutation rate variation within pentamers") +
 
  theme_classic()+
  theme(axis.text.y=element_text(size=15),axis.text.x=element_text(size=15),axis.title=element_text(size=17),
legend.text=element_text(size=17),plot.title = element_text(hjust = .5))
    

p
ggsave(filename="C:/Vova/study/NGM/Figures/Figure1/Decile_ratio_denovo.svg", plot = p,device = "svg", path = NULL,
       scale = 1, width = 160, height = 90, units = c("mm"),
       dpi = 300, limitsize = TRUE)

data<-read.table(file="maternal_denovo", header = F, sep="\t")

data$V2 <- factor(data$V2, levels = unique(data$V2))
p<-ggplot(data,aes(x=V1, y=V3*(10^8/6232),fill=V2))+ geom_bar(stat="identity",  width = 0.8,position = "dodge")+#geom_smooth()+
  #position_dodge2(width = 0.8,padding = 0.3, reverse = TRUE)+
  #geom_smooth(span = 0.1,method = 'loess',n=2000,size=2.1, se=FALSE)+
  # ylab("mutation rate") + 
  scale_fill_manual(name = "", values=(c("red","pink","black","grey")))+ 
                    #labels=c("fraction of clusters within maternal regions","fraction of the the maternal regions"))+
  ylab("Mut. rate per site x 10^-8") + 
  xlab("mutation type") + theme_bw()+#xlim(c(genelb,generb))+
  # geom_segment(aes(x = genel, y = 0, xend = gener, yend = 0,size=1.5, color='Sig5'))+guides(size=FALSE)+
  #ylim(c(0,1.05))+
  #scale_color_manual(values=c('Sig9','Sig9'))+
  #ggtitle(chr1)+
  #scale_y_continuous(trans='log2')+
  theme(plot.title = element_text(hjust = 0.5))+
  theme( axis.text.y=element_text(size=15),axis.text.x=element_text(size=15),axis.title=element_text(size=17),
         legend.text=element_text(size=14))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + theme(legend.position=c(.75,.9))#theme(legend.justification=c(1,1), legend.position=c(1,1))

p
ggsave(filename="C:/Vova/study/NGM/Figures/Supplementary_Fig//maternal_regions_denovo.svg", plot = p,device = "svg", path = NULL,
       scale = 1, width = 270, height = 120, units = c("mm"),
       dpi = 300, limitsize = TRUE)

