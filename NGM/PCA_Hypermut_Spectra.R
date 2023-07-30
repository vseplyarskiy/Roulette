require(ggplot2)
require(reshape2)
setwd("C:/Vova/study/R_scripts_data/R_data/NGM/")
data<-read.table(file="for_PCA", header = F, sep=" ")
data<-read.table(file="for_PCA_oe", header = F, sep=" ")


res.pca <- prcomp(data[,c(-1,-8)], scale = TRUE)
summary(res.pca)
as.data.frame(res.pca[["x"]])->coods


p<-ggplot(coods,aes(x=PC1,y=PC2,color=as.factor(data[,1])))+geom_point(asize=5.9,alpha=0.45)+
  #geom_point(data5,mapping=aes(y=rep(-5,length(data5[,2])),x=(rowSums(data5[,c(2,3)])/2)),size=4.4,alpha=0.8)+
  #geom_line(aes(y=SMA(df$RT_der,100)),color="red",size=0.8, se=FALSE,alpha=0.9)+
  #geom_smooth(span = 0.07,method = 'loess',n=100,size=1, se=FALSE)+
  #ylab("Rare SNV in 100 nuc. window") + 
  #xlab("Coordinate, Chr. 2") +
  theme_classic()+#xlim(c(genelb,generb))+
  scale_color_viridis_d(name = "")+
  #ylim(c(-10.5,yma))+
  # scale_x_continuous(breaks=c(20000000, 22500000, 25000000), limits = c(20000000, 25600000))+
  
  #scale_color_manual(values=c('Sig9','Sig9'))+
  # ggtitle(chr1)+theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme( axis.text.y=element_text(size=15),axis.text.x=element_text(size=15),axis.title=element_text(size=17),
         legend.text=element_text(size=14))+
  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border.y = element_blank()
  )+
  ggtitle("Mutational spectra of hotspots")
p

ggsave(filename="C:/Vova/study/NGM/Figures/Figure2/PCA_spectra_PC1_2_NORM.svg", plot = p,device = "svg", path = NULL,
       scale = 1, width = 150, height = 90, units = c("mm"),
       dpi = 300, limitsize = TRUE)

