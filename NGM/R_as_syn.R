require(ggplot2)
library(DescTools)
require(PropCIs)
setwd("C:/Vova/study/R_scripts_data/R_data/NGM/")
data<-read.table(file="rate_asRT2", header = F, sep=" ")
#asym<-abs(log(data[,12]))
asym<-abs(log(data[,2]/data[,7]))

data<-cbind(data,asym)
data<-cbind(data,(data[,c(2,3,7,8)]/apply(data[,c(2,7)], 1, FUN = min)))
newdata <- data[order(-data[,13]),]
#newdata <- data[order(-abs(log(data[,2]/data[,7]))),]
#newdata[c(1:5),13]
res1<-cbind(newdata[,c(1,4,14,13)],'lagging','Predicted')
res2<-cbind(newdata[,c(1,5,15,13)],'lagging','Observed')
res3<-cbind(newdata[,c(1,9,16,13)],'leading','Predicted')
res4<-cbind(newdata[,c(1,10,17,13)],'leading','Observed')
colnames(res1)<-c('context','mut','norm','asym','strand','type')
colnames(res2)<-c('context','mut','norm','asym','strand','type')
colnames(res3)<-c('context','mut','norm','asym','strand','type')
colnames(res4)<-c('context','mut','norm','asym','strand','type')
rbind(res1,res2,res3,res4)->res
res[,2]<-as.integer(res[,2])
res <- res[order(-res[,4]),]
res<-res[!(res[,1] =='GCT>G'),]
CI<-PoissonCI(res[,2], n = 1, conf.level = 0.95, sides = c("two.sided"))
CI<-CI*(res[,3]/res[,2])

ress<-res[c(1:16),]
ress$type <- factor(ress$type, levels = unique(ress$type))
CI1<-CI[c(1:16),]
p<-ggplot(ress, aes( x=context, y=norm,color=type,shape=strand)) + 
  geom_point(position=position_dodge2(.1),size=.2) +
  scale_color_manual( name="",values = c("red","black")) +
  # scale_alpha_manual(name="strand",values = c(1,.3))+
  ylab("normalised mutation rate")+
  #ylab("rate")+
  theme_classic() +
  geom_pointrange( aes( ymin=CI1[,2], ymax=CI1[,3]),stat="identity",
                   position=position_dodge2(.1))+
  #theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1))+
  ggtitle("Replicational asymmetry")+theme(plot.title = element_text(hjust = 0.5))+
  #theme(legend.position = "none") + 
  theme(plot.title = element_text(hjust = 0.5))+
  theme( axis.text.y=element_text(size=15),axis.text.x=element_text(size=15),axis.title=element_text(size=17),
         legend.text=element_text(size=14))+
  theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())# + theme(legend.position=c(.75,.9))#theme(legend.justification=c(1,1), legend.position=c(1,1))

p
ggsave(filename="C:/Vova/study/NGM/Figures/Supplementary_Fig/Replicational_assymmetry.svg", plot = p,device = "svg", path = NULL,
       scale = 1, width = 160, height = 90, units = c("mm"),
       dpi = 300, limitsize = TRUE)
#xlab("Regions")+  facet_wrap( ~ V2, nrow = 1)+
#theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1))+
#scale_y_continuous(trans='log2', breaks =c(1,2,4),limits = c(0.9,4))
