require(ggplot2)
require(epitools)

setwd("C:/Vova/study/R_scripts_data/R_data/NGM//")

data<-read.table(file="extended_10kb_C_G_4",sep=" ")

dataGCobs<-mean(data[data$V2 == 'G_C' & data$V4 =='observed',3],na.rm = T)
dataGCfull<-mean(data[data$V2 == 'G_C' & data$V4 =='full_model',3],na.rm = T)
dataGCtri<-mean(data[data$V2 == 'G_C' & data$V4 =='trinucleotide',3],na.rm = T)




data<-read.table(file="extended_10kb_C_G_8",sep=" ")
dataGC<-data[data$V2 == 'G_C',]
is.na(dataGC[,3])<-0
dataGC[dataGC$V4 == 'trinucleotide',3]<-dataGC[dataGC$V4 =='trinucleotide',3]*(dataGCobs/dataGCtri)
dataGC[dataGC$V4 == 'full_model',3]<-dataGC[dataGC$V4 =='full_model',3]*(dataGCobs/dataGCfull)


p<-ggplot(dataGC,aes(x=V1, y=V3,color=as.factor(V4)))+geom_point(size=.4,alpha=0.3)+#geom_point(size=0.8,alpha=1/3)+#geom_smooth()+
  geom_smooth(span = 0.07,method = 'loess',n=100,size=1, se=FALSE)+
  ylab("SNVs in 10kb window") + 
  xlab("position,chr8") + 
  theme_bw()+xlim(c(1000000,30000000))+
  labs(color = "")+
  theme_bw()+
  #annotate("rect", xmin = x1, xmax = x2, ymin = y1, ymax = y2,alpha = .5,fill=str)+
  #geom_rect(data6, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=str), color="black", alpha=0.5)
  #geom_segment(aes(x = genel, y = -.2, xend = gener, yend = -.2,size=1.5, color='IC5'))+guides(size=FALSE)+
  #ylim(c(-0.1*yma,yma))+
  #scale_x_continuous(breaks=c(20000000, 22500000, 25000000), limits = c(20000000, 25600000))+
  #scale_x_continuous(breaks=c(75000000, 79000000, 83000000),limits = c(74000000,84000000))+
  #scale_color_manual(values=c('Sig9','Sig9'))+
  #ggtitle("FHIT")+theme(plot.title = element_text(hjust = 0.5))+
  theme( axis.text.y=element_text(size=12),axis.text.x=element_text(size=12),axis.title=element_text(size=15),
         legend.text=element_text(size=10))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())#+ theme(legend.position = "none")#+scale_x_continuous(breaks=c(20000000,22500000,25000000))#+
#geom_point(data=data1,colour="black", size=3.7, alpha=0.05)+geom_point(data=data3, colour="red",size=3.7,alpha=0.1)
p

x<-data[data$V2 == 'G_C' & data$V4 =='trinucleotide',]



