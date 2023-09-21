require(ggplot2)
require(reshape2)
require(TTR)

setwd("C:/Vova/study/R_scripts_data/R_data/NGM/")
#data<-read.table(file="distance_low_FDR", header = F, sep=" ")
data<-read.table(file="hash_100_2", header = T, sep=" ")

data5<-read.table(file="IGK_coordinates", header = F, sep=" ")

chr=1
genel=2938099
gener=6994972

genelb=88830510 #FHIT 59749310 #CSMD1  0938099 
generb=88959370 #FHIT  61251450 #CSMD1  6994972


df1<-cbind(data[,c(1,26)],
           rep('Expected',length(data[,c(26)])))
names(df1)<-(c('cord','rate','Signature'))  
df2<-cbind(data[,c(1,27)],
           rep('Observed',length(data[,c(27)])))
names(df2)<-(c('cord','rate','Signature'))    


df<-rbind(df1,df2)

#chr1<-paste("Chromosome", chr)

n<-length(df$rate)
yma=sort(df$rate)[n]
ymi=0-0.4*yma
#ymi=sort(df$rate)[1]
size=yma-ymi

x1<-data5[,2]
x2<-data5[,3]
n1<-length(x1)
y1<-rep(-10,n1)
y2<-rep(0,n1)

data6<-as.data.frame(cbind(x1,x2,y1,y2))  
#data6$str<-data5$strand[data5$chr == as.character(chr)]


#colors(as.factor(str))
p<-ggplot(df)+geom_point(aes(x=cord, y=rate,color=as.factor(Signature)),size=1.9,alpha=0.25)+
  geom_point(data5,mapping=aes(y=rep(-5,length(data5[,2])),x=(rowSums(data5[,c(2,3)])/2)),size=4.4,alpha=0.8)+
  #geom_line(aes(y=SMA(df$RT_der,100)),color="red",size=0.8, se=FALSE,alpha=0.9)+
  #geom_smooth(span = 0.07,method = 'loess',n=100,size=1, se=FALSE)+
  ylab("Rare SNV in 100 nuc. window") + 
  xlab("Coordinate, Chr. 2") + theme_classic()+xlim(c(genelb,generb))+
  scale_color_manual(name = "",values=(c("orange","cyan")))+
    ylim(c(-10.5,yma))+
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
  ggtitle("Imunoglobins are mutational hotspots")
  p

  ggsave(filename="C:/Vova/study/NGM/Figures/Figure2/IGK.svg", plot = p,device = "svg", path = NULL,
         scale = 1, width = 150, height = 90, units = c("mm"),
         dpi = 300, limitsize = TRUE)
#+  geom_point(data=data3,colour="blue", size=1.1,shape="star")



