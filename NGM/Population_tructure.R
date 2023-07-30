library(ggplot2)
setwd("C:/Vova/study/R_scripts_data/R_data/NGM/")



data<-read.table(file="pop_vector", sep=" ", header = T);
#data<-data[data[,1] != 'GSE30267.ELL2.hct116_starved',]
p<-ggplot(data,aes(x = rate,color=mut,y=data[,6]+data[,7]))+#scale_x_discrete(limits = positions)+#geom_bar()
  geom_point()+ #scale_size(range = c(1, 1))+
  theme_bw()+
 # xlim(c(0.4,1.2))
 # ylab("Density") +
#  scale_fill_discrete(name="Mutation\nType") +
  #scale_x_log10(breaks=x_breaks)
  #scale_x_continuous(trans = 'log2', breaks =  c(1,2.5,5,10,25,50,100,200),lim=c(1,200))+
  #scale_x_continuous(trans='log2') +
  #scale_y_continuous(trans='log2')+
  xlab("Mutation rate") +
  ylab("fraction of SNPs in more than two populations")+
  theme(legend.text = element_text(size = 15),legend.title = element_text( size=15, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(  axis.text=element_text(size = 15),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())#+theme(legend.position = "none")
p

data<-read.table(file="rate_filter", sep=" ", header = F);
#data<-data[data[,1] != 'GSE30267.ELL2.hct116_starved',]
data1<-data[(data[,1] =='T>C' | data[,1] =='C>T' | data[,2]<0.5),]
p<-ggplot(data1,aes(x = V2,color=V1,y=V3))+#scale_x_discrete(limits = positions)+#geom_bar()
  geom_point()+ #scale_size(range = c(1, 1))+
  theme_bw()+
  #xlim(c(.7,1.7))
# ylab("Density") +
    scale_color_discrete(name="Mutation\nType") +
  #scale_x_log10(breaks=x_breaks)
  #scale_x_continuous(trans = 'log2', breaks =  c(1,2.5,5,10,25,50,100,200),lim=c(1,200))+
  #scale_x_continuous(trans='log2') +
  #scale_y_continuous(trans='log2')+
  xlab("Predicted mutation rate") +
  ylab("Observed mutation rate")+
  theme(legend.text = element_text(size = 15),legend.title = element_text( size=15, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(  axis.text=element_text(size = 15),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())#+theme(legend.position = "none")
p