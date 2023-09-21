library(ggplot2)
setwd("C:/Vova/study/R_scripts_data/R_data/NGM/")



data<-read.table(file="RNU_ind_qc_cat", sep=" ");
#data<-data[data[,1] != 'GSE30267.ELL2.hct116_starved',]
p<-ggplot(data,aes(x = V3,color=V2))+#scale_x_discrete(limits = positions)+#geom_bar()
  geom_density()+ #scale_size(range = c(1, 1))+
  theme_classic()+
  
  ylab("Density") +
  scale_fill_discrete(name="Mutation\nType") +
  #scale_x_log10(breaks=x_breaks)
  #scale_x_continuous(trans = 'log2', breaks =  c(1,2.5,5,10,25,50,100,200),lim=c(1,200))+
  #scale_x_continuous(trans='log2') +
  #scale_y_continuous(trans='log2')+
  xlab("AS_SOR per SNV") +
  theme(legend.text = element_text(size = 15),legend.title = element_text( size=15, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(  axis.text=element_text(size = 15),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())#+theme(legend.position = "none")
p