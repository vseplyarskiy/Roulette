library(ggplot2)
setwd("C:/Vova/study/R_scripts_data/R_data/NGM/")


data<-read.table(file="rate_chipseq_pol", sep=" ");
data<-data[data[,1] != 'GSE30267.ELL2.hct116_starved',]
p<-ggplot(data,aes(x = V2))+#scale_x_discrete(limits = positions)+#geom_bar()
  geom_histogram( position="dodge2",bins = 150)+ #scale_size(range = c(1, 1))+
  theme_bw()+
  
  ylab("Density") +
  scale_fill_discrete(name="Mutation\nType") +
  #scale_x_log10(breaks=x_breaks)
  #scale_x_continuous(trans = 'log2', breaks =  c(1,2.5,5,10,25,50,100,200),lim=c(1,200))+
  #scale_x_continuous(trans='log2') +
  #scale_y_continuous(trans='log2')+
  xlab("Observed/Expected") +
  theme(legend.text = element_text(size = 15),legend.title = element_text( size=15, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(  axis.text=element_text(size = 15),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())#+theme(legend.position = "none")
p




data<-read.table(file="rate_chipseq_den", sep=" ");
data<-data[c(875,876),]
p<-ggplot(data,aes(x = V1,fill=V6,y=V3))+#scale_x_discrete(limits = positions)+#geom_bar()
  geom_bar( position="dodge2",stat="identity")+ #scale_size(range = c(1, 1))+
  theme_bw()+
  
  ylab("# de novo mutations") +
  scale_fill_discrete(name="") +
  #scale_x_log10(breaks=x_breaks)
  #scale_x_continuous(trans = 'log2', breaks =  c(1,2.5,5,10,25,50,100,200),lim=c(1,200))+
  xlab("TFBS") +
  theme(legend.text = element_text(size = 15),legend.title = element_text( size=15, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(  axis.text=element_text(size = 15),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
           panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())#+theme(legend.position = "none")
p


data<-read.table(file="rate_ench", sep=" ");
data1<-read.table(file="rate_ench_cntrl", sep=" ");

data2<-rbind(data,data1)

data<-data[c(5:10),]
p<-ggplot(data2,aes(x = V1,color=V6,y=V2))+#scale_x_discrete(limits = positions)+#geom_bar()
  geom_point( position=position_dodge(width=0.2))+ #scale_size(range = c(1, 1))+
  theme_minimal()+
  
   xlab("Commpartment") +
  ylab("observed/expected") +
  ylim(c(0.8,1.2))+
  scale_color_discrete(name="") +
  #scale_x_log10(breaks=x_breaks)
  #scale_x_continuous(trans = 'log2', breaks =  c(1,2.5,5,10,25,50,100,200),lim=c(1,200))+

  theme(axis.text.x = element_text(angle = 90, vjust = 0.35,size=10))+
  theme(legend.text = element_text(size = 15),legend.title = element_text( size=15, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(  axis.text=element_text(size = 15),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())#+theme(legend.position = "none")
p


data<-read.table(file="rate_mut_ench", sep=" ");


p<-ggplot(data,aes(x = V6,color=V7,y=V2))+#scale_x_discrete(limits = positions)+#geom_bar()
 # geom_line()+
  geom_point(size=0.4)+
  geom_smooth(n=95,span=0.01,se=F)+
  
  #scale_size(range = c(1, 1))+
  theme_classic()+
  
  xlab("Distance from CTCF+promoter") +
  ylab("observed/expected") +
  ylim(c(0.9,1.8))+
  scale_color_discrete(name="") +
  #scale_x_log10(breaks=x_breaks)
  #scale_x_continuous(trans = 'log2', breaks =  c(1,2.5,5,10,25,50,100,200),lim=c(1,200))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme( axis.text.y=element_text(size=15),axis.text.x=element_text(size=15),axis.title=element_text(size=17),
         legend.text=element_text(size=14))+
  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border.y = element_blank()
  )+
  ggtitle("Mutagenic CTCF") 

p

ggsave(filename="C:/Vova/study/NGM/Figures/Figure2/CTCF.svg", plot = p,device = "svg", path = NULL,
       scale = 1, width = 150, height = 90, units = c("mm"),
       dpi = 300, limitsize = TRUE)



data<-read.table(file="rate_mut_RNU", sep=" ");

library(stringr)

p<-ggplot(data,aes(x = V1,color=V2,y=V4,
                   alpha = str_match(V3, "transcribed|non-transcribed")[,1]))+
   #geom_point(size=1,position=position_dodge(width = 0.1))+
  geom_line(size=1)+
   xlim(c(-505,505))+
  #geom_smooth(n=95,span=0.01,se=F)+
  
  #scale_size(range = c(1, 1))+
  theme_classic()+
  
  xlab("RNU center") +
  ylab("observed/expected") +
 # ylim(c(0.8,10))+
  scale_color_discrete(name="mutation") +
  scale_alpha_manual(name="strand",values = c(1,0.5))+
  #scale_x_log10(breaks=x_breaks)
  annotate("rect", xmin = -80, xmax = 80, ymin = .7, ymax = 2.2,alpha = .8)+
    annotate("text", x = 0, y = 1.5,  
             label = "RNU" , color="orange",size=6,fontface="bold")+ 
  theme(plot.title = element_text(hjust = 0.5))+
  theme( axis.text.y=element_text(size=15),axis.text.x=element_text(size=15),axis.title=element_text(size=17),
         legend.text=element_text(size=14))+
  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border.y = element_blank()
  )+
  ggtitle("Hypermutable RNU") 

p

ggsave(filename="C:/Vova/study/NGM/Figures/Figure2/RNU.svg", plot = p,device = "svg", path = NULL,
       scale = 1, width = 150, height = 90, units = c("mm"),
       dpi = 300, limitsize = TRUE)



data<-read.table(file="rate_mut_tRNA", sep=" ");

library(stringr)
library(ggplot2)

p<-ggplot(data,aes(x = V1,color=V2,y=V4,
                   alpha = str_match(V3, "transcribed|non-transcribed")[,1]))+
  geom_line(size=1)+
  xlim(c(-105,105))+
  #geom_smooth(n=95,span=0.01,se=F)+
  
  #scale_size(range = c(1, 1))+
  theme_classic()+
  
  xlab("tRNA center") +
  ylab("observed/expected") +
  ylim(c(0.6,12))+
  scale_color_discrete(name="mutation") +
  scale_alpha_manual(name="strand",values = c(1,0.5))+
  #scale_x_log10(breaks=x_breaks)
  annotate("rect", xmin = -35, xmax = 35, ymin = .65, ymax = 1.7,alpha = .8)+
  annotate("text", x = 0, y = 1.2,  
           label = "tRNA" , color="orange",size=6,fontface="bold")+ 
  theme(plot.title = element_text(hjust = 0.5))+
  theme( axis.text.y=element_text(size=15),axis.text.x=element_text(size=15),axis.title=element_text(size=17),
         legend.text=element_text(size=14))+
  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border.y = element_blank()
  )+
  ggtitle("Hypermutable tRNA") 

p

ggsave(filename="C:/Vova/study/NGM/Figures/Figure2/tRNA.svg", plot = p,device = "svg", path = NULL,
       scale = 1, width = 150, height = 90, units = c("mm"),
       dpi = 300, limitsize = TRUE)

