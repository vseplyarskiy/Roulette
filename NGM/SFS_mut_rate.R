library(ggplot2)
setwd("C:/Vova/study/R_scripts_data/R_data/NGM/")

data<-read.table(file="rate_evan_rec", sep=" ");
#data<-read.table(file="rate_evan", sep=" ");

#data1<-data[data[,1]>=.5,]

p<-ggplot(data,aes(x = V3,y=V2))+#scale_x_discrete(limits = positions)+#geom_bar()
  geom_bar(stat = 'identity', position="dodge2",aes(fill=V1)) +
  scale_fill_gradient(low = "blue", high = "red", name="mutation rate")+
   #scale_y_sqrt()+
  #scale_y_log10()+
  scale_y_continuous(trans='log10')+
  facet_wrap(~V6,nrow = 4,
             labeller = as_labeller(c(`1`="recombination bin 1",
                                        `2`="recombination bin 2",
                                         `3`="recombination bin 3",
                                        `4`="recombination bin 4")))+

  scale_x_continuous(breaks=c(1,2,3,4,5,6),
                     labels=c( "0","(0,0.5e-5]",
"(0.5e-5,0.5e-4]",  "(0.5e-4,0.5e-3]",
 "(0.5e-3,0.5e-2]",  "(0.5e-2,0.2]"))+
xlab("frequency bin")+ylab("proportion of possible mutations")+
  theme_minimal()+
  theme(legend.text = element_text(size = 12),axis.title =element_text(size = 16),
        axis.text  =element_text(size = 12),
        legend.title = element_text( size=12, face="bold"))+
  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border.y = element_blank()
  )
p


ggsave(filename="C:/Vova/study/NGM/Presentation/Figures/Reccurence_SFS.svg", plot = p,device = "svg", path = NULL,
       scale = 1, width = 350, height = 100, units = c("mm"),
       dpi = 300, limitsize = TRUE)





data<-read.table(file="rate_evan3", sep=" ");
p<-ggplot(data,aes(x = V3,y=V2))+#scale_x_discrete(limits = positions)+#geom_bar()
  geom_bar(stat = 'identity', position="dodge2",aes(fill=V1)) +
  scale_fill_gradient(low = "blue", high = "red")+
  scale_y_sqrt()+
  scale_x_continuous(breaks=c(3,4,5,6),
                     labels=c( 
                               "(0.5e-5,0.5e-4]",  "(0.5e-4,0.5e-3]",
                               "(0.5e-3,0.5e-2]",  "(0.5e-2,0.2]"))+
  xlab("frequency bin")+ylab("proportion of possible mutations")+
  theme_minimal()
p




  scale_fill_gradient2(low='orange', mid='snow', high='blue')

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
