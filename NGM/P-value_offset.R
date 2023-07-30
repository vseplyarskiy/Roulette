library(ggplot2)
setwd("C:/Vova/study/R_scripts_data/R_data/NGM/")



data<-read.table(file="top_genes_FP_syn",  sep=" ",header = F) ###http://pklab.med.harvard.edu/ruslan/spacemut/tracks_update/
#data[,c(7,8)]<-data[,c(7,8)]+0.0000001
data[,c(9,10,11)]<-data[,c(9,10,11)]+0.00000001


                 ggplot(data,aes(x = V9,y=V10))+#scale_x_discrete(limits = positions)+#geom_bar()
                   geom_point()+ #scale_size(range = c(1, 1))+
                   theme_bw()+
                   
                   ylab("Carlson p-value") +
                   xlab("Roullete p-value") +
                   
                    scale_y_continuous(trans='log2')+ 
                    scale_x_continuous(trans = 'log2')+
                    geom_abline(intercept = 0, slope = 1, color="grey")
                   
                   scale_fill_discrete(name="Mutation\nType") +
                   
                   #scale_x_continuous(trans='log2') +
                   
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