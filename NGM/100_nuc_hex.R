library(ggplot2)

setwd("C:/Vova/study/R_scripts_data/R_data/NGM/")
#data<-read.table(file="distance_low_FDR", header = F, sep=" ")
data<-read.table(file="hash", header = T, sep=" ")

data[data[,1]==14419000,]
p1<-ggplot(data, aes(x=data[,26], y=data[,27])) + geom_hex() + scale_fill_viridis(trans = "log10")+
  ylab("observed")+xlab("expected")+
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=1)+
  geom_abline(intercept = 70, slope = 0, color="gray", 
              linetype="dashed",size=1.5)+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme( axis.text.y=element_text(size=15),axis.text.x=element_text(size=15),axis.title=element_text(size=17),
         legend.text=element_text(size=14))+
  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border.y = element_blank()
  )+
ggtitle("Rare SNVs in 100 nucleotide window, Chr1") 

p1


ggsave(filename="C:/Vova/study/NGM/Figures/Figure2/hex100.svg", plot = p1,device = "svg", path = NULL,
       scale = 1, width = 150, height = 90, units = c("mm"),
       dpi = 300, limitsize = TRUE)