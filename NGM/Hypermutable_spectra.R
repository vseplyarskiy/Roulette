library(ggplot2)
library(jtools) # Load jtools

library(ggplot2)
setwd("C:/Vova/study/R_scripts_data/R_data/NGM/")


data<-read.table(file="rate_spectra_v2", sep=" ");
data[data[,1]=='other',1]<-'Remaining genome'
library(stringr)
library(DescTools)
data$V4/data$V5
CI<-PoissonCI(data[,4], n = 1, conf.level = 0.95, sides = c("two.sided"),
              method = c("exact"))
sum(CI[,3])

p<-ggplot(data,aes(x = V2,fill=V2,y=V3))+
                  # alpha = str_match(V3, "transcribed|non-transcribed")[,1])

  #geom_point(size=1,position=position_dodge(width = 0.1))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_y_continuous(trans='log2',limits=c(0.94, 14))+ 
  geom_errorbar(aes(ymin=CI[,2]/data$V5, ymax=CI[,3]/data$V5), width=.2,
                position=position_dodge(.9))+
  facet_wrap(~V1,strip.position="top",scales = "free")+
  
  
  theme_classic()+
  
  xlab("") +
  ylab("observed/expected") +
  #ylim(c(0,18))+
  #scale_color_discrete(name="mutation") +
  #scale_alpha_manual(name="strand",values = c(1,0.5))+
  #scale_x_log10(breaks=x_breaks)
  #annotate("rect", xmin = -80, xmax = 80, ymin = .7, ymax = 2.2,alpha = .8)+
  #annotate("text", x = 0, y = 1.5,  
   #        label = "RNU" , color="orange",size=6,fontface="bold")+ 
  theme(plot.title = element_text(hjust = 0.5,size=22))+
  theme( axis.text.y=element_text(size=15),axis.text.x=element_text(size=15),axis.title=element_text(size=17),
         legend.text=element_text(size=14),
         strip.text.x = element_text(size = 20)
  )+
  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border.y = element_blank()
  )+ theme(legend.position="none")
  
  ggtitle("Mutational spectra") 

p

ggsave(filename="C:/Vova/study/NGM/Figures/Supplementary_Fig/RNU_facet_spectra.svg", plot = p,device = "svg", path = NULL,
       scale = 1, width = 240, height = 140, units = c("mm"),
       dpi = 300, limitsize = TRUE)