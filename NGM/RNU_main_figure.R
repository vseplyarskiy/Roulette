library(ggplot2)
setwd("C:/Vova/study/R_scripts_data/R_data/NGM/")


data<-read.table(file="rate_mut_RNU", sep=" ");

library(stringr)
library(DescTools)

CI<-PoissonCI(data[,5], n = 1, conf.level = 0.95, sides = c("two.sided"),
              method = c("exact"))
sum(CI[,3])

p<-ggplot(data,aes(x = V1,color=V2,y=V4,
                   alpha = str_match(V3, "transcribed|non-transcribed")[,1])
)+
  #geom_point(size=1,position=position_dodge(width = 0.1))+
  geom_line(size=1)+
  xlim(c(-505,505))+
  #geom_ribbon(aes(ymin=CI[,2]/data$V6, ymax=CI[,3]/data$V6, group=data$V3), alpha=0.1)+
  #geom_smooth(n=95,span=0.01,se=F)+
  
  #scale_size(range = c(1, 1))+
  theme_classic()+
  
  xlab("RNU center") +
  ylab("observed/expected") +
  ylim(c(0,18))+
  scale_color_discrete(name="mutation") +
  scale_alpha_manual(name="strand",values = c(1,0.5))+
  #scale_x_log10(breaks=x_breaks)
  annotate("rect", xmin = -80, xmax = 80, ymin = .7, ymax = 2.2,alpha = .8)+
  annotate("text", x = 0, y = 1.5,  
           label = "RNU" , color="orange",size=6,fontface="bold")+ 
  theme(plot.title = element_text(hjust = 0.5))+
  theme( axis.text.y=element_text(size=15),axis.text.x=element_text(size=15),axis.title=element_text(size=17),
         legend.text=element_text(size=14),
         strip.text.x = element_text(size = 20)
  )+
  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border.y = element_blank()
  )
  #facet_wrap(~V2,strip.position="top",scales = "free")+
  
  #ggtitle("Hypermutable RNU") 

p

#ggsave(filename="C:/Vova/study/NGM/Figures/Supplementary_Fig/RNU_facet_errorbars.svg", plot = p,device = "svg", path = NULL,
#       scale = 1, width = 480, height = 280, units = c("mm"),
#       dpi = 300, limitsize = TRUE)

ggsave(filename="C:/Vova/study/NGM/Figures/Figure2/RNU.svg", plot = p,device = "svg", path = NULL,
       scale = 1, width = 150, height = 90, units = c("mm"),
       dpi = 300, limitsize = TRUE)

