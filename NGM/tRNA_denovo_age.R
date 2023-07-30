library(ggplot2)
library(jtools) # Load jtools


setwd("C:/Vova/study/R_scripts_data/R_data/NGM");
data<-read.table(file="tRNA_age")
sum(data[data[,1]=='tRNA',3])
sum(data[data[,1]=='Imunoglobulin',3])
sum(data[data[,1]=='RNU',3])

summary(glm(V3 ~ V2, family=poisson(link = "identity"),data[data[,1]=='tRNA',]))
summary(glm(V3 ~ V2, family=poisson(),data[data[,1]=='RNU',]))
RtRNA<-(glm(V3 ~ V2, family=poisson(),data[data[,1]=='RNU',]))
tRNAage<-effect_plot(RtRNA, pred = V2,plot.points=F, interval = TRUE,  
            x.label="Age of Father", y.label = "#mutaions in RNU")
fi<-paste("C:/Vova/study/NGM/Figures/Supplementary_Fig/","RNUage",".svg",sep ="")
#fi<-paste("C:/Vova/Vova/study/Signature germline/Figures/SFigures/Rev1/Signature_positive_axis/legend.svg",sep ="")

ggsave(filename=fi, plot = tRNAage,device = "svg", path = NULL,
       scale = 1, width = 215, height = 85, units = c("mm"),
       dpi = 300, limitsize = TRUE)

20*0.0006687