require(ggplot2)

library(scales)

setwd("C:/Vova/Vova/study/R_scripts_data/R_data/NGM//")
data<-read.table(file="clin_var", sep=" ");

ggplot(data, aes(x = V3, y=1,color=V2)) +
  geom_point()+
  #scale_x_continuous(trans='log2')+ 
  xlab("Experimental assay")+
  #ylab("T-asymmetry")+
  #stat_function(fun = age_ef_rev, colour = "black")+
  theme_bw()+
  theme(axis.line =element_line(colour = "black"),
    #    axis.text.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

setwd("C:/Vova/Vova/study/R_scripts_data/R_data/NGM//")
data<-read.table(file="Grep_per_nucl", sep=" ");
ggplot(data, aes(x = V5, y=V6)) +
  geom_point(size=1.5, alpha=0.3)+
  #scale_x_continuous(trans='log2')+ 
  xlab("ADRB2")+
  ylab("MC4R")+
  #stat_function(fun = age_ef_rev, colour = "black")+
  theme_bw()+
  theme(axis.line =element_line(colour = "black"),
        #    axis.text.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


setwd("C:/Vova/Vova/study/R_scripts_data/R_data/NGM//")
data<-read.table(file="Grep_per_nucl", sep=" ");
ggplot(data, aes(x = V5)) +
  geom_density()+
  #scale_x_continuous(trans='log2')+ 
  xlab("GERRP,ADRB2")+
  #xlab("GERRP, MC4R")+
  ylab("density")+
  #stat_function(fun = age_ef_rev, colour = "black")+
  theme_bw()+
  theme(axis.line =element_line(colour = "black"),
        #    axis.text.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


cor.test(data$V5,data$V6,method="spearman")