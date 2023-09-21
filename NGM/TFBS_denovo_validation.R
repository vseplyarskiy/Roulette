
require(ggplot2)
library(DescTools)
setwd("C:/Vova/study/R_scripts_data/R_data/NGM/")
data<-read.table(file="rate_mut_denovo_TFBS", header = F, sep=" ")
data[data[,3] !='CG>A'& data[,3] !='CG>T' & data[,3] !='CG>G',]->data
data[data[2]==0,2]<-'0-100 nucleotides'
data[data[2]==1,2]<-'101-200 nucleotides'
norm<-sum(data[,6])/sum(data[,5])
CI<-PoissonCI(data$V5, n = 1, conf.level = 0.95, sides = c("two.sided"))
CI<-CI*norm
#data<-read.table(file="enr_norm", header = F, sep=" ")
#val<-cbind(sapply(1:nrow(data),function(i){rateratio.test(c(data[i,3],data[i,4]),c(c(data[i,5],data[i,6])))$conf.int[c(1,2)]}))
#val<-t(val)
data$V1 <- factor(data$V1, levels = unique(data$V1))

p1<-ggplot(data, aes( x=V3, y=(V5/data$V6)*norm,fill=V1)) + 
  geom_bar(stat="identity",width = 0.8,position=position_dodge(.9)) +
  scale_fill_manual( name="tissue",values = c("red","yellow","violet", "green")) +
  ylab("observed/expected")+
  #ylab("rate")+
  theme_classic() +
  geom_errorbar( aes( ymin=CI[,2]/data$V6, ymax=CI[,3]/data$V6),stat="identity",
                 width = 0.5,position=position_dodge(.9))+
  #theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1))+
  ggtitle("SNV rate at TFBS")+theme(plot.title = element_text(hjust = 0.5))+
  #theme(legend.position = "none") + 
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank())+
  #xlab("Regions")+  
  facet_wrap( ~ V2, nrow = 1)#+
  #theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1))+
  #scale_y_continuous(breaks =c(1,2,4),limits = c(0.9,4))
ggsave(filename="C:/Vova/study/NGM/Figures/Supplementary_Fig/TFBS_denovo.svg", plot = p1,device = "svg", path = NULL,
       scale = 1, width = 150, height = 90, units = c("mm"),
       dpi = 300, limitsize = TRUE)
p1
