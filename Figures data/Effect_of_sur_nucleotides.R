setwd("C:/Vova/study/R_scripts_data/R_data/NGM/")
data<-read.table(file="differnt_effect", header = F, sep=" ")
names(data)[c(6,7)]<-c('nucleotide','compartment')
data$V5<-factor(data$V5, levels = unique(data$V5))
p<-ggplot(data, aes( x=V5, y=V4,color=nucleotide,shape=compartment)) + 
  geom_point(position=position_dodge2(.1),size=1.2) +
 # scale_color_manual( name="",values = c("red","black")) +
  # scale_alpha_manual(name="strand",values = c(1,.3))+
  ylab("Effect on mutation rate")+
  xlab("nucleotide position")+
  theme_classic() +
  #geom_pointrange( aes( ymin=CI1[,2], ymax=CI1[,3]),stat="identity",position=position_dodge2(.1))+
  #theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1))+
  ggtitle("Effect of surronding nucleotides\n
          on rate of AGCAC>A mutations")+theme(plot.title = element_text(hjust = 0.5))+
  #theme(legend.position = "none") + 
  theme(plot.title = element_text(hjust = 0.5))+
  theme( axis.text.y=element_text(size=15),axis.text.x=element_text(size=15),axis.title=element_text(size=17),
         legend.text=element_text(size=14))+
  #theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())# + theme(legend.position=c(.75,.9))#theme(legend.justification=c(1,1), legend.position=c(1,1))

p
ggsave(filename="C:/Vova/study/NGM/Figures/Supplementary_Fig/Effect_of_sur_nucleotides.svg", plot = p,device = "svg", path = NULL,
       scale = 1, width = 160, height = 90, units = c("mm"),
       dpi = 300, limitsize = TRUE)