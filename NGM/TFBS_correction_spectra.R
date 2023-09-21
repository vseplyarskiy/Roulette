library(ggplot2)

setwd("C:/Vova/study/R_scripts_data/R_data/NGM/")
data<-read.table(file="Fig_TFBS_improvment_v2", header = F, sep="\t")
nam<-data[1,]
data<-data[-1,]
max(data[,6])
for (numb in 2:5){
  #numb=3
  n2=numb*2;
  n3=numb*2+1;
  
  val=as.numeric(data[,n2]);
  name<-nam[n3]
  #val=data$V10
  compon<-paste(name)
  
  yma=max(as.numeric(data[,10]))*1.2;
  ymi=min(as.numeric(data[,10]))
  
  y1=c(rep(yma*0.9,6),rep(yma*0.93,6))
  y2=c(rep(yma*0.96,6),rep(yma*0.92,6))
  y01=rep(yma*0.99,6)
  y02=rep(yma,6)
  ttext=c("C>A","C>G","C>T","T>A","T>C","T>G")
  
  x01<-c("ACA>A","ACA>G","ACA>T","ATA>A","ATA>C","ATA>G")
  x02<-c("TCT>A","TCT>G","TCT>T","TTT>A","TTT>C","TTT>G")
  x001<-c("CCT>A","CCT>G","CCT>T","CTT>A","CTT>C","CTT>G")
  
  x1<-rep(x01,2)
  x2<-rep(x02,2)
  
  
  
  
  
  data$V1 <- factor(data$V1, levels = unique(data$V1))
  data$V2 <- factor(data$V2, levels = unique(data$V2))
  cols <- c("blue","black","red","grey","darkgreen","pink")
  cols3=NULL
  cols1=NULL
  for (i in 1:6){
    for (b in 1:1){
      for (cc in 1:2){
        if (cc ==2){
          cols1<-c(cols1,adjustcolor(cols[i],0.4))
          cols3<-c(cols3,adjustcolor(cols[i],0.4))
        }
        else {cols1<-c(cols1,cols[i])}
      }}}
  cols4<-c(cols,cols3)
  cols2=cols1
  #cols2[seq(2,2*length(cols1),2)] = adjustcolor(cols1,0.3)
  
  p<-ggplot(data, aes(x=V1, y=val*1,fill=V2)) +
    geom_bar(stat="identity",  width = 1,position = "dodge2")+
    theme_classic()+#ggtitle("rep")+#ylim(c(0,0.031))+
    #scale_fill_discrete(name="Mutation\ntype")+
    #scale_fill_manual("legend",values = cols)+ 
    scale_fill_manual(values=cols2, name="mutation\ntype")+
    
    theme(title=element_text(size=25))+ylab("obs/exp")+xlab("")+
    theme( axis.text.y=element_text(size=25),axis.text.x=element_text(size=4),axis.title.y =element_text(size=25))+theme(axis.text.x = element_text(angle = 90, vjust = 0.35,size=8))+
    #theme(legend.text = element_text(size = 9))+
    # ylim(c(ymi,yma))+
    annotate("rect", xmin = x1, xmax = x2, ymin = y1, ymax = y2, fill=cols4)+
    annotate("text", x = x001, y = (y01+y02)/2,  color="black",label=ttext, size=6)+
    
    theme(legend.key.width = unit(0.5, "cm"),legend.key.height = unit(0.5, "cm"))+
    ggtitle(compon)+theme(plot.title = element_text(hjust = 0.5))+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())+
    scale_y_continuous(trans = 'log2',breaks = c(0.9,1.2,1.5))+
    theme(legend.position = "none")#+
  #guides(fill=guide_legend(ncol=3))+ theme(legend.title.align=0.5,legend.title =  element_text(size=25))+
  #scale_fill_manual(name="", labels=c("b"="remaning genome", "r"="high intensity"))+
  
  #theme(legend.position = "none")+
  #theme(axis.title.x=element_blank(),
  #     axis.text.x=element_blank(),
  #    axis.ticks.x=element_blank(),axis.line.x=element_blank())
  p
  fi<-paste("C:/Vova/study/NGM/Figures/Supplementary_Fig/",name,".svg",sep ="")
  #fi<-paste("C:/Vova/Vova/study/Signature germline/Figures/SFigures/Rev1/Signature_positive_axis/legend.svg",sep ="")
  
  ggsave(filename=fi, plot = p,device = "svg", path = NULL,
         scale = 1, width = 215, height = 85, units = c("mm"),
         dpi = 300, limitsize = TRUE)
}

