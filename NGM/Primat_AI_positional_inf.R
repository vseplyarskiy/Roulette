library(ggplot2)

setwd("C:/Vova/Vova/study/R_scripts_data/R_data/NGM//")
data<-read.table(file="Primat_AI_positional_inf",header = FALSE,  sep=" ")
dim(data)
nro<-nrow(data)
nco<-ncol(data)

dens=NULL
for (i in 3:641){
  if(sum(as.numeric((is.na(data[,i]))))>100)
    next
  print (i)
  
  for (j in i:641){
    if(sum(as.numeric((is.na(data[,j]))))>100)
      next
      if(i==j)
      next
    #print (i)
    #print (j)
    
    res<-cor.test(data[,i],data[,j], method="pearson")
    rs<-as.character(res$estimate)
    dens<-rbind(dens,rs)
  }}
head(dens1)
dens1<-as.double(dens)
dens1<-as.data.frame(dens1)
ggplot(dens1,aes(x=dens)) +
  #geom_bar(stat="identity",  width = 1,position = "dodge2")+
  geom_density()+
  theme_bw()
hist(dens1$dens1,main="PrimatAI",xlab = "Pairwise positional correlation")
sum(as.numeric((is.na(data[,j]))))

data<-read.table(file="EVAN_Tri",header = FALSE,  sep=" ")
cor(data$V3,data$V6,method = "spearman")

data<-read.table(file="Evan_func",header = FALSE,  sep=" ")
(cor.test(data$V3[data$V4>0.01],data$V4[data$V4>0.005],method = "spearman"))
length(data$V3[data$V4<0.01 & data$V3>1]))&& data$V3[data$V4>0.01]
length(data$V4[data$V4>0.01])
chisq.test(cbind(c(31,0),c(157,55)))
55/212
ggplot(data,aes(x=V3 ,y=V4))   +
  #geom_bar(stat="identity",  width = 1,position = "dodge2")+
  geom_point(color="black")+
  theme_bw()+  #geom_smooth(method = "lm", se = FALSE)+
  #scale_y_continuous(trans='log2')+
  theme(title=element_text(size=12))+ylab("Selection constraint")+xlab("Intolerance(High = Tolernat)")+
  theme( axis.text.y=element_text(size=15),axis.text.x=element_text(size=15))+
  #ylim(c(0,maxy))+xlim(c(0,maxx))+
  ggtitle('')+theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())