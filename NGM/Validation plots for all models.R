library(dplyr)
library(tidyr)
library(pROC)
require(rms)
library(ggplot2)
setwd("C:/Vova/study/R_scripts_data/R_data/NGM/")
data<-read.table(file="AACGT_T_gen_4_met5_cl", sep="\t");
data<-na.omit(data)

#data[,13]<-data[,13]+0.01
data[,c(13)]<-as.numeric(data[,c(13)])
data[,c(11)]<-as.numeric(data[,c(11)])*1
data<-na.omit(data)

data[,c(4:11,13)]<-log(data[,c(4:11,13)])
#data[,c(2)]<-as.factor(data[,c(2)])
data[,c(12)]<-as.factor(data[,c(12)])
data<-na.omit(data)
data<-data[data$V13>-10 & data$V13<10,]
trc<-sample(length(data$V2),as.integer(length(data$V2))*0.5)

data1<-data[trc,]
data2<-data[-trc,]
data2e<-exp(data2[,c(4:11,13)])

mylogit1<-glm(V2~(V4+V5+V6+V7+V8+V9+V10+V11+V12+V13), data = data1, family = "binomial")
mylogit2<-glm(V2~(V4+V5+V6+V7+V8+V9+V10+V11+V12+V13)^2, data = data1, family = "binomial")
#mylogit1lr<-lrm(V2~(V4+V5+V6+V7+V8+V9+V10+V11+V12+V13), data = data1)
summary(mylogit1)
summary(mylogit2)

pr1<-predict(mylogit1, newdata = data2, type="response")
pr2<-predict(mylogit2, newdata = data2, type="response")
pr3<-mean(data2$V2)*data2e$V4*data2e$V5*data2e$V6*data2e$V7*data2e$V8*data2e$V10*data2e$V11*data2e$V13



pr3[pr3>=1]<-0.995

r1<-roc(data2[,2], pr1)$auc
r2<-roc(data2[,2], pr2)$auc
r3<-roc(data2[,2], pr3)$auc

print (c(r1,r2,r3))

LL1<-sum(-log(abs((1-data2[,2])-pr1)))
LL2<-sum(-log(abs((1-data2[,2])-pr2)))
LL3<-sum(-log(abs((1-data2[,2])-pr3)))
LL4<-sum(-log(abs((1-data2[,2])-mean(data2[,2]))))

print (c(LL1,LL2,LL3,LL4))

setwd("C:/Vova/study/R_scripts_data/R_data/NGM/")
data<-read.table(file="l", sep=" ");

r1<-roc(data[,6], data[,3])$auc

LL1<-sum(-log(abs((1-data2[,2])-pr1)))

hist_val<-function(bin,val,pred){
out<-NULL;
out<-as.data.frame(out)
pr<-pred;
#bin=100;
  for (i in 1:bin){
    #i=1; val=data[,6];pred<-data[,3];pr<-pred;
    
  vb=(i-1)/bin
  vt=(i)/bin
  obs<-mean(val[pr>quantile(pred,vb) & pr<=quantile(pred,vt)])
  expc<-mean(pr[pr>quantile(pred,vb) & pr<=quantile(pred,vt)])
  out<-rbind(out,c(i,obs,expc))
  #out<-rbind(out,c(i,expc,'predicted'))
  #out<-rbind(out,c(i,expc,'predicted'))
  }
  colnames(out)<-c('bin','observed','expected')
return(out)
}

hh<-hist_val(10,data[,6],data[,4])
hh1<-hist_val(5,data2[,2],pr1)

data<-read.table(file="rate_denovo", sep=" ");
#data<-cbind(data,data[2]-data[1],data[,2]/data[,1])
#colnames(data)[5]<-c('diff')
data1<-data[data[,5]=='FULL_model',]
data1<-data[data[,3]>=20,]

sum(data[,4][data[,5]>0.05])
sum(data[,4][data[,2]/data[,1]>1.5])
data[,1][data[,2]/data[,1]>1.5]
sum(data[,4])
  colnames(hh1)<-c('bin','observed','expected')
p<-ggplot(data1,aes(x = V1,y=V2,size = V4))+#scale_x_discrete(limits = positions)+#geom_bar()
      geom_point()+ scale_size(range = c(1, 1))+
    theme_bw()+
  #xlim(c(0,5000))+ylim(c(0,5000))+
  geom_abline(intercept = 0, slope = 1, color="grey", 
              linetype="dashed", size=.5)+
      ylab("Trinucl mutations") +
    scale_color_discrete(name="Mutation\nType") + 
  scale_x_continuous(trans='log2') +scale_y_continuous(trans='log2')+
  xlab("Predicted mutation rate") +
  theme(legend.text = element_text(size = 15),legend.title = element_text( size=15, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(  axis.text=element_text(size = 15),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          #axis.text.y=element_blank(),
          #axis.ticks.x=element_blank(),
          #axis.ticks.y=element_blank(),
          #axis.line.x = ,
          #axis.line.y = black,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          #panel.border = element_blank(),
          panel.background = element_blank())+theme(legend.position = "none")
p
ggplot(data,aes(x = abs(V2-1),fill=V5))+#scale_x_discrete(limits = positions)+#geom_bar()
  #geom_histogram(bins = 100, alpha(0.4))+
  #scale_x_continuous(trans='log2')+
  geom_density(alpha=.2)+
  theme_bw()

quantile(abs(data$V2[data$V5=='FULL_model']-1),.995)
quantile(abs(data$V2[data$V5=='Penata_model']-1),0.995)


median(data[,2])
pr<-data[,3]
tr1=0.99
test<-data[pr>=quantile(pr,.99),]
mean(test[,6],)
mean(pr[pr>=quantile(pr,.99)])

tr2=0.03
test<-data[pr<=quantile(pr,tr2),]
mean(test[,6],)
mean(pr[pr<=quantile(pr,tr2)])