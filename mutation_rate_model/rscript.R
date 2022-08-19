#setwd("C:/Vova/study/R_scripts_data/R_data/NGM/")
library(pROC)

  args <- commandArgs(trailingOnly = TRUE)
  fn <- args[1]

fnf<-paste("/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/for_LASSO/context/", fn,sep = "")
data<-read.table(file=fnf, sep="\t");
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
#dim(data2)


#summary
mylogit1<-glm(V2~(V4+V5+V6+V7+V8+V9+V10+V11+V12+V13), data = data1, family = "binomial")
#summary
mylogit2<-glm(V2~(V4+V5+V6+V7+V8+V9+V10+V11+V12+V13)^2, data = data1, family = "binomial")
#mylogit1lr<-lrm(V2~(V4+V5+V6+V7+V8+V9+V10+V11+V12+V13), data = data1)

pr1<-predict(mylogit1, newdata = data2, type="response")
pr2<-predict(mylogit2, newdata = data2, type="response")
pr3<-mean(data2$V2)*data2e$V4*data2e$V5*data2e$V6*data2e$V7*data2e$V8*data2e$V10*data2e$V11*data2e$V13

pr3[pr3>=1]<-0.999

r1<-roc(data2[,2], pr1)$auc
r2<-roc(data2[,2], pr2)$auc
r3<-roc(data2[,2], pr3)$auc



LL1<-sum(-log(abs((1-data2[,2])-pr1)))
LL2<-sum(-log(abs((1-data2[,2])-pr2)))
LL3<-sum(-log(abs((1-data2[,2])-pr3)))
LL4<-sum(-log(abs((1-data2[,2])-mean(data2[,2]))))
#print (c(LL1,LL2,LL3,LL4));
#print (c(r1,r2,r3));

auc_m<-max(r1,r2,r3)
out<-mean(data[,2])
if (auc_m==r1 && LL1<LL4){out<-mylogit1$coefficients}
if (auc_m==r2 && LL2<LL4){out<-mylogit2$coefficients}
if (auc_m==r3 && LL3<LL4){out<-'independ';out<-c(mean(data[,2]),out)}

fnf1<-paste("/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/for_LASSO/Logit/", fn,".model",sep = "")
fnf2<-paste("/net/home/vseplyarsky/GNOMAD_Tufts/mut_model/for_LASSO/Logit/", fn,".est",sep = "")

out<-c(fn,out);
#write(out,file=fnf1,sep="\t",row.names=FALSE,quote = FALSE,col.names = FALSE);
if (length(data[,2])<1000){out<-'too_short'}
AU<-rbind(c(r1,r2,r3,0.5),c(LL1,LL2,LL3,LL4));


write.table(t(out),file=fnf1,sep="\t",row.names=FALSE,quote = FALSE,col.names = FALSE);
write.table((AU),file=fnf2,sep="\t",row.names=FALSE,quote = FALSE,col.names = FALSE);

