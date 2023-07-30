library(dplyr)
library(tidyr)
library(pROC)
require(rms)

setwd("C:/Vova/study/R_scripts_data/R_data/NGM/")
fn="test_AACGT_T_inter_met1"
data<-read.table(file="test_AATTT_A_inter_met", sep="\t"); ### example of the table
data<-na.omit(data)

#data[,13]<-data[,13]+0.01
#data[,c(13)]<-as.numeric(data[,c(13)])
#data[,c(11)]<-as.numeric(data[,c(11)])*1
#data<-na.omit(data)

#### table structure: V1 - chromosome and coordinate; V2 has values "1" if rare SNV is absent and "0" if present
#### table structure: V3 - context name; V4-V11 mean effect of nucleotides in positions 4-8 to the left and to the right from the site 
#### table structure: V12 - direction of replication fork 1 - is to the left and 7 is to the right; V13 effect of local mutation rate 


data[,c(4:11,13)]<-log(data[,c(4:11,13)]) ### these columns reflect odds ratio for mutation to happen
### in columns 4-11 due to neighboring nucleotides and 13 due to regional mutation rate
### we expected these factors to have multiplicative effect on mutation rate, thus performed log transformation befor logit.
#data[,c(2)]<-as.factor(data[,c(2)])
data[,c(12)]<-as.factor(data[,c(12)])
data<-na.omit(data)
data<-data[data$V13>-10 & data$V13<10,]  ### exclude regions with unrealistic effect of local mutation rate these regions corresponds to less than 0.01% of the data
trc<-sample(length(data$V2),as.integer(length(data$V2))*0.5) ### select training and testing cohorts 

data1<-data[trc,]
data2<-data[-trc,]

mylogit1<-glm(V2~(V4+V5+V6+V7+V8+V9+V10+V11+V12+V13), data = data1, family = "binomial") ### logistic regression no intercations
mylogit2<-glm(V2~(V4+V5+V6+V7+V8+V9+V10+V11+V12+V13)^2, data = data1, family = "binomial") ### logistic regression pairwise intercations
#mylogit1lr<-lrm(V2~(V4+V5+V6+V7+V8+V9+V10+V11+V12+V13), data = data1)
summary(mylogit1)
summary(mylogit2)

#### predicting mutation probabilities for the test data [logistic regression with no (pr1) or with pairwise interactions (pr2)]
pr1<-predict(mylogit1, newdata = data2, type="response")
pr2<-predict(mylogit2, newdata = data2, type="response")

###### testing just naive multiplicative model
data2e<-exp(data2[,c(4:11,13)]) 
pr3<-mean(data2$V2)*data2e$V4*data2e$V5*data2e$V6*data2e$V7*data2e$V8*data2e$V10*data2e$V11*data2e$V13
pr3[pr3>=1]<-0.995 #### regularising probabilities higher than 1, from multiplicative approach

## Calculate receiver operating characteristic curve for testing set for three models

r1<-roc(data2[,2], pr1)$auc
r2<-roc(data2[,2], pr2)$auc
r3<-roc(data2[,2], pr3)$auc


## Calculate log likelihood for testing set for three models

LL1<-sum(-log(abs((1-data2[,2])-pr1)))
LL2<-sum(-log(abs((1-data2[,2])-pr2)))
LL3<-sum(-log(abs((1-data2[,2])-pr3)))
LL4<-sum(-log(abs((1-data2[,2])-mean(data2[,2]))))

## selecting output with highest ROC value and log likelihood better than just intercept will produce 
auc_m<-max(r1,r2,r3)
out<-mean(data[,2])
if (auc_m==r1 && LL1<LL4){out<-mylogit1$coefficients}
if (auc_m==r2 && LL2<LL4){out<-mylogit2$coefficients}
if (auc_m==r3 && LL3<LL4){out<-'independ'}

#### Save "out"

###Decile plot; Only for visualisation
DD_p<-NULL;

for (ii in 1:10){
  qq=ii/10;
  qqb=qq-0.1;
  pr<-pr2
  test<-data2[pr>=quantile(pr,qqb) & pr<quantile(pr,qq),]
  tq<-(1-mean(test[,2],))
  pq<-(1-mean(pr[pr>=quantile(pr,qqb) & pr<quantile(pr,qq)]))
  
  DD_p<-rbind(DD_p,cbind(tq,ii,'observed'),cbind(pq,ii,'expected'))
}
DD_p<-as.data.frame(DD_p)
DD_p[,1]<-as.numeric(DD_p[,1])
DD_p[,2]<-as.numeric(DD_p[,2])


### Plot for deciles of predicted mutation rate;
ggplot(DD_p, aes(x=ii, y = tq, color=V3)) +
  geom_point()+ylab("SNV rate")+xlab("Decile")+
  theme_classic()