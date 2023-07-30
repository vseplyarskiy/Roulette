
library(randomForest)
setwd("C:/Vova/Vova/study/R_scripts_data/R_data/NGM/")




MC4R<-read.table(file="full_tm_MC4R_full", header = T, sep=" ")

evan<-read.table(file="KL_lhdiff.csv", header = T, sep=",")
#effect$mutAAfreq<-log(effect$mutAAfreq+0.003)
effect<-read.table(file="full_table_KL_ev_phyl", header = T, sep=" ")
fae<-effect$avarege_cAMP
fae[effect$avarege_cAMP<.8]<-'D'
fae[effect$avarege_cAMP>1.7]<-'A'
fae[effect$avarege_cAMP>0.8 & effect$avarege_cAMP<1.7]<-'N'
sum(fae=='A' )
effect<-cbind(effect[fae=='A' | fae=='D' | fae=='N',],fae[fae=='A' | fae=='D' | fae=='N'])
names(effect)[25]<-'fae'
names(effect)[11]
rf1 <- randomForest(
  y=effect$avarege_cAMP[effect$avarege_cAMP<.9] ,x=effect[effect$avarege_cAMP<.9,c(6,7,11,23,24)],
  #x=effect[effect$avarege_cAMP>0.7,c(8,9,11,23,24,14,17,20)],#x=effect[,c(6,7,8,9,11,23,24)],
  data=effect
)
rf1$confusion
summary(rf1)
y_pred = as.vector(predict(rf1, x=MC4R[,c(8,9,10)],newdata = MC4R))
ypred1<-y_pred
ypred1<-droplevels(ypred1)
ypred1[y_pred=='A']<-'Activating'
ypred1[y_pred=='N']<-'Neutral'
ypred1[y_pred=='D']<-'Lower activity'
cbind(MC4R$avarege_cAMP,ypred1)
cbind(MC4R$avarege_cAMP,y_pred)

y_pred = predict(rf1, x = MC4R$wtAAfreq+ MC4R$mutAAfreq+ MC4R$PIscore+ MC4R$KL+ MC4R$lhr,data=MC4R)




rf1 <- randomForest(
  effect$fae ~effect$PIscore+effect$KL+effect$lhr,
  data=effect
)
rf1$confusion


rf2 <- randomForest(
  evan$factor ~evan$lh.diff+evan$primat+evan$entropy+evan$KL,mtry=4,
  data=evan
)
cbind()
factor[evan$activity<0.7]<-'D'
factor[evan$activity>1.2]<-'A'
factor[evan$activity>0.7 & evan$activity<1.2]<-'N'
factor<-evan$activity

evan<-cbind(evan,factor)

rf2$rsq
rf2$confusion

rf1$rsq
T_as_mat<-round(importance(rf1), 2)