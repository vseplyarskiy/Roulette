
library(randomForest)
library(ggplot2)
setwd("C:/Vova/study/R_scripts_data/R_data/NGM/")




MC4R<-read.table(file="full_tm_MC4R_full", header = T, sep=" ")

evan<-read.table(file="frequency_pat", header = F, sep=" ")
evan2<-read.table(file="b2ar_forrand_gemme", header = F, sep=" ")

evan1<-read.table(file="b2ar_forrand_gemme_all", header = F, sep=" ")

#effect$mutAAfreq<-log(effect$mutAAfreq+0.003)
rf1<-randomForest(y=evan2$V5,x=evan2[,c(1,6,8,10,12,14)],ntree=500)
rf1<-randomForest(y=evan2$V5,x=evan2[,c(1,7,9,11,13,15)])
rf1<-randomForest(y=evan2$V5,x=evan1[,c(1,8)],ntree=500)

cor.test(abs(evan1[,6]**1.9),abs(evan2[,6]**2),method="pearson")
length(evan2[,8])
evan2[,6]<-(abs(evan2[,6]))**0.1
evan2[,7]<-(abs(evan2[,7]))**0.1
evan2[,8]<-((evan2[,8]))**1/3
summary(lmv<-lm((V5)~ V1+V7+V8,data=evan2))







rf1<-randomForest(y=evan1$avarege_cAMP,x=evan1[,c(1,6,7)])
summary(lmv<-lm((evan$V5)~#evan$V1*(abs(evan$V7))+
             #evan$V1*sqrt(abs(evan$V11))+
           #  evan$V1*(abs(evan$V8))+
          #   evan$V1*(abs(evan$V9))+
             evan$V1*(abs(evan$V6+evan$V7+evan$V8+evan$V9+evan$V10))
             #evan$V1*(abs(evan$V10))
             ))

evan2[,6]<-sqrt(abs(evan2[,6]))
evan2[,7]<-sqrt(abs(evan2[,7]))
evan2[,8]<-sqrt(abs(evan2[,8]))
evan2[,9]<-sqrt(abs(evan2[,9]))
evan2<-read.table(file="b2ar_forrand", header = F, sep=" ")

evan2[,6]<-(abs(evan2[,6]))**0.5
evan2[,7]<-(abs(evan2[,7]))**0.5
evan2[,8]<-(abs(evan2[,8]))**0.5
evan2[,9]<-(abs(evan2[,9]))**0.5
summary(lmv<-lm(V5~V1*V9+V7+V6+V8,data=evan2))


evan3<-read.table(file="GRPR", header = F, sep=" ")
evan3[,7]<-evan3[,6]
evan3[,5]<-evan3[,10]


evan3[,6]<-(abs(evan3[,6]))**0.5
evan3[,7]<-(abs(evan3[,7]))**0.5
evan3[,8]<-(abs(evan3[,8]))**0.5
evan3[,9]<-(abs(evan3[,9]))**0.5

evan4<-read.table(file="MC4R_psic", header = F, sep=" ")
evan4[,7]<-evan4[,6]
evan4[,7]<-(abs(evan4[,7]))**0.5
evan4[,8]<-(abs(evan4[,8]))**0.5
evan4[,9]<-(abs(evan4[,9]))**0.5

pr4<-as.data.frame(predict(lmv,evan4))
evan4_pr<-cbind(evan4,pr4)
names(evan4_pr)[10]<-'pred'

pr3<-as.data.frame(predict(lmv,evan3))
evan3_pr<-cbind(evan3,pr3)
names(evan3_pr)[11]<-'pred'


ggplot(evan4_pr,aes(y=V5,x=pred))   +
  geom_point() +
  geom_smooth(method='lm', formula= y~x)+
  theme_bw()+
  xlab("prediction LM")+ylab("cAMP")+
  ggtitle("MC4R, freq")+theme(plot.title = element_text(hjust = 0.5))

ggplot(evan3_pr,aes(y=V10,x=pred))   +
  geom_point() +
  geom_smooth()+
  theme_bw()+
  xlab("prediction LM")+ylab("emax")+
  ggtitle("GRPR, Freq")+theme(plot.title = element_text(hjust = 0.5))
pr<-as.data.frame(predict(lmv,evan3))
dim(evan3)
evan3_pr<-cbind(evan3,pr)
names(evan3_pr)[9]<-'pred'


pr<-as.data.frame(predict(lmv,evan2))
evan2_pr<-cbind(evan2,pr)
names(evan2_pr)[10]<-'pred'

ggplot(evan3_pr,aes(x=V5,y=pred))   +
  geom_violin()+geom_point() +  theme_bw()+
  ylab("prediction LM")+xlab("GRPR class")

ggplot(evan3_pr,aes(y=V10,x=pred))   +
  geom_point() +
  geom_smooth()+
  theme_bw()+
  xlab("prediction LM")+ylab("emax")+
  ggtitle("GRPR")+theme(plot.title = element_text(hjust = 0.5))


ggplot(evan2_pr,aes(y=V5,x=pred))   +
  geom_point(alpha=0.2) + ylim(c(0,3))+
  geom_smooth()+
  theme_bw()+
  xlab("prediction LM")+ylab("e50")+
  ggtitle("B2AR")+theme(plot.title = element_text(hjust = 0.5))

+  #geom_smooth(method = "lm", se = FALSE)+
  theme(title=element_text(size=12))+ylab("CpG>CpA mutations in a window")+
  xlab("CpG>TpG mutations in a window")+
  theme( axis.text.y=element_text(size=15),axis.text.x=element_text(size=15))+
  #ylim(c(0,3.3))+xlim(c(0,3.3))+geom_abline(slope = 1,color="grey")+
  #ggtitle(compon)+theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

p

                  #evan$V1*(abs(evan$V7))+

           #       evan2$V1*(abs(evan2$V8))+
          #        evan2$V1*(abs(evan2$V9))+
            #      evan2$V1*(abs(evan2$V10))
)

evan[,6]<-sqrt(abs(evan[,6]))
evan[,7]<-sqrt(abs(evan[,7]))

pr<-as.data.frame(predict(lmv,evan))
names(pr)<-'act'
ggplot(pr,aes(x=act))+geom_density()

evan5<-read.table(file="frequency_exp", header = F, sep=" ")
evan5[,6]<-sqrt(abs(evan5[,6]))
evan5[,7]<-sqrt(abs(evan5[,7]))

pr<-as.data.frame(predict(lmv,evan5))
names(pr)<-'act'
ggplot(pr,aes(x=act))+geom_density()
summary(lm(evan5$V5~pr[,1]))

ggplot(evan5,aes(y=V5,x=pr[,1]))+geom_point()



pr
summary(lmv<-lm((evan1$avarege_cAMP)~#evan$V1*(abs(evan$V7))+
                  #evan$V1*sqrt(abs(evan$V11))+
                  #  evan$V1*(abs(evan$V8))+
                  #   evan$V1*(abs(evan$V9))+
                  evan1$TM+sqrt(abs(log(evan1$mutAAfreq+0.01/(evan1$mutAAfreq+evan1$wtAAfreq+0.1))))#*sqrt(abs(evan2$V7))
                #evan$V1*(abs(evan$V10))
))



ggplot(evan,aes(y=V5,x=lmv$fitted.values))+geom_point()+geom_smooth(method="lm")+
  xlab("lnear model on the sum of log freq. of all clusters")+ylab("activity B2AR")+
  theme_classic()

       
       
lmv$fitted.values
glm.D93 <- glm(evan$V5 ~ evan[,1], family = gaussian(link = "identity"))


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