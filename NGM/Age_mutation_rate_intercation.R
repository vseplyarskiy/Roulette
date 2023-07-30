require(ggplot2)

setwd("C:/Vova/Vova/study/R_scripts_data/R_data/NGM/")

C_Tb1<-read.table(file="C_A_table", header = F, sep=" ")
C_Tb2<-read.table(file="C_G_table", header = F, sep=" ")
C_Tb3<-read.table(file="C_T_table", header = F, sep=" ")
C_Tb4<-read.table(file="CG_T_table", header = F, sep=" ")
C_Tb5<-read.table(file="T_A_table", header = F, sep=" ")
C_Tb6<-read.table(file="T_C_table", header = F, sep=" ")
C_Tb7<-read.table(file="T_G_table", header = F, sep=" ")

C_Tb<-cbind(C_Tb4,rowSums(cbind(C_Tb1$V1,C_Tb2$V1,C_Tb3$V1,C_Tb5$V1,C_Tb6$V1,C_Tb7$V1)))
names(C_Tb)[7]<-'V7'
#sum(C_Tb$V1)cbinzzzz
#sum(C_Tb$V1[C_Tb$V2=='E1' | C_Tb$V4==10001,])/
(C_Tb[C_Tb$V2=='E1' & C_Tb$V4==10001,3]/C_Tb[C_Tb$V2=='E15' & C_Tb$V4==10001,3])/(C_Tb1[C_Tb$V2=='E1' & C_Tb$V4==10001,3]/C_Tb1[C_Tb$V2=='E15' & C_Tb$V4==10001,3])
summary(m1 <- glm(cbind(V1, V3) ~ V5*V2, family=binomial(link = "identity"), data=C_Tb, control = list(maxit = 70)))
summary(m1 <- glm(cbind(V1, V3) ~ V5*V2, family=binomial(link = "log"), data=C_Tb, control = list(maxit = 70)))

summary(m1 <- glm(cbind(V1, V3) ~ V5*V2, family=binomial(link = "log"), data=C_Tb, control = list(maxit = 70)))
sum (C_Tb4$V1[C_Tb$V2=='E14'])/sum (C_Tb1$V3[C_Tb$V2=='E2'])


to_bar_plot<-function(mat, est,tr1,tr2){
  v1<-sum (mat$V1[mat$V2==est & mat$V5<=tr1],na.rm = T)
  v2<-sum (mat$V1[mat$V2==est & mat$V5>tr1 & mat$V5<tr2],na.rm = T)
  v3<-sum (mat$V1[mat$V2==est & mat$V5>=tr2],na.rm = T)
  s1<-sum (mat$V3[mat$V2==est & mat$V5<=tr1],na.rm = T)
  s2<-sum (mat$V3[mat$V2==est & mat$V5>tr1 & mat$V5<tr2],na.rm = T)
  s3<-sum (mat$V3[mat$V2==est & mat$V5>=tr2],na.rm = T)
  out<-rbind(c(1,v1,s1),c(2,v2,s2),c(3,v3,s3))
  return(out)
}
E5v<-to_bar_plot(C_Tb6,'E4',25,35)
E15v<-to_bar_plot(C_Tb6,'E15',25,35)

mode(E5v)<-'numeric'
mode(E15v)<-'numeric'
#E5v<-as.data.frame(E5v)
#E15v<-as.data.frame(E15v)

fig<-rbind(cbind(E5v,rep('genes',nrow(E5v))),cbind(E15v,rep('heterochrom',nrow(E15v))))

fig<-cbind(fig,as.numeric(fig[,1])/as.numeric(fig[,2]))
fig1<-as.data.frame(fig)
fig1$V5<-as.numeric(fig[,2])/as.numeric(fig[,3])
library(binom)

bc<-binom.confint(as.numeric(fig[,2]), as.numeric(fig[,3]), conf.level = 0.95, methods = "asymptotic")

ggplot(fig1,aes(x=V1, fill=V4,y=V5))+ #geom_violin(trim=FALSE, fill="gray", adjust=2)+
   geom_bar(position="dodge", stat="identity")+
  geom_errorbar(aes(ymin=bc[,5], ymax=bc[,6]),stat="identity", position=position_dodge())+
  
  #scale_color_manual(name = "Mutation type", values =colorConverter(data$V2))+
  ylab("") + 
  xlab("age bin") + theme_bw()+#xlim(c(genelb,generb))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme( axis.text.y=element_text(size=9),#axis.title.x=element_blank(), axis.text.x=element_blank(),
         #axis.ticks.x=element_blank(),
         axis.title=element_text(size=9),
         legend.text=element_text(size=9))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())#+ theme(legend.position = "none")


sum (C_Tb4$V1[C_Tb$V2=='E15' & C_Tb$V5<=22],na.rm = T)
sum (C_Tb4$V1[C_Tb$V2=='E15' & C_Tb$V5>22 & C_Tb$V5<40],na.rm = T)
sum (C_Tb4$V1[C_Tb$V2=='E15' & C_Tb$V5>=40],na.rm = T)

sum (C_Tb4$V1[C_Tb$V2=='E1' & C_Tb$V5<=22],na.rm = T)
sum (C_Tb4$V1[C_Tb$V2=='E1' & C_Tb$V5>22 & C_Tb$V5<40],na.rm = T)
sum (C_Tb4$V1[C_Tb$V2=='E1' & C_Tb$V5>=40],na.rm = T)




sum (C_Tb4$V1)
summary(m1 <- glm(cbind(non.C_G, phase) ~ age+rate, family=binomial(link = "log"), data = data))
