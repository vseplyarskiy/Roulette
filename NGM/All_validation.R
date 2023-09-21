library(ggplot2)
setwd("C:/Vova/study/R_scripts_data/R_data/NGM/")

data<-read.table(file="cont_distr_density", sep=" ");
p<-ggplot(data,aes(x = V3,fill=V2))+#scale_x_discrete(limits = positions)+#geom_bar()
  geom_histogram( position="dodge2")+ #scale_size(range = c(1, 1))+
  theme_bw()+
 
  ylab("Density") +
  scale_fill_discrete(name="Mutation\nType") +
  #scale_x_log10(breaks=x_breaks)
  scale_x_continuous(trans = 'log2',
                     breaks =  c(1,2.5,5,10,25,50,100,200),lim=c(1,200))+
  #scale_x_continuous(trans='log2') +
  #scale_y_continuous(trans='log2')+
  xlab("Fold change") +
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
          panel.background = element_blank())#+theme(legend.position = "none")
p

data<-read.table(file="cont_distr_AATTT_A", sep=" ");
p<-ggplot(data[c(-73,-74),],aes(x = -log(V2),y=-log(V3)))+#scale_x_discrete(limits = positions)+#geom_bar()
  geom_point()+ #scale_size(range = c(1, 1))+
  scale_x_continuous(breaks =  c(0,1,2),limits = c(0,2.7))+
  scale_y_continuous(breaks =  c(0,1,2),limits = c(0,2.7))+
  
  theme_bw()+
  
  ylab("Observed") +
  xlab("Predicted") +
  
 # scale_fill_discrete(name="Mutation\nType") +
  #scale_x_log10(breaks=x_breaks)
  #scale_x_continuous(trans = 'log2',                  breaks =  c(1,2.5,5,10,25,50,100,200),lim=c(1,200))+
  #scale_x_continuous(trans='log2') +
  #scale_y_continuous(trans='log2')+
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
          panel.background = element_blank())#+theme(legend.position = "none")
p



library(reshape2)

data1<-read.table(file="spatial_validation", sep=" ");
colnames(data1)<-c('position','observed','our model','Carlson et al')
df <- melt(data1,id=1)  #the function melt reshapes it from wide to long
#df[,2] <- c('observed','our model','Carlson et al') #add a rowid identifying variable
p<-ggplot(df,aes(x = position,y=value,color = variable))+#scale_x_discrete(limits = positions)+#geom_bar()
  geom_point(size=0.3,alpha=0.3)+
  geom_smooth(span = 0.03,method = 'loess',n=100,size=.8,alpha=0.7, se=FALSE)+
  xlim(c(5000000,31000000))+
  ylim(c(0,0.25))+
  
  #scale_size(range = c(1, 1))+
  theme_bw()+
  #xlim(c(0,5000))+ylim(c(0,5000))+
  #geom_abline(intercept = 0, slope = 0, color="grey", linetype="dashed", size=.5)+
  ylab("C>G frequency") +
  scale_color_discrete(name="Model") + 
  #scale_x_continuous(trans='log2') +scale_y_continuous(trans='log2')+
  xlab("Position, chr8") +
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
          panel.background = element_blank())#+theme(legend.position = "none")
p

data1<-read.table(file="Effect_of_R_as_GTCGT_T_inter_met5", sep=" ");
colnames(data1)[1:3]<-c('direction','observed','expected')
df <- melt(data1[,1:3],id=1)  #the function melt reshapes it from wide to long

p<-ggplot(df,aes(x = direction,y=-log(value),fill = variable))+#scale_x_discrete(limits = positions)+#geom_bar()
  geom_bar(stat = "identity",position = "dodge2")+
   
  #scale_size(range = c(1, 1))+
  theme_bw()+
  #xlim(c(0,5000))+ylim(c(0,5000))+
  #geom_abline(intercept = 0, slope = 0, color="grey", linetype="dashed", size=.5)+
  ylab("CpG frequency per site") +
  scale_fill_discrete(name="") + 
  #scale_x_continuous(trans='log2') +scale_y_continuous(trans='log2')+
  xlab("Replication fork direction") +
  theme(legend.text = element_text(size = 15),legend.title = element_text( size=15, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(  axis.text=element_text(size = 15),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          #axis.ticks.y=element_blank(),
          #axis.line.x = ,
          #axis.line.y = black,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          #panel.border = element_blank(),
          panel.background = element_blank())#+theme(legend.position = "none")
p
