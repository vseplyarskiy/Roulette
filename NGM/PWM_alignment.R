library(Biostrings)
library(Rtools)
require(ggplot2)
require(ggseqlogo)
library(installr)




setwd("C:/Vova/study/R_scripts_data/R_data/NGM/")
seq<-read.table(file="LoopTM3_al_Loop_cls3", header = F, sep=" ")
seq<-read.table(file="TM_TMall_v1s3", header = F, sep=" ")

sequences <-as.matrix(seq[,1])# c("ERLRELETKIGSLQEDLC","CERLRELETKIGSLQEDL","ERLRECLETKIGSLQEDL")
string_set <- AAStringSet(sequences)
#PSSM <- consensusMatrix(string_set, as.prob = F)
cpwm<-consensusMatrix(sequences)
cpwm<-(cpwm)/colSums(cpwm)
ggseqlogo( cpwm[,c(1:10)], seq_type='aa')+ggtitle("TM7")+theme(plot.title = element_text(hjust = 0.5))


for (tm in 1:7){
  for (cl in 1:6){
    if (cl == 5){
      next
    }
  #  tm=3;cl=6;
    fi<-paste("C:/Vova/study/R_scripts_data/R_data/NGM/freq_TM",tm,"_",cl,sep ="")
    data<-read.table(file=fi, header = F, sep=" ")
  row.names(data)<-data[,1]
  data<-as.matrix(data[,-1])
  tt<-paste("TM",tm," cluster",cl)
  p<-ggseqlogo( data, seq_type='aa')+ggtitle(tt)+theme(plot.title = element_text(hjust = 0.5))
  fi<-paste("C:/Vova/study/R_scripts_data/R_data/NGM/Entropy_plots/freq_TM",tm,"_",cl,".svg",sep ="")
  p
  ggsave(filename=fi, plot = p,device = "svg", path = NULL,
         scale = 1, width = 215, height = 85, units = c("mm"),
         dpi = 300, limitsize = TRUE)
  }
}


