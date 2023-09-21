
library(ggplot2)
library(car)

setwd("C:/Vova/Vova/study/R_scripts_data/R_data/Signatures_germline//")

data1<-read.table(file="mut_old_30kb_as.rt_bi",header = TRUE,  sep=" ")
data<-log2(data1)
data$RT_as<-data1$RT_as
#cormat <- round(cor(data,use="complete.obs", method = "spearman"),4)
cormat1<-round(cor(data,use="complete.obs", method = "spearman"),4)

nam<-as.vector(names(data1))

d_dir=data[FALSE,FALSE];

x <- nam[99]

pred <- lm(log2(data1[, x]) ~ RT_as, data = data1)
summary(pred)
#effect_plot(pred, pred = RT_as, interval = TRUE, plot.points = TRUE)
ggplot(data1,aes(y=log2(data1[, x]),x=RT_as))+stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=4, direction=-1)  +xlim(c(-.8,.8))+ylab("Asymetry")+
  xlab("Direction of Transcription")+
  theme_bw()+geom_smooth(method="lm", color="black")

all_anova <- function(y, data1) {
  x <- nam[y]
  pred <- lm(log2((data1[, x])) ~ RT_as, data = data1)
  #Anova(pred, type = "II")
  #sp <- Anova(pred, type = "II")
  #as.vector(sp$`Sum Sq`)
  val <- pred$coefficients[-1]
  write.table(
    val,
    file = x,
    append = FALSE,
    quote = FALSE,
    sep = " ",
    eol = "\n",
    na = "NA",
    dec = ".",
    row.names = FALSE,
    col.names = TRUE,
    qmethod = c("escape", "double"),
    fileEncoding = ""
  )
}
setwd("C:/Vova/Vova/study/R_scripts_data/R_data/NGM/mut_model/R-as")

for (y in 3:98) {
  all_anova(y, data1)
}
qplot(data1$RT_as)
quantile(data1$RT_as,0.25)
var(data1$RT_as)
