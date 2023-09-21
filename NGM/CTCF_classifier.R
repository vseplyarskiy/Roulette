library(randomForest)
setwd("C:/Vova/study/R_scripts_data/R_data/NGM/")
library(pROC)
library(epitools)


data<-read.table(file="class_CTCF_only_CTCFPLS", header = F, sep=" ")
data1<-data[,c(-14)]
data1[,1] <-as.factor(data[,14]) #ifelse(df$quality>5,1,0)

index <- sample(1:nrow(data1),size = 0.8*nrow(data1))
train <- data1[index,]
test <- data1[-index,]
rf_model <- randomForest(V1 ~., data = train)
rf_prediction <- predict(rf_model, test, type = "prob")
ROC_rf <- roc(test$V1, rf_prediction[,2])


data2<-read.table(file="class_PLS_CTCFPLS", header = F, sep=" ")
data3<-data2[,c(-14)]
data3[,1] <-as.factor(data2[,14]) #ifelse(df$quality>5,1,0)

index <- sample(1:nrow(data3),size = 0.8*nrow(data3))
train <- data3[index,]
test <- data3[-index,]
rf_model <- randomForest(V1 ~., data = train)
rf_prediction <- predict(rf_model, test, type = "prob")
ROC_rf2 <- roc(test$V1, rf_prediction[,2])



plot(ROC_rf, col = "green", main = "CTCF classification")
lines(ROC_rf2, col = "red")
ROC_rf_auc <- auc(ROC_rf)
ROC_rf2_auc <- auc(ROC_rf2)
# Area Under Curve (AUC) for each ROC curve (higher -> better)

# plot ROC curves
plot(ROC_rf, col = "green", main = "ROC For Random Forest (GREEN) vs Logistic Regression (RED)")
lines(ROC_lr, col = "red")
