# Random Forest
library(randomForest)
library(ROCR)

# setwd("C:\Users\dinhtuanphan\Google Drive\R\SweatProfile")
df <- read.csv("dC_and_dothers_cutoff_074_woTskin.csv", header = TRUE)

set.seed(123)

df.rf_temp <- df[, -c(1)]
df.rf <- df.rf_temp[, -c(2:3)]

head(df.rf)
rf_model <- randomForest(Dehydration ~ ., df.rf, ntree=500, importantce = T)#, nodesize = 2)
varImpPlot(rf_model, main = 'Variable Importance', cex.main = 0.9)

# roc curve using OOB estimations
OOB.pred <- rf_model$votes[ ,2]
pred.obj <- prediction(OOB.pred, df.rf$Dehydration)
ROC.perf <- performance(pred.obj, "tpr","fpr")

# calculate AUC
perf.AUC <- performance(pred.obj, "auc"); AUC <- perf.AUC@y.values[[1]]
plot(ROC.perf, main = "ROC Curve from OOB Measurements", cex.main = 0.9)
abline(0, 1)
text(0.8, 0.2, paste("AUC = ", format(AUC, digits = 3, scientific = F)))

rf_model
write.csv(rf_model$importance, file="rf_importance.csv")







