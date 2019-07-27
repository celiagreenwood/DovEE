
########################################################
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(caret)
library(DT)
library(MLmetrics)
library(ggplot2)
library(pROC)

#########################################3
# Datasets are called train.x, train.y, test.x and test.y
# train.x and train.y are used for model building
# test.x and test.y are used for model testing
# train.x and test.x contain all the features.
# test.x and test.y contain a variable called "Category" which indicates
#   whether the patient has a cancer or a Benign condition
# We created also variables symptom1 and symptom2 to facilitate 
#   assessments of performance
# caret package is used.
##############################################################

symtoms1 <- "Benign"
symtoms2 <- c("ECH","OCHG","ECL","LMP","OCG","OCLG")

cctrl1 <- trainControl(method = "LOOCV",
                       classProbs = TRUE, 
                       summaryFunction = twoClassSummary, 
                       seeds = seeds)


##################################################################
#select features using stepwise selection
load(paste0(PATH.new.save  ,"/name_fts.RData"),verbose=TRUE )
##################################################################

wsrf.oc.ec <- train(
  train.x, train.y, 
  method = "wsrf", 
  trControl = cctrl1,
  ntree = 13,
  parallel = FALSE,
  metric='ROC',
  preProc = c("center", "scale"  ) )

tmp <- predict(wsrf.oc.ec, test.x, type = "prob")

train.pred <- predict(wsrf.oc.ec, train.x, type="prob")

roc_obj <-roc(as.factor(test.y),tmp[,which(colnames(tmp)=="symtoms1")], levels=c("symtoms1","symtoms2"))

auc.val <- auc(roc_obj)

print(auc.val)

train.prob <- data.frame(
  
  "patientID" = data.oc.ec$patientID,
  "Diagnosis" =  data.oc.ec$Category,
  Model.prob.cancer = round( train.pred$symtoms2 ,3),
  Model.prob.benign = round( train.pred$symtoms1,3),
  
  Datatype="train"
)

test.prob <- data.frame(
  
  "patientID" = sind.data.test$patientID,
  "Diagnosis" =  sind.data.test$Category,
  Model.prob.cancer = round( tmp$symtoms2,3),
  Model.prob.benign = round( tmp$symtoms1,3),
  
  Datatype="test"
)

total.probs <- rbind(test.prob, train.prob)


write.table(total.probs,"manuscriptS10.txt",sep="\t",row.names=FALSE)

############################################################################
