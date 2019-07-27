##########################################################3
# This code chooses an optimal feature subset.
# Results are saved and used in the primary model building
##########################################################33
library(caret)
library(pROC)

# the dataset is called data.all.
# Category2 contains the diagnosis information.

####################################################################
cv_pred <- function(data.all, name.fts, n_folds = 3 ){
  
  set.seed(998)
  cv.folds <- createFolds(data.all$Category2, k = n_folds, list = TRUE, returnTrain = FALSE)
  
  #the probability matrix that are saved
  prob.fts <- data.frame(patientID=data.all$patientID, Category2 = data.all$Category2,
                         probs = rep(0,nrow(data.all)) )
  #for each fold we use it as test cases and train model on other folds 
  for(i_folds in 1:n_folds){
    
    #using selected features, and selected cv-folds
    train.x <- data.all[-cv.folds[[i_folds]],match(name.fts, colnames(data.all))]
    test.x <- data.all[cv.folds[[i_folds]],match(name.fts, colnames(data.all))]
    train.y <- data.all$Category2[-cv.folds[[i_folds]]]
    test.y <- data.all$Category2[cv.folds[[i_folds]]]
    
    seeds <- vector(mode = "list", length = nrow(data.oc.ec) + 1)
    seeds <- lapply(seeds, function(x) 1:20)
    cctrl1 <- trainControl(method = "LOOCV",
                           classProbs = TRUE, 
                           summaryFunction = twoClassSummary, 
                           seeds = seeds)
    wsrf.oc.ec <- caret::train(
      train.x, train.y, 
      method = "wsrf", 
      trControl = cctrl1,
      ntree = 10,
      parallel = FALSE,
      metric = "ROC", 
      preProc = c("center", "scale"  ) )  
    
    
    tmp <- predict(wsrf.oc.ec, test.x, type = "prob")
    prob.fts$probs[cv.folds[[i_folds]]] <- tmp$symtoms2
  }
  
  ######now collect the probability to be cancer. Then calculate AUC
  roc_obj <- roc(data.all$Category2=="symtoms2", prob.fts$probs)
  auc.val <- pROC::auc(roc_obj)
  
  ##output auc
  list(auc =auc.val, FN = sum( (prob.fts$probs<0.5)*(data.all$Category2=="symtoms2")), 
       prob=prob.fts$probs, roc.obj = roc_obj  )
}

##################################################################
##################################################################

auc.fts <- crit.auc <- 0.6

#start initial value with one features in top 10
name.fts <- all.fts[sample( c(1:10), size=1, replace=FALSE   )]

while(auc.fts >= crit.auc){
  
  crit.auc <- auc.fts
  ct.features <- name.fts
  exist.fts <- all.fts[is.na(match(all.fts,ct.features))]
  
  id.fts <- sample( c(1:length( exist.fts )), size=1, replace=FALSE   )
    #the complete features we use for model fitting
    name.fts <- c(ct.features, exist.fts[id.fts])
    #the auc for this feature on training cases
    cv.fit <- cv_pred(data.oc.ec, name.fts, n_folds = 10 )
    
    auc.fts <- cv.fit$auc
    FN.fts  <- cv.fit$FN
 
    if(auc.fts>crit.auc){
      #change line and output
      cat(sprintf("\n"), file= "optimal_auc.txt", append=TRUE )
      cat(auc.fts, file= "optimal_auc.txt", append=TRUE, sep = "\t")
      
      #change line and output
      cat(sprintf("\n"), file= "optimal_fts.txt", append=TRUE )
      cat(id.fts, file= "optimal_fts.txt", append=TRUE, sep = "\t")
      
      #change line and output false positive
      cat(sprintf("\n"), file= "optimal_FN.txt", append=TRUE )
      cat(FN.fts, file= "optimal_FN.txt", append=TRUE, sep = "\t")
    }

    #if not better, switch back
    else{
      name.fts <- ct.features 
    }

}

cat(name.fts, file= "optimal.txt", append=TRUE, sep = "\t")
#save(name.fts, file=paste0(PATH.save,"/name_fts.RData" ) )





