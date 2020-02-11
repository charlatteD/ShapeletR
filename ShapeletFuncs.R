#-------------------------------------------------------------------------------------------------------------#
library(rdist)
library(dtwclust)
library(randomForest)
library(TSPred)
library(pROC)
library(matrixStats)
library(qlcMatrix)
library(TSEntropies)
library(dplyr)
library(plyr)
library(zoo)
library(gbm)
library(rpart)
library(rpart.plot)
library(lme4)
library(ggplot2)
#-------------------------------------------------------------------------------------------------------------#
#set the color
set.seed(222)
colorlist = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
color1=sample(colorlist, 10)
#-------------------------------------------------------------------------------------------------------------#
#Insert shape to simulated data
getDataWithShapeV2<-function(ExampleShape1,
                             data,
                             perc,
                             period1,
                             period2,
                             prob1,
                             prob2,
                             ChangeMagnitude=TRUE){
  
  lExampleShape<-length(ExampleShape1) #Length of the shape
  
  #(Example) Extract the characteristic of the time-series 
  lrow<-nrow(data) 
  lcol<-ncol(data)
  lcol2<-lcol-lExampleShape
  
  totalObsWShape<-floor(lrow*perc)
  
  #randomly select n*perc columns
  randcol1 <- sample(1:period1,lrow*perc*0.1,replace=TRUE) #Morning time
  randcol2 <- sample(period2:lcol2,lrow*perc*0.9,replace=TRUE) #After/night time, start from 75 to avoid morning time   
  randcol<-c(randcol1,randcol2)#randomly select n*perc columns
  
  #randomly select n*perc columns (equal to the number of random columns)
  randrow <- sample(1:lrow,length(randcol),replace=FALSE) #randomly select n*perc rows 
  
  #Insert the example shape to the location (a:Control the height of the shape)
  for(i in 1:length(randcol)){
    a<-ifelse(ChangeMagnitude,runif(1,min=0.01,max=3),1)
    data[randrow[i],randcol[i]:c(randcol[i]+lExampleShape-1)]<-a*jitter(ExampleShape1)
  }
  data<-rbind(data[randrow,],data[-randrow,])#Re-arrange the order of the T-S 
  
  #Assign the value of Y based on the probability of P(y=1)
  Y<-rep(0,nrow(data))
  prob<-c(rep(prob1,length(randcol)),rep(prob2,(nrow(data)-length(randcol))))
  Y<-rbinom(nrow(data),1,prob)
  
  data1<-data[which(Y==1),]
  data0<-data[which(Y==0),]
  
  return(list(data,Y,prob,nrow(data),length(randrow),randcol,data1,data0))
}
#Normalize the data
minmax<-function(input){
  normalized = (input-min(input))/(max(input)-min(input))
  return(normalized)
}
#-------------------------------------------------------------------------------------------------------------#
#Randomly generate a priori shapelets of different scale from the a priori shapes library
#inc: Increment factor, lshp: shapelet lengths, func: shape function 
getshplist<-function(inc,lshp,func){
  lseq<-seq(0,1,by=inc)
  expgrid<-expand.grid(lseq,lseq)
  expgrid2<-expgrid[expgrid[,2]>expgrid[,1],]
  
  mat<-matrix(NA,nrow=dim(expgrid2)[1],ncol=lshp)
  for(i in 1:dim(expgrid2)[1]){
    mat[i,]<-func(lshp-1,expgrid2[i,1],expgrid2[i,2])[,2]
  }
  return(mat)
}
#-------------------------------------------------------------------------------------------------------------#
#Create a sliding window dataset (lshp: shapelet lengths)
#input: Original time series, lshp: shapelet lengths
GetSlidingWindow<-function(input,lshp){
  inputlist <- split(input, seq(nrow(input))) #split the data
  inputlist2<-lapply(inputlist,as.numeric) #Change from data frame to vector
  mwlist<-lapply(inputlist2,slidingWindows,swSize=lshp) #Create sliding window
  return(mwlist)
}
#-------------------------------------------------------------------------------------------------------------#
#Slope function (for feature-based initialization method)
slopefun<-function(vector1){
  lvector1<-length(vector1)
  slope1<-(vector1[lvector1]-vector1[1])/(lvector1-1)
  return(abs(slope1))
}

#Lag function (for feature-based initialization method)
#ar<-function(vector1){
#  a<-summary(lm(vector1~lag(vector1)))$coefficient[2,1]
#  return(a)
#}
#-------------------------------------------------------------------------------------------------------------#
# Select the sliding window with the top characteristics and return initialized shapelets (for stratified analysis)
#Input: Input (TrainX), lshp: shapelet lengths
GetUniqueSubsequence<-function(input,lshp){
  
  #Transpose time series: row: observations, column: time stamps 
  input_transpose<-t(input)
  input_transpose<-as.data.frame(input_transpose)
  
  #Rolling mean 
  mean1<-rollapply(input_transpose, lshp, mean)
  loc_vector<-which(mean1==max(mean1),arr.ind=T)
  loc_vector<-rbind(loc_vector,which(mean1==min(mean1),arr.ind=T))
  #Rolling median 
  median1<-rollapply(input_transpose, lshp, median)
  loc_vector<-rbind(loc_vector,which(median1==min(median1),arr.ind=T))
  loc_vector<-rbind(loc_vector,which(median1==max(median1),arr.ind=T))
  #SD 
  SD1<-rollapply(input_transpose, lshp, sd)
  loc_vector<-rbind(loc_vector,which(SD1==max(SD1),arr.ind=T))
  #Max
  max1<-rollapply(input_transpose, lshp, max)
  loc_vector<-rbind(loc_vector,which(max1==max(max1),arr.ind=T))
  #Min
  min1<-rollapply(input_transpose, lshp, min)
  loc_vector<-rbind(loc_vector,which(min1==min(min1),arr.ind=T))
  #Slope fun
  slopefun1<-rollapply(input_transpose, lshp, function(x) slopefun(x))
  loc_vector<-rbind(loc_vector,which(slopefun1==max(slopefun1),arr.ind=T))
  #AR(1)
  #ar1<-rollapply(input_transpose, lshp, function(x) ar(x))
  #loc_vector<-rbind(loc_vector,which(ar1==max(ar1),arr.ind=T))
  
  #Return the shapes with the unique characteristics 
  Cand_shp<-matrix(NA,nrow=dim(loc_vector)[1],ncol=lshp)
  for(i in 1:dim(loc_vector)[1]){
    Cand_shp[i,]<-input_transpose[loc_vector[i,1]:(loc_vector[i,1]+lshp-1),loc_vector[i,2]]
  }
  
  Cand_shp<-as.data.frame(Cand_shp)
  #Get rid of the redundant rows
  Cand_shp2<- Cand_shp %>% 
    distinct 
  
  return(Cand_shp2)
}
#-------------------------------------------------------------------------------------------------------------#
#Apply K mean clustering to the dataset after the sliding window transformation 
#mwlist: List of the moving window by observation, nshp: number of shapelets
GetKmean<-function(mwlist,nshp){
  mwdf<-do.call(rbind,mwlist)
  
  #Randomly selected 10,000 rows with replaces
  selectedRow<-sample(1:nrow(mwdf),10000,replace=TRUE)
  mwdfreduced<-mwdf[selectedRow,]
  
  km1<-kmeans(mwdfreduced, centers=nshp, nstart = 100, iter.max = 30, algorithm="MacQueen")
  center<-km1$centers
  
  return(center)
}
#-------------------------------------------------------------------------------------------------------------#
#Calculate the Euclidean distance (shp: shp matrix)
#mwlist: List of the moving window by observation, shp: candidate shapelets
GetDistanceVectorsE<-function(mwlist,shp){
  dist2shp<-lapply(mwlist,cdist,shp) #Calculate distance between the list to the shapelet 
  
  mindf<-matrix(NA,nrow=length(mwlist),ncol=dim(shp)[1])
  minlocdf<-matrix(NA,nrow=length(mwlist),ncol=dim(shp)[1])
  for(i in 1:dim(shp)[1]){
    mindf[,i]<-as.numeric(sapply(dist2shp, function(x) min(x[,i]))) #Select the min distance
    minlocdf[,i]<-as.numeric(sapply(dist2shp, function(x) which.min(x[,i]))) #Select the location of the min distance
  }
  return(list(mindf,minlocdf))
}
#-------------------------------------------------------------------------------------------------------------#
#Calculate the DTW distance (shp: shp matrix)
GetDistanceVectorsDTW<-function(mwlist,shp){
  dist2shp<-lapply(mwlist,dtw_lb,y=shp,window.size=dim(shp)[2]/2) #Calculate distance between the list to the shapelet 
  
  mindf<-matrix(NA,nrow=length(mwlist),ncol=dim(shp)[1])
  minlocdf<-matrix(NA,nrow=length(mwlist),ncol=dim(shp)[1])
  for(i in 1:dim(shp)[1]){
    mindf[,i]<-as.numeric(sapply(dist2shp, function(x) min(x[,i]))) #Select the min distance
    minlocdf[,i]<-as.numeric(sapply(dist2shp, function(x) which.min(x[,i]))) #Select the location of the min distance
  }
  return(list(mindf,minlocdf))
}
#-------------------------------------------------------------------------------------------------------------#
#Calculate the DTW distance and return predicted groups 
getpredY2<-function(xtest,xtrain,ytrain,ytest){
  miny<-rep(NA,dim(xtest)[1])
  loc<-dtw_lb(xtrain,xtest,window.size=15)
  miny<-ytrain[apply(loc,2,which.min)]
  
  a1<-multiclass.roc(as.numeric(ytest),as.numeric(miny),quiet=T)
  confmatrix<-table(ytest, miny)
  accuracy <- sum(diag(confmatrix)) / sum(confmatrix)
  tab1<-c(a1$auc,accuracy)
  tab1
  return(miny)
}
#-------------------------------------------------------------------------------------------------------------#
Shapelet_Kmeans<-function(TrainX,TrainY,TestX,TestY,lshp,ncentroid,bag.fraction1,train.fraction1){
  
  #Return number of shapelet
  #If lshp <- c(10,15), meaning two shapelets, one with length of 10 and another with length of 15
  nshp<-length(lshp)
  
  #Get an empty list for shapelets, distance matrix (Train & Test), location of the minimum distance (Train & Test)
  shplists <- list()
  inputDistTrain <- list()
  inputLocTrain <- list()
  inputDistTest <- list()
  inputlocTest <- list()
  
  for(i in 1:nshp){
    #Get sliding windown
    SWtrain<-GetSlidingWindow(minmax(TrainX),lshp[i])
    SWtest<-GetSlidingWindow(minmax(TestX),lshp[i])
    
    #Get n centroids (shapelet (candidate))
    shplists[[i]]<-GetKmean(SWtrain,ncentroid)
    
    #Get the minimum distance and the location of minimum distance between the data and shapelet (candidate)
    inputDTrain1<-GetDistanceVectorsE(SWtrain,shplists[[i]])
    inputDistTrain[[i]]<-inputDTrain1[[1]]
    inputLocTrain[[i]]<-inputDTrain1[[2]]
    inputDTest1<-GetDistanceVectorsE(SWtest,shplists[[i]])
    inputDistTest[[i]]<-inputDTest1[[1]]
    inputlocTest[[i]]<-inputDTest1[[2]]
    
    #For rbind.fill function
    shplists[[i]]<-as.data.frame(shplists[[i]])
  }
  #rbind across list (rbind.fill: for data frames with different lengths of the row)
  shplists1<-do.call(rbind.fill,shplists)
  
  inputDistTrain1<-do.call(cbind,inputDistTrain)
  inputLocTrain1<-do.call(cbind,inputLocTrain)
  inputDistTest1<-do.call(cbind,inputDistTest)
  inputlocTest1<-do.call(cbind,inputlocTest)
  
  #Generate the outputs 
  output1<-ModelPredictionBV2(inputDistTrain1,TrainY,inputDistTest1,TestY,bag.fraction1,train.fraction1)
  
  #AUC/ACC
  AUC_ACC<-output1[[1]]
  
  #gbm.model1 
  gbm.model1<-output1[[3]]
  
  #gbm.best.iter1 
  gbm.best.iter1<-output1[[4]]
  
  return(list(AUC_ACC=AUC_ACC,shplists1=shplists1,
              gbm.model1=gbm.model1,gbm.best.iter1=gbm.best.iter1,
              inputDistTrain=inputDistTrain1,inputDistTest=inputDistTest1,
              inputLocTrain=inputLocTrain1,inputlocTest=inputlocTest1))
}

#Shapelet_Kmeans(TrainX,TrainY,TestX,TestY,c(10,15),10)
#-------------------------------------------------------------------------------------------------------------#
#Proposed shape (output)
Shapelet_ProposedShapelet<-function(TrainX,TrainY,TestX,TestY,MinMaxDiff,lshp,shpfunc,bag.fraction1,train.fraction1){
  #Return number of shapelet
  #If lshp <- c(10,15), meaning two shapelets, one with length of 10 and another with length of 15
  nshp<-length(lshp)
  #Return number of proposed a priori functions
  #If shpfunc <- c("shp1","shp2"), meaning two a priori functions
  nshpfunc<-length(shpfunc)
  #Get the combination of number of shapelet and number of proposed functions for evaluations 
  nshp_nfunc_combine <- expand.grid(x = 1:nshp, y = 1:nshpfunc)
  
  #Get an empty list for shapelets, distance matrix (Train & Test), location of the minimum distance (Train & Test)
  shplists <- list()
  inputDistTrain <- list()
  inputLocTrain <- list()
  inputDistTest <- list()
  inputlocTest <- list()
  
  for(i in 1:dim(nshp_nfunc_combine)[1]){ 
    #For each combination 
    lshp_loc <- nshp_nfunc_combine[i,1] 
    lshpfuc_loc <- nshp_nfunc_combine[i,2]
    
    lshp1<-lshp[lshp_loc]#Get number of shapelets
    shpfunc1<-shpfunc[lshpfuc_loc]#Get propposed shapelet function
    shplists[[i]]<-getshplist(MinMaxDiff,lshp1,get(shpfunc1)) #Get propposed shapelets
    
    #Get sliding windown
    SWtrain<-GetSlidingWindow(minmax(TrainX),lshp1)
    SWtest<-GetSlidingWindow(minmax(TestX),lshp1)
    
    #Get the minimum distance and the location of minimum distance
    inputDTrain1<-GetDistanceVectorsE(SWtrain,shplists[[i]])
    inputDistTrain[[i]]<-inputDTrain1[[1]]
    inputLocTrain[[i]]<-inputDTrain1[[2]]
    inputDTest1<-GetDistanceVectorsE(SWtest,shplists[[i]])
    inputDistTest[[i]]<-inputDTest1[[1]]
    inputlocTest[[i]]<-inputDTest1[[2]]
    
    #For rbind.fill function
    shplists[[i]]<-as.data.frame(shplists[[i]])
  }
  #rbind across list (rbind.fill: for data frames with different lengths of the row)
  shplists1<-do.call(rbind.fill,shplists)
  
  inputDistTrain1<-do.call(cbind,inputDistTrain)
  inputLocTrain1<-do.call(cbind,inputLocTrain)
  inputDistTest1<-do.call(cbind,inputDistTest)
  inputlocTest1<-do.call(cbind,inputlocTest)
  
  #Generate the outputs 
  output1<-ModelPredictionBV2(inputDistTrain1,TrainY,inputDistTest1,TestY,bag.fraction1,train.fraction1)
  
  #AUC/ACC
  AUC_ACC<-output1[[1]]
  
  #gbm.model1 
  gbm.model1<-output1[[3]]
  
  #gbm.best.iter1 
  gbm.best.iter1<-output1[[4]]
  
  return(list(AUC_ACC=AUC_ACC,shplists1=shplists1,
              gbm.model1=gbm.model1,gbm.best.iter1=gbm.best.iter1,
              inputDistTrain=inputDistTrain1,inputDistTest=inputDistTest1,
              inputLocTrain=inputLocTrain1,inputlocTest=inputlocTest1))
}
#Shapelet_ProposedShapelet(TrainX,TrainY,TestX,TestY,0.3,c(10,30),c("shp1","shp5"))
#-------------------------------------------------------------------------------------------------------------#
#Shapes with unique features (output)
Shapelet_UniqueFeatures<-function(TrainX,TrainY,TestX,TestY,lshp,bag.fraction1,train.fraction1){
  
  temp_position<-data.frame(TrainY,1:length(unlist(TrainY)))
  colnames(temp_position)<-c("TrainY","Position")
  uniqueY<-unlist(unique(TrainY))
  nshp<-length(lshp)
  nshp_nuniqueY_combine <- expand.grid(x = 1:nshp, y = 1:length(uniqueY))
  
  shplists <- list()
  
  inputDistTrain <- list()
  inputLocTrain <- list()
  inputDistTest <- list()
  inputlocTest <- list()
  
  TrainX<-minmax(TrainX)
  TestX<-minmax(TestX)
  
  for(i in 1:dim(nshp_nuniqueY_combine)[1]){
    
    lshp_loc <- nshp_nuniqueY_combine[i,1]
    lY_loc <- nshp_nuniqueY_combine[i,2]
    
    #Get propposed shapelets
    lshp1<-lshp[lshp_loc]
    uniqueY1<-uniqueY[lY_loc]
    
    Pos1<-temp_position$Position[temp_position$TrainY==uniqueY1]
    TrainX1<-TrainX[Pos1,]
    TestX1<-TestX[Pos1,]
    shplists[[i]]<-GetUniqueSubsequence(TrainX1,lshp1)
    
    #Get sliding windown
    SWtrain<-GetSlidingWindow(TrainX,lshp1)
    SWtest<-GetSlidingWindow(TestX,lshp1)
    
    #Get the minimum distance and the location of minimum distance
    inputDTrain1<-GetDistanceVectorsE(SWtrain,shplists[[i]])
    inputDistTrain[[i]]<-inputDTrain1[[1]]
    inputLocTrain[[i]]<-inputDTrain1[[2]]
    inputDTest1<-GetDistanceVectorsE(SWtest,shplists[[i]])
    inputDistTest[[i]]<-inputDTest1[[1]]
    inputlocTest[[i]]<-inputDTest1[[2]]
    
    #For rbind.fill function
    shplists[[i]]<-as.data.frame(shplists[[i]])
  }
  #rbind across list (rbind.fill: for data frames with different lengths of the row)
  shplists1<-do.call(rbind.fill,shplists)
  
  inputDistTrain1<-do.call(cbind,inputDistTrain)
  inputLocTrain1<-do.call(cbind,inputLocTrain)
  inputDistTest1<-do.call(cbind,inputDistTest)
  inputlocTest1<-do.call(cbind,inputlocTest)
  
  #Generate the outputs 
  output1<-ModelPredictionBV2(inputDistTrain1,TrainY,inputDistTest1,TestY,bag.fraction1,train.fraction1)
  
  #AUC/ACC
  AUC_ACC<-output1[[1]]
  
  #gbm.model1 
  gbm.model1<-output1[[3]]
  
  #gbm.best.iter1 
  gbm.best.iter1<-output1[[4]]
  
  return(list(AUC_ACC=AUC_ACC,shplists1=shplists1,
              gbm.model1=gbm.model1,gbm.best.iter1=gbm.best.iter1,
              inputDistTrain=inputDistTrain1,inputDistTest=inputDistTest1,
              inputLocTrain=inputLocTrain1,inputlocTest=inputlocTest1))
}
#-------------------------------------------------------------------------------------------------------------#
#Fit into a model (Binary response, inverse distance matrix)
ModelPredictionBV2<-function(inputDistTrain,outcomeTrain,inputDistTest,outcomeTest,bag.fraction1,train.fraction1,interdep=1,randomState=1){
  
  #Take the inverse of the shapelet distance 
  inputDistTrain_inv<-1/(inputDistTrain+0.1)
  inputDistTest_inv<-1/(inputDistTest+0.1)
  
  #Combine the outcome and distance data and assignt the variable name of the outcome to "y"
  #Training set
  finaldtTrain<-data.frame(outcomeTrain,inputDistTrain_inv)
  colnames(finaldtTrain)[1] <- "y" 
  
  #Test set
  finaldtTest<-data.frame(outcomeTest,inputDistTest_inv)
  colnames(finaldtTest)[1] <- "y" 
  
  set.seed(randomState)
  #GBM 
  gbm1<-gbm(y ~ ., data = finaldtTrain,
            distribution = "bernoulli", n.trees = 1000, 
            shrinkage = 0.01,interaction.depth = interdep, 
            bag.fraction = bag.fraction1, train.fraction = train.fraction1,  
            n.minobsinnode = 10, cv.folds = 3,n.cores = 3)  
  best.iter <- gbm.perf(gbm1, method = "cv", plot.it = FALSE)
  
  #GBM (CV:AUC)
  cv_prediction<-exp(gbm1$cv.fitted)/(1+exp(gbm1$cv.fitted))
  a.gbmcv<-roc(unlist(outcomeTrain[1:length(cv_prediction)]),as.numeric(cv_prediction),quiet=T)
  a.gbmcv.out<-coords(a.gbmcv, "best", ret=c("threshold", "specificity", "sensitivity", "accuracy",
                                             "precision", "recall"), best.method=c("youden"), transpose = TRUE)
  a.gbmcv.out.f1<-2*a.gbmcv.out["precision"]*a.gbmcv.out["recall"]/(a.gbmcv.out["precision"]+a.gbmcv.out["recall"])
  
  #GBM (Train: AUC,accuracy)
  PredTrainY_gbm1<-predict(gbm1,newdata=finaldtTrain, n.trees = best.iter,type="response")
  a.gbm1<-roc(unlist(outcomeTrain),as.numeric(PredTrainY_gbm1),quiet=T)
  a.gbm1.out<-coords(a.gbm1, "best", ret=c("threshold", "specificity", "sensitivity", "accuracy",
                                           "precision", "recall"), best.method=c("youden"), transpose = TRUE)
  a.gbm1.out.f1<-2*a.gbm1.out["precision"]*a.gbm1.out["recall"]/(a.gbm1.out["precision"]+a.gbm1.out["recall"])
  #GBM (Test: AUC,accuracy)
  PredTestY_gbm1<-predict(gbm1,newdata=finaldtTest, n.trees = best.iter,type="response")
  a.gbm2<-roc(as.numeric(unlist(outcomeTest)),as.numeric(PredTestY_gbm1),quiet=T)
  a.gbm2.out<-coords(a.gbm2, "best", ret=c("threshold", "specificity", "sensitivity", "accuracy",
                                           "precision", "recall"), best.method=c("youden"), transpose = TRUE)
  a.gbm2.out.f1<-2*a.gbm2.out["precision"]*a.gbm2.out["recall"]/(a.gbm2.out["precision"]+a.gbm2.out["recall"])
  
  #Factorize the outcome for random forest 
  finaldtTrain$y<-as.factor(finaldtTrain$y)
  #Fill in missing value
  finaldtTrain.roughfix <- na.roughfix(finaldtTrain)
  finaldtTest.roughfix <- na.roughfix(finaldtTest)
  #Random forest 
  set.seed(randomState)
  model1 <- randomForest(y ~ ., data = finaldtTrain.roughfix, importance = TRUE)
  
  #random forest (Train: AUC,accuracy)
  PredTrainY1<-predict(model1,newdata=finaldtTrain.roughfix,type="prob")
  a1<-roc(unlist(outcomeTrain),as.numeric(PredTrainY1[,2]),quiet=T)
  a1.out<-coords(a1, "best", ret=c("threshold", "specificity", "sensitivity", "accuracy",
                                   "precision", "recall"), best.method=c("youden"), transpose = TRUE)
  
  #random forest (Test: AUC,accuracy)
  PredTestY1<-predict(model1,newdata=finaldtTest.roughfix,type="prob")
  a2<-roc(unlist(outcomeTest),as.numeric(PredTestY1[,2]),quiet=T)
  a2.out<-coords(a2, "best", ret=c("threshold", "specificity", "sensitivity", "accuracy",
                                   "precision", "recall"), best.method=c("youden"), transpose = TRUE)
  #random forest (Variable importance list)
  implist<-as.data.frame(importance(model1))
  implist<-implist[c("MeanDecreaseAccuracy","MeanDecreaseGini")]
  
  #Logistic regression 
  model1 <- glm(y~ ., data = finaldtTrain, family = binomial)
  
  #Logistic (Train: AUC,accuracy)
  PredTrainY2<-predict(model1,finaldtTrain,type="response")
  a3<-roc(unlist(outcomeTrain),as.numeric(PredTrainY2),quiet=T)
  a3.out<-coords(a3, "best", ret=c("threshold", "specificity", "sensitivity", "accuracy",
                                   "precision", "recall"), best.method=c("youden"), transpose = TRUE)
  
  #Logistic (Train: AUC,accuracy)
  PredTestY2<-predict(model1,newdata=finaldtTest,type="response")
  a4<-roc(unlist(outcomeTest),as.numeric(PredTestY2),quiet=T)
  a4.out<-coords(a4, "best", ret=c("threshold", "specificity", "sensitivity", "accuracy",
                                   "precision", "recall"), best.method=c("youden"), transpose = TRUE)
  
  #Combine the output
  tab1<-data.frame(round(a.gbmcv$auc,4),round(a.gbmcv.out["accuracy"],4),round(a.gbmcv.out["precision"],4),round(a.gbmcv.out["recall"],4),round(a.gbmcv.out.f1,4),
                   round(a.gbm1$auc,4),round(a.gbm1.out["accuracy"],4),round(a.gbm1.out["precision"],4),round(a.gbm1.out["recall"],4),round(a.gbm1.out.f1,4),
                   round(a.gbm2$auc,4),round(a.gbm2.out["accuracy"],4),round(a.gbm2.out["precision"],4),round(a.gbm2.out["recall"],4),round(a.gbm2.out.f1,4))
  colnames(tab1)<-c("GBM:CVAUC","GBM:CVACC","GBM:CVPrecision","GBM:CVRecall","GBM:CVF1",
                    "GBM:trainAUC","GBM:trainACC","GBM:trainPrecision","GBM:trainRecall","GBM:trainF1",
                    "GBM:testAUC","GBM:testACC","GBM:testPrecision","GBM:testRecall","GBM:testF1")
  
  tab2_RFLR<-data.frame(round(a1$auc,2),round(a1.out["accuracy"],2),round(a2$auc,2),round(a2.out["accuracy"],2),
                        round(a3$auc,2),round(a3.out["accuracy"],2),round(a4$auc,2),round(a4.out["accuracy"],2))
  colnames(tab2_RFLR)<-c("RF:trainAUC","RF:trainACC","RF:testAUC","RF:testACC",
                         "LR:trainAUC","LR:trainACC","LR:testAUC","LR:testACC")
  return(list(tab1,tab2_RFLR,gbm1,best.iter))
}

