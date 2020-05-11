#-------------------------------------------------------------------------------------------------------------#
rm(list = ls())
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
#Import a priori shapes
setwd("PATH")
source("ShapeletFuncs.R")
#-------------------------------------------------------------------------------------------------------------#
####################################################################################################
#(Simulated data) Data analysis 
####################################################################################################
#-------------------------------------------------------------------------------------------------------------#
#Data step (Simulated data)
#-------------------------------------------------------------------------------------------------------------#
#"True" shape
ExampleShape<-c(0,10,10,9.5,9.5,9,9,7,6,5.5)
#Training data
set.seed(123)
datacommon_sim<-as.data.frame(t(replicate(100, rnorm(100))))
sim1<-getDataWithShapeV2(ExampleShape,datacommon_sim,0.7,50,51,0.95,0.05,F)
TrainX<-sim1[[1]]
TrainY<-sim1[[2]]
#Test data
set.seed(222)
datacommon_sim<-as.data.frame(t(replicate(100, rnorm(100))))
sim2<-getDataWithShapeV2(ExampleShape,datacommon_sim,0.7,50,51,0.95,0.05,F)
TestX<-sim2[[1]]
TestY<-sim2[[2]]

#Kmean clustering 
set.seed(123)
start_time <- Sys.time()
Kmean_output<-Shapelet_Kmeans(TrainX,TrainY,TestX,TestY,c(10,15),10,1,1)
kmean_time<-Sys.time()-start_time

#Unique shapes (stratified)
set.seed(123)
start_time <- Sys.time()
UniqueFeatures_output<-Shapelet_UniqueFeatures(TrainX,TrainY,TestX,TestY,c(10,15),1,1)
Unique_time<-Sys.time()-start_time

#Plot the three shapelets
PlotShapePartialDependencetemp1<-function(gbm.model,best.iter,topvalue,shplists,inputDistTrain,input,ylimc,outcome,shploclist,ylab1){
  varlist<-summary(gbm.model,n.trees = best.iter, plotit = FALSE)[,1]
  varlist2<-as.numeric(substring(varlist, first=2, last = 1000000L))
  topgbm<-varlist2[1:topvalue]
  
  #Transform the input for plotting
  minmaxInput<-minmax(input)
  
  #Plot the shapelet
  plot(as.numeric(shplists[varlist2[1],]),type="l",xlab="", ylab="",ylim=c(0,1),main="Shapelet")
  
  #Plot the best matching location 
  inputDist_minloc<-which.min(inputDistTrain[,varlist2[1]])
  plot(1:length(as.numeric(input[inputDist_minloc,])), as.numeric(minmaxInput[inputDist_minloc,]),
       type="l",ylab="",lwd=2,ylim=c(0,1),xlab="",col="grey",main="Shapelet on its\n best match location")
  #,main=paste0("Outcome = ",outcome[varlist2[i],]))
  lines(c(shploclist[inputDist_minloc,varlist2[1]]):(shploclist[inputDist_minloc,varlist2[1]]+length(shplists[varlist2[1],])-1),shplists[varlist2[1],], 
        col="black",lwd=3)
  
  #Plot the minimum distance vs. predicted value
  p1<-plot(gbm.model, i.var = varlist2[1], n.trees = best.iter,return.grid=T)
  plot(p1[,1],p1[,2], type="l",xlab="Inverse minimum distance", ylab="log(OR)",ylim=ylimc,main="Partial dependence plot")
  
}
PlotShapePartialDependencetemp2<-function(gbm.model,best.iter,topvalue,shplists,inputDistTrain,input,ylimc,outcome,shploclist,ylab1){
  varlist<-summary(gbm.model,n.trees = best.iter, plotit = FALSE)[,1]
  varlist2<-as.numeric(substring(varlist, first=2, last = 1000000L))
  topgbm<-varlist2[1:topvalue]
  
  #Transform the input for plotting
  minmaxInput<-minmax(input)
  #Plot the shapelet
  plot(as.numeric(shplists[varlist2[1],]),type="l",xlab="", ylab="",ylim=c(0,1),main="")
  
  #Plot the best matching location 
  inputDist_minloc<-which.min(inputDistTrain[,varlist2[1]])
  plot(1:length(as.numeric(input[inputDist_minloc,])), as.numeric(minmaxInput[inputDist_minloc,]),
       type="l",ylab="",lwd=2,ylim=c(0,1),xlab="",col="grey",main="")
  #,main=paste0("Outcome = ",outcome[varlist2[i],]))
  lines(c(shploclist[inputDist_minloc,varlist2[1]]):(shploclist[inputDist_minloc,varlist2[1]]+length(shplists[varlist2[1],])-1),shplists[varlist2[1],], 
        col="black",lwd=3)
  
  #Plot the minimum distance vs. predicted value
  p1<-plot(gbm.model, i.var = varlist2[1], n.trees = best.iter,return.grid=T)
  plot(p1[,1],p1[,2], type="l",xlab="", ylab="",ylim=ylimc,main="")
}
pdf("Sim_all_partial.pdf",width=8,height=8)
par(mfrow=c(2,3),mar=c(4,4,5,1)+ 0.1, mgp=c(3,.7,0), tck=-.01)
PlotShapePartialDependencetemp1(Kmean_output$gbm.model1,Kmean_output$gbm.best.iter1,1,Kmean_output$shplists1,Kmean_output$inputDistTrain,TrainX,TrainY,Kmean_output$inputLocTrain,ylimc=c(-0.5,2.5))
PlotShapePartialDependencetemp2(UniqueFeatures_output$gbm.model1,UniqueFeatures_output$gbm.best.iter1,1,UniqueFeatures_output$shplists1,UniqueFeatures_output$inputDistTrain,TrainX,TrainY,UniqueFeatures_output$inputLocTrain,ylimc=c(-0.5,2.5))
mtext("K-means initialization", side = 3, line = -5, outer = TRUE)
mtext("Extreme statistics initialization", side = 3, line = -32, outer = TRUE)
dev.off()
