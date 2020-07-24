#-------------------------------------------------------------------------------------------------------------#
rm(list = ls())
#-------------------------------------------------------------------------------------------------------------#
library(rdist)
library(TSPred)
library(pROC)
library(dplyr)
library(plyr)
library(zoo)
library(gbm)
library(pdp)
#-------------------------------------------------------------------------------------------------------------#
#Import a priori shapes
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
Kmean_output<-Shapelet_Kmeans(TrainX,TrainY,TestX,TestY,c(10,15),10)

#Unique shapes (stratified)
set.seed(123)
UniqueFeatures_output<-Shapelet_UniqueFeatures(TrainX,TrainY,TestX,TestY,c(10,15))

#Plot the three shapelets
PlotShapePartialDependencetemp1<-function(gbm.model,
                                          shplists,
                                          inputDistTrain,
                                          input,
                                          ylimc,
                                          outcome,
                                          shploclist){
  
  inputDistTrain = data.frame(inputDistTrain)
  
  #Get the most important predictor
  varlist2 <- xgb.importance(model=gbm.model)$Feature
  first <- varlist2[1]
  firstval<- as.numeric(substring(first,2,3))
  
  #Transform the input for plotting
  minmaxInput<-minmax(input)
  
  #Plot the shapelet
  plot(as.numeric(shplists[firstval,]),type="l",xlab="", ylab="",ylim=c(0,1),main="Shapelet")
  
  #Plot the best matching location 
  inputDist_minloc<-which.min(inputDistTrain[,firstval])
  plot(1:length(as.numeric(input[inputDist_minloc,])), as.numeric(minmaxInput[inputDist_minloc,]),
       type="l",ylab="",lwd=2,ylim=c(0,1),xlab="",col="grey",main="Shapelet on its\n best match location")
  lines(c(shploclist[inputDist_minloc,firstval]):(shploclist[inputDist_minloc,firstval]+length(shplists[firstval,])-1),shplists[firstval,], 
        col="black",lwd=3)
  
  #Plot the minimum distance vs. predicted value
  p1<-pdp::partial(gbm.model, train=1/(inputDistTrain+0.1), pred.var=first)
  plot(p1[,1],p1[,2], type="l",xlab="Inverse minimum distance", ylab="log(OR)",main="Partial dependence plot")
  
}
pdf("Sim_all_partial.pdf",width=8,height=8)
par(mfrow=c(2,3),mar=c(4,4,5,1)+ 0.1, mgp=c(3,.7,0), tck=-.01)
PlotShapePartialDependencetemp1(Kmean_output$gbm1,
                                Kmean_output$shplists1,
                                Kmean_output$inputDistTrain+0.1,
                                TrainX,
                                ylimc=c(-2.0,-1.7),
                                TrainY,
                                Kmean_output$inputLocTrain)

PlotShapePartialDependencetemp1(UniqueFeatures_output$gbm1,
                                UniqueFeatures_output$shplists1,
                                UniqueFeatures_output$inputDistTrain,
                                TrainX,
                                ylimc=c(-0.5,2.5),
                                TrainY,
                                UniqueFeatures_output$inputLocTrain)
mtext("K-means initialization", side = 3, line = -5, outer = TRUE)
mtext("Extreme statistics initialization", side = 3, line = -32, outer = TRUE)
dev.off()
