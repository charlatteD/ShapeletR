#Shape 1 a-b-a
shp1<-function(lshp,minshp,maxshp){
  midpoint<-floor(lshp/2)
  
  slope1<-(maxshp-minshp)/(midpoint-0)
  intercept1<-minshp-0*slope1
  slope2<-(minshp-maxshp)/(lshp-midpoint)
  intercept2<-minshp-lshp*slope2
  
  x<-c(0:lshp)
  y1st<-x[0:midpoint]*slope1+intercept1
  y2nd<-x[(midpoint+1):(lshp+1)]*slope2+intercept2
  
  y<-c(y1st,y2nd)
  
  outdata<-data.frame(x,y)
  return(outdata)
}

#Shape 2 b-a-b
shp2<-function(lshp,minshp,maxshp){
  
  midpoint<-floor(lshp/2)
  
  slope1<-(minshp-maxshp)/(lshp-midpoint)
  intercept1<-maxshp-0*slope1
  slope2<-(maxshp-minshp)/(midpoint-0)
  intercept2<-maxshp-lshp*slope2
  
  x<-c(0:lshp)
  y1st<-x[0:midpoint]*slope1+intercept1
  y2nd<-x[(midpoint+1):(lshp+1)]*slope2+intercept2
  
  y<-c(y1st,y2nd)
  
  outdata<-data.frame(x,y)
  return(outdata)
}

#Shape 3 a-b-b-b 
shp3<-function(lshp,minshp,maxshp){
  bpoint1<-floor(lshp/4)
  
  slope1<-(maxshp-minshp)/(bpoint1-0)
  intercept1<-minshp-0*slope1
  
  x<-c(0:lshp)
  y1st<-x[0:(bpoint1+1)]*slope1+intercept1
  y2nd<-rep(maxshp,length(x[(bpoint1+1):(lshp)]))
  
  y<-c(y1st,y2nd)
  
  outdata<-data.frame(x,y)
  return(outdata)
}

#Shape 4 b-b-a
shp4<-function(lshp,minshp,maxshp){
  bpoint1<-floor(3*lshp/4)
  
  slope1<-(minshp-maxshp)/(lshp-bpoint1)
  intercept1<-minshp-lshp*slope1
  
  x<-c(0:lshp)
  y1st<-rep(maxshp,length(x[0:(bpoint1)]))
  y2nd<-x[(bpoint1+1):(lshp+1)]*slope1+intercept1
  
  y<-c(y1st,y2nd)
  
  outdata<-data.frame(x,y)
  return(outdata)
}

#Shape 5 a-a-a-a
shp5<-function(lshp,minshp,maxshp){
  
  x<-c(0:lshp)
  y<-rep(maxshp,length(x))
  
  outdata<-data.frame(x,y)
  return(outdata)
}

#Shape 6 a-b-b-a
shp6<-function(lshp,minshp,maxshp){
  
  bpoint1<-floor(lshp/4)
  bpoint2<-floor(lshp*3/4)
  
  slope1<-(maxshp-minshp)/(bpoint1-0)
  intercept1<-minshp-0*slope1
  slope2<-(minshp-maxshp)/(lshp-bpoint2)
  intercept2<-minshp-lshp*slope2
  
  x<-c(0:lshp)
  y1st<-x[0:bpoint1]*slope1+intercept1
  y2nd<-rep(maxshp,length(x[bpoint1:bpoint2]))
  y3rd<-x[(bpoint2+2):(lshp+1)]*slope2+intercept2
  
  y<-c(y1st,y2nd,y3rd)
  
  outdata<-data.frame(x,y)
  return(outdata)
}

#Shape 7 b-a-a-b
shp7<-function(lshp,minshp,maxshp){
  
  bpoint1<-floor(lshp/4)
  bpoint2<-floor(lshp*3/4)
  
  slope1<-(minshp-maxshp)/(bpoint1-0)
  intercept1<-minshp-bpoint1*slope1
  slope2<-(maxshp-minshp)/(lshp-bpoint2)
  intercept2<-minshp-bpoint2*slope2
  
  x<-c(0:lshp)
  y1st<-x[0:bpoint1]*slope1+intercept1
  y2nd<-rep(minshp,length(x[bpoint1:bpoint2]))
  y3rd<-x[(bpoint2+2):(lshp+1)]*slope2+intercept2
  
  y<-c(y1st,y2nd,y3rd)
  
  outdata<-data.frame(x,y)
  return(outdata)
}

#Shape 8 b-c-a-b
shp8<-function(lshp,minshp,maxshp){
  bpoint1<-floor(lshp/4)
  bpoint2<-floor(3*lshp/4)
  
  slope1<-((maxshp+minshp)/2-minshp)/(bpoint1-0)
  intercept1<-(maxshp+minshp)/2-0*slope1
  slope2<-(minshp-maxshp)/(bpoint2-(bpoint1+1))
  intercept2<-maxshp-(bpoint1+1)*slope2
  slope3<-((maxshp+minshp)/2-minshp)/(lshp-(bpoint2+1))
  intercept3<-(maxshp+minshp)/2-lshp*slope3
  
  x<-c(0:lshp)
  y1st<-x[0:bpoint1+1]*slope1+intercept1
  y2nd<-x[(bpoint1+2):(bpoint2+1)]*slope2+intercept2
  y3rd<-x[(bpoint2+2):(lshp+1)]*slope3+intercept3
  
  y<-c(y1st,y2nd,y3rd)
  
  outdata<-data.frame(x,y)
  return(outdata)
}

#Shape 9 b-a-c-b
shp9<-function(lshp,minshp,maxshp){
  bpoint1<-floor(lshp/4)
  bpoint2<-floor(3*lshp/4)
  
  slope1<-(minshp-(maxshp+minshp)/2)/(bpoint1-0)
  intercept1<-(maxshp+minshp)/2-0*slope1
  slope2<-(maxshp-minshp)/(bpoint2-(bpoint1+1))
  intercept2<-maxshp-(bpoint2)*slope2
  slope3<-((maxshp+minshp)/2-maxshp)/(lshp-(bpoint2+1))
  intercept3<-(maxshp+minshp)/2-lshp*slope3
  
  x<-c(0:lshp)
  y1st<-x[0:bpoint1+1]*slope1+intercept1
  y2nd<-x[(bpoint1+2):(bpoint2+1)]*slope2+intercept2
  y3rd<-x[(bpoint2+2):(lshp+1)]*slope3+intercept3
  
  y<-c(y1st,y2nd,y3rd)
  
  outdata<-data.frame(x,y)
  return(outdata)
}

#Shape 10 b-a-a
shp10<-function(lshp,minshp,maxshp){
  bpoint1<-floor(lshp/4)
  bpoint2<-floor(3*lshp/4)
  
  slope2<-(minshp-maxshp)/(bpoint2-(bpoint1+1))
  intercept2<-maxshp-(bpoint1+1)*slope2
  
  x<-c(0:lshp)
  y1st<-rep(maxshp,length(x[0:bpoint1+1]))
  y2nd<-x[(bpoint1+2):(bpoint2+1)]*slope2+intercept2
  y3rd<-rep(minshp,length(x[0:bpoint1+1]))
  
  y<-c(y1st,y2nd,y3rd)
  
  outdata<-data.frame(x,y)
  return(outdata)
}

#Shape 11 a-a-b
shp11<-function(lshp,minshp,maxshp){
  bpoint1<-floor(lshp/4)
  bpoint2<-floor(3*lshp/4)
  
  slope2<-(maxshp-minshp)/(bpoint2-(bpoint1+1))
  intercept2<-maxshp-(bpoint2)*slope2
  
  x<-c(0:lshp)
  y1st<-rep(minshp,length(x[0:bpoint1+1]))
  y2nd<-x[(bpoint1+2):(bpoint2+1)]*slope2+intercept2
  y3rd<-rep(maxshp,length(x[0:bpoint1+1]))
  
  y<-c(y1st,y2nd,y3rd)
  
  outdata<-data.frame(x,y)
  return(outdata)
}