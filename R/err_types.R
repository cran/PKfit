err_types<-function(pick,par,factor){

sim_value<-NULL

if(pick==1){
  sim_value<-par+rnorm(1,mean=0,sd=factor)
  while(sim_value<=0){
  sim_value<-par+rnorm(1,mean=0,sd=factor)}
}
if(pick==2){
  sim_value<-par+runif(1,min=-factor,max=factor)
  while(sim_value<=0){
  sim_value<-par+runif(1,min=-factor,max=factor)}
}
if(pick==3){
  sim_value<-par*rnorm(1,mean=0,sd=factor)+par
  while(sim_value<=0){
  sim_value<-par*rnorm(1,mean=0,sd=factor)+par}
}
if(pick==4){
  sim_value<-par*runif(1,min=-factor,max=factor)+par
  while(sim_value<=0){
  sim_value<-par*runif(1,min=-factor,max=factor)+par}
}
return(sim_value)
}
