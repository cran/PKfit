mcsim.demo<-function(){
cat("\n\n")
options(warn=-1)
Subject<-1

PKtime<-c(0.25,0.5,1,2,4,6,8,10,12,16,20,24)

Dose<-250

par1<-1.1
par2<-0.14
par3<-22.9   

re<-1600

factor1<-0.02
factor2<-0.02
factor3<-0.1

defun <- function(time, y, parms) { 
     dy1dt <- -parms["ka"] * y[1]
     dy2dt <-  parms["ka"] * y[1]/parms["Vd"] - parms["kel"] * y[2]
     list(c(dy1dt,dy2dt)) 
} 

C1.lsoda<-list()
for (j in 1:re){
  ka<-par1+runif(1,min=-factor1,max=factor1)
  while(ka<=0){
     ka<-par1+runif(1,min=-factor1,max=factor1)}
  kel<-par2+runif(1,min=-factor2,max=factor2)
  while(kel<=0){
     kel<-par2+runif(1,min=-factor2,max=factor2)}   
  Vd<-par3+runif(1,min=-factor3,max=factor3)
  while(Vd<=0){
     Vd<-par3+runif(1,min=-factor3,max=factor3)}
  time<-PKtime
  parms<-c(ka=ka,kel=kel,Vd=Vd)  
  XX<-data.frame(lsoda(c(Dose,0),c(0,time), defun, parms))  
  C1.lsoda[[j]]<-data.frame(XX[2:(length(time)+1),1],XX[2:(length(time)+1),3])
  colnames(C1.lsoda[[j]])<-list("time","concentration") 
} 
     
for( i in 1:Subject) 
conc<-rowMeans(as.data.frame(lapply(C1.lsoda,"[","concentration")))
C1.lsoda<-do.call("rbind",C1.lsoda)
rownames(C1.lsoda)<-seq(nrow(C1.lsoda))
conc<-ifelse(conc<=0, 0, conc)
PKindex<-data.frame(i,time,conc)
windows(record=TRUE)
par(mfrow=c(2,2), ask = FALSE)
x<-C1.lsoda$time
y<-ifelse(C1.lsoda$concentration<=0, 0, C1.lsoda$concentration)
  

plot(C1.lsoda,type="p",main="Monte-Carlo Simulation",xlab="Time (hr)",ylab="Drug Plasma Conc.")
lines(PKindex$time,PKindex$conc,type="l",lty=1,col="firebrick3",lwd="2")
mtext("Linear",side=3,cex=0.8)
  
plot(x,y,log="y",type='p',main="Monte-Carlo Simulation",xlab="Time (hr)",ylab="Drug Plasma Conc.")
lines(PKindex$time,PKindex$conc,log="y",type="l",lty=1,col="firebrick3",lwd="2")
mtext("Semi-log",side=3,cex=0.8) 
  
len<-length(time)
sub<-len*re
a<-seq(0,sub,by=len)
AA<-a[2:(length(a)-1)]
colnames(PKindex)<-list("Subject","Time","Conc")

for(i in rev(AA))
  C1.lsoda<-C1.lsoda[append(1:(sub+(length(AA)-1)),NA,after=i),]   
 
plot(C1.lsoda,type="l",main="Monte-Carlo Simulation",xlab="Time (hr)",ylab="Drug Plasma Conc.")  
lines(PKindex$time,PKindex$conc,type="l",lty=1,col="firebrick3",lwd="2")
mtext("Linear plot",side=3,cex=0.8)
  
x<-C1.lsoda$time
y<-C1.lsoda$concentration
  
plot(x,y,log="y",type='l',main="Monte-Carlo Simulation",xlab="Time (hr)",ylab="Drug Plasma Conc.")
lines(PKindex$time,PKindex$conc,log="y",type="l",lty=1,col="firebrick3",lwd="2")
mtext("Semi-log plot",side=3,cex=0.8) 
cat("\n\n")
show(PKindex)
cat("\n\n")
}