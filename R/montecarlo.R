montecarlo<-function(C1.lsoda,time1,i,re)
{
  #options(warn=-1)
  conc<-rowMeans(as.data.frame(lapply(C1.lsoda,"[","concentration")))
  C1.lsoda<-do.call("rbind",C1.lsoda)
  rownames(C1.lsoda)<-seq(nrow(C1.lsoda))
  
  conc<-ifelse(conc<=0,0,conc)
  
  PKindex<-data.frame(i,time1,conc)
##  dev.new()
  par(mfrow=c(2,2))
  x<-C1.lsoda$time
  #y<-C1.lsoda$concentration
  
  y<-ifelse(C1.lsoda$concentration<=0,0,C1.lsoda$concentration)
  
  main<-paste(c("Subject:-", i),collapse=" ")
  plot(C1.lsoda,type="p",main=main,xlab="Time after dosing",ylab="Drug plasma conc.")
  lines(PKindex$time1,PKindex$conc,type="l",lty=1,col="firebrick3",lwd="2")
  mtext("Linear",side=3,cex=0.8)
  
  plot(x,y,log="y",type='p',main=main,xlab="Time after dosinge",ylab="Drug plasma conc.")
  lines(PKindex$time1,PKindex$conc,log="y",type="l",lty=1,col="firebrick3",lwd="2")
  mtext("Semi-log",side=3,cex=0.8) 
  
  len<-length(time1)
  sub<-len*re
  a<-seq(0,sub,by=len)
  AA<-a[2:(length(a)-1)]
  colnames(PKindex)<-list("Subject","time","conc")
  cat("\n");show(PKindex)
  
  for(i in rev(AA))
      C1.lsoda<-C1.lsoda[append(1:(sub+(length(AA)-1)),NA,after=i),]   
  
  plot(C1.lsoda,type="l",main=main,xlab="Time after dosing",ylab="Drug plasma conc.")  
  lines(PKindex$time,PKindex$conc,type="l",lty=1,col="firebrick3",lwd="2")
  mtext("Linear",side=3,cex=0.8)
  
  x<-C1.lsoda$time
  y<-C1.lsoda$conc
  
  plot(x,y,log="y",type='l',main=main,xlab="Time after dosing",ylab="Drug plasma conc.")
  lines(PKindex$time,PKindex$conc,log="y",type="l",lty=1,col="firebrick3",lwd="2")
  mtext("Semi-log",side=3,cex=0.8) 
  
  return(PKindex)
}