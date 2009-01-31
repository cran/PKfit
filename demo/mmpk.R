library(PKfit)
library(stats4)
options(warn=-1)
PKindex<-data.frame(Subject=c(1),time=c(0.25,0.5,0.75,1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10,12,16,20,24),
conc=c(1.03,2.05,3.06,4.06,6.07,8.07,10.08,12.08,11.86,11.64,11.42,11.21,10.77,10.33,9.89,9.46,9.02,8.15,6.41,4.69,3.00))
Dose<-300
Tinf<-3

defun<-function(time, y, parms) { 
        if(time<=Tinf)  
           dCpdt <- (Dose/Tinf)/parms["Vd"] -
           (parms["Vm"]/parms["Vd"])*y[1]/(parms["Km"]/parms["Vd"]+y[1]) 
         else
           dCpdt <- -(parms["Vm"]/parms["Vd"])*y[1]/(parms["Km"]/parms["Vd"]+y[1])
         list(dCpdt)
}
modfun<-function(time,Vm,Km,Vd) { 
        out <- lsoda(0,c(0,time),defun,parms=c(Vm=Vm,Km=Km,Vd=Vd),rtol=1e-5,atol=1e-5)
        out[-1,2] 
} 

objfun <- function(par) {
        out <- modfun(PKindex$time, par[1], par[2], par[3])
        gift <- which( PKindex$conc != 0 )
        sum((PKindex$conc[gift]-out[gift])^2)
}        

gen <- genoud(objfun,nvars=3,max=FALSE,pop.size=30,max.generations=20,
              wait.generations=10,
              starting.value=c(40,8,12),BFGS=FALSE,
              print.level=0,boundary.enforcement=0,
              Domains=matrix(c(0,0,0,100,100,100),3,2),
              MemoryMatrix=TRUE)
namegen<-c("Vm","Km","Vd")
outgen<-c(gen$par[1],gen$par[2],gen$par[3])
F<-objfun(gen$par)      
opt<-optim(c(gen$par[1],gen$par[2],gen$par[3]),objfun,method="Nelder-Mead")  
nameopt<-c("Vm","Km","Vd")
outopt<-c(opt$par[1],opt$par[2],opt$par[3])     
fm<-nls(conc~modfun(time,Vm,Km,Vd),data=PKindex,
        start=list(Vm=opt$par[1],Km=opt$par[2],Vd=opt$par[3]),trace=TRUE,
        nls.control(tol=1))

x<-PKindex$time
y<-PKindex$conc
cal<-predict(fm,list(time=x))
wei<-ifelse(y==0.0, 0, y-cal)
add<-function(time,conc){
     auc<-0 ; aumc<-0
     for(i in 2:length(time)) {
     auc[i]<-1/2*(time[i]-time[i-1])*(conc[i]+conc[i-1])
     auc[i]<-auc[i]+auc[i-1]
     aumc[i]<-1/2*(time[i]-time[i-1])*(conc[i]*time[i]+conc[i-1]*time[i-1])
     aumc[i]<-aumc[i]+aumc[i-1]
     }
     return(list(auc=auc,aumc=aumc))
  }
add1 <- add(x,y)
AUC<-c(NA,add1$auc[-1])
AUMC<-c(NA,add1$aumc[-1])
output<-data.frame(x,y,cal,wei,AUC,AUMC)
colnames(output)<-list("Time","Observed","Calculated","Wt. Residuals","AUC","AUMC")

windows(record=TRUE)

par(mfrow=c(2,2))

plot(y~x,data=PKindex,type='p',main="time - drug conc. curve", 
     xlab="Time", ylab="concentration",pch=15,col="black",bty="l",
     font.lab=2,cex.lab=1,cex.axis=1,cex.main=1) 
lines(x,predict(fm,list(time=x)),type="l",lty=1,
      col="firebrick3",lwd="2")
mtext("Linear",side=3,cex=0.88)
    
plot(x,y,log="y",type='p',main="time - drug conc. curve",
     xlab="Time", ylab="concentration",pch=15,col="black",bty="l",
     font.lab=2,cex.lab=1,cex.axis=1,cex.main=1) 
lines(x,predict(fm,list(time=x)),type="l",lty=1,
      col="firebrick3",lwd="2")
mtext("Semi-log",side=3,cex=0.88)
     
plot(x,wei,pch=15,col="blue",bty="l",xlab="Time",
     ylab="Weighted Residual",main="Residual Plots",cex.lab=1,
     cex.axis=1,cex.main=1,font.lab=2)
abline(h=0,lwd=2,col="black",lty=2)
    
plot(cal,wei,pch=15,col="blue",bty="l",xlab="Calc Cp(i)",
     ylab="Weighted Residual",main="Residual Plots",cex.lab=1,
     cex.axis=1,cex.main=1,font.lab=2)
abline(h=0,lwd=2,col="black",lty=2)    

output

AIC(fm)
    
logLik(fm) 

BIC(fm)
    
print(summary(fm))    

       