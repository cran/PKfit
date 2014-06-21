iv.bolus.demo<-function(){
cat("\n\n")
options(warn=-1)
modfun<-NULL
conc<-NULL

PKindex<-data.frame(Subject=c(1),time=c(1,2,3,4,6,10,12),
                    conc=c(14.94,13.73,10.55,8.16,5.21,3.19,2.62))
Dose<-500
defun<- function(time, y, parms) { 
      dCpdt <- -parms["kel"] * y[1] 
      list(dCpdt) 
} 
    
modfun <<- function(time,kel, Vd) {  
      out <- lsoda(Dose/Vd,c(0,time),defun,parms=c(kel=kel,Vd=Vd),
                   rtol=1e-6,atol=1e-6) 
      out[-1,2] 
}

objfun <- function(par) {
        out <- modfun(PKindex$time, par[1], par[2])
        gift <- which( PKindex$conc != 0 )
        ### sum((PKindex$conc[gift]-out[gift])^2)
        sum(((PKindex$conc[gift]-out[gift])/PKindex$conc[gift])^2)
}        
cat("- weighting scheme: 1/Cp^2\n")
cat("- model selection: a one-compartment, iv bolus pk model with\n  1st-ordered elim.\n\n")
### gen<-genoud(objfun,nvars=2,max=FALSE,pop.size=30,
###             max.generations=20,wait.generations=10,
###             starting.values=c(0.13,20),
###             BFGS=FALSE,print.level=0,boundary.enforcement=2,
###             Domains=matrix(c(0.01,0.01,100,100),2,2),
###             MemoryMatrix=TRUE)  
### namegen<-c("kel","Vd")
### outgen<-c(gen$par[1],gen$par[2]) 
     
opt<-optim(c(0.13,20),objfun,method="Nelder-Mead",control=list(maxit=5000))  
nameopt<-c("kel","Vd")
outopt<-c(opt$par[1],opt$par[2])
cat("<< PK parameters obtained from Nelder-Mead Simplex algorithm >>\n")
print(data.frame(Parameter=nameopt,Value=outopt))
cat("\n\n")

  if(opt$par[1]<0) {opt$par[1]<-0.01}
  if(opt$par[2]<0) {opt$par[2]<-0.01}

fm<-nlsLM(conc ~ modfun(time, kel, Vd),data=PKindex,start=list(kel=opt$par[1],Vd=opt$par[2]),
         control=nls.lm.control(maxiter=500),weights=(1/conc^2)) ### lower of Vd should not be zero due to Dose/Vd. --YJ
        
coef<-data.frame(coef(fm)["kel"])
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
add1<-add(x,y)
AUC<-c(NA,add1$auc[-1])
AUMC<-c(NA,add1$aumc[-1])
output<-data.frame(x,y,cal,wei,AUC,AUMC)
colnames(output)<-list("Time","Observed","Calculated","Wt. Residuals","AUC","AUMC")

auc.infinity<-y[length(y)]/coef[1,1]
auc<-AUC[length(y)]+auc.infinity

aumc.infinity<-(x[length(x)]*y[length(y)])/coef[1,1]+y[length(y)]/((coef[1,1])^2)
aumc<-AUMC[length(y)]+aumc.infinity

### windows(record=TRUE)
dev.new()

par(mfrow=c(2,2), ask = FALSE)

plot(y~x,data=PKindex,type='p',main="Drug Plasma Conc. vs. Time Curve", 
     xlab="Time", ylab="Drug Plasma Conc.",pch=15,col="black",bty="l",
     font.lab=2,cex.lab=1,cex.axis=1,cex.main=1) 
lines(x,predict(fm,list(time=x)),type="l",lty=1,
      col="firebrick3",lwd="2")
mtext("Linear",side=3,cex=0.88)
    
plot(x,y,log="y",type='p',main="Drug Plasma Conc. vs. Time Curve",
     xlab="Time", ylab="Drug Plasma Conc.",pch=15,col="black",bty="l",
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
cat("<< Residual sum-of-squares and PK parameter values with nlsLM >>\n\n")
show(output)
aicllsbc(fm)
cat("\n\n")
}
       