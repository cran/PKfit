mmpk.demo<-function(){
cat("\n\n")
modfun<-NULL
conc<-NULL

options(warn=-1)
PKindex<-data.frame(Subject=c(1),time=c(0,6,12,24,36,48,72,96,144),
conc=c(7.88,7.10,6.79,5.27,4.52,3.43,1.97,1.01,0.23))
Dose<-300
par<-data.frame(Parameter=c("Vm","Km","Vd"),Initial=c(17,5,10))
description_version()
cat("- This is an iv-bolus drug exhibiting with Mechaelis-Menten elimination -\n")
cat(" PK and drug plasma samples were collected and assayed after dosing.\n")
cat(" Data was from: Fitting Models to Biological Data using Linear and \n")
cat(" Nonlinear Regression, A pratical guide to curve fitting. GraphPad \n")
cat(" PRISM, version 4.0, Harvey Motulsky & Arthur Christopoulos. 2ed. ed.,\n")
cat(" 2003, pp. 74-77. (Note: we do correct Vm here with Vd.)\n\n")
cat(" 1. Dose =",Dose,"\n")
cat(" 2. Input data:-\n")
show(PKindex);cat("\n\n")
cat(" 3. Initial values for parameters:-\n")
show(par);cat("\n")

defun<-function(time, y, parms) { 
         dCpdt <- -parms["Vm"]*y[1]/(parms["Km"]+y[1])
         list(dCpdt)
}
modfun<<-function(time,Vm,Km,Vd) { 
        out <- lsoda(Dose/Vd,c(0,time),defun,parms=c(Vm=Vm,Km=Km,Vd=Vd),rtol=1e-6,atol=1e-6)
        out[-1,2] 
} 

objfun <- function(par) {
        out <- modfun(PKindex$time, par[1], par[2], par[3])
        gift <- which( PKindex$conc != 0 )
        ### sum((PKindex$conc[gift]-out[gift])^2)
        sum(((PKindex$conc[gift]-out[gift])/PKindex$conc[gift])^2)
}        
cat("- weighting scheme: equal weight\n\n")
cat("- model selection: a one-compartment, iv bolus pk model with\n  M-M elim.\n\n")
### gen <- genoud(objfun,nvars=3,max=FALSE,pop.size=30,max.generations=20,
###               wait.generations=10,
###               starting.values=c(40,8,12),BFGS=FALSE,
###               print.level=0,boundary.enforcement=0,
###               Domains=matrix(c(0,0,0,100,100,100),3,2),
###               MemoryMatrix=TRUE)
### namegen<-c("Vm","Km","Vd")
### outgen<-c(gen$par[1],gen$par[2],gen$par[3])
### cat("<< PK parameters obtained from genetic algorithm >>\n\n")
### print(data.frame(Parameter=namegen,Value=outgen))  
### F<-objfun(gen$par)

cat(" running optimx() right now...\n\n")
opt<-optimx(c(par[1,2],par[2,2],par[3,2]),objfun,method="Nelder-Mead",control=list(maxit=5000))
nameopt<-c("Vm","Km","Vd")
outopt<-c(opt$p1,opt$p2,opt$p3) 
cat("<< PK parameters obtained from Nelder-Mead Simplex algorithm >>\n\n")
print(data.frame(Parameter=nameopt,Value=outopt))
cat("\n\n")
if(opt$p1<0) {opt$p1<-0.001}
if(opt$p2<0) {opt$p2<-0.001}
if(opt$p3<0) {opt$p3<-0.001}

fm<-nlsLM(conc ~ modfun(time,Vm,Km,Vd), data=PKindex,start=list(Vm=opt$p1,Km=opt$p2,Vd=opt$p3),
         control=nls.lm.control(maxiter=500,maxfev=5000,factor=100),weights=(1/conc^0))  ### set 'lower=c(...)' may cause crashed.  --YJ

### tried to use self-starter function for nls()
### fm<-nls(conc~SSmicmen(time,Vm,Km),data=PKindex)  ### no Vd?  so SSmicmen() is not useful for this.
### print(summary(fm))

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

### windows(record=TRUE)
dev.new()

par(mfrow=c(2,2), ask = FALSE)

plot(y~x,data=PKindex,type='p',main="Phenytoin Plasma Conc. vs. Time Curve", 
     xlab="Time after dosing (hr)", ylab="Phenytoin Plasma Conc. (mg/L)",pch=15,col="black",bty="l",
     font.lab=2,cex.lab=1,cex.axis=1,cex.main=1) 
lines(x,predict(fm,list(time=x)),type="l",lty=1,
     col="firebrick3",lwd="2")
mtext("Linear plot",side=3,cex=0.88)
    
plot(x,y,log="y",type='p',main="Phenytoin Plasma Conc. vs. Time Curve",
     xlab="Time after dosing (hr)", ylab="Phenytoin Plasma Conc. (mg/L)",pch=15,col="black",bty="l",
     font.lab=2,cex.lab=1,cex.axis=1,cex.main=1) 
lines(x,predict(fm,list(time=x)),type="l",lty=1,
     col="firebrick3",lwd="2")
mtext("Semi-log plot",side=3,cex=0.88)
     
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
       