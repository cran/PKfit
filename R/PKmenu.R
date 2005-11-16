#PKfit main menu
PKmenu <- function() 
{
  cat("\n")
  file.menu <- c("Normal Fitting", 
                 "Simulation", 
                 "Quit")
  pick <- menu(file.menu, title = "<< PK startup menu >>")  
  if (pick == 1){
     cat("\n\n")
     nor.fit()
  } 
  else if (pick == 2){
     cat("\n\n") 
     PK.sim()
  } 
  else if (pick == 3){
     cat("\nQuitting menu !!\n\n")     
     
  }  
}

#Normal fitting
nor.fit <- function()
{
  file.menu <- c("Data Manipulation",
                 "Selection of PK Model",
                 "Go Back One Upper Level")
  pick <- menu(file.menu, title = "<< Normal fitting >>")
  if (pick == 1){
     cat("\n\n")      
     data.manipulate()
  } 
  else if (pick == 2){
     cat("\n\n") 
     PK.fit()
  }      
  else if (pick == 3){
     cat("\n\n") 
     PKmenu()
  }  
}

#Data manipulation
data.manipulate <- function()
{
  file.menu <- c("Load Data Files", 
                 "Key in Data", 
                 "Go Back One Upper Level")  
  pick <- menu(file.menu, title = "<< Data edit >>")
  
  if (pick == 1){
     cat("\n\n<< Enter data file name(.csv) >>\n")
     cat("<< The data should consist of subject number, time, concnetration >>\n")
     cat("<< The file should be put in the working directory >>\n")
     PK.file<-scan(what=character(),strip.white=TRUE,nlines=1,quiet=TRUE)
     cnames<-c("Subject", "time", "conc")
     PKindex<-read.csv(PK.file,header=TRUE,sep=",",row.names=NULL,col.names=cnames)
     PKindex<-edit(PKindex)
     show(PKindex)
     cat("\n\n")
     save(PKindex,file="PK.RData")
     return(nor.fit())   
  } 
  else if (pick == 2){
     cat("\n\n") 
     PKindex<-data.frame(Subject=c(1),time=c(0),conc=c(0))
     PKindex<-edit(PKindex)
     show(PKindex)     
     cat("\n\n")
     save(PKindex,file="PK.RData")
     return(nor.fit())
  } 
  else if (pick == 3){
     cat("\n\n") 
     nor.fit()
  } 
}

#PK model option
PK.fit <- function()
{
  file.menu <- c("One Compartment PK Model", 
                 "Two Compartment PK Model",
                 "Macroconstant Exponential Functions",
                 "Go Back One Upper Level")
  pick <- menu(file.menu, title = "<< List of PK models >>")
  if (pick== 1){
     cat("\n\n")  
     one.list()
  }     
  else if (pick == 2){
     cat("\n\n")
     two.list()
  }
  else if (pick == 3){
     cat("\n\n")
     macro()
  }
  else if (pick == 4){
     cat("\n\n")
     nor.fit()
  }        
}

#------------Normal fitting menu----------------
one.list <- function()
{
  file.menu <- c("IV (Bolus, Infusion)", 
                 "Non IV route",
                 "Go Back One Upper Level")
  pick <- menu(file.menu, title = "<< Administration >>")
  if (pick == 1){
     cat("\n\n")  
     iv.route()
  }
  else if (pick == 2){
     cat("\n\n")
     noniv.route()
  }
  else if (pick == 3){
     cat("\n\n")
     PK.fit()
  }              
}

two.list <- function()
{
  file.menu <- c("IV bolus single dose", 
                 "IV infusion single dose",
                 "Extravascular single dose first-order absorption without lag time",
                 "Go Back One Upper Level")
  pick <- menu(file.menu, title = "<< Administration >>")
  if (pick == 1){
     cat("\n\n")  
     fbolus2()
  }
  else if (pick == 2){
     cat("\n\n")
     finfu2()
  }
  else if (pick == 3){
     cat("\n\n")
     ffirst2()
  }
  else if (pick == 4){
     cat("\n\n")
     PK.fit()
  }
}

PK.sim <- function()
{
  file.menu <- c("One compartment model: IV (Bolus, Infusion, Intermediate infusion)", 
                 "One compartment model: non IV route",
                 "Two compartment model",
                 "Macroconstant exponential functions",
                 "Go Back One Upper Level")
  pick <- menu(file.menu, title = "<< Selection of PK Model >>")
  if (pick == 1){
     cat("\n\n")  
     sone.iv.route()
  }
  else if (pick == 2){
     cat("\n\n")
     sone.noniv.route()
  }
  else if (pick == 3){
     cat("\n\n")
     stwo.all()
  }      
  else if (pick == 4){
     cat("\n\n")
     smacro()
  }
  else if (pick == 5){
     cat("\n\n")
     PK.sim()
  }
}

iv.route <- function()
{
  file.menu <- c("IV Bolus Single Dose", 
                 "IV Bolus Single Dose Nonlinear Elimination", 
                 "IV Infusion Single Dose",
                 "IV Infusion Single Dose Nonlinear Elimination",
                 "Go back one upper level")
  pick <- menu(file.menu, title = "<< IV (Bolus, Infusion, Intermediate infusion) >>")
  if (pick ==1){
     cat("\n\n")
     fbolus1()
  }
  else if (pick == 2){
     cat("\n\n") 
     fbolus.mm()
  }
  else if (pick == 3){
     cat("\n\n") 
     finfu1()
  }
  else if (pick == 4){
     cat("\n\n") 
     finfu.mm()
  }
  else if (pick == 5){
     cat("\n\n") 
     one.list()
  }
}

noniv.route <- function()
{
  file.menu <- c("Single Dose First-Order Absorption with Lag Time",
                 "Single Dose First-Order Absorption without Lag Time",
                 "Single Dose Zero-Order Absorption without Lag Time",
                 "Single Dose First-Order Absorption with Lag Time Nonlinear Elimination",
                 "Single Dose First-Order Absorption without Lag Time Nonlinear Elimination",
                 "Single Dose Zero-Order Absorption without Lag Time Nonlinear Elimination",
                 "Go back one upper level")
  pick <- menu(file.menu, title = "<< Non IV route >>")
  if (pick == 1){
     cat("\n\n") 
     ffirst.lag()
  }
  else if (pick == 2){
     cat("\n\n") 
     ffirst.nolag()
  }
  else if (pick == 3){
     cat("\n\n") 
     fzero.nolag()
  }
  else if (pick == 4){
     cat("\n\n") 
     ffirst.lagm()
  }
  else if (pick == 5){
     cat("\n\n") 
     ffirst.nolagm()
  }
  else if (pick == 6){
     cat("\n\n") 
     fzero.nolagm()
  }
  else if (pick == 7){
     cat("\n\n") 
     one.list()
  }
}

macro <- function()
{
  file.menu <- c("One Exponential Term", 
                 "Two Exponential Terms",
                 "Three Exponential Terms",
                 "Go Back One Upper Level")
  pick <- menu(file.menu, title = "<< Administration >>")
  if (pick ==1){
     cat("\n\n")
     fmacro.one()
  }
  else if (pick == 2){
     cat("\n\n") 
     fmacro.two()
  }
  else if (pick == 3){
     cat("\n\n") 
     fmacro.three()
  }
  else if (pick == 4){
     cat("\n\n") 
     PK.fit()
  }
}

#------------Simulation menu----------------
sone.iv.route <- function()
{
  file.menu <- c("IV Bolus Single Dose", 
                 "IV Bolus Single Dose Nonlinear Elimination", 
                 "IV Infusion Single Dose",
                 "IV Infusion Single Dose Nonlinear Elimination",
                 "Go back one upper level")
  pick <- menu(file.menu, title = "<< IV (Bolus, Infusion, Intermediate infusion) >>")
  if (pick ==1 ){
     cat("\n\n")
     sbolus1()
  }
  else if (pick == 2){
     cat("\n\n") 
     sbolus.mm()
  }
  else if (pick == 3){
     cat("\n\n") 
     sinfu1()
  }
  else if (pick == 4){
     cat("\n\n") 
     sinfu.mm()
  }
  else if (pick == 5){
     cat("\n\n") 
     PK.sim()
  }
}

sone.noniv.route <- function()
{
  file.menu <- c("Single Dose First-Order Absorption with Lag Time",
                 "Single Dose First-Order Absorption without Lag Time",
                 "Single Dose Zero-Order Absorption without Lag Time",
                 "Single Dose First-Order Absorption with Lag Time Nonlinear Elimination",
                 "Single Dose First-Order Absorption without Lag Time Nonlinear Elimination",
                 "Single Dose Zero-Order Absorption without Lag Time Nonlinear Elimination",
                 "Go back one upper level")
  pick <- menu(file.menu, title = "<< Non IV route >>")
  if (pick == 1){
     cat("\n\n") 
     sfirst.lag()
  }
  else if (pick == 2){
     cat("\n\n") 
     sfirst.nolag()
  }
  else if (pick == 3){
     cat("\n\n") 
     szero.nolag()
  }
  else if (pick == 4){
     cat("\n\n") 
     sfirst.lagm()
  }
  else if (pick == 5){
     cat("\n\n") 
     sfirst.nolagm()
  }
  else if (pick == 6){
     cat("\n\n") 
     szero.nolagm()
  }
  else if (pick == 7){
     cat("\n\n") 
     PK.sim()
  }
}

smacro <- function()
{
  file.menu <- c("One exponential term", 
                 "Two exponential terms", 
                 "Three exponential terms",
                 "Go back one upper level")
  pick <- menu(file.menu, title = "<< Macroconstant >>")
  if (pick ==1){
     cat("\n\n")
     smacro.one()
  }
  else if (pick == 2){
     cat("\n\n") 
     smacro.two()
  }
  else if (pick == 3){
     cat("\n\n") 
     smacro.three()
  }
  else if (pick == 4){
     cat("\n\n") 
     PK.sim()
  }
}

stwo.all <- function()
{
  file.menu <- c("IV bolus single dose", 
                 "IV infusion single dose",
                 "Extravascular single dose first-order absorption without lag time",
                 "Go Back One Upper Level")
  pick <- menu(file.menu, title = "<< Two compartment >>")
  if (pick ==1){
     cat("\n\n")
     sbolus2()
  }
  else if (pick == 2){
     cat("\n\n") 
     sinfu2()
  }
  else if (pick == 3){
     cat("\n\n") 
     sfirst2()
  }
  else if (pick == 4){
     cat("\n\n") 
     PK.sim()
  }  
}

#---------------plot for linear----------------
plotting.lin <- function (PKindex, fm, i, pick, coef, xaxis, yaxis)
{               
  j<-1:length(PKindex$time[PKindex$Subject==i])
  x<-PKindex$time[PKindex$Subject==i]
  y<-PKindex$conc[PKindex$Subject==i]
        
 #Calculated concentration
  cal<-predict(fm,list(time=x))
    
 #Weighted residuals
  if (pick ==1){
     wei<-ifelse(y[j]==0.0, 0, y[j]-cal[j])
  }
  else if (pick ==2){
     wei<-ifelse(y[j]==0.0, 0, sqrt(1/(y[j]))*(y[j]-cal[j]))
  }
  else if (pick ==3){  
     wei<-ifelse(y[j]==0.0, 0, sqrt(1/((y[j])^2))*(y[j]-cal[j]))
  }
  
 #calculate AUC and AUMC   
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
  add<-add(x,y)
  AUC<-add$auc
  AUMC<-add$aumc
        
 #Output   
  cat("<< Output >>\n")  
  output<-data.frame(x,y,cal,wei,AUC,AUMC)
  colnames(output)<-list("time","Observed","Calculated","Wtd Residuals","AUC","AUMC")
  show(output)  
  
 #AUC (0 to infinity)              
  cat("\n<< AUC (0 to infinity) computed by trapezoidal rule >>\n\n")
  auc.infinity<-y[length(y)]/coef[1,1]
  auc<-AUC[length(y)]+auc.infinity
  show(auc) 
  
 #AUMC (0 to infinity) 
  cat("\n<< AUMC (0 to infinity) computed by trapezoidal rule >>\n\n")
  aumc.infinity<-(x[length(x)]*y[length(y)])/coef[1,1]+x[length(x)]/((coef[1,1])^2)
  aumc<-AUMC[length(y)]+aumc.infinity
  show(aumc)   
  cat("\n")
         
  aicllsbc(fm)
  cat("\n\n")
    
 #Divide plot console into four parts
  windows()
  par(mfrow=c(2,2))
  main<-paste(c("Subject", i ,"plot"),collapse=" ")
     
 #Linear plot
  plot(y~x,data=PKindex,type='p',main=main, 
  xlab=xaxis, ylab=yaxis,pch=15,col="black",bty="l",
  font.lab=2,cex.lab=1,cex.axis=1,cex.main=1) 
  lines(x,predict(fm,list(time=x)),type="l",lty=1,
  col="firebrick3",lwd="2")
  mtext("Linear",side=3,cex=0.88)
    
 #Semi-log plot
  plot(x,y,log="y",type='p',main=main,
  xlab=xaxis, ylab=yaxis,pch=15,col="black",bty="l",
  font.lab=2,cex.lab=1,cex.axis=1,cex.main=1) 
  lines(x,predict(fm,list(time=x)),type="l",lty=1,
  col="firebrick3",lwd="2")
  mtext("Semi-log",side=3,cex=0.88)
     
 #Residual plot, time vs weighted residual
  plot(x,wei,pch=15,col="blue",bty="l",xlab=xaxis,
  ylab="Weighted Residual",main="Residual Plots",cex.lab=1,
  cex.axis=1,cex.main=1,font.lab=2)
  abline(h=0,lwd=2,col="black",lty=2)
    
 #Residual plot, calculated concentration vs weigthed residual
  plot(cal,wei,pch=15,col="blue",bty="l",xlab="Calc Cp(i)",
  ylab="Weighted Residual",main="Residual Plots",cex.lab=1,
  cex.axis=1,cex.main=1,font.lab=2)
  abline(h=0,lwd=2,col="black",lty=2)          
}    

#---------------plot for nonlinear----------------
plotting.non <- function (PKindex, fm, i, pick, xaxis, yaxis)
{  
  options(warn=-1)
              
  j<-1:length(PKindex$time[PKindex$Subject==i])
  x<-PKindex$time[PKindex$Subject==i]
  y<-PKindex$conc[PKindex$Subject==i]
   
 #Calculated concentration
  cal<-predict(fm,list(time=x))
    
 #Weighted residuals
  if (pick ==1){
     wei<-ifelse(y[j]==0.0, 0, y[j]-cal[j])
  }
  else if (pick ==2){
     wei<-ifelse(y[j]==0.0, 0, sqrt(1/(y[j]))*(y[j]-cal[j]))
  }
  else if (pick ==3){  
     wei<-ifelse(y[j]==0.0, 0, sqrt(1/((y[j])^2))*(y[j]-cal[j]))
  }
                   
 #calculate AUC and AUMC       
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
  add<-add(x,y)
  AUC<-add$auc
  AUMC<-add$aumc
              
 #Output  
  cat("<< Output >>\n\n")      
  output<-data.frame(x,y,cal,wei,AUC,AUMC)
  colnames(output)<-list("time","Observed","Calculated","Wtd Residuals","AUC","AUMC")
  show(output)  
  cat("\n") 
        
  aicllsbc(fm)
  cat("\n\n")
    
 #Divide plot console into four parts
  windows()
  par(mfrow=c(2,2))
  main<-paste(c("Subject", i ,"plot"),collapse=" ")
     
 #Linear plot
  plot(y~x,data=PKindex,type='p',main=main, 
  xlab="Time", ylab="Concentration",pch=15,col="black",bty="l",
  font.lab=2,cex.lab=1,cex.axis=1,cex.main=1) 
  lines(x,predict(fm,list(time=x)),type="l",lty=1,
  col="firebrick3",lwd="2")
  mtext("Linear",side=3,cex=0.88)
    
 #Semi-log plot
  plot(x,y,log="y",type='p',main=main,
  xlab="Time", ylab="Concentration",pch=15,col="black",bty="l",
  font.lab=2,cex.lab=1,cex.axis=1,cex.main=1) 
  lines(x,predict(fm,list(time=x)),type="l",lty=1,
  col="firebrick3",lwd="2")
  mtext("Semi-log",side=3,cex=0.88)
     
 #Residual plot, time vs weighted residual----- 
  plot(x,wei,pch=15,col="blue",bty="l",xlab="Time",
  ylab="Weighted Residual",main="Residual Plots",cex.lab=1,
  cex.axis=1,cex.main=1,font.lab=2)
  abline(h=0,lwd=2,col="black",lty=2)
    
 #Residual plot, calcukated concentration vs weigthed residual-----
  plot(cal,wei,pch=15,col="blue",bty="l",xlab="Calc Cp(i)",
  ylab="Weighted Residual",main="Residual Plots",cex.lab=1,
  cex.axis=1,cex.main=1,font.lab=2)
  abline(h=0,lwd=2,col="black",lty=2)
}    

#----------------plot for simulation----------------
plotting.sim <- function(i,x,y)
{
  main<-paste(c("Subject", i ,"plot"),collapse=" ")
  options(warn=-1)
  windows()
  par(mfrow=c(2,2))
   
 #linear plot
  plot(y~x,type='p',main=main, 
  xlab="Time",ylab="Concentration",pch=15,col="black",bty="l",
  font.lab=2,cex.lab=1,cex.axis=1,cex.main=1) 
   text("Linear",side=3,cex=0.88)

 #semi-log plot
  plot(x,y,log="y",type='p',main=main,
  xlab="Time",ylab="Concentration",pch=15,col="black",bty="l",
  font.lab=2,cex.lab=1,cex.axis=1,cex.main=1) 
  mtext("Semi-log",side=3,cex=0.88) 
  
  cat("\n\n")
}

#Estimate model fitting: AIC, Log likelihood, SBC
aicllsbc <- function(fm)
{   
  cat("\n") 
  cat("<< Akaike's Information Criterion (AIC) >>\n\n")
  show(AIC(fm))
    
  cat("\n<< Log likelihood >>\n\n")
  show(logLik(fm))
    
  cat("\n<< Schwarz's Bayesian Criterion (SBC) >>\n\n")
  show(BIC(fm))
  cat("\n")     
    
 #Summary the results of nls
  print(summary(fm))      
}  

#for iv bolus 
sbolus1.out<-function(PKindex,kel,Vd,defun,par1,par2,Dose,i)
{
  time<-PKindex$time
  parms<-c(kel=kel,Vd=Vd) 
  C1.lsoda<-data.frame(lsoda(Dose/Vd,c(0,time),defun,parms)) 
  cat("\n")
  sim1<-kel
  sim2<-Vd
  sim<-matrix(c(sim1,sim2,par1,par2),2,2)
  dimnames(sim)<-list(c("kel","Vd"),c("Value","Original"))
  show(sim)
  cat("\n\n<< Output >>\n\n")
  good<-ifelse(C1.lsoda[2:(length(time)+1),2]<=0,0,C1.lsoda[2:(length(time)+1),2])
  PKdata<-data.frame(i,C1.lsoda[2:(length(time)+1),1],good)
  colnames(PKdata)<-list("Subject","time","concentration")
  show(PKdata)
  x<-C1.lsoda[2:(length(time)+1),1]
  y<-good
  plotting.sim(i,x,y)
}

#for iv bolus nonlinear elimination 
sbolus.mm.out<-function(PKindex,Vm,Km,Vd,defun,par1,par2,par3,Dose,i)
{ 
  time<-PKindex$time
  parms<-c(Vm=Vm,Km=Km,Vd=Vd)  
  C1.lsoda<-data.frame(lsoda(Dose/Vd,c(0,time),defun,parms))     
  cat("\n")
  sim1<-Vm
  sim2<-Km
  sim3<-Vd
  sim<-matrix(c(sim1,sim2,sim3,par1,par2,par3),3,2)
  dimnames(sim)<-list(c("Vm","Km","Vd"),c("Value","Original"))
  show(sim)
  cat("\n\n<< Output >>\n\n")  
  good<-ifelse(C1.lsoda[2:(length(time)+1),2]<=0,0,C1.lsoda[2:(length(time)+1),2])
  PKdata<-data.frame(i,C1.lsoda[2:(length(time)+1),1],good)
  colnames(PKdata)<-list("Subject","time","concentration")
  show(PKdata)  
  x<-C1.lsoda[2:(length(time)+1),1]
  y<-good
  plotting.sim(i,x,y)
}

#for iv infusion 
sinfu1.out<-function(PKindex,kel,Vd,defun,par1,par2,Dose,i)
{
  time<-PKindex$time
  parms<-c(kel=kel,Vd=Vd) 
  C1.lsoda<-data.frame(lsoda(0,c(0,time),defun,parms)) 
  cat("\n")
  sim1<-kel
  sim2<-Vd
  sim<-matrix(c(sim1,sim2,par1,par2),2,2)
  dimnames(sim)<-list(c("kel","Vd"),c("Value","Original"))
  show(sim)
  cat("\n\n<< Output >>\n\n")
  good<-ifelse(C1.lsoda[2:(length(time)+1),2]<=0,0,C1.lsoda[2:(length(time)+1),2])
  PKdata<-data.frame(i,C1.lsoda[2:(length(time)+1),1],good)
  colnames(PKdata)<-list("Subject","time","concentration")
  show(PKdata)
  x<-C1.lsoda[2:(length(time)+1),1]
  y<-good  
  plotting.sim(i,x,y)
}

#for iv infusion nonlinear elimination
sinfu.mm.out<-function(PKindex,Vm,Km,Vd,defun,par1,par2,par3,Dose,i)
{
  time<-PKindex$time
  parms<-c(Vm=Vm,Km=Km,Vd=Vd)  
  C1.lsoda<-data.frame(lsoda(0,c(0,time),defun,parms))    
  cat("\n")
  sim1<-Vm
  sim2<-Km
  sim3<-Vd
  sim<-matrix(c(sim1,sim2,sim3,par1,par2,par3),3,2)
  dimnames(sim)<-list(c("Vm","Km","Vd"),c("Value","Original"))
  show(sim)    
  cat("\n\n<< Output >>\n\n")   
  good<-ifelse(C1.lsoda[2:(length(time)+1),2]<=0,0,C1.lsoda[2:(length(time)+1),2])   
  PKdata<-data.frame(i,C1.lsoda[2:(length(time)+1),1],good)
  colnames(PKdata)<-list("Subject","time","concentration")
  show(PKdata) 
  x<-C1.lsoda[2:(length(time)+1),1]
  y<-good  
  plotting.sim(i,x,y)
}

#for extravascular first order absorption with/without lag time
sfirst1.out<-function(PKindex,ka,kel,Vd,defun,par1,par2,par3,Dose,i)  
{   
  time<-PKindex$time
  parms<-c(ka=ka,kel=kel,Vd=Vd)  
  C1.lsoda<-data.frame(lsoda(c(Dose,0),c(0,time), defun, parms))  
  cat("\n\n")
  sim1<-ka
  sim2<-kel
  sim3<-Vd
  sim<-matrix(c(sim1,sim2,sim3,par1,par2,par3),3,2)
  dimnames(sim)<-list(c("ka","kel","Vd"),c("Value","Original"))
  show(sim)
  cat("\n\n<< Output >>\n\n") 
  good<-ifelse(C1.lsoda[2:(length(time)+1),3]<=1e-5, 0, C1.lsoda[2:(length(time)+1),3])   
  PKdata<-data.frame(i,C1.lsoda[2:(length(time)+1),1],good)
  colnames(PKdata)<-list("Subject","time","concentration")
  show(PKdata)
  x<-C1.lsoda[2:(length(time)+1),1]
  y<-good 
  plotting.sim(i,x,y)
}  

#for extravascular first order absorption with/without lag time nonlinear elimination
sfirst.mm.out<-function(PKindex,ka,Vm,Km,Vd,defun,par1,par2,par3,par4,Dose,i)  
{       
  time<-PKindex$time
  parms<-c(ka=ka,Vm=Vm,Km=Km,Vd=Vd) 
  C1.lsoda<-data.frame(lsoda(c(Dose,0),c(0,time), defun, parms)) 
  cat("\n\n")
  sim1<-ka
  sim2<-Vm
  sim3<-Km
  sim4<-Vd
  sim<-matrix(c(sim1,sim2,sim3,sim4,par1,par2,par3,par4),4,2)
  dimnames(sim)<-list(c("ka","Vm","Km","Vd"),c("Value","Original"))
  show(sim)
  cat("\n\n<< Output >>\n\n") 
  good<-ifelse(C1.lsoda[2:(length(time)+1),3]<=1e-5, 0, C1.lsoda[2:(length(time)+1),3])   
  PKdata<-data.frame(i,C1.lsoda[2:(length(time)+1),1],good)
  colnames(PKdata)<-list("Subject","time","concentration")
  show(PKdata)
  x<-C1.lsoda[2:(length(time)+1),1]
  y<-good
  plotting.sim(i,x,y) 
}

#for extravascular zero order absorption without lag time
szero.out<-function(PKindex,Tabs,kel,Vd,defun,par1,par2,par3,Dose,i)    
{
  time<-PKindex$time
  parms<-c(Tabs=Tabs,kel=kel,Vd=Vd)  
  C1.lsoda<-data.frame(lsoda(0, c(0,time), defun, parms))   
  cat("\n")
  sim1<-Tabs
  sim2<-kel
  sim3<-Vd
  sim<-matrix(c(sim1,sim2,sim3,par1,par2,par3),3,2)
  dimnames(sim)<-list(c("Tabs","kel","Vd"),c("Value","Original"))
  show(sim)  
  cat("\n\n<< Output >>\n\n")   
  good<-ifelse(C1.lsoda[2:(length(time)+1),2]<=0,0,C1.lsoda[2:(length(time)+1),2])  
  PKdata<-data.frame(i,C1.lsoda[2:(length(time)+1),1],good)
  colnames(PKdata)<-list("Subject","time","concentration")
  show(PKdata)
  x<-C1.lsoda[2:(length(time)+1),1]
  y<-C1.lsoda[2:(length(time)+1),2]
  plotting.sim(i,x,y)
}

#for extravascular zero order absorption without lag time nonlinear elimination
szero.mm.out<-function(PKindex,Tabs,Vm,Km,Vd,defun,par1,par2,par3,par4,Dose,i)
{                      
  time<-PKindex$time
  parms<-c(Tabs=Tabs,Vm=Vm,Km=Km,Vd=Vd) 
  C1.lsoda<-data.frame(lsoda(0,c(0,time), defun, parms))   
  cat("\n\n")
  sim1<-Tabs
  sim2<-Vm
  sim3<-Km
  sim4<-Vd
  sim<-matrix(c(sim1,sim2,sim3,sim4,par1,par2,par3,par4),4,2)
  dimnames(sim)<-list(c("Tabs","Vm","Km","Vd"),c("Value","Original"))
  show(sim)
  cat("\n\n<< Output >>\n\n") 
  good<-ifelse(C1.lsoda[2:(length(time)+1),2]<=0,0,C1.lsoda[2:(length(time)+1),2]) 
  PKdata<-data.frame(i,C1.lsoda[2:(length(time)+1),1],good)
  colnames(PKdata)<-list("Subject","time","concentration")
  show(PKdata) 
  x<-C1.lsoda[2:(length(time)+1),1]
  y<-C1.lsoda[2:(length(time)+1),2]   
  plotting.sim(i,x,y)
}

#for one exponential term
smacro.one.out<-function(PKindex,A,a,defun,par1,par2,Dose,i)
{         
  time<-PKindex$time
  defun<- A*exp(-a*time)
  output<-data.frame(time,defun) 
  cat("\n")
  sim1<-A
  sim2<-a
  sim<-matrix(c(sim1,sim2,par1,par2),2,2)
  dimnames(sim)<-list(c("A","a"),c("Value","Original"))
  show(sim)
  cat("\n\n<< Output >>\n\n")
  PKdata<-data.frame(i,output)
  colnames(PKdata)<-list("Subject","time","concentration")
  show(PKdata)
  x<-PKdata[,2]
  y<-PKdata[,3]
  plotting.sim(i,x,y)
}

#for two exponential term
smacro.two.out<-function(PKindex,A,a,B,b,defun,par1,par2,par3,par4,Dose,i)      
{   
  time<-PKindex$time
  defun<- A*exp(-a*time)+B*exp(-b*time)
  output<-data.frame(time,defun) 
  cat("\n")
  sim1<-A
  sim2<-a
  sim3<-B
  sim4<-b
  sim<-matrix(c(sim1,sim2,sim3,sim4,par1,par2,par3,par4),4,2)
  dimnames(sim)<-list(c("A","a","B","b"),c("Value","Original"))
  show(sim)
  cat("\n\n<< Output >>\n\n")
  PKdata<-data.frame(i,output)
  colnames(PKdata)<-list("Subject","time","concentration")
  show(PKdata)
  x<-PKdata[,2]
  y<-PKdata[,3]
  plotting.sim(i,x,y)
}

#for three exponential term
smacro.three.out<-function(PKindex,A,a,B,b,C,c,defun,par1,par2,par3,par4,par5,par6,Dose,i)      
{         
  time<-PKindex$time
  defun<- A*exp(-a*time)+B*exp(-b*time)+C*exp(-c*time)
  output<-data.frame(time,defun) 
  cat("\n")
  sim1<-A
  sim2<-a
  sim3<-B
  sim4<-b
  sim5<-C
  sim6<-c
  sim<-matrix(c(sim1,sim2,sim3,sim4,sim5,sim6,par1,par2,par3,par4,par5,par6),6,2)
  dimnames(sim)<-list(c("A","a","B","b","C","c"),c("Value","Original"))
  show(sim)
  cat("\n\n<< Output >>\n\n")
  PKdata<-data.frame(i,output)
  colnames(PKdata)<-list("Subject","time","concentration")
  show(PKdata)
  x<-PKdata[,2]
  y<-PKdata[,3]
  plotting.sim(i,x,y)
}

#for two compartment iv bolus
sbolus2.out<-function(PKindex,kel,k12,k21,Vd,defun,par1,par2,par3,par4,Dose,i)
{ 
  time<-PKindex$time
  parms<-c(kel=kel,k12=k12,k21=k21,Vd=Vd) 
  C1.lsoda<-data.frame(lsoda(c(Dose/Vd,0),c(0,time),defun,parms)) 
  cat("\n")
  sim1<-kel
  sim2<-k12
  sim3<-k21
  sim4<-Vd
  sim<-matrix(c(sim1,sim2,sim3,sim4,par1,par2,par3,par4),4,2)
  dimnames(sim)<-list(c("kel","k12","k21","Vd"),c("Value","Original"))
  show(sim)
  cat("<<\n\nOutput\n\n>>")
  good<-ifelse(C1.lsoda[2:(length(time)+1),2]<=0,0,C1.lsoda[2:(length(time)+1),2])
  PKdata<-data.frame(i,C1.lsoda[2:(length(time)+1),1],good)
  colnames(PKdata)<-list("Subject","time","concentration")
  show(PKdata)
  x<-C1.lsoda[2:(length(time)+1),1]
  y<-C1.lsoda[2:(length(time)+1),2]
  plotting.sim(i,x,y)
}

#for two compartment iv infusion
sinfu2.out<-function(PKindex,kel,k12,k21,Vd,defun,par1,par2,par3,par4,Dose,i)
{       
  time<-PKindex$time
  parms<-c(kel=kel,k12=k12,k21=k21,Vd=Vd) 
  C1.lsoda<-data.frame(lsoda(c(0,0),c(0,time),defun,parms)) 
  cat("\n")
  sim1<-kel
  sim2<-k12
  sim3<-k21
  sim4<-Vd
  sim<-matrix(c(sim1,sim2,sim3,sim4,par1,par2,par3,par4),4,2)
  dimnames(sim)<-list(c("kel","k12","k21","Vd"),c("Value","Original"))
  show(sim)
  cat("\n\n<< Output >>\n\n")
  good<-ifelse(C1.lsoda[2:(length(time)+1),2]<=0,0,C1.lsoda[2:(length(time)+1),2]) 
  PKdata<-data.frame(i,C1.lsoda[2:(length(time)+1),1],good)
  colnames(PKdata)<-list("Subject","time","concentration")
  show(PKdata)
  x<-C1.lsoda[2:(length(time)+1),1]
  y<-C1.lsoda[2:(length(time)+1),2]
  plotting.sim(i,x,y)
}

#for two compartment extravascular first order absorption
sfirst2.out<-function(PKindex,ka,kel,k12,k21,Vd,defun,par1,par2,par3,par4,par5,Dose,i)
{                   
  time<-PKindex$time
  parms<-c(ka=ka,kel=kel,k12=k12,k21=k21,Vd=Vd) 
  C1.lsoda<-data.frame(lsoda(c(Dose,0,0),c(0,time),defun,parms)) 
  cat("\n")
  sim1<-ka
  sim2<-kel
  sim3<-k12
  sim4<-k21
  sim5<-Vd
  sim<-matrix(c(sim1,sim2,sim3,sim4,sim5,par1,par2,par3,par4,par5),5,2)
  dimnames(sim)<-list(c("ka","kel","k12","k21","Vd"),c("Value","Original"))
  show(sim)
  cat("\n\n<< Output >>\n\n")
  good<-ifelse(C1.lsoda[2:(length(time)+1),2]<=0,0,C1.lsoda[2:(length(time)+1),3])
  PKdata<-data.frame(i,C1.lsoda[2:(length(time)+1),1],good)
  colnames(PKdata)<-list("Subject","time","concentration")
  show(PKdata)
  x<-C1.lsoda[2:(length(time)+1),1]
  y<-C1.lsoda[2:(length(time)+1),3]
  plotting.sim(i,x,y)
}
