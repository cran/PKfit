#Normal fitting
#One compartment PK model extravascular single dose zero-order absorption without lag time
fzero.nolag <- function(PKindex)
{ 
  options(warn=-1)
   
 #lsoda is belong to odesolve package
  library(odesolve)
   
 #BIC function is belong to stats4 package
  library(stats4)
   
 #genoud function is belong to rgenoud package   
  library(rgenoud)
   
 #Input dose and initial value for Tabs, kel, Vd
  cat("Enter Dose Value\n")
  Dose<-scan(nlines=1,quiet=TRUE)
  cat("\n")
   
  par<-data.frame(Parameter=c("Tabs","kel","Vd"),Initial=c(0))
  par<-edit(par)
  show(par)
  cat("\n")
   
 #User-supplied function
  defun<- function(time, y, parms) { 
     if(time<=parms["Tabs"]) 
       dCpdt <- (Dose/parms["Tabs"])/parms["Vd"] - parms["kel"] * y[1]
     else
       dCpdt <- - parms["kel"] * y[1]
     list(dCpdt) 
  } 
   
  modfun <- function(time,Tabs,kel,Vd) { 
     out <- lsoda(0,c(0,time),defun,parms=c(Tabs=Tabs,kel=kel,Vd=Vd),
                  rtol=1e-8,atol=1e-8) 
     out[-1,2]
  } 
   
 #Selection of weighting schemes-----
  file.menu <- c("equal weight", 
                 "1/Cp",
                 "1/Cp^2")          
  pick <- menu(file.menu, title = "<< Weighting Schemes >>")
  if (pick ==1){
     cat("\nEnter the title of x-axis(time)\n\n")
     xaxis<-readline()
     if (substr(xaxis, 1, 1) == "")  xaxis<-"Time"  else xaxis<-xaxis
     cat("\nEnter the title of y-axis(Cp)\n\n")
     yaxis<-readline()
     if (substr(yaxis, 1, 1) == "")  yaxis<-"Concentration"  else yaxis<-yaxis
     for( i in 1:length(unique(PKindex$Subject)))  {
     cat("\n\n               << Subject",i,">>\n\n" )  
     objfun <- function(par) { 
     out<-modfun(PKindex$time[PKindex$Subject==i],par[1],par[2],par[3])
     gift<- which( PKindex$conc[PKindex$Subject==i] != 0 )
     sum((PKindex$conc[PKindex$Subject==i][gift] - out[gift])^2)}
     gen<-genoud(objfun,nvars=3,max=FALSE,pop.size=30,max.generations=20,wait.generations=10,
          starting.value=c(par[1,2],par[2,2],par[3,2]),BFGS=FALSE,print.level=0,boundary.enforcement=0,
          Domains=matrix(c(1,0.01,1,50,1,100),3,2),MemoryMatrix=TRUE,lexical=FALSE)
     cat("<< The value of parameter fitted by genetic algorithm >>\n\n")   
     namegen<-c("Tabs","kel","Vd")
     outgen<-c(gen$par[1],gen$par[2],gen$par[3])
     print(data.frame(Parameter=namegen,Value=outgen))  
     F<-objfun(gen$par)
     opt<-optim(c(gen$par[1],gen$par[2],gen$par[3]),objfun, method="Nelder-Mead")
     cat("\n<< The value of parameter fitted by Nelder-Mead Simplex slgorithm >>\n\n")   
     nameopt<-c("Tabs","kel","Vd")
     outopt<-c(opt$par[1],opt$par[2],opt$par[3])
     print(data.frame(Parameter=nameopt,Value=outopt))  
     cat("\n<< Residual sum-of-squares and parameter values fitted by nls >>\n\n")
     fm<-nls(conc~modfun(time,Tabs,kel,Vd),data=PKindex, algorithm="default",subset=Subject==i,
         start=list(Tabs=opt$par[1],kel=opt$par[2],Vd=opt$par[3]),trace=TRUE,
         nls.control(tol=1))
     cat("\n")
     coef<-data.frame(coef(fm)["kel"])
     plotting.lin(PKindex, fm, i, pick, coef, xaxis, yaxis)
     }           
  }
  else if (pick == 2){
     cat("\nEnter the title of x-axis(time)\n\n")
     xaxis<-readline()
     if (substr(xaxis, 1, 1) == "")  xaxis<-"Time"  else xaxis<-xaxis
     cat("\nEnter the title of y-axis(Cp)\n\n")
     yaxis<-readline()
     if (substr(yaxis, 1, 1) == "")  yaxis<-"Concentration"  else yaxis<-yaxis
     for( i in 1:length(unique(PKindex$Subject)))  {
     cat("\n\n               << Subject",i,">>\n\n" )  
     objfun <- function(par) { 
     out<-modfun(PKindex$time[PKindex$Subject==i],par[1],par[2],par[3])
     gift<- which( PKindex$conc[PKindex$Subject==i] != 0 )
     sum((PKindex$conc[PKindex$Subject==i][gift]-out[gift])^2/PKindex$conc[PKindex$Subject==i][gift])}
     gen<-genoud(objfun,nvars=3,max=FALSE,pop.size=30,max.generations=20,wait.generations=10,
          starting.value=c(par[1,2],par[2,2],par[3,2]),BFGS=FALSE,print.level=0,boundary.enforcement=0,
          Domains=matrix(c(1,0.01,1,50,1,100),3,2),MemoryMatrix=TRUE,lexical=FALSE)
     cat("<< The value of parameter fitted by genetic algorithm >>\n\n")   
     namegen<-c("Tabs","kel","Vd")
     outgen<-c(gen$par[1],gen$par[2],gen$par[3])
     print(data.frame(Parameter=namegen,Value=outgen))  
     F<-objfun(gen$par)
     opt<-optim(c(gen$par[1],gen$par[2],gen$par[3]),objfun, method="Nelder-Mead")
     cat("\n<< The value of parameter fitted by Nelder-Mead Simplex slgorithm >>\n\n")   
     nameopt<-c("Tabs","kel","Vd")
     outopt<-c(opt$par[1],opt$par[2],opt$par[3])
     print(data.frame(Parameter=nameopt,Value=outopt))
     cat("\n<< Residual sum-of-squares and parameter values fitted by nls >>\n\n")
     fm<-nls(conc~modfun(time,Tabs,kel,Vd),data=PKindex, algorithm="default",subset=Subject==i,
         start=list(Tabs=opt$par[1],kel=opt$par[2],Vd=opt$par[3]),trace=TRUE,
         nls.control(tol=1))
     cat("\n")
     coef<-data.frame(coef(fm)["kel"])
     plotting.lin(PKindex, fm, i, pick, coef, xaxis, yaxis)
     }
  }           
  else if (pick == 3){
     cat("\nEnter the title of x-axis(time)\n\n")
     xaxis<-readline()
     if (substr(xaxis, 1, 1) == "")  xaxis<-"Time"  else xaxis<-xaxis
     cat("\nEnter the title of y-axis(Cp)\n\n")
     yaxis<-readline()
     if (substr(yaxis, 1, 1) == "")  yaxis<-"Concentration"  else yaxis<-yaxis
     for( i in 1:length(unique(PKindex$Subject)))  {
     cat("\n\n               << Subject",i,">>\n\n" )  
     objfun <- function(par) { 
     out<-modfun(PKindex$time[PKindex$Subject==i],par[1],par[2],par[3])
     gift<- which( PKindex$conc[PKindex$Subject==i] != 0 ) 
     sum(((PKindex$conc[PKindex$Subject==i][gift] - out[gift])/PKindex$conc[PKindex$Subject==i][gift])^2)}
     gen<-genoud(objfun,nvars=3,max=FALSE,pop.size=30,max.generations=20,wait.generations=10,
          starting.value=c(par[1,2],par[2,2],par[3,2]),BFGS=FALSE,print.level=0,boundary.enforcement=0,
          Domains=matrix(c(1,0.01,1,50,1,100),3,2),MemoryMatrix=TRUE,lexical=FALSE)
     cat("<< The value of parameter fitted by genetic algorithm >>\n\n")   
     namegen<-c("Tabs","kel","Vd")
     outgen<-c(gen$par[1],gen$par[2],gen$par[3])
     print(data.frame(Parameter=namegen,Value=outgen)) 
     F<-objfun(gen$par)
     opt<-optim(c(gen$par[1],gen$par[2],gen$par[3]),objfun, method="Nelder-Mead")
     cat("\n<< The value of parameter fitted by Nelder-Mead Simplex slgorithm >>\n\n")   
     nameopt<-c("Tabs","kel","Vd")
     outopt<-c(opt$par[1],opt$par[2],opt$par[3])
     print(data.frame(Parameter=nameopt,Value=outopt))
     cat("\n<< Residual sum-of-squares and parameter values fitted by nls >>\n\n")
     fm<-nls(conc~modfun(time,Tabs,kel,Vd),data=PKindex, algorithm="default",subset=Subject==i,
         start=list(Tabs=opt$par[1],kel=opt$par[2],Vd=opt$par[3]),trace=TRUE,
         nls.control(tol=1))
     cat("\n")
     coef<-data.frame(coef(fm)["kel"])
     plotting.lin(PKindex, fm, i, pick, coef, xaxis, yaxis)
     }
  }      
  cat("\n")  
}
    
#One compartment PK model extravascular single dose zero-order absorption without lag time nonlinear elimination
fzero.nolagm <- function(PKindex)
{
  options(warn=-1)
   
 #lsoda is belong to odesolve package
  library(odesolve)
   
 #BIC function is belong to stats4 package
  library(stats4)
 #genoud function is belong to rgenoud package 
  library(rgenoud)
   
 #Input dose and initial value for Tabs, Vm, Km, Vd
  cat("Enter Dose Value\n")
  Dose<-scan(nlines=1,quiet=TRUE)
  cat("\n")
   
  par<-data.frame(Parameter=c("Tabs","Vm","Km","Vd"),Initial=c(0))
  par<-edit(par)
  show(par)
  cat("\n")
   
 #User-supplied function-----
  defun<- function(time, y, parms) { 
     if(time<=parms["Tabs"]) 
       dCpdt <- (Dose/parms["Tabs"])/parms["Vd"]-(parms["Vm"]/parms["Vd"])*y[1]/(parms["Km"]/parms["Vd"]+y[1])            
     else
       dCpdt <- -(parms["Vm"]/parms["Vd"])*y[1]/(parms["Km"]/parms["Vd"]+y[1])            
     list(dCpdt) 
  }   
   
  modfun <- function(time,Tabs,Vm,Km,Vd) { 
     out <- lsoda(0,c(0,time),defun,parms=c(Tabs=Tabs,Vm=Vm,Km=Km,Vd=Vd),
                  rtol=1e-6,atol=1e-8) 
     out[-1,2]
  } 
   
 #Selection of weighting schemes-----
  file.menu <- c("equal weight", 
                 "1/Cp",
                 "1/Cp^2")           
  pick <- menu(file.menu, title = "<< Weighting Schemes >>")
  if (pick ==1){
     cat("\nEnter the title of x-axis(time)\n\n")
     xaxis<-readline()
     if (substr(xaxis, 1, 1) == "")  xaxis<-"Time"  else xaxis<-xaxis
     cat("\nEnter the title of y-axis(Cp)\n\n")
     yaxis<-readline()
     if (substr(yaxis, 1, 1) == "")  yaxis<-"Concentration"  else yaxis<-yaxis
     for( i in 1:length(unique(PKindex$Subject)))  {
     cat("\n\n               << Subject",i,">>\n\n" ) 
     objfun <- function(par) { 
     out<-modfun(PKindex$time[PKindex$Subject==i],par[1],par[2],par[3],par[4])
     gift<- which( PKindex$conc[PKindex$Subject==i] != 0 )
     sum((PKindex$conc[PKindex$Subject==i][gift] - out[gift])^2)}
     gen<-genoud(objfun,nvars=4,max=FALSE,pop.size=30,max.generations=20,wait.generations=10,
          starting.value=c(par[1,2],par[2,2],par[3,2],par[4,2]),BFGS=FALSE,print.level=0,boundary.enforcement=0,
          Domains=matrix(c(1,1,1,1,50,100,100,100),4,2),MemoryMatrix=TRUE,lexical=FALSE)
     cat("<< The value of parameter fitted by genetic algorithm >>\n\n")   
     namegen<-c("Tabs","Vm","Km","Vd")
     outgen<-c(gen$par[1],gen$par[2],gen$par[3],gen$par[4])
     print(data.frame(Parameter=namegen,Value=outgen)) 
     F<-objfun(gen$par)
     opt<-optim(c(gen$par[1],gen$par[2],gen$par[3],gen$par[4]),objfun, method="Nelder-Mead")
     cat("\n<< The value of parameter fitted by Nelder-Mead Simplex slgorithm >>\n\n")   
     nameopt<-c("Tabs","Vm","Km","Vd")
     outopt<-c(opt$par[1],opt$par[2],opt$par[3],opt$par[4])
     print(data.frame(Parameter=nameopt,Value=outopt))
     cat("\n<< Residual sum-of-squares and parameter values fitted by nls >>\n\n")
     fm<-nls(conc~modfun(time,Tabs,Vm,Km,Vd),data=PKindex, algorithm="default",subset=Subject==i,
         start=list(Tabs=opt$par[1],Vm=opt$par[2],Km=opt$par[3],Vd=opt$par[4]),trace=TRUE,
         nls.control(maxiter=50,tol=10,minFactor=1/1024))
     cat("\n")
     plotting.non(PKindex, fm, i, pick,xaxis,yaxis)
     }           
  }
  else if (pick == 2){
     cat("\nEnter the title of x-axis(time)\n\n")
     xaxis<-readline()
     if (substr(xaxis, 1, 1) == "")  xaxis<-"Time"  else xaxis<-xaxis
     cat("\nEnter the title of y-axis(Cp)\n\n")
     yaxis<-readline()
     if (substr(yaxis, 1, 1) == "")  yaxis<-"Concentration"  else yaxis<-yaxis
     for( i in 1:length(unique(PKindex$Subject)))  {
     cat("\n\n               << Subject",i,">>\n\n" ) 
     objfun <- function(par) { 
     out<-modfun(PKindex$time[PKindex$Subject==i],par[1],par[2],par[3],par[4])
     gift<- which( PKindex$conc[PKindex$Subject==i] != 0 )
     sum((PKindex$conc[PKindex$Subject==i][gift]-out[gift])^2/PKindex$conc[PKindex$Subject==i][gift])} 
     gen<-genoud(objfun,nvars=4,max=FALSE,pop.size=30,max.generations=20,wait.generations=10,
          starting.value=c(par[1,2],par[2,2],par[3,2],par[4,2]),BFGS=FALSE,print.level=0,boundary.enforcement=0,
          Domains=matrix(c(1,1,1,1,50,100,100,100),4,2),MemoryMatrix=TRUE,lexical=FALSE)
     cat("<< The value of parameter fitted by genetic algorithm >>\n\n")   
     namegen<-c("Tabs","Vm","Km","Vd")
     outgen<-c(gen$par[1],gen$par[2],gen$par[3],gen$par[4])
     print(data.frame(Parameter=namegen,Value=outgen))  
     F<-objfun(gen$par)
     opt<-optim(c(gen$par[1],gen$par[2],gen$par[3],gen$par[4]),objfun, method="Nelder-Mead")
     cat("\n<< The value of parameter fitted by Nelder-Mead Simplex slgorithm >>\n\n")   
     nameopt<-c("Tabs","Vm","Km","Vd")
     outopt<-c(opt$par[1],opt$par[2],opt$par[3],opt$par[4])
     print(data.frame(Parameter=nameopt,Value=outopt))
     cat("\n<< Residual sum-of-squares and parameter values fitted by nls >>\n\n")
     fm<-nls(conc~modfun(time,Tabs,Vm,Km,Vd),data=PKindex, algorithm="default",subset=Subject==i,
         start=list(Tabs=opt$par[1],Vm=opt$par[2],Km=opt$par[3],Vd=opt$par[4]),trace=TRUE,
         nls.control(maxiter=50,tol=10,minFactor=1/1024))
     cat("\n")
     plotting.non(PKindex, fm, i, pick,xaxis,yaxis)
     }
  }           
  else if (pick == 3){
     cat("\nEnter the title of x-axis(time)\n\n")
     xaxis<-readline()
     if (substr(xaxis, 1, 1) == "")  xaxis<-"Time"  else xaxis<-xaxis
     cat("\nEnter the title of y-axis(Cp)\n\n")
     yaxis<-readline()
     if (substr(yaxis, 1, 1) == "")  yaxis<-"Concentration"  else yaxis<-yaxis
     for( i in 1:length(unique(PKindex$Subject)))  {
     cat("\n\n               << Subject",i,">>\n\n" ) 
     objfun <- function(par) { 
     out<-modfun(PKindex$time[PKindex$Subject==i],par[1],par[2],par[3],par[4])
     gift<- which( PKindex$conc[PKindex$Subject==i] != 0 ) 
     sum(((PKindex$conc[PKindex$Subject==i][gift] - out[gift])/PKindex$conc[PKindex$Subject==i][gift])^2)}
     gen<-genoud(objfun,nvars=4,max=FALSE,pop.size=30,max.generations=20,wait.generations=10,
          starting.value=c(par[1,2],par[2,2],par[3,2],par[4,2]),BFGS=FALSE,print.level=0,boundary.enforcement=0,
          Domains=matrix(c(1,1,1,1,50,100,100,100),4,2),MemoryMatrix=TRUE,lexical=FALSE)
     cat("<< The value of parameter fitted by genetic algorithm >>\n\n")   
     namegen<-c("Tabs","Vm","Km","Vd")
     outgen<-c(gen$par[1],gen$par[2],gen$par[3],gen$par[4])
     print(data.frame(Parameter=namegen,Value=outgen))  
     F<-objfun(gen$par)
     opt<-optim(c(gen$par[1],gen$par[2],gen$par[3],gen$par[4]),objfun, method="Nelder-Mead")
     cat("\n<< The value of parameter fitted by Nelder-Mead Simplex slgorithm >>\n\n")   
     nameopt<-c("Tabs","Vm","Km","Vd")
     outopt<-c(opt$par[1],opt$par[2],opt$par[3],opt$par[4])
     print(data.frame(Parameter=nameopt,Value=outopt))
     cat("\n<< Residual sum-of-squares and parameter values fitted by nls >>\n\n")
     fm<-nls(conc~modfun(time,Tabs,Vm,Km,Vd),data=PKindex, algorithm="default",subset=Subject==i,
         start=list(Tabs=opt$par[1],Vm=opt$par[2],Km=opt$par[3],Vd=opt$par[4]),trace=TRUE,
         nls.control(maxiter=50,tol=10,minFactor=1/1024))
     cat("\n")
     plotting.non(PKindex, fm, i, pick,xaxis,yaxis)
     }
  }
  cat("\n")  
}

#simulation
#One compartment PK model extravascular single dose zero-order absorption without lag time
szero.nolag <- function()
{
  cat("How many subject do you want?\n")
  Subject<-scan(nlines=1,quiet=TRUE)
  PKtime<-data.frame(time=c(0))
  PKtime<-edit(PKtime) 
  cat("\n")
  show(PKtime)
  cat("\nEnter dose value\n")
  Dose<-scan(nlines=1,quiet=TRUE)
  cat("\n")

  par<-data.frame(Parameter=c("Tabs","kel","Vd"),Initial=c(0))
  par<-edit(par)
  cat("\n")
 
  par1<-par[1,2]
  par2<-par[2,2]    
  par3<-par[3,2]
    
  library(odesolve)
  defun<- function(time, y, parms) { 
     if(time<=parms["Tabs"]) 
        dCpdt <- (Dose/parms["Tabs"])/parms["Vd"] - parms["kel"] * y[1]
     else
        dCpdt <- - parms["kel"] * y[1]
     list(dCpdt) 
  } 
  file.menu <- c("no error",
                 "error=normal error", 
                 "error=uniform error",
                 "error=normal error*true value",
                 "error=uniform error*true value")
  pick <- menu(file.menu, title = "<< Error type >>")
  if (pick==1){
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" )         
     Tabs<-par1
     kel<-par2
     Vd<-par3     
     szero.out(PKtime,Tabs,kel,Vd,defun,par1,par2,par3,Dose,i)    
     }      
  }   
  else if (pick == 2){
     cat("\n\nEnter error factor for Tabs\n")
     factor1<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for kel\n")
     factor2<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for Vd\n")
     factor3<-scan(nlines=1,quiet=TRUE)
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" ) 
     Tabs<-par1+rnorm(1,mean=0,sd=factor1)
     while(Tabs<=0){
           Tabs<-par1+rnorm(1,mean=0,sd=factor1)}
     kel<-par2+rnorm(1,mean=0,sd=factor2)
     while(kel<=0){
           kel<-par2+rnorm(1,mean=0,sd=factor2)}
     Vd<-par3+rnorm(1,mean=0,sd=factor3)
     while(Vd<=0){
           Vd<-par3+rnorm(1,mean=0,sd=factor3)}
     szero.out(PKtime,Tabs,kel,Vd,defun,par1,par2,par3,Dose,i)         
     }       
  }
  else if (pick == 3){
     cat("\n\nEnter error factor for Tabs\n")
     factor1<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for kel\n")
     factor2<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for Vd\n")
     factor3<-scan(nlines=1,quiet=TRUE)
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" ) 
     Tabs<-par1+runif(1,min=-factor1,max=factor1)
     while(Tabs<=0){
           Tabs<-par1+runif(1,min=-factor1,max=factor1)}
     kel<-par2+runif(1,min=-factor2,max=factor2)
     while(kel<=0){
           kel<-par2+runif(1,min=-factor2,max=factor2)}  
     Vd<-par3+runif(1,min=-factor3,max=factor3)
     while(Vd<=0){
           Vd<-par3+runif(1,min=-factor3,max=factor3)}
     szero.out(PKtime,Tabs,kel,Vd,defun,par1,par2,par3,Dose,i)      
     }       
  }
  else if (pick == 4){
     cat("\n\nEnter error factor for Tabs\n")
     factor1<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for kel\n")
     factor2<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for Vd\n")
     factor3<-scan(nlines=1,quiet=TRUE)
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" ) 
     Tabs<-par1*rnorm(1,mean=0,sd=factor1)+par1
     while(Tabs<=0){
           Tabs<-par1*rnorm(1,mean=0,sd=factor1)+par1}
     kel<-par2*rnorm(1,mean=0,sd=factor2)+par2
     while(kel<=0){
           kel<-par2*rnorm(1,mean=0,sd=factor2)+par2}
     Vd<-par3*rnorm(1,mean=0,sd=factor3)+par3
     while(Vd<=0){
           Vd<-par3*rnorm(1,mean=0,sd=factor3)+par3}
     szero.out(PKtime,Tabs,kel,Vd,defun,par1,par2,par3,Dose,i)     
     }       
  }
  else if (pick == 5){
     cat("\n\nEnter error factor for Tabs\n")
     factor1<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for kel\n")
     factor2<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for Vd\n")
     factor3<-scan(nlines=1,quiet=TRUE)
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" ) 
     Tabs<-par1*runif(1,min=-factor1,max=factor1)+par1
     while(Tabs<=0){
           Tabs<-par1*runif(1,min=-factor1,max=factor1)+par1}
     kel<-par2*runif(1,min=-factor2,max=factor2)+par2
     while(kel<=0){
           kel<-par2*runif(1,min=-factor2,max=factor2)+par2}
     Vd<-par3*runif(1,min=-factor3,max=factor3)+par3
     while(Vd<=0){
           Vd<-par3*runif(1,min=-factor3,max=factor3)+par3}
     szero.out(PKtime,Tabs,kel,Vd,defun,par1,par2,par3,Dose,i)        
     }       
  }
}

#One compartment PK model extravascular single dose zero-order absorption without lag time nonlinear elimination
szero.nolagm <- function()
{
  cat("How many subject do you want?\n")
  Subject<-scan(nlines=1,quiet=TRUE)
  PKtime<-data.frame(time=c(0))
  PKtime<-edit(PKtime) 
  cat("\n")
  show(PKtime)
  cat("\nEnter dose value\n")
  Dose<-scan(nlines=1,quiet=TRUE)
  cat("\n")

  par<-data.frame(Parameter=c("Tabs","Vm","Km","Vd"),Initial=c(0))
  par<-edit(par)
  cat("\n")
 
  par1<-par[1,2]
  par2<-par[2,2]    
  par3<-par[3,2]
  par4<-par[4,2]   
    
  library(odesolve)
  defun<- function(time, y, parms) { 
     if(time<=parms["Tabs"]) 
       dCpdt <- (Dose/parms["Tabs"])/parms["Vd"]-(parms["Vm"]/parms["Vd"])*y[1]/(parms["Km"]/parms["Vd"]+y[1])            
     else
       dCpdt <- -(parms["Vm"]/parms["Vd"])*y[1]/(parms["Km"]/parms["Vd"]+y[1])            
     list(dCpdt) 
  }   

  file.menu <- c("no error",
                 "error=normal error", 
                 "error=uniform error",
                 "error=normal error*true value",
                 "error=uniform error*true value")
  pick <- menu(file.menu, title = "<< Error type >>")
  if (pick==1){
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" )         
     Tabs<-par1
     Vm<-par2
     Km<-par3
     Vd<-par4        
     szero.mm.out(PKtime,Tabs,Vm,Km,Vd,defun,par1,par2,par3,par4,Dose,i)
     }  
  }   
  else if (pick == 2){
     cat("\n\nEnter error factor for Tabs\n")
     factor1<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for Vm\n")
     factor2<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for Km\n")
     factor3<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for Vd\n")
     factor4<-scan(nlines=1,quiet=TRUE)
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" ) 
     Tabs<-par1+rnorm(1,mean=0,sd=factor1)
     while(Tabs<=0){
           Tabs<-par1+rnorm(1,mean=0,sd=factor1)}
     Vm<-par2+rnorm(1,mean=0,sd=factor2)
     while(Vm<=0){
           Vm<-par2+rnorm(1,mean=0,sd=factor2)}
     Km<-par3+rnorm(1,mean=0,sd=factor3)
     while(Km<=0){
           Km<-par3+rnorm(1,mean=0,sd=factor3)}
     Vd<-par4+rnorm(1,mean=0,sd=factor4)
     while(Vd<=0){
           Vd<-par4+rnorm(1,mean=0,sd=factor4)}
     szero.mm.out(PKtime,Tabs,Vm,Km,Vd,defun,par1,par2,par3,par4,Dose,i)     
     }       
  }
  else if (pick == 3){
     cat("\n\nEnter error factor for Tabs\n")
     factor1<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for Vm\n")
     factor2<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for Km\n")
     factor3<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for Vd\n")
     factor4<-scan(nlines=1,quiet=TRUE)
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" ) 
     Tabs<-par1+runif(1,min=-factor1,max=factor1)
     while(Tabs<=0){
           Tabs<-par1+runif(1,min=-factor1,max=factor1)}
     Vm<-par2+runif(1,min=-factor2,max=factor2)
     while(Vm<=0){
           Vm<-par2+runif(1,min=-factor2,max=factor2)}  
     Km<-par3+runif(1,min=-factor3,max=factor3)
     while(Km<=0){
           Km<-par3+runif(1,min=-factor3,max=factor3)}  
     Vd<-par4+runif(1,min=-factor4,max=factor4)
     while(Vd<=0){
           Vd<-par4+runif(1,min=-factor4,max=factor4)}
     szero.mm.out(PKtime,Tabs,Vm,Km,Vd,defun,par1,par2,par3,par4,Dose,i)          
     }       
  }
  else if (pick == 4){
     cat("\n\nEnter error factor for Tabs\n")
     factor1<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for Vm\n")
     factor2<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for Km\n")
     factor3<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for Vd\n")
     factor4<-scan(nlines=1,quiet=TRUE)
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" ) 
     Tabs<-par1*rnorm(1,mean=0,sd=factor1)+par1
     while(Tabs<=0){
           Tabs<-par1*rnorm(1,mean=0,sd=factor1)+par1}
     Vm<-par2*rnorm(1,mean=0,sd=factor2)+par2
     while(Vm<=0){
           Vm<-par2*rnorm(1,mean=0,sd=factor2)+par2}
     Km<-par3*rnorm(1,mean=0,sd=factor3)+par3
     while(Km<=0){
           Km<-par3*rnorm(1,mean=0,sd=factor3)+par3}
     Vd<-par4*rnorm(1,mean=0,sd=factor4)+par4
     while(Vd<=0){
           Vd<-par4*rnorm(1,mean=0,sd=factor4)+par4}
     szero.mm.out(PKtime,Tabs,Vm,Km,Vd,defun,par1,par2,par3,par4,Dose,i)       
     }       
  }
  else if (pick == 5){
     cat("\n\nEnter error factor for Tabs\n")
     factor1<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for Vm\n")
     factor2<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for Km\n")
     factor3<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for Vd\n")
     factor4<-scan(nlines=1,quiet=TRUE)
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" ) 
     Tabs<-par1*runif(1,min=-factor1,max=factor1)+par1
     while(Tabs<=0){
           Tabs<-par1*runif(1,min=-factor1,max=factor1)+par1}
     Vm<-par2*runif(1,min=-factor2,max=factor2)+par2
     while(Vm<=0){
           Vm<-par2*runif(1,min=-factor2,max=factor2)+par2}
     Km<-par3*runif(1,min=-factor3,max=factor3)+par3
     while(Km<=0){
           Km<-par3*runif(1,min=-factor3,max=factor3)+par3}
     Vd<-par4*runif(1,min=-factor4,max=factor4)+par4
     while(Vd<=0){
           Vd<-par4*runif(1,min=-factor4,max=factor4)+par4}
     szero.mm.out(PKtime,Tabs,Vm,Km,Vd,defun,par1,par2,par3,par4,Dose,i)         
     }       
  }
}   