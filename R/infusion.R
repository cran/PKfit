#Normal fitting
#One compartment PK model iv infusion single dose
finfu1 <- function()
{
  load("PK.RData")
  PKindex
  options(warn=-1)
    
 #lsoda is belong to odesolve package
  library(odesolve) 
    
 #BIC is belong to stats4 package
  library(stats4)
  
 #genoud is belong to rgenoud package 
  library(rgenoud)
        
 #Input dose and Tinf and initial value for kel and Vd
  cat("Enter Dose value\n")
  Dose<-scan(nlines=1,quiet=TRUE)

  cat("\nEnter infusion duration\n")
  Tinf<-scan(nlines=1,quiet=TRUE)
  cat("\n")
    
  par<-data.frame(Parameter=c("kel","Vd"),Initial=c(0))
  par<-edit(par)
  show(par)
  cat("\n")
    
 #User-supplied function
  defun<- function(time,y,parms) { 
     if(time<=Tinf)         
       dCpdt <- (Dose/Tinf)/parms["Vd"]-parms["kel"]*y[1]
     else
       dCpdt <- -parms["kel"]*y[1]
     list(dCpdt) 
  } 
   
  modfun <- function(time,kel,Vd) { 
     out<-lsoda(0,c(0,time),defun, parms=c(kel=kel,Vd=Vd),
                rtol=1e-5,atol=1e-5)
     out[-1,2] 
  } 
     
 #Select weighting schemes
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
     out <- modfun(PKindex$time[PKindex$Subject==i], par[1], par[2])
     gift<- which( PKindex$conc[PKindex$Subject==i] != 0 )
     sum((PKindex$conc[PKindex$Subject==i][gift]-out[gift])^2)}
     gen<-genoud(objfun,nvars=2,max=FALSE,pop.size=30,max.generations=20,wait.generations=10,
          starting.value=c(par[1,2],par[2,2]),BFGS=FALSE,print.level=0,boundary.enforcement=0,
          Domains=matrix(c(0.01,0.01,100,100),2,2),MemoryMatrix=TRUE,lexical=FALSE)       
     cat("<< The value of parameter fitted by genetic algorithm >>\n\n")   
     namegen<-c("kel","Vd")
     outgen<-c(gen$par[1],gen$par[2])
     print(data.frame(Parameter=namegen,Value=outgen))  
     F<-objfun(gen$par)         
     opt<-optim(c(gen$par[1],gen$par[2]),objfun,method="Nelder-Mead")  
     cat("\n<< The value of parameter fitted by Nelder-Mead Simplex slgorithm >>\n\n")   
     nameopt<-c("kel","Vd")
     outopt<-c(opt$par[1],opt$par[2])
     print(data.frame(Parameter=nameopt,Value=outopt))
     cat("\n<< Residual sum-of-squares and parameter values fitted by nls >>\n\n")
     fm<-nls(conc ~ modfun(time, kel, Vd), data=PKindex,
         start=list(kel=opt$par[1],Vd=opt$par[2]),trace=TRUE,subset=Subject==i,
         nls.control(tol=1))
     cat("\n")
     coef<-data.frame(coef(fm)["kel"])
     plotting.lin(PKindex, fm, i, pick, coef, xaxis, yaxis)
     }           
  }
  else if ( pick == 2 ){
     cat("\nEnter the title of x-axis(time)\n\n")
     xaxis<-readline()
     if (substr(xaxis, 1, 1) == "")  xaxis<-"Time"  else xaxis<-xaxis
     cat("\nEnter the title of y-axis(Cp)\n\n")
     yaxis<-readline()
     if (substr(yaxis, 1, 1) == "")  yaxis<-"Concentration"  else yaxis<-yaxis
     for( i in 1:length(unique(PKindex$Subject)))  {
     cat("\n\n               << Subject",i,">>\n\n" )  
     objfun <- function(par) { 
     out <- modfun(PKindex$time[PKindex$Subject==i], par[1], par[2])
     gift<- which( PKindex$conc[PKindex$Subject==i] != 0 )
     sum((PKindex$conc[PKindex$Subject==i][gift]-out[gift])^2/PKindex$conc[gift])}
     gen<-genoud(objfun,nvars=2,max=FALSE,pop.size=30,max.generations=20,wait.generations=10,
          starting.value=c(par[1,2],par[2,2]),BFGS=FALSE,print.level=0,boundary.enforcement=0,
          Domains=matrix(c(0.01,0.01,100,100),2,2),MemoryMatrix=TRUE,lexical=FALSE)       
     cat("<< The value of parameter fitted by genetic algorithm >>\n\n")   
     namegen<-c("kel","Vd")
     outgen<-c(gen$par[1],gen$par[2])
     print(data.frame(Parameter=namegen,Value=outgen))  
     F<-objfun(gen$par)         
     opt<-optim(c(gen$par[1],gen$par[2]),objfun,method="Nelder-Mead")  
     cat("\n<< The value of parameter fitted by Nelder-Mead Simplex slgorithm >>\n\n")   
     nameopt<-c("kel","Vd")
     outopt<-c(opt$par[1],opt$par[2])
     print(data.frame(Parameter=nameopt,Value=outopt))
     cat("\n<< Residual sum-of-squares and parameter values fitted by nls >>\n\n")
     fm<-nls(conc ~ modfun(time, kel, Vd), data=PKindex,
         start=list(kel=opt$par[1],Vd=opt$par[2]),trace=TRUE,subset=Subject==i,
         nls.control(tol=1))
     cat("\n")
     coef<-data.frame(coef(fm)["kel"])
     plotting.lin(PKindex, fm, i, pick, coef, xaxis, yaxis)
     }
  }           
  else if ( pick == 3 ){
     cat("\nEnter the title of x-axis(time)\n\n")
     xaxis<-readline()
     if (substr(xaxis, 1, 1) == "")  xaxis<-"Time"  else xaxis<-xaxis
     cat("\nEnter the title of y-axis(Cp)\n\n")
     yaxis<-readline()
     if (substr(yaxis, 1, 1) == "")  yaxis<-"Concentration"  else yaxis<-yaxis
     for( i in 1:length(unique(PKindex$Subject)))  {
     cat("\n\n               << Subject",i,">>\n\n" )  
     objfun <- function(par) { 
     out <- modfun(PKindex$time[PKindex$Subject==i], par[1], par[2])
     gift<- which( PKindex$conc[PKindex$Subject==i] != 0 ) 
     sum(((PKindex$conc[PKindex$Subject==i][gift] - out[gift])/PKindex$conc[gift])^2)}
     gen<-genoud(objfun,nvars=2,max=FALSE,pop.size=30,max.generations=20,wait.generations=10,
          starting.value=c(par[1,2],par[2,2]),BFGS=FALSE,print.level=0,boundary.enforcement=0,
          Domains=matrix(c(0.01,0.01,100,100),2,2),MemoryMatrix=TRUE,lexical=FALSE)       
     cat("<< The value of parameter fitted by genetic algorithm >>\n\n")   
     namegen<-c("kel","Vd")
     outgen<-c(gen$par[1],gen$par[2])
     print(data.frame(Parameter=namegen,Value=outgen))  
     F<-objfun(gen$par)         
     opt<-optim(c(gen$par[1],gen$par[2]),objfun,method="Nelder-Mead")  
     cat("\n<< The value of parameter fitted by Nelder-Mead Simplex slgorithm >>\n\n")   
     nameopt<-c("kel","Vd")
     outopt<-c(opt$par[1],opt$par[2])
     print(data.frame(Parameter=nameopt,Value=outopt))
     cat("\n<< Residual sum-of-squares and parameter values fitted by nls >>\n\n")
     fm<-nls(conc ~ modfun(time, kel, Vd), data=PKindex,
         start=list(kel=opt$par[1],Vd=opt$par[2]),trace=TRUE,subset=Subject==i,
         nls.control(tol=1))
     cat("\n")
     coef<-data.frame(coef(fm)["kel"])
     plotting.lin(PKindex, fm, i, pick, coef, xaxis, yaxis)
     }
  }
  cat("\n")   
}

#One compartment PK model iv infusion single dose michaelis-menten elimination
finfu.mm <- function()
{
  load("PK.Rdata")
  options(warn=-1)
   
 #lsoda is belong to odesolve package
  library(odesolve) 
   
 #BIC function is belong to stats4 package
  library(stats4)
  
 #genoud is belong to rgenoud package  
  library(rgenoud)
   
 #Input dose and Tinf and initial value for Vm, Km and Vd
  cat("Enter Dose value\n")
  Dose<-scan(nlines=1,quiet=TRUE)

  cat("\nEnter infusion duration\n")
  Tinf<-scan(nlines=1,quiet=TRUE)
  cat("\n")
   
  par<-data.frame(Parameter=c("Vm","Km","Vd"),Initial=c(0))
  par<-edit(par)
  show(par)
  cat("\n")
       
 #User-supplied function-----
  defun<-function(time, y, parms) { 
     if(time<=Tinf)  
       dCpdt <- (Dose/Tinf)/parms["Vd"]-(parms["Vm"]/parms["Vd"])*y[1]/(parms["Km"]/parms["Vd"]+y[1]) 
     else
       dCpdt <- -(parms["Vm"]/parms["Vd"])*y[1]/(parms["Km"]/parms["Vd"]+y[1])
     list(dCpdt)
  }
    
  modfun<-function(time,Vm,Km,Vd) { 
     out <- lsoda(0,c(0,time),defun,parms=c(Vm=Vm,Km=Km,Vd=Vd),rtol=1e-5,atol=1e-5)
     out[-1,2] 
  } 
    
 #Select weighting schemes-----
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
     out <- modfun(PKindex$time[PKindex$Subject==i], par[1], par[2],par[3])
     gift<- which( PKindex$conc[PKindex$Subject==i] != 0 )
     sum((PKindex$conc[PKindex$Subject==i][gift] - out[gift])^2)}
     gen<-genoud(objfun,nvars=3,max=FALSE,pop.size=30,max.generations=20,wait.generations=10,
          starting.value=c(par[1,2],par[2,2],par[3,2]),BFGS=FALSE,print.level=0,boundary.enforcement=0,
          Domains=matrix(c(0,0,0,100,100,100),3,2),MemoryMatrix=TRUE,lexical=FALSE) 
     cat("<< The value of parameter fitted by genetic algorithm >>\n\n")   
     namegen<-c("Vm","Km","Vd")
     outgen<-c(gen$par[1],gen$par[2],gen$par[3])
     print(data.frame(Parameter=namegen,Value=outgen)) 
     F<-objfun(gen$par)
     opt<-optim(c(gen$par[1],gen$par[2],gen$par[3]),objfun, method="Nelder-Mead")
     cat("\n<< The value of parameter fitted by Nelder-Mead Simplex slgorithm >>\n\n")   
     nameopt<-c("Vm","Km","Vd")
     outopt<-c(opt$par[1],opt$par[2],opt$par[3])
     print(data.frame(Parameter=nameopt,Value=outopt)) 
     cat("\n<< Residual sum-of-squares and parameter values fitted by nls >>\n\n")
     fm<-nls(conc~modfun(time,Vm,Km,Vd),data=PKindex,subset=Subject==i,
         start=list(Vm=opt$par[1],Km=opt$par[2],Vd=opt$par[3]),trace=TRUE,
         nls.control(tol=1))
     cat("\n")
     plotting.non(PKindex, fm, i, pick, xaxis, yaxis)
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
     out <- modfun(PKindex$time[PKindex$Subject==i], par[1], par[2],par[3])
     gift<- which( PKindex$conc[PKindex$Subject==i] != 0 )
     sum((PKindex$conc[PKindex$Subject==i][gift]-out[gift])^2/PKindex$conc[PKindex$Subject==i][gift])}
     gen<-genoud(objfun,nvars=3,max=FALSE,pop.size=30,max.generations=20,wait.generations=10,
          starting.value=c(par[1,2],par[2,2],par[3,2]),BFGS=FALSE,print.level=0,boundary.enforcement=0,
          Domains=matrix(c(0,0,0,100,100,100),3,2),MemoryMatrix=TRUE,lexical=FALSE) 
     cat("<< The value of parameter fitted by genetic algorithm >>\n\n")   
     namegen<-c("Vm","Km","Vd")
     outgen<-c(gen$par[1],gen$par[2],gen$par[3])
     print(data.frame(Parameter=namegen,Value=outgen)) 
     F<-objfun(gen$par)
     opt<-optim(c(gen$par[1],gen$par[2],gen$par[3]),objfun, method="Nelder-Mead")
     cat("\n<< The value of parameter fitted by Nelder-Mead Simplex slgorithm >>\n\n")   
     nameopt<-c("Vm","Km","Vd")
     outopt<-c(opt$par[1],opt$par[2],opt$par[3])
     print(data.frame(Parameter=nameopt,Value=outopt)) 
     cat("\n<< Residual sum-of-squares and parameter values fitted by nls >>\n\n")
     fm<-nls(conc~modfun(time,Vm,Km,Vd),data=PKindex,subset=Subject==i,
         start=list(Vm=opt$par[1],Km=opt$par[2],Vd=opt$par[3]),trace=TRUE,
         nls.control(tol=1))
     cat("\n")
     plotting.non(PKindex, fm, i, pick, xaxis, yaxis)
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
     out <- modfun(PKindex$time[PKindex$Subject==i], par[1], par[2],par[3])
     gift<- which( PKindex$conc[PKindex$Subject==i] != 0 ) 
     sum(((PKindex$conc[PKindex$Subject==i][gift] - out[gift])/PKindex$conc[PKindex$Subject==i][gift])^2)}
     gen<-genoud(objfun,nvars=3,max=FALSE,pop.size=30,max.generations=20,wait.generations=10,
          starting.value=c(par[1,2],par[2,2],par[3,2]),BFGS=FALSE,print.level=0,boundary.enforcement=0,
          Domains=matrix(c(0,0,0,100,100,100),3,2),MemoryMatrix=TRUE,lexical=FALSE) 
     cat("<< The value of parameter fitted by genetic algorithm >>\n\n")   
     namegen<-c("Vm","Km","Vd")
     outgen<-c(gen$par[1],gen$par[2],gen$par[3])
     print(data.frame(Parameter=namegen,Value=outgen)) 
     F<-objfun(gen$par)
     opt<-optim(c(gen$par[1],gen$par[2],gen$par[3]),objfun, method="Nelder-Mead")
     cat("\n<< The value of parameter fitted by Nelder-Mead Simplex slgorithm >>\n\n")   
     nameopt<-c("Vm","Km","Vd")
     outopt<-c(opt$par[1],opt$par[2],opt$par[3])
     print(data.frame(Parameter=nameopt,Value=outopt)) 
     cat("\n<< Residual sum-of-squares and parameter values fitted by nls >>\n\n")
     fm<-nls(conc~modfun(time,Vm,Km,Vd),data=PKindex,subset=Subject==i,
         start=list(Vm=opt$par[1],Km=opt$par[2],Vd=opt$par[3]),trace=TRUE,
         nls.control(tol=1))
     cat("\n")
     plotting.non(PKindex, fm, i, pick, xaxis, yaxis)
     }
  }
  cat("\n")
}    

#Simulation
#One compartment PK model iv infusion single dose
sinfu1 <- function()
{
  cat("How many subject do you want?\n")
  Subject<-scan(nlines=1,quiet=TRUE)
  PKindex<-data.frame(time=c(0))
  PKindex<-edit(PKindex) 
  cat("\n")
  show(PKindex)
  cat("\nEnter dose value\n")
  Dose<-scan(nlines=1,quiet=TRUE)
  cat("\nEnter value for infusion time\n")
  Tinf<-scan(nlines=1,quiet=TRUE)
  cat("\n")

  par<-data.frame(Parameter=c("kel","Vd"),Initial=c(0))
  par<-edit(par)
  cat("\n")
 
  par1<-par[1,2]
  par2<-par[2,2]
    
  library(odesolve)    
  defun<- function(time,y,parms) { 
     if(time<=Tinf)         
       dCpdt <- (Dose/Tinf)/parms["Vd"]-parms["kel"]*y[1]
      else
       dCpdt <- -parms["kel"]*y[1]
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
     kel<-par1
     Vd<-par2   
     sinfu1.out(PKindex,kel,Vd,defun,par1,par2,Dose,i) 
     }       
  }   
  else if (pick == 2){
     cat("\n\nEnter error factor for kel\n")
     factor1<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for Vd\n")
     factor2<-scan(nlines=1,quiet=TRUE)
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" )  
     kel<-par1+rnorm(1,mean=0,sd=factor1)
     while(kel<=0){
           kel<-par1+rnorm(1,mean=0,sd=factor1)}
     Vd<-par2+rnorm(1,mean=0,sd=factor2)
     while(Vd<=0){
           Vd<-par2+rnorm(1,mean=0,sd=factor2)}
     sinfu1.out(PKindex,kel,Vd,defun,par1,par2,Dose,i)       
     }      
  }
  else if (pick == 3){
     cat("\n\nEnter error factor for kel\n")
     factor1<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for Vd\n")
     factor2<-scan(nlines=1,quiet=TRUE)
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" )
     kel<-par1+runif(1,min=-factor1,max=factor1)
     while(kel<=0){
           kel<-par1+runif(1,min=-factor1,max=factor1)}
     Vd<-par2+runif(1,min=-factor2,max=factor2)
     while(Vd<=0){
           Vd<-par2+runif(1,min=-factor2,max=factor2)}
     sinfu1.out(PKindex,kel,Vd,defun,par1,par2,Dose,i)     
     }      
  }
  else if ( pick == 4 ){
     cat("\n\nEnter error factor for kel\n")
     factor1<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for Vd\n")
     factor2<-scan(nlines=1,quiet=TRUE)
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" )
     kel<-par1*rnorm(1,mean=0,sd=factor1)+par1
     while(kel<=0){
           kel<-par1*rnorm(1,mean=0,sd=factor1)+par1}
     Vd<-par2*rnorm(1,mean=0,sd=factor2)+par2
     while(Vd<=0){
           Vd<-par2*rnorm(1,mean=0,sd=factor2)+par2}
     sinfu1.out(PKindex,kel,Vd,defun,par1,par2,Dose,i)        
     }      
  }
  else if ( pick == 5 ){
     cat("\n\nEnter error factor for kel\n")
     factor1<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for Vd\n")
     factor2<-scan(nlines=1,quiet=TRUE)
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" )
     kel<-par1*runif(1,min=-factor1,max=factor1)+par1
     while(kel<=0){
           kel<-par1*runif(1,min=-factor1,max=factor1)+par1}
     Vd<-par2*runif(1,min=-factor2,max=factor2)+par2
     while(Vd<=0){
           Vd<-par2*runif(1,min=-factor2,max=factor2)+par2}
     sinfu1.out(PKindex,kel,Vd,defun,par1,par2,Dose,i)      
     }      
  }
}

#One compartment PK model iv infusion single dose michaelis-menten elimination
sinfu.mm <- function()
{
  cat("How many subject do you want?\n")
  Subject<-scan(nlines=1,quiet=TRUE)
  PKindex<-data.frame(time=c(0))
  PKindex<-edit(PKindex) 
  cat("\n")
  show(PKindex)
  cat("\nEnter dose value\n")
  Dose<-scan(nlines=1,quiet=TRUE)
  cat("\nEnter value for infusion time\n")
  Tinf<-scan(nlines=1,quiet=TRUE)
  cat("\n")

  par<-data.frame(Parameter=c("Vm","Km","Vd"),Initial=c(0))
  par<-edit(par)
  cat("\n")
 
  par1<-par[1,2]
  par2<-par[2,2]
  par3<-par[3,2]
   
  library(odesolve)
  defun<- function(time, y, parms) { 
     if(time<=Tinf)  
       dCpdt <- (Dose/Tinf)/parms["Vd"]-(parms["Vm"]/parms["Vd"])*y[1]/(parms["Km"]/parms["Vd"]+y[1]) 
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
     Vm<-par1
     Km<-par2
     Vd<-par3       
     sinfu.mm.out(PKindex,Vm,Km,Vd,defun,par1,par2,par3,Dose,i)
     }   
  }   
  else if (pick == 2){
     cat("\n\nEnter error factor for Vm\n")
     factor1<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for Km\n")
     factor2<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for Vd\n")
     factor3<-scan(nlines=1,quiet=TRUE)
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" ) 
     Vm<-par1+rnorm(1,mean=0,sd=factor1)
     while(Vm<=0){
           Vm<-par1+rnorm(1,mean=0,sd=factor1)}
     Km<-par2+rnorm(1,mean=0,sd=factor2)
     while(Km<=0){
           Km<-par2+rnorm(1,mean=0,sd=factor2)}
     Vd<-par3+rnorm(1,mean=0,sd=factor3)
     while(Vd<=0){
           Vd<-par3+rnorm(1,mean=0,sd=factor3)}
     sinfu.mm.out(PKindex,Vm,Km,Vd,defun,par1,par2,par3,Dose,i)     
     }       
  }
  else if (pick == 3){
     cat("\n\nEnter error factor for Vm\n")
     factor1<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for Km\n")
     factor2<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for Vd\n")
     factor3<-scan(nlines=1,quiet=TRUE)
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" ) 
     Vm<-par1+runif(1,min=-factor1,max=factor1)
     while(Vm<=0){
           Vm<-par1+runif(1,min=-factor1,max=factor1)}
     Km<-par2+runif(1,min=-factor2,max=factor2)
     while(Km<=0){
           Km<-par2+runif(1,min=-factor2,max=factor2)}
     Vd<-par3+runif(1,min=-factor3,max=factor3)
     while(Vd<=0){
           Vd<-par3+runif(1,min=-factor3,max=factor3)}
     sinfu.mm.out(PKindex,Vm,Km,Vd,defun,par1,par2,par3,Dose,i)     
     }       
  }
  else if (pick == 4){
     cat("\n\nEnter error factor for Vm\n")
     factor1<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for Km\n")
     factor2<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for Vd\n")
     factor3<-scan(nlines=1,quiet=TRUE)
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" ) 
     Vm<-par1*rnorm(1,mean=0,sd=factor1)+par1
     while(Vm<=0){
           Vm<-par1*rnorm(1,mean=0,sd=factor1)+par1}
     Km<-par2*rnorm(1,mean=0,sd=factor2)+par2
     while(Km<=0){
           Km<-par2*rnorm(1,mean=0,sd=factor2)+par2}
     Vd<-par3*rnorm(1,mean=0,sd=factor3)+par3
     while(Vd<=0){
           Vd<-par3*rnorm(1,mean=0,sd=factor3)+par3}
     sinfu.mm.out(PKindex,Vm,Km,Vd,defun,par1,par2,par3,Dose,i)      
     }       
  }
  else if (pick == 5){
     cat("\n\nEnter error factor for Vm\n")
     factor1<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for Km\n")
     factor2<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for Vd\n")
     factor3<-scan(nlines=1,quiet=TRUE)
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" ) 
     Vm<-par1*runif(1,min=-factor1,max=factor1)+par1
     while(Vm<=0){
           Vm<-par1*runif(1,min=-factor1,max=factor1)+par1}
     Km<-par2*runif(1,min=-factor2,max=factor2)+par2
     while(Km<=0){
           Km<-par2*runif(1,min=-factor2,max=factor2)+par2}
     Vd<-par3*runif(1,min=-factor3,max=factor3)+par3
     while(Vd<=0){
           Vd<-par3*runif(1,min=-factor3,max=factor3)+par3}
     sinfu.mm.out(PKindex,Vm,Km,Vd,defun,par1,par2,par3,Dose,i)       
     }       
  }
}      
