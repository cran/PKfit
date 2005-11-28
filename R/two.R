#Normal fitting
#Two compartment PK model iv bolus single dose
fbolus2 <- function(PKindex)
{
  options(warn=-1)
   
 #lsoda is belong to odesolve package
  library(odesolve)
    
 #BIC function is belong to stats4 package
  library(stats4)  
   
 #genoud function is belong to rgenoud package
  library(rgenoud)
   
 #Input dose and initial value for kel,k12,k21,Vd
  cat("Enter Dose value \n")
  Dose<-scan(nlines=1,quiet=TRUE)
  cat("\n")
   
  par<-data.frame(Parameter=c("kel","k12","k21","Vd"),Initial=c(0))
  par<-edit(par)
  show(par)
  cat("\n")
  
 #User-supplied function------
  defun<- function(time, y, parms) { 
     dCp1dt <- -parms["kel"]*y[1]-parms["k12"]*y[1]+parms["k21"]*y[2] 
     dCp2dt <-  parms["k12"]*y[1]-parms["k21"]*y[2]
     list(c(dCp1dt,dCp2dt)) 
  } 
    
  modfun <- function(time,kel,k12,k21, Vd) { 
     out <- lsoda(c(Dose/Vd,0),c(0,time),defun,parms=c(kel=kel,k12=k12,k21=k21,Vd=Vd),
                  rtol=1e-5,atol=1e-5) 
     out[-1,2] 
  }        
 
 #Select weighting schemes------ 
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
     out <- modfun(PKindex$time[PKindex$Subject==i], par[1], par[2], par[3], par[4])
     sum((PKindex$conc[PKindex$Subject==i] - out)^2)}
     gen<-genoud(objfun,nvars=4,max=FALSE,pop.size=30,max.generations=20,wait.generations=10,
          starting.value=c(par[1,2],par[2,2],par[3,2],par[4,2]),BFGS=FALSE,print.level=0,boundary.enforcement=2,
          Domains=matrix(c(0.01,0.01,0.01,1,1,10,1,100),4,2))    
     cat("<< The value of parameter fitted by genetic algorithm >>\n\n")   
     namegen<-c("kel","k12","k21","Vd")
     outgen<-c(gen$par[1],gen$par[2],gen$par[3],gen$par[4])
     print(data.frame(Parameter=namegen,Value=outgen))  
     F<-objfun(gen$par)
     opt <- optim(c(gen$par[1],gen$par[2],gen$par[3],gen$par[4]),objfun,method="Nelder-Mead")
     cat("\n<< The value of parameter fitted by Nelder-Mead Simplex slgorithm >>\n\n")   
     nameopt<-c("kel","k12","k21","Vd")
     outopt<-c(opt$par[1],opt$par[2],opt$par[3],opt$par[4])
     print(data.frame(Parameter=nameopt,Value=outopt))
     cat("\n<< Residual sum-of-squares and parameter values fitted by nls >>\n\n")
     fm<-nls(conc ~ modfun(time,kel,k12,k21,Vd), data=PKindex,subset=Subject==i,
         start=list(kel=opt$par[1],k12=opt$par[2],k21=opt$par[3],Vd=opt$par[4]),trace=TRUE,
         nls.control(maxiter=50,tol=1,minFactor=1/1024))
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
     out <- modfun(PKindex$time[PKindex$Subject==i], par[1], par[2], par[3], par[4])
     sum((PKindex$conc[PKindex$Subject==i] - out)^2/PKindex$conc[PKindex$Subject==i])}
     gen<-genoud(objfun,nvars=4,max=FALSE,pop.size=30,max.generations=20,wait.generations=10,
          starting.value=c(par[1,2],par[2,2],par[3,2],par[4,2]),BFGS=FALSE,print.level=0,boundary.enforcement=2,
          Domains=matrix(c(0.01,0.01,0.01,1,1,10,1,100),4,2))    
     cat("<< The value of parameter fitted by genetic algorithm >>\n\n")   
     namegen<-c("kel","k12","k21","Vd")
     outgen<-c(gen$par[1],gen$par[2],gen$par[3],gen$par[4])
     print(data.frame(Parameter=namegen,Value=outgen)) 
     F<-objfun(gen$par)
     opt <- optim(c(gen$par[1],gen$par[2],gen$par[3],gen$par[4]),objfun,method="Nelder-Mead")
     cat("\n<< The value of parameter fitted by Nelder-Mead Simplex slgorithm >>\n\n")   
     nameopt<-c("kel","k12","k21","Vd")
     outopt<-c(opt$par[1],opt$par[2],opt$par[3],opt$par[4])
     print(data.frame(Parameter=nameopt,Value=outopt))
     cat("\n<< Residual sum-of-squares and parameter values fitted by nls >>\n\n")
     fm<-nls(conc ~ modfun(time,kel,k12,k21,Vd), data=PKindex,subset=Subject==i,
         start=list(kel=opt$par[1],k12=opt$par[2],k21=opt$par[3],Vd=opt$par[4]),trace=TRUE,
         nls.control(maxiter=50,tol=1,minFactor=1/1024))
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
     out <- modfun(PKindex$time[PKindex$Subject==i], par[1], par[2], par[3], par[4])
     sum(((PKindex$conc[PKindex$Subject==i] - out)/(PKindex$conc[PKindex$Subject==i]))^2)}
     gen<-genoud(objfun,nvars=4,max=FALSE,pop.size=30,max.generations=20,wait.generations=10,
          starting.value=c(par[1,2],par[2,2],par[3,2],par[4,2]),BFGS=FALSE,print.level=0,boundary.enforcement=2,
          Domains=matrix(c(0.01,0.01,0.01,1,1,10,1,100),4,2))    
     cat("<< The value of parameter fitted by genetic algorithm >>\n\n")   
     namegen<-c("kel","k12","k21","Vd")
     outgen<-c(gen$par[1],gen$par[2],gen$par[3],gen$par[4])
     print(data.frame(Parameter=namegen,Value=outgen)) 
     F<-objfun(gen$par)
     opt <- optim(c(gen$par[1],gen$par[2],gen$par[3],gen$par[4]),objfun,method="Nelder-Mead")
     cat("\n<< The value of parameter fitted by Nelder-Mead Simplex slgorithm >>\n\n")   
     nameopt<-c("kel","k12","k21","Vd")
     outopt<-c(opt$par[1],opt$par[2],opt$par[3],opt$par[4])
     print(data.frame(Parameter=nameopt,Value=outopt))
     cat("\n<< Residual sum-of-squares and parameter values fitted by nls >>\n\n")
     fm<-nls(conc ~ modfun(time,kel,k12,k21,Vd), data=PKindex,subset=Subject==i,
         start=list(kel=opt$par[1],k12=opt$par[2],k21=opt$par[3],Vd=opt$par[4]),trace=TRUE,
         nls.control(maxiter=50,tol=1,minFactor=1/1024))
     cat("\n")
     coef<-data.frame(coef(fm)["kel"])     
     plotting.lin(PKindex, fm, i, pick, coef, xaxis, yaxis)
     }
  }
  cat("\n")
}
    
#Two compartment PK model iv infusion single dose    
finfu2 <- function(PKindex)
{
  options(warn=-1)
    
 #lsoda is belong to odesolve package
  library(odesolve) 
    
 #BIC is belong to stats4 package
  library(stats4)
    
 #genoud function is belong to rgenoud package
  library(rgenoud)
    
 #Input dose and Tinf and initial value for kel,k12,k21,Vd
  cat("Enter Dose value\n")
  Dose<-scan(nlines=1,quiet=TRUE)

  cat("\nEnter infusion duration\n")
  Tinf<-scan(nlines=1,quiet=TRUE)
  cat("\n")
    
  par<-data.frame(Parameter=c("kel","k12","k21","Vd"),Initial=c(0))
  par<-edit(par)
  show(par)
  cat("\n")
    
 #User-supplied function
  defun<- function(time, y, parms) { 
  if(time<=Tinf) {
     dCp1dt<- (Dose/Tinf)/parms["Vd"]+parms["k21"]*y[2]-parms["kel"]*y[1]-parms["k12"]*y[1]
     dCp2dt<- parms["k12"]*y[1]-parms["k21"]*y[2]
  }
  else {       
     dCp1dt <- -parms["kel"]*y[1]-parms["k12"]*y[1]+parms["k21"]*y[2] 
     dCp2dt <- parms["k12"]*y[1]-parms["k21"]*y[2]
  }
     list(c(dCp1dt,dCp2dt)) 
  } 
      
  modfun <- function(time,kel,k12,k21,Vd) { 
     out <- lsoda(c(0,0),c(0,time),defun,parms=c(kel=kel,k12=k12,k21=k21,Vd=Vd),
                  rtol=1e-8,atol=1e-8) 
     out[-1,2] 
  }
     
 #Select weighting schemes-----
  file.menu <- c("equal weight", 
                 "1/Cp",
                 "1/Cp^2")           
  pick <- menu(file.menu, title = "<< Weighting Schemes >>")
  if (pick ==1 ){
     cat("\nEnter the title of x-axis(time)\n\n")
     xaxis<-readline()
     if (substr(xaxis, 1, 1) == "")  xaxis<-"Time"  else xaxis<-xaxis
     cat("\nEnter the title of y-axis(Cp)\n\n")
     yaxis<-readline()
     if (substr(yaxis, 1, 1) == "")  yaxis<-"Concentration"  else yaxis<-yaxis
     for( i in 1:length(unique(PKindex$Subject)))  {
     cat("\n\n               << Subject",i,">>\n\n" )  
     objfun <- function(par) { 
     out <- modfun(PKindex$time[PKindex$Subject==i], par[1], par[2], par[3], par[4])
     gift<- which( PKindex$conc[PKindex$Subject==i] != 0 )
     sum((PKindex$conc[PKindex$Subject==i][gift]-out[gift])^2)}
     gen<-genoud(objfun,nvars=4,max=FALSE,pop.size=30,max.generations=20,wait.generations=10,
          starting.value=c(par[1,2],par[2,2],par[3,2],par[4,2]),BFGS=FALSE,print.level=0,boundary.enforcement=0,
          Domains=matrix(c(0.01,0.01,0.01,1,1,10,1,100),4,2),MemoryMatrix=TRUE,lexical=FALSE)    
     cat("<< The value of parameter fitted by genetic algorithm >>\n\n")   
     namegen<-c("kel","k12","k21","Vd")
     outgen<-c(gen$par[1],gen$par[2],gen$par[3],gen$par[4])
     print(data.frame(Parameter=namegen,Value=outgen)) 
     F<-objfun(gen$par)
     opt <- optim(c(gen$par[1],gen$par[2],gen$par[3],gen$par[4]),objfun,method="Nelder-Mead")
     cat("\n<< The value of parameter fitted by Nelder-Mead Simplex slgorithm >>\n\n")   
     nameopt<-c("kel","k12","k21","Vd")
     outopt<-c(opt$par[1],opt$par[2],opt$par[3],opt$par[4])
     print(data.frame(Parameter=nameopt,Value=outopt))
     cat("\n<< Residual sum-of-squares and parameter values fitted by nls >>\n\n")
     fm<-nls(conc ~ modfun(time,kel,k12,k21,Vd), data=PKindex,subset=Subject==i,
         start=list(kel=opt$par[1],k12=opt$par[2],k21=opt$par[3],Vd=opt$par[4]),trace=TRUE,
         nls.control(maxiter=50,tol=1,minFactor=1/1024))
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
     out <- modfun(PKindex$time[PKindex$Subject==i], par[1], par[2], par[3], par[4])
     gift<- which( PKindex$conc[PKindex$Subject==i] != 0 )
     sum((PKindex$conc[PKindex$Subject==i][gift]-out[gift])^2/PKindex$conc[PKindex$Subject==i][gift])}
     gen<-genoud(objfun,nvars=4,max=FALSE,pop.size=30,max.generations=20,wait.generations=10,
          starting.value=c(par[1,2],par[2,2],par[3,2],par[4,2]),BFGS=FALSE,print.level=0,boundary.enforcement=0,
          Domains=matrix(c(0.01,0.01,0.01,1,1,10,1,100),4,2),MemoryMatrix=TRUE,lexical=FALSE)    
     cat("<< The value of parameter fitted by genetic algorithm >>\n\n")   
     namegen<-c("kel","k12","k21","Vd")
     outgen<-c(gen$par[1],gen$par[2],gen$par[3],gen$par[4])
     print(data.frame(Parameter=namegen,Value=outgen)) 
     F<-objfun(gen$par)
     opt <- optim(c(gen$par[1],gen$par[2],gen$par[3],gen$par[4]),objfun,method="Nelder-Mead")
     cat("\n<< The value of parameter fitted by Nelder-Mead Simplex slgorithm >>\n\n")   
     nameopt<-c("kel","k12","k21","Vd")
     outopt<-c(opt$par[1],opt$par[2],opt$par[3],opt$par[4])
     print(data.frame(Parameter=nameopt,Value=outopt))
     cat("\n<< Residual sum-of-squares and parameter values fitted by nls >>\n\n")
     fm<-nls(conc ~ modfun(time,kel,k12,k21,Vd), data=PKindex,subset=Subject==i,
         start=list(kel=opt$par[1],k12=opt$par[2],k21=opt$par[3],Vd=opt$par[4]),trace=TRUE,
         nls.control(maxiter=50,tol=1,minFactor=1/1024))
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
     out <- modfun(PKindex$time[PKindex$Subject==i], par[1], par[2], par[3], par[4])
     gift<- which( PKindex$conc[PKindex$Subject==i] != 0 ) 
     sum(((PKindex$conc[PKindex$Subject==i][gift] - out[gift])/PKindex$conc[PKindex$Subject==i][gift])^2)}
     gen<-genoud(objfun,nvars=4,max=FALSE,pop.size=30,max.generations=20,wait.generations=10,
          starting.value=c(par[1,2],par[2,2],par[3,2],par[4,2]),BFGS=FALSE,print.level=0,boundary.enforcement=0,
          Domains=matrix(c(0.01,0.01,0.01,1,1,10,1,100),4,2),MemoryMatrix=TRUE,lexical=FALSE)    
     cat("<< The value of parameter fitted by genetic algorithm >>\n\n")   
     namegen<-c("kel","k12","k21","Vd")
     outgen<-c(gen$par[1],gen$par[2],gen$par[3],gen$par[4])
     print(data.frame(Parameter=namegen,Value=outgen)) 
     F<-objfun(gen$par)
     opt <- optim(c(gen$par[1],gen$par[2],gen$par[3],gen$par[4]),objfun,method="Nelder-Mead")
     cat("\n<< The value of parameter fitted by Nelder-Mead Simplex slgorithm >>\n\n")   
     nameopt<-c("kel","k12","k21","Vd")
     outopt<-c(opt$par[1],opt$par[2],opt$par[3],opt$par[4])
     print(data.frame(Parameter=nameopt,Value=outopt))
     cat("\n<< Residual sum-of-squares and parameter values fitted by nls >>\n\n")
     fm<-nls(conc ~ modfun(time,kel,k12,k21,Vd), data=PKindex,subset=Subject==i,
         start=list(kel=opt$par[1],k12=opt$par[2],k21=opt$par[3],Vd=opt$par[4]),trace=TRUE,
         nls.control(maxiter=50,tol=1,minFactor=1/1024))
     cat("\n")
     coef<-data.frame(coef(fm)["kel"])     
     plotting.lin(PKindex, fm, i, pick, coef, xaxis, yaxis)
     }
  }
  cat("\n")    
}

#Two compartment PK model extravascular single dose first order absorption 
ffirst2 <- function(PKindex)
{
  options(warn=-1)
    
 #lsoda is belong to odesolve package-----
  library(odesolve)
    
 #BIC function is belong to stats4 package-----
  library(stats4)
    
 #genoud function is belong to rgenoud package
  library(rgenoud)
       
 #Input dose and initial value for ka, kel, k12, k21 and  Vd-----
  cat("Enter Dose value\n")
  Dose<-scan(nlines=1,quiet=TRUE)
  cat("\n")
    
  par<-data.frame(Parameter=c("ka","kel","k12","k21","Vd"),Initial=c(0))
  par<-edit(par)
  show(par)
  cat("\n")
    
 #User-supplied function----
  defun<- function(time, y, parms) { 
     dCp1dt <- -parms["ka"]*y[1]
     dCp2dt <- parms["ka"]*y[1]/parms["Vd"]-parms["kel"]*y[2]+parms["k21"]*y[3]-parms["k12"]*y[2]
     dCp3dt <- parms["k12"]*y[2]-parms["k21"]*y[3]
     list(c(dCp1dt,dCp2dt,dCp3dt)) 
  } 
    
  modfun <- function(time,ka,kel,k12,k21,Vd) { 
     out <- lsoda(c(Dose,0,0),c(0,time),defun,parms=c(ka=ka,kel=kel,k12=k12,k21=k21,Vd=Vd),
                  rtol=1e-6,atol=1e-8) 
     out[-1,3] 
  }
     
 #selection of weighting schemes
  file.menu <- c("equal weight", 
                 "1/Cp",
                 "1/Cp^2")           
  pick <- menu(file.menu, title = "<< Weighting Schemes >>")
  if (pick ==1 ){
     cat("\nEnter the title of x-axis(time)\n\n")
     xaxis<-readline()
     if (substr(xaxis, 1, 1) == "")  xaxis<-"Time"  else xaxis<-xaxis
     cat("\nEnter the title of y-axis(Cp)\n\n")
     yaxis<-readline()
     if (substr(yaxis, 1, 1) == "")  yaxis<-"Concentration"  else yaxis<-yaxis
     for( i in 1:length(unique(PKindex$Subject)))  {
     cat("\n\n               << Subject",i,">>\n\n" )  
     objfun <- function(par) { 
     out <- modfun(PKindex$time[PKindex$Subject==i],par[1],par[2],par[3],par[4],par[5])
     gift<- which( PKindex$conc[PKindex$Subject==i] != 0 )
     sum((PKindex$conc[PKindex$Subject==i][gift]-out[gift])^2)}
     gen<-genoud(objfun,nvars=5,max=FALSE,pop.size=30,max.generations=20,wait.generations=10,
          starting.value=c(par[1,2],par[2,2],par[3,2],par[4,2],par[5,2]),BFGS=FALSE,print.level=0,boundary.enforcement=0,
          Domains=matrix(c(0.01,0.01,0.01,0.01,1,10,1,10,1,100),5,2),MemoryMatrix=TRUE,lexical=FALSE) 
     cat("<< The value of parameter fitted by genetic algorithm >>\n\n")   
     namegen<-c("ka","kel","k12","k21","Vd")
     outgen<-c(gen$par[1],gen$par[2],gen$par[3],gen$par[4],gen$par[5])
     print(data.frame(Parameter=namegen,Value=outgen)) 
     F<-objfun(gen$par)
     opt<-optim(c(gen$par[1],gen$par[2],gen$par[3],gen$par[4],gen$par[5]),objfun, method="Nelder-Mead")
     cat("\n<< The value of parameter fitted by Nelder-Mead Simplex slgorithm >>\n\n")   
     nameopt<-c("ka","kel","k12","k21","Vd")
     outopt<-c(opt$par[1],opt$par[2],opt$par[3],opt$par[4],opt$par[5])
     print(data.frame(Parameter=nameopt,Value=outopt))
     cat("\n<< Residual sun-of-squares and parameter value of nls >>\n\n")
     fm<-nls(conc ~ modfun(time,ka,kel,k12,k21,Vd), data=PKindex,subset=Subject==i,
         start=list(ka=opt$par[1],kel=opt$par[2],k12=opt$par[3],k21=opt$par[4],Vd=opt$par[5]),trace=TRUE,
         nls.control(maxiter=50,tol=1e-2,minFactor=1/1024))
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
     out <- modfun(PKindex$time[PKindex$Subject==i],par[1],par[2],par[3],par[4],par[5])
     gift<- which( PKindex$conc[PKindex$Subject==i] != 0 )
     sum((PKindex$conc[PKindex$Subject==i][gift]-out[gift])^2/PKindex$conc[PKindex$Subject==i][gift])}
     gen<-genoud(objfun,nvars=5,max=FALSE,pop.size=30,max.generations=20,wait.generations=10,
          starting.value=c(par[1,2],par[2,2],par[3,2],par[4,2],par[5,2]),BFGS=FALSE,print.level=0,boundary.enforcement=0,
          Domains=matrix(c(0.01,0.01,0.01,0.01,1,10,1,10,1,100),5,2),MemoryMatrix=TRUE,lexical=FALSE) 
     cat("<< The value of parameter fitted by genetic algorithm >>\n\n")   
     namegen<-c("ka","kel","k12","k21","Vd")
     outgen<-c(gen$par[1],gen$par[2],gen$par[3],gen$par[4],gen$par[5])
     print(data.frame(Parameter=namegen,Value=outgen))
     F<-objfun(gen$par)
     opt<-optim(c(gen$par[1],gen$par[2],gen$par[3],gen$par[4],gen$par[5]),objfun, method="Nelder-Mead")
     cat("\n<< The value of parameter fitted by Nelder-Mead Simplex slgorithm >>\n\n")   
     nameopt<-c("ka","kel","k12","k21","Vd")
     outopt<-c(opt$par[1],opt$par[2],opt$par[3],opt$par[4],opt$par[5])
     print(data.frame(Parameter=nameopt,Value=outopt))
     cat("\n<< Residual sun-of-squares and parameter value of nls >>\n\n")
     fm<-nls(conc ~ modfun(time,ka,kel,k12,k21,Vd), data=PKindex,subset=Subject==i,
         start=list(ka=opt$par[1],kel=opt$par[2],k12=opt$par[3],k21=opt$par[4],Vd=opt$par[5]),trace=TRUE,
         nls.control(maxiter=50,tol=1e-2,minFactor=1/1024))
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
     out <- modfun(PKindex$time[PKindex$Subject==i],par[1],par[2],par[3],par[4],par[5])
     gift<- which( PKindex$conc[PKindex$Subject==i] != 0 ) 
     sum(((PKindex$conc[PKindex$Subject==i][gift]-out[gift])/PKindex$conc[PKindex$Subject==i][gift])^2)}
     gen<-genoud(objfun,nvars=5,max=FALSE,pop.size=30,max.generations=20,wait.generations=10,
          starting.value=c(par[1,2],par[2,2],par[3,2],par[4,2],par[5,2]),BFGS=FALSE,print.level=0,boundary.enforcement=0,
          Domains=matrix(c(0.01,0.01,0.01,0.01,1,10,1,10,1,100),5,2),MemoryMatrix=TRUE,lexical=FALSE) 
     cat("<< The value of parameter fitted by genetic algorithm >>\n\n")   
     namegen<-c("ka","kel","k12","k21","Vd")
     outgen<-c(gen$par[1],gen$par[2],gen$par[3],gen$par[4],gen$par[5])
     print(data.frame(Parameter=namegen,Value=outgen)) 
     F<-objfun(gen$par)
     opt<-optim(c(gen$par[1],gen$par[2],gen$par[3],gen$par[4],gen$par[5]),objfun, method="Nelder-Mead")
     cat("\n<< The value of parameter fitted by Nelder-Mead Simplex slgorithm >>\n\n")   
     nameopt<-c("ka","kel","k12","k21","Vd")
     outopt<-c(opt$par[1],opt$par[2],opt$par[3],opt$par[4],opt$par[5])
     print(data.frame(Parameter=nameopt,Value=outopt))
     cat("\n<< Residual sun-of-squares and parameter value of nls >>\n\n")
     fm<-nls(conc ~ modfun(time,ka,kel,k12,k21,Vd), data=PKindex,subset=Subject==i,
         start=list(ka=opt$par[1],kel=opt$par[2],k12=opt$par[3],k21=opt$par[4],Vd=opt$par[5]),trace=TRUE,
         nls.control(maxiter=50,tol=1e-2,minFactor=1/1024))
     cat("\n")
     coef<-data.frame(coef(fm)["kel"])     
     plotting.lin(PKindex, fm, i, pick, coef, xaxis, yaxis)
     }
  }
   cat("\n")
}

#Simulation
#Two compartment PK model iv bolus single dose
sbolus2 <- function()
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

  par<-data.frame(Parameter=c("kel","k12","k21","Vd"),Initial=c(0))
  par<-edit(par)
  cat("\n")
  
  par1<-par[1,2]
  par2<-par[2,2]
  par3<-par[3,2]
  par4<-par[4,2]

  library(odesolve)
  defun<- function(time, y, parms) { 
     dCp1dt <- -parms["kel"]*y[1]-parms["k12"]*y[1]+parms["k21"]*y[2] 
     dCp2dt <-  parms["k12"]*y[1]-parms["k21"]*y[2]
     list(c(dCp1dt,dCp2dt)) 
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
     cat("\n\n") 
     kel<-par1
     k12<-par2
     k21<-par3
     Vd<-par4   
     sbolus2.out(PKtime,kel,k12,k21,Vd,defun,par1,par2,par3,par4,Dose,i)
     }   
  }   
  else if (pick == 2){
     cat("\n\nEnter error factor for kel\n")
     factor1<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for k12\n")
     factor2<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for k21\n")
     factor3<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for Vd\n")
     factor4<-scan(nlines=1,quiet=TRUE)
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" )  
     kel<-par1+rnorm(1,mean=0,sd=factor1)
     while(kel<=0){
           kel<-par1+rnorm(1,mean=0,sd=factor1)} 
     k12<-par2+rnorm(1,mean=0,sd=factor2)
     while(k12<=0){
           k12<-par2+rnorm(1,mean=0,sd=factor2)}
     k21<-par3+rnorm(1,mean=0,sd=factor3)
     while(k21<=0){
           k21<-par3+rnorm(1,mean=0,sd=factor3)}        
     Vd<-par4+rnorm(1,mean=0,sd=factor4)
     while(Vd<=0){
           Vd<-par4+rnorm(1,mean=0,sd=factor4)}
     sbolus2.out(PKtime,kel,k12,k21,Vd,defun,par1,par2,par3,par4,Dose,i)
     }       
  }
  else if (pick == 3){
     cat("\n\nEnter error factor for kel\n")
     factor1<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for k12\n")
     factor2<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for k21\n")
     factor3<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for Vd\n")
     factor4<-scan(nlines=1,quiet=TRUE)
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" )  
     kel<-par1+runif(1,min=-factor1,max=factor1)
     while(kel<=0){
           kel<-par1+runif(1,min=-factor1,max=factor1)}
     k12<-par2+runif(1,min=-factor2,max=factor2)
     while(k12<=0){
           k12<-par2+runif(1,min=-factor2,max=factor2)}
     k21<-par3+runif(1,min=-factor3,max=factor3)
     while(k21<=0){
           k21<-par3+runif(1,min=-factor3,max=factor3)}
     Vd<-par4+runif(1,min=-factor4,max=factor4)
     while(Vd<=0){
           Vd<-par4+runif(1,min=-factor4,max=factor4)}
     sbolus2.out(PKtime,kel,k12,k21,Vd,defun,par1,par2,par3,par4,Dose,i)      
     }      
  }
  else if (pick == 4){
     cat("\n\nEnter error factor for kel\n")
     factor1<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for k12\n")
     factor2<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for k21\n")
     factor3<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for Vd\n")
     factor4<-scan(nlines=1,quiet=TRUE)
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" )  
     kel<-par1*rnorm(1,mean=0,sd=factor1)+par1
     while(kel<=0){
           kel<-par1*rnorm(1,mean=0,sd=factor1)+par1}
     k12<-par2*rnorm(1,mean=0,sd=factor2)+par2
     while(k12<=0){
           k12<-par2*rnorm(1,mean=0,sd=factor2)+par2}
     k21<-par3*rnorm(1,mean=0,sd=factor3)+par3
     while(k21<=0){
           k21<-par3*rnorm(1,mean=0,sd=factor3)+par3}
     Vd<-par4*rnorm(1,mean=0,sd=factor4)+par4
     while(Vd<=0){
           Vd<-par4*rnorm(1,mean=0,sd=factor4)+par4}
     sbolus2.out(PKtime,kel,k12,k21,Vd,defun,par1,par2,par3,par4,Dose,i)      
     }       
  }
  else if (pick == 5){
     cat("\n\nEnter error factor for kel\n")
     factor1<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for k12\n")
     factor2<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for k21\n")
     factor3<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for Vd\n")
     factor4<-scan(nlines=1,quiet=TRUE)
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" ) 
     kel<-par1*runif(1,min=-factor1,max=factor1)+par1
     while(kel<=0){
           kel<-par1*runif(1,min=-factor1,max=factor1)+par1}
     k12<-par2*runif(1,min=-factor2,max=factor2)+par2
     while(k12<=0){
           k12<-par2*runif(1,min=-factor2,max=factor2)+par2}
     k21<-par3*runif(1,min=-factor3,max=factor3)+par3
     while(k21<=0){
           k21<-par3*runif(1,min=-factor3,max=factor3)+par3}
     Vd<-par4*runif(1,min=-factor4,max=factor4)+par4
     while(Vd<=0){
           Vd<-par4*runif(1,min=-factor4,max=factor4)+par4}
     sbolus2.out(PKtime,kel,k12,k21,Vd,defun,par1,par2,par3,par4,Dose,i)         
     }       
  }
}      

#Two compartment PK model iv infusion single dose  
sinfu2 <- function()
{
  cat("How many subject do you want?\n")
  Subject<-scan(nlines=1,quiet=TRUE)
  PKtime<-data.frame(time=c(0))
  PKtime<-edit(PKtime) 
  cat("\n")
  show(PKtime)  
  cat("\nEnter dose value\n")
  Dose<-scan(nlines=1,quiet=TRUE)
  cat("\nEnter value for infusion time\n")
  Tinf<-scan(nlines=1,quiet=TRUE)
  cat("\n")
  
  par<-data.frame(Parameter=c("kel","k12","k21","Vd"),Initial=c(0))
  par<-edit(par)
  cat("\n")

  par1<-par[1,2]
  par2<-par[2,2]
  par3<-par[3,2]
  par4<-par[4,2]

  library(odesolve)
  defun<- function(time, y, parms) { 
       if(time<=Tinf) {
           dCp1dt<- (Dose/Tinf)/parms["Vd"]+parms["k21"]*y[2]-parms["kel"]*y[1]-parms["k12"]*y[1]
           dCp2dt<- parms["k12"]*y[1]-parms["k21"]*y[2]
         }
       else {       
           dCp1dt <- -parms["kel"]*y[1]-parms["k12"]*y[1]+parms["k21"]*y[2] 
           dCp2dt <- parms["k12"]*y[1]-parms["k21"]*y[2]
         }
       list(c(dCp1dt,dCp2dt)) 
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
     k12<-par2
     k21<-par3
     Vd<-par4   
     sinfu2.out(PKtime,kel,k12,k21,Vd,defun,par1,par2,par3,par4,Dose,i)  
     }   
  }   
  else if (pick == 2){
     cat("\n\nEnter error factor for kel\n")
     factor1<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for k12\n")
     factor2<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for k21\n")
     factor3<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for Vd\n")
     factor4<-scan(nlines=1,quiet=TRUE)
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" ) 
     kel<-par1+rnorm(1,mean=0,sd=factor1)
     while(kel<=0){
           kel<-par1+rnorm(1,mean=0,sd=factor1)} 
     k12<-par2+rnorm(1,mean=0,sd=factor2)
     while(k12<=0){
           k12<-par2+rnorm(1,mean=0,sd=factor2)}
     k21<-par3+rnorm(1,mean=0,sd=factor3)
     while(k21<=0){
           k21<-par3+rnorm(1,mean=0,sd=factor3)}
     Vd<-par4+rnorm(1,mean=0,sd=factor4)
     while(Vd<=0){
           Vd<-par4+rnorm(1,mean=0,sd=factor4)}
     sinfu2.out(PKtime,kel,k12,k21,Vd,defun,par1,par2,par3,par4,Dose,i)      
     }       
  }
  else if (pick == 3){
     cat("\n\nEnter error factor for kel\n")
     factor1<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for k12\n")
     factor2<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for k21\n")
     factor3<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for Vd\n")
     factor4<-scan(nlines=1,quiet=TRUE)
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" ) 
     kel<-par1+runif(1,min=-factor1,max=factor1)
     while(kel<=0){
           kel<-par1+runif(1,min=-factor1,max=factor1)}
     k12<-par2+runif(1,min=-factor2,max=factor2)
     while(k12<=0){
           k12<-par2+runif(1,min=-factor2,max=factor2)}
     k21<-par3+runif(1,min=-factor3,max=factor3)
     while(k21<=0){
           k21<-par3+runif(1,min=-factor3,max=factor3)}
     Vd<-par4+runif(1,min=-factor4,max=factor4)
     while(Vd<=0){
           Vd<-par4+runif(1,min=-factor4,max=factor4)}
     sinfu2.out(PKtime,kel,k12,k21,Vd,defun,par1,par2,par3,par4,Dose,i)       
     }       
  }
  else if (pick == 4){
     cat("\n\nEnter error factor for kel\n")
     factor1<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for k12\n")
     factor2<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for k21\n")
     factor3<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for Vd\n")
     factor4<-scan(nlines=1,quiet=TRUE)
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" ) 
     kel<-par1*rnorm(1,mean=0,sd=factor1)+par1
     while(kel<=0){
           kel<-par1*rnorm(1,mean=0,sd=factor1)+par1}
     k12<-par2*rnorm(1,mean=0,sd=factor2)+par2
     while(k12<=0){
           k12<-par2*rnorm(1,mean=0,sd=factor2)+par2}
     k21<-par3*rnorm(1,mean=0,sd=factor3)+par3
     while(k21<=0){
           k21<-par3*rnorm(1,mean=0,sd=factor3)+par3}
     Vd<-par4*rnorm(1,mean=0,sd=factor4)+par4
     while(Vd<=0){
           Vd<-par4*rnorm(1,mean=0,sd=factor4)+par4}
     sinfu2.out(PKtime,kel,k12,k21,Vd,defun,par1,par2,par3,par4,Dose,i)        
     }       
  }
  else if (pick == 5){
     cat("\n\nEnter error factor for kel\n")
     factor1<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for k12\n")
     factor2<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for k21\n")
     factor3<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for Vd\n")
     factor4<-scan(nlines=1,quiet=TRUE)
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" ) 
     kel<-par1*runif(1,min=-factor1,max=factor1)+par1
     while(kel<=0){
           kel<-par1*runif(1,min=-factor1,max=factor1)+par1}
     k12<-par2*runif(1,min=-factor2,max=factor2)+par2
     while(k12<=0){
           k12<-par2*runif(1,min=-factor2,max=factor2)+par2}
     k21<-par3*runif(1,min=-factor3,max=factor3)+par3
     while(k21<=0){
           k21<-par3*runif(1,min=-factor3,max=factor3)+par3}
     Vd<-par4*runif(1,min=-factor4,max=factor4)+par4
     while(Vd<=0){
           Vd<-par4*runif(1,min=-factor4,max=factor4)+par4}
     sinfu2.out(PKtime,kel,k12,k21,Vd,defun,par1,par2,par3,par4,Dose,i)          
     }       
  }
}

#Two compartment PK model extravascular single dose first order absorption
sfirst2 <- function()
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
   
  par<-data.frame(Parameter=c("ka","kel","k12","k21","Vd"),Initial=c(0))
  par<-edit(par)
  cat("\n")

  par1<-par[1,2]
  par2<-par[2,2]
  par3<-par[3,2]
  par4<-par[4,2]
  par5<-par[5,2]

  library(odesolve)
  defun<- function(time, y, parms) { 
       dCp1dt <- -parms["ka"]*y[1]
       dCp2dt <- parms["ka"]*y[1]/parms["Vd"]-parms["kel"]*y[2]+parms["k21"]*y[3]-parms["k12"]*y[2]
       dCp3dt <- parms["k12"]*y[2]-parms["k21"]*y[3]
       list(c(dCp1dt,dCp2dt,dCp3dt)) 
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
     ka<-par1
     kel<-par2
     k12<-par3
     k21<-par4
     Vd<-par5   
     sfirst2.out(PKtime,ka,kel,k12,k21,Vd,defun,par1,par2,par3,par4,par5,Dose,i)  
     }   
  }   
  else if (pick == 2){
     cat("\n\nEnter error factor for ka\n")
     factor1<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for kel\n")
     factor2<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for k12\n")
     factor3<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for k21\n")
     factor4<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for Vd\n")
     factor5<-scan(nlines=1,quiet=TRUE)
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" ) 
     ka<-par1+rnorm(1,mean=0,sd=factor1)
     while(ka<=0){
           ka<-par1+rnorm(1,mean=0,sd=factor1)}          
     kel<-par2+rnorm(1,mean=0,sd=factor2)
     while(kel<=0){
           kel<-par2+rnorm(1,mean=0,sd=factor2)} 
     k12<-par3+rnorm(1,mean=0,sd=factor3)
     while(k12<=0){
           k12<-par3+rnorm(1,mean=0,sd=factor3)}
     k21<-par4+rnorm(1,mean=0,sd=factor4)
     while(k21<=0){
           k21<-par4+rnorm(1,mean=0,sd=factor4)}
     Vd<-par5+rnorm(1,mean=0,sd=factor5)
     while(Vd<=0){
           Vd<-par5+rnorm(1,mean=0,sd=factor5)}
     sfirst2.out(PKtime,ka,kel,k12,k21,Vd,defun,par1,par2,par3,par4,par5,Dose,i)      
     }       
  }
  else if (pick == 3){
     cat("\n\nEnter error factor for ka\n")
     factor1<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for kel\n")
     factor2<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for k12\n")
     factor3<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for k21\n")
     factor4<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for Vd\n")
     factor5<-scan(nlines=1,quiet=TRUE)
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" ) 
     ka<-par1+runif(1,min=-factor1,max=factor1)
     while(ka<=0){
           ka<-par1+runif(1,min=-factor1,max=factor1)}      
     kel<-par2+runif(1,min=-factor2,max=factor2)
     while(kel<=0){
           kel<-par2+runif(1,min=-factor2,max=factor2)}
     k12<-par3+runif(1,min=-factor3,max=factor3)
     while(k12<=0){
           k12<-par3+runif(1,min=-factor3,max=factor3)}
     k21<-par4+runif(1,min=-factor4,max=factor4)
     while(k21<=0){
           k21<-par4+runif(1,min=-factor4,max=factor4)}
     Vd<-par5+runif(1,min=-factor5,max=factor5)
     while(Vd<=0){
           Vd<-par5+runif(1,min=-factor5,max=factor5)}
     sfirst2.out(PKtime,ka,kel,k12,k21,Vd,defun,par1,par2,par3,par4,par5,Dose,i)       
     }       
  }
  else if (pick == 4){
     cat("\n\nEnter error factor for ka\n")
     factor1<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for kel\n")
     factor2<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for k12\n")
     factor3<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for k21\n")
     factor4<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for Vd\n")
     factor5<-scan(nlines=1,quiet=TRUE)
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" ) 
     ka<-par1*rnorm(1,mean=0,sd=factor1)+par1
     while(ka<=0){
           ka<-par1*rnorm(1,mean=0,sd=factor1)+par1}         
     kel<-par2*rnorm(1,mean=0,sd=factor2)+par2
     while(kel<=0){
           kel<-par2*rnorm(1,mean=0,sd=factor2)+par2}
     k12<-par3*rnorm(1,mean=0,sd=factor3)+par3
     while(k12<=0){
           k12<-par3*rnorm(1,mean=0,sd=factor3)+par3}
     k21<-par4*rnorm(1,mean=0,sd=factor4)+par4
     while(k21<=0){
           k21<-par4*rnorm(1,mean=0,sd=factor4)+par4}
     Vd<-par5*rnorm(1,mean=0,sd=factor5)+par5
     while(Vd<=0){
           Vd<-par5*rnorm(1,mean=0,sd=factor5)+par5}
     sfirst2.out(PKtime,ka,kel,k12,k21,Vd,defun,par1,par2,par3,par4,par5,Dose,i)       
     }       
  }
  else if (pick == 5){
     cat("\n\nEnter error factor for ka\n")
     factor1<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for kel\n")
     factor2<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for k12\n")
     factor3<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for k21\n")
     factor4<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for Vd\n")
     factor5<-scan(nlines=1,quiet=TRUE)
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" ) 
     ka<-par1*runif(1,min=-factor1,max=factor1)+par1
     while(ka<=0){
           ka<-par1*runif(1,min=-factor1,max=factor1)+par1}
     kel<-par2*runif(1,min=-factor2,max=factor2)+par2
     while(kel<=0){
           kel<-par2*runif(1,min=-factor2,max=factor2)+par2}
     k12<-par3*runif(1,min=-factor3,max=factor3)+par3
     while(k12<=0){
           k12<-par3*runif(1,min=-factor3,max=factor3)+par3}
     k21<-par4*runif(1,min=-factor4,max=factor4)+par4
     while(k21<=0){
           k21<-par4*runif(1,min=-factor4,max=factor4)+par4}
     Vd<-par5*runif(1,min=-factor5,max=factor5)+par5
     while(Vd<=0){
           Vd<-par5*runif(1,min=-factor5,max=factor5)+par5}
     sfirst2.out(PKtime,ka,kel,k12,k21,Vd,defun,par1,par2,par3,par4,par5,Dose,i)       
     }       
  }
}
