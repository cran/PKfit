#Normal fitting
#One exponential term
fmacro.one <- function(PKindex)
{
  options(warn=-1)
    
 #lsoda is belong to odesolve package-----
  library(odesolve)
    
 #BIC function is belong to stats4 package-----
  library(stats4)  
    
 #genoud function is belong to rgenoud package 
  library(rgenoud)
    
  par<-data.frame(Parameter=c("A","a"),Initial=c(0))
  par<-edit(par)
  show(par)
  cat("\n")
    
  defun<-function(time,A,a){
     pred<- A*exp(-a*time)
  }
    
  file.menu <-c("equal weight",
                "1/Cp",
                "1/Cp^2")
  pick <- menu(file.menu, title="<< weighting schemes >>")
  if (pick==1){
     cat("\nEnter the title of x-axis(time)\n\n")
     xaxis<-readline()
     if (substr(xaxis, 1, 1) == "")  xaxis<-"Time"  else xaxis<-xaxis
     cat("\nEnter the title of y-axis(Cp)\n\n")
     yaxis<-readline()
     if (substr(yaxis, 1, 1) == "")  yaxis<-"Concentration"  else yaxis<-yaxis
     for( i in 1:length(unique(PKindex$Subject)))  {
     cat("\n\n               << Subject",i,">>\n\n" )  
     objfun<-function(par) {
     out<-defun(PKindex$conc[PKindex$Subject==i],par[1],par[2])
     gift<- which( PKindex$conc[PKindex$Subject==i] != 0 )
     sum((PKindex$conc[PKindex$Subject==i][gift]-out[gift])^2)}
     gen<-genoud(objfun,nvars=2,max=FALSE,pop.size=30,max.generations=20,wait.generations=10,
          starting.value=c(par[1,2],par[2,2]),BFGS=FALSE,print.level=0,boundary.enforcement=0,
          Domains=matrix(c(0,0,100,10),2,2),MemoryMatrix=TRUE,lexical=FALSE)
     cat("<< The value of parameter fitted by genetic algorithm >>\n\n")   
     namegen<-c("A","a")
     outgen<-c(gen$par[1],gen$par[2])
     print(data.frame(Parameter=namegen,Value=outgen)) 
     F<-objfun(gen$par)
     #opt<-optim(c(gen$par[1],gen$par[2]),objfun,method="Nelder-Mead")       
     #cat("\n<< The value of parameter fitted by Nelder-Mead Simplex slgorithm >>\n\n")   
     #nameopt<-c("A","a")
     #outopt<-c(opt$par[1],opt$par[2])
     #print(data.frame(Parameter=nameopt,Value=outopt))
     cat("\n<< Residuals sum-of-squares and parameter values fitted by nls >>\n\n")
     fm<-nls(conc~defun(time,A,a),data=PKindex, algorithm="default",subset=Subject==i,
         start=list(A=gen$par[1],a=gen$par[2]),trace=TRUE,
         nls.control(maxiter=500,tol=1e-2,minFactor=1/1024))
     cat("\n")   
     coef<-data.frame(coef(fm)["a"])      
     plotting.lin(PKindex, fm, i, pick, coef, xaxis, yaxis)
     }
  }   
  else if (pick==2){
     cat("\nEnter the title of x-axis(time)\n\n")
     xaxis<-readline()
     if (substr(xaxis, 1, 1) == "")  xaxis<-"Time"  else xaxis<-xaxis
     cat("\nEnter the title of y-axis(Cp)\n\n")
     yaxis<-readline()
     if (substr(yaxis, 1, 1) == "")  yaxis<-"Concentration"  else yaxis<-yaxis
     for( i in 1:length(unique(PKindex$Subject)))  {
     cat("\n\n               << Subject",i,">>\n\n" )  
     objfun<-function(par) {
     out<-defun(PKindex$conc[PKindex$Subject==i],par[1],par[2])
     gift<- which( PKindex$conc[PKindex$Subject==i] != 0 )
     sum((PKindex$conc[PKindex$Subject==i][gift]-out[gift])^2/PKindex$conc[PKindex$Subject==i][gift])}
     gen<-genoud(objfun,nvars=2,max=FALSE,pop.size=30,max.generations=20,wait.generations=10,
          starting.value=c(par[1,2],par[2,2]),BFGS=FALSE,print.level=0,boundary.enforcement=0,
          Domains=matrix(c(0,0,100,10),2,2),MemoryMatrix=TRUE,lexical=FALSE)
     cat("<< The value of parameter fitted by genetic algorithm >>\n\n")   
     namegen<-c("A","a")
     outgen<-c(gen$par[1],gen$par[2])
     print(data.frame(Parameter=namegen,Value=outgen)) 
     F<-objfun(gen$par)
     #opt<-optim(c(gen$par[1],gen$par[2]),objfun,method="Nelder-Mead")       
     #cat("\n<< The value of parameter fitted by Nelder-Mead Simplex slgorithm >>\n\n")   
     #nameopt<-c("A","a")
     #outopt<-c(opt$par[1],opt$par[2])
     #print(data.frame(Parameter=nameopt,Value=outopt))
     cat("\n<< Residuals sum-of-squares and parameter values fitted by nls >>\n\n")
     fm<-nls(conc~defun(time,A,a),data=PKindex, algorithm="default",subset=Subject==i,
         start=list(A=gen$par[1],a=gen$par[2]),trace=TRUE,
         nls.control(maxiter=500,tol=1e-2,minFactor=1/1024))
     cat("\n")   
     coef<-data.frame(coef(fm)["a"])      
     plotting.lin(PKindex, fm, i, pick, coef, xaxis, yaxis)
     }
  }        
  else if (pick==3){
     cat("\nEnter the title of x-axis(time)\n\n")
     xaxis<-readline()
     if (substr(xaxis, 1, 1) == "")  xaxis<-"Time"  else xaxis<-xaxis
     cat("\nEnter the title of y-axis(Cp)\n\n")
     yaxis<-readline()
     if (substr(yaxis, 1, 1) == "")  yaxis<-"Concentration"  else yaxis<-yaxis
     for( i in 1:length(unique(PKindex$Subject)))  {
     cat("\n\n               << Subject",i,">>\n\n" )  
     objfun<-function(par) {
     out<-defun(PKindex$conc[PKindex$Subject==i],par[1],par[2])
     gift<- which( PKindex$conc[PKindex$Subject==i] != 0 )
     sum(((PKindex$conc[PKindex$Subject==i][gift]-out[gift])/PKindex$conc[PKindex$Subject==i][gift])^2)}
     gen<-genoud(objfun,nvars=2,max=FALSE,pop.size=30,max.generations=20,wait.generations=10,
          starting.value=c(par[1,2],par[2,2]),BFGS=FALSE,print.level=0,boundary.enforcement=0,
          Domains=matrix(c(0,0,100,10),2,2),MemoryMatrix=TRUE,lexical=FALSE)
     cat("<< The value of parameter fitted by genetic algorithm >>\n\n")   
     namegen<-c("A","a")
     outgen<-c(gen$par[1],gen$par[2])
     print(data.frame(Parameter=namegen,Value=outgen)) 
     F<-objfun(gen$par)
     #opt<-optim(c(gen$par[1],gen$par[2]),objfun,method="Nelder-Mead")       
     #cat("\n<< The value of parameter fitted by Nelder-Mead Simplex slgorithm >>\n\n")   
     #nameopt<-c("A","a")
     #outopt<-c(opt$par[1],opt$par[2])
     #print(data.frame(Parameter=nameopt,Value=outopt))
     cat("\n<< Residuals sum-of-squares and parameter values fitted by nls >>\n\n")
     fm<-nls(conc~defun(time,A,a),data=PKindex, algorithm="default",subset=Subject==i,
         start=list(A=gen$par[1],a=gen$par[2]),trace=TRUE,
         nls.control(maxiter=500,tol=1e-2,minFactor=1/1024))
     cat("\n")   
     coef<-data.frame(coef(fm)["a"])      
     plotting.lin(PKindex, fm, i, pick, coef, xaxis, yaxis)
     }
  }    
  cat("\n")  
}

#Two exponential term
fmacro.two <- function(PKindex)
{
  options(warn=-1)
    
 #lsoda is belong to odesolve package
  library(odesolve)
    
 #BIC function is belong to stats4 package
  library(stats4)  
    
 #genoud function is belong to rgenoud package 
  library(rgenoud)
    
  par<-data.frame(Parameter=c("A","a","B","b"),Initial=c(0))
  par<-edit(par)
  show(par)
  cat("\n")
      
  defun<-function(time,A,a,B,b){
     pred<-A*exp(-a*time)+B*exp(-b*time)
  }    
    
  file.menu <-c("equal weight",
                "1/Cp",
                "1/Cp^2")
  pick <- menu(file.menu, title="<< weighting schemes >>")
  if (pick==1){
     cat("\nEnter the title of x-axis(time)\n\n")
     xaxis<-readline()
     if (substr(xaxis, 1, 1) == "")  xaxis<-"Time"  else xaxis<-xaxis
     cat("\nEnter the title of y-axis(Cp)\n\n")
     yaxis<-readline()
     if (substr(yaxis, 1, 1) == "")  yaxis<-"Concentration"  else yaxis<-yaxis
     for( i in 1:length(unique(PKindex$Subject)))  {
     cat("\n\n               << Subject",i,">>\n\n" )  
     objfun<-function(par) {
     out<-defun(PKindex$time[PKindex$Subject==i],par[1],par[2],par[3],par[4])
     gift<- which( PKindex$conc[PKindex$Subject==i] != 0 )
     sum((PKindex$conc[PKindex$Subject==i][gift]-out[gift])^2)}
     gen<-genoud(objfun,nvars=4,max=FALSE,pop.size=30,max.generations=20,wait.generations=10,
          starting.value=c(par[1,2],par[2,2],par[3,2],par[4,2]),BFGS=FALSE,print.level=0,boundary.enforcement=0,
          Domains=matrix(c(1,0.01,0.1,0.01,100,10,50,1),4,2),MemoryMatrix=TRUE,lexical=FALSE)      
     cat("<< The value of parameter fitted by genetic algorithm >>\n\n")   
     namegen<-c("A","a","B","b")
     outgen<-c(gen$par[1],gen$par[2],gen$par[3],gen$par[4])
     print(data.frame(Parameter=namegen,Value=outgen)) 
     F<-objfun(gen$par)
     opt<-optim(c(gen$par[1],gen$par[2],gen$par[3],gen$par[4]),objfun,method="Nelder-Mead")          
     cat("\n<< The value of parameter fitted by Nelder-Mead Simplex slgorithm >>\n\n")   
     nameopt<-c("A","a","B","b")
     outopt<-c(opt$par[1],opt$par[2],opt$par[3],opt$par[4])
     print(data.frame(Parameter=nameopt,Value=outopt))
     cat("\n<< Residuals sum-of-squares and parameter values fitted by nls >>\n\n")
     fm<-nls(conc~defun(time,A,a,B,b),data=PKindex, algorithm="default",subset=Subject==i,
         start=list(A=opt$par[1], a=opt$par[2], B=opt$par[3], b=opt$par[4]),trace=TRUE,
         nls.control(maxiter=500,tol=1e-2,minFactor=1/1024))
     cat("\n")        
     coef<-data.frame(coef(fm)["b"])     
     plotting.lin(PKindex, fm, i, pick, coef, xaxis, yaxis)
     }
  }   
  else if (pick==2){
     cat("\nEnter the title of x-axis(time)\n\n")
     xaxis<-readline()
     if (substr(xaxis, 1, 1) == "")  xaxis<-"Time"  else xaxis<-xaxis
     cat("\nEnter the title of y-axis(Cp)\n\n")
     yaxis<-readline()
     if (substr(yaxis, 1, 1) == "")  yaxis<-"Concentration"  else yaxis<-yaxis
     for( i in 1:length(unique(PKindex$Subject)))  {
     cat("\n\n               << Subject",i,">>\n\n" )  
     objfun<-function(par) {
     out<-defun(PKindex$time[PKindex$Subject==i],par[1],par[2],par[3],par[4])
     gift<- which( PKindex$conc[PKindex$Subject==i] != 0 )
     sum((PKindex$conc[PKindex$Subject==i][gift]-out[gift])^2/PKindex$conc[PKindex$Subject==i][gift])}
     gen<-genoud(objfun,nvars=4,max=FALSE,pop.size=30,max.generations=20,wait.generations=10,
          starting.value=c(par[1,2],par[2,2],par[3,2],par[4,2]),BFGS=FALSE,print.level=0,boundary.enforcement=0,
          Domains=matrix(c(1,0.01,0.1,0.01,100,10,50,1),4,2),MemoryMatrix=TRUE,lexical=FALSE)      
     cat("<< The value of parameter fitted by genetic algorithm >>\n\n")   
     namegen<-c("A","a","B","b")
     outgen<-c(gen$par[1],gen$par[2],gen$par[3],gen$par[4])
     print(data.frame(Parameter=namegen,Value=outgen)) 
     F<-objfun(gen$par)
     opt<-optim(c(gen$par[1],gen$par[2],gen$par[3],gen$par[4]),objfun,method="Nelder-Mead")          
     cat("\n<< The value of parameter fitted by Nelder-Mead Simplex slgorithm >>\n\n")   
     nameopt<-c("A","a","B","b")
     outopt<-c(opt$par[1],opt$par[2],opt$par[3],opt$par[4])
     print(data.frame(Parameter=nameopt,Value=outopt))
     cat("\n<< Residuals sum-of-squares and parameter values fitted by nls >>\n\n")
     fm<-nls(conc~defun(time,A,a,B,b),data=PKindex, algorithm="default",subset=Subject==i,
         start=list(A=opt$par[1], a=opt$par[2], B=opt$par[3], b=opt$par[4]),trace=TRUE,
         nls.control(maxiter=500,tol=1e-2,minFactor=1/1024))
     cat("\n")        
     coef<-data.frame(coef(fm)["b"])     
     plotting.lin(PKindex, fm, i, pick, coef, xaxis, yaxis)
     }
  }        
  else if (pick==3){
     cat("\nEnter the title of x-axis(time)\n\n")
     xaxis<-readline()
     if (substr(xaxis, 1, 1) == "")  xaxis<-"Time"  else xaxis<-xaxis
     cat("\nEnter the title of y-axis(Cp)\n\n")
     yaxis<-readline()
     if (substr(yaxis, 1, 1) == "")  yaxis<-"Concentration"  else yaxis<-yaxis
     for( i in 1:length(unique(PKindex$Subject)))  {
     cat("\n\n               << Subject",i,">>\n\n" )  
     objfun<-function(par) {
     out<-defun(PKindex$time[PKindex$Subject==i],par[1],par[2],par[3],par[4])
     gift<- which( PKindex$conc[PKindex$Subject==i] != 0 )
     sum(((PKindex$conc[PKindex$Subject==i][gift]-out[gift])/PKindex$conc[PKindex$Subject==i][gift])^2)}
     gen<-genoud(objfun,nvars=4,max=FALSE,pop.size=30,max.generations=20,wait.generations=10,
          starting.value=c(par[1,2],par[2,2],par[3,2],par[4,2]),BFGS=FALSE,print.level=0,boundary.enforcement=0,
          Domains=matrix(c(1,0.01,0.1,0.01,100,10,50,1),4,2),MemoryMatrix=TRUE,lexical=FALSE)      
     cat("<< The value of parameter fitted by genetic algorithm >>\n\n")   
     namegen<-c("A","a","B","b")
     outgen<-c(gen$par[1],gen$par[2],gen$par[3],gen$par[4])
     print(data.frame(Parameter=namegen,Value=outgen)) 
     F<-objfun(gen$par)
     opt<-optim(c(gen$par[1],gen$par[2],gen$par[3],gen$par[4]),objfun,method="Nelder-Mead")          
     cat("\n<< The value of parameter fitted by Nelder-Mead Simplex slgorithm >>\n\n")   
     nameopt<-c("A","a","B","b")
     outopt<-c(opt$par[1],opt$par[2],opt$par[3],opt$par[4])
     print(data.frame(Parameter=nameopt,Value=outopt))
     cat("\n<< Residuals sum-of-squares and parameter values fitted by nls >>\n\n")
     fm<-nls(conc~defun(time,A,a,B,b),data=PKindex, algorithm="default",subset=Subject==i,
         start=list(A=opt$par[1], a=opt$par[2], B=opt$par[3], b=opt$par[4]),trace=TRUE,
         nls.control(maxiter=500,tol=1e-2,minFactor=1/1024))
     cat("\n")        
     coef<-data.frame(coef(fm)["b"])     
     plotting.lin(PKindex, fm, i, pick, coef, xaxis, yaxis)
     }
  }     
  cat("\n")   
}

#Three exponential term
fmacro.three <- function(PKindex)
{
  options(warn=-1)
    
 #lsoda is belong to odesolve package-----
  library(odesolve)
    
 #BIC function is belong to stats4 package-----
  library(stats4)     
    
 #genoud function is belong to rgenoud package  
  library(rgenoud)
   
  par<-data.frame(Parameter=c("A","a","B","b","C","c"),Initial=c(0))
  par<-edit(par)
  show(par)
  cat("\n")
     
  defun<-function(time,A,a,B,b,C,c){
     pred<- A*exp(-a*time)+B*exp(-b*time)+C*exp(-c*time)
  }    
    
  file.menu <-c("equal weight",
                "1/Cp",
                "1/Cp^2")
  pick <- menu(file.menu, title="<< weighting schemes >>")
  if (pick==1){
     cat("\nEnter the title of x-axis(time)\n\n")
     xaxis<-readline()
     if (substr(xaxis, 1, 1) == "")  xaxis<-"Time"  else xaxis<-xaxis
     cat("\nEnter the title of y-axis(Cp)\n\n")
     yaxis<-readline()
     if (substr(yaxis, 1, 1) == "")  yaxis<-"Concentration"  else yaxis<-yaxis
     for( i in 1:length(unique(PKindex$Subject)))  {
     cat("\n\n               << Subject",i,">>\n\n" )  
     objfun<-function(par) {
     out<-defun(PKindex$time[PKindex$Subject==i],par[1],par[2],par[3],par[4],par[5],par[6])
     gift<- which( PKindex$conc[PKindex$Subject==i] != 0 )
     sum((PKindex$conc[PKindex$Subject==i][gift]-out[PKindex$Subject==i][gift])^2)}
     gen<-genoud(objfun,nvars=6,max=FALSE,pop.size=30,max.generations=20,wait.generations=10,
          starting.value=c(par[1,2],par[2,2],par[3,2],par[4,2],par[5,2],par[6,2]),BFGS=FALSE,print.level=0,boundary.enforcement=0,
          Domains=matrix(c(1,0.1,1,0.01,0.1,0.01,100,10,50,10,10,1),6,2),MemoryMatrix=TRUE,lexical=FALSE)
     cat("<< The value of parameter fitted by genetic algorithm >>\n\n")   
     namegen<-c("A","a","B","b","C","c")
     outgen<-c(gen$par[1],gen$par[2],gen$par[3],gen$par[4],gen$par[5],gen$par[6])
     print(data.frame(Parameter=namegen,Value=outgen)) 
     F<-objfun(gen$par)
     opt<-optim(c(gen$par[1],gen$par[2],gen$par[3],gen$par[4],gen$par[5],gen$par[6]),objfun,method="Nelder-Mead")       
     cat("\n<< The value of parameter fitted by Nelder-Mead Simplex slgorithm >>\n\n")   
     nameopt<-c("A","a","B","b","C","c")
     outopt<-c(opt$par[1],opt$par[2],opt$par[3],opt$par[4],opt$par[5],opt$par[6])
     print(data.frame(Parameter=nameopt,Value=outopt)) 
     cat("\n<< Residuals sum-of-squares and parameter values fitted by nls >>\n\n")
     fm<-nls(conc~defun(time,A,a,B,b,C,c),data=PKindex, algorithm="default",subset=Subject==i,
         start=list(A=opt$par[1],a=opt$par[2],B=opt$par[3],b=opt$par[4],C=opt$par[5],c=opt$par[6]),
         trace=TRUE,nls.control(maxiter=500,tol=1,minFactor=1/2048))
     cat("\n")
     coef<-data.frame(coef(fm)["c"])     
     plotting.lin(PKindex, fm, i, pick, coef, xaxis, yaxis)
     }
  }   
  else if (pick==2){
     cat("\nEnter the title of x-axis(time)\n\n")
     xaxis<-readline()
     if (substr(xaxis, 1, 1) == "")  xaxis<-"Time"  else xaxis<-xaxis
     cat("\nEnter the title of y-axis(Cp)\n\n")
     yaxis<-readline()
     if (substr(yaxis, 1, 1) == "")  yaxis<-"Concentration"  else yaxis<-yaxis
     for( i in 1:length(unique(PKindex$Subject)))  {
     cat("\n\n               << Subject",i,">>\n\n" )  
     objfun<-function(par) {
     out<-defun(PKindex$time[PKindex$Subject==i],par[1],par[2],par[3],par[4],par[5],par[6])
     gift<- which( PKindex$conc[PKindex$Subject==i] != 0 )
     sum((PKindex$conc[PKindex$Subject==i][gift]-out[PKindex$Subject==i][gift])^2/PKindex$conc[PKindex$Subject==i][gift])}
     gen<-genoud(objfun,nvars=6,max=FALSE,pop.size=30,max.generations=20,wait.generations=10,
          starting.value=c(par[1,2],par[2,2],par[3,2],par[4,2],par[5,2],par[6,2]),BFGS=FALSE,print.level=0,boundary.enforcement=0,
          Domains=matrix(c(1,0.1,1,0.01,0.1,0.01,100,10,50,10,10,1),6,2),MemoryMatrix=TRUE,lexical=FALSE)
     cat("<< The value of parameter fitted by genetic algorithm >>\n\n")   
     namegen<-c("A","a","B","b","C","c")
     outgen<-c(gen$par[1],gen$par[2],gen$par[3],gen$par[4],gen$par[5],gen$par[6])
     print(data.frame(Parameter=namegen,Value=outgen)) 
     F<-objfun(gen$par)
     opt<-optim(c(gen$par[1],gen$par[2],gen$par[3],gen$par[4],gen$par[5],gen$par[6]),objfun,method="Nelder-Mead")       
     cat("\n<< The value of parameter fitted by Nelder-Mead Simplex slgorithm >>\n\n")   
     nameopt<-c("A","a","B","b","C","c")
     outopt<-c(opt$par[1],opt$par[2],opt$par[3],opt$par[4],opt$par[5],opt$par[6])
     print(data.frame(Parameter=nameopt,Value=outopt))
     cat("\n<< Residuals sum-of-squares and parameter values fitted by nls >>\n\n")
     fm<-nls(conc~defun(time,A,a,B,b,C,c),data=PKindex, algorithm="default",subset=Subject==i,
         start=list(A=opt$par[1],a=opt$par[2],B=opt$par[3],b=opt$par[4],C=opt$par[5],c=opt$par[6]),
         trace=TRUE,nls.control(maxiter=500,tol=1,minFactor=1/2048))
     cat("\n")
     coef<-data.frame(coef(fm)["c"])     
     plotting.lin(PKindex, fm, i, pick, coef, xaxis, yaxis)
     }
  }        
  else if (pick==3){
     cat("\nEnter the title of x-axis(time)\n\n")
     xaxis<-readline()
     if (substr(xaxis, 1, 1) == "")  xaxis<-"Time"  else xaxis<-xaxis
     cat("\nEnter the title of y-axis(Cp)\n\n")
     yaxis<-readline()
     if (substr(yaxis, 1, 1) == "")  yaxis<-"Concentration"  else yaxis<-yaxis
     for( i in 1:length(unique(PKindex$Subject)))  {
     cat("\n\n               << Subject",i,">>\n\n" )  
     objfun<-function(par) {
     out<-defun(PKindex$time[PKindex$Subject==i],par[1],par[2],par[3],par[4],par[5],par[6])
     gift<- which( PKindex$conc[PKindex$Subject==i] != 0 )
     sum(((PKindex$conc[PKindex$Subject==i][gift]-out[PKindex$Subject==i][gift])/PKindex$conc[PKindex$Subject==i][gift])^2)}
     gen<-genoud(objfun,nvars=6,max=FALSE,pop.size=30,max.generations=20,wait.generations=10,
          starting.value=c(par[1,2],par[2,2],par[3,2],par[4,2],par[5,2],par[6,2]),BFGS=FALSE,print.level=0,boundary.enforcement=0,
          Domains=matrix(c(1,0.1,1,0.01,0.1,0.01,100,10,50,10,10,1),6,2),MemoryMatrix=TRUE,lexical=FALSE)
     cat("<< The value of parameter fitted by genetic algorithm >>\n\n")   
     namegen<-c("A","a","B","b","C","c")
     outgen<-c(gen$par[1],gen$par[2],gen$par[3],gen$par[4],gen$par[5],gen$par[6])
     print(data.frame(Parameter=namegen,Value=outgen))  
     F<-objfun(gen$par)
     opt<-optim(c(gen$par[1],gen$par[2],gen$par[3],gen$par[4],gen$par[5],gen$par[6]),objfun,method="Nelder-Mead")       
     cat("\n<< The value of parameter fitted by Nelder-Mead Simplex slgorithm >>\n\n")   
     nameopt<-c("A","a","B","b","C","c")
     outopt<-c(opt$par[1],opt$par[2],opt$par[3],opt$par[4],opt$par[5],opt$par[6])
     print(data.frame(Parameter=nameopt,Value=outopt))
     cat("\n<< Residuals sum-of-squares and parameter values fitted by nls >>\n\n")
     fm<-nls(conc~defun(time,A,a,B,b,C,c),data=PKindex, algorithm="default",subset=Subject==i,
         start=list(A=opt$par[1],a=opt$par[2],B=opt$par[3],b=opt$par[4],C=opt$par[5],c=opt$par[6]),
         trace=TRUE,nls.control(maxiter=500,tol=1,minFactor=1/2048))
     cat("\n")
     coef<-data.frame(coef(fm)["c"])     
     plotting.lin(PKindex, fm, i, pick, coef, xaxis, yaxis)
     }
  }     
  cat("\n")       
}

  
#Simulation
#One exponential term
smacro.one <-function()
{
  cat("How many subject do you want?\n")
  Subject<-scan(nlines=1,quiet=TRUE)
  PKtime<-data.frame(time=c(0))
  PKtime<-edit(PKtime) 
  cat("\n")
  show(PKtime)

  par<-data.frame(Parameter=c("A","a"),Initial=c(0))
  par<-edit(par)
  cat("\n")
  
  par1<-par[1,2]
  par2<-par[2,2]
 
  file.menu <- c("no error",
                 "error=normal error", 
                 "error=uniform error",
                 "error=normal error*true value",
                 "error=uniform error*true value")
  pick <- menu(file.menu, title = "<< Error type >>")
  if (pick==1){
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" ) 
     A<-par1
     a<-par2     
     smacro.one.out(PKtime,A,a,defun,par1,par2,Dose,i)
     }  
  }   
  else if (pick == 2){
     cat("\n\nEnter error factor for A\n")
     factor1<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for a\n")
     factor2<-scan(nlines=1,quiet=TRUE)
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" ) 
     A<-par1+rnorm(1,mean=0,sd=factor1)
     while(A<=0){
           A<-par1+rnorm(1,mean=0,sd=factor1)}
     a<-par2+rnorm(1,mean=0,sd=factor2)
     while(a<=0){
           a<-par2+rnorm(1,mean=0,sd=factor2)}
     smacro.one.out(PKtime,A,a,defun,par1,par2,Dose,i)  
     }    
  }
  else if (pick == 3){
     cat("\n\nEnter error factor for A\n")
     factor1<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for a\n")
     factor2<-scan(nlines=1,quiet=TRUE)
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" )
     A<-par1+runif(1,min=-factor1,max=factor1)
     while(A<=0){
           A<-par1+runif(1,min=-factor1,max=factor1)}
     a<-par2+runif(1,min=-factor2,max=factor2)
     while(a<=0){
           a<-par2+runif(1,min=-factor2,max=factor2)}
     smacro.one.out(PKtime,A,a,defun,par1,par2,Dose,i)       
     }       
  }
  else if (pick == 4){
     cat("\n\nEnter error factor for A\n")
     factor1<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for a\n")
     factor2<-scan(nlines=1,quiet=TRUE)
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" )
     A<-par1*rnorm(1,mean=0,sd=factor1)+par1
     while(A<=0){
           A<-par1*rnorm(1,mean=0,sd=factor1)+par1}
     a<-par2*rnorm(1,mean=0,sd=factor2)+par2
     while(a<=0){
           a<-par2*rnorm(1,mean=0,sd=factor2)+par2}
     smacro.one.out(PKtime,A,a,defun,par1,par2,Dose,i)      
     }       
  }
  else if (pick == 5){
     cat("\n\nEnter error factor for A\n")
     factor1<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for a\n")
     factor2<-scan(nlines=1,quiet=TRUE)
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" )
     A<-par1*runif(1,min=-factor1,max=factor1)+par1
     while(A<=0){
           A<-par1*runif(1,min=-factor1,max=factor1)+par1}
     a<-par2*runif(1,min=-factor2,max=factor2)+par2
     while(a<=0){
           a<-par2*runif(1,min=-factor2,max=factor2)+par2}
     smacro.one.out(PKtime,A,a,defun,par1,par2,Dose,i)      
     }       
  }
}      

#Two exponential term
smacro.two <-function()
{
  cat("How many subject do you want?\n")
  Subject<-scan(nlines=1,quiet=TRUE)
  PKtime<-data.frame(time=c(0))
  PKtime<-edit(PKtime) 
  cat("\n")
  show(PKtime)

  par<-data.frame(Parameter=c("A","a","B","b"),Initial=c(0))
  par<-edit(par)
  cat("\n")
  
  par1<-par[1,2]
  par2<-par[2,2]
  par3<-par[3,2]
  par4<-par[4,2]  
 
  file.menu <- c("no error",
                 "error=normal error", 
                 "error=uniform error",
                 "error=normal error*true value",
                 "error=uniform error*true value")
  pick <- menu(file.menu, title = "<< Error type >>")
  if (pick==1){
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" ) 
     A<-par1
     a<-par2   
     B<-par3
     b<-par4  
     smacro.two.out(PKtime,A,a,B,b,defun,par1,par2,par3,par4,Dose,i)
     }  
  }   
  else if (pick == 2){
     cat("\n\nEnter error factor for A\n")
     factor1<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for a\n")
     factor2<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for B\n")
     factor3<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for b\n")
     factor4<-scan(nlines=1,quiet=TRUE)
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" )
     A<-par1+rnorm(1,mean=0,sd=factor1)
     while(A<=0){
           A<-par1+rnorm(1,mean=0,sd=factor1)}
     a<-par2+rnorm(1,mean=0,sd=factor2)
     while(a<=0){
           a<-par2+rnorm(1,mean=0,sd=factor2)}
     B<-par3+rnorm(1,mean=0,sd=factor3)
     while(B<=0){
           B<-par3+rnorm(1,mean=0,sd=factor3)}
     b<-par4+rnorm(1,mean=0,sd=factor4)
     while(b<=0){
           b<-par4+rnorm(1,mean=0,sd=factor4)}
     smacro.two.out(PKtime,A,a,B,b,defun,par1,par2,par3,par4,Dose,i)     
     }       
  }
  else if (pick == 3){
     cat("\n\nEnter error factor for A\n")
     factor1<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for a\n")
     factor2<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for B\n")
     factor3<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for b\n")
     factor4<-scan(nlines=1,quiet=TRUE)
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" )
     A<-par1+runif(1,min=-factor1,max=factor1)
     while(A<=0){
           A<-par1+runif(1,min=-factor1,max=factor1)}
     a<-par2+runif(1,min=-factor2,max=factor2)
     while(a<=0){
           a<-par2+runif(1,min=-factor2,max=factor2)}
     B<-par3+runif(1,min=-factor3,max=factor3)
     while(B<=0){
           B<-par1+runif(1,min=-factor3,max=factor3)}
     b<-par4+runif(1,min=-factor4,max=factor4)
     while(b<=0){
           b<-par4+runif(1,min=-factor4,max=factor4)}
     smacro.two.out(PKtime,A,a,B,b,defun,par1,par2,par3,par4,Dose,i)      
     }       
  }
  else if (pick == 4){
     cat("\n\nEnter error factor for A\n")
     factor1<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for a\n")
     factor2<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for B\n")
     factor3<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for b\n")
     factor4<-scan(nlines=1,quiet=TRUE)
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" )
     A<-par1*rnorm(1,mean=0,sd=factor1)+par1
     while(A<=0){
           A<-par1*rnorm(1,mean=0,sd=factor1)+par1}
     a<-par2*rnorm(1,mean=0,sd=factor2)+par2
     while(a<=0){
           a<-par2*rnorm(1,mean=0,sd=factor2)+par2}
     B<-par3*rnorm(1,mean=0,sd=factor3)+par3
     while(B<=0){
           B<-par3*rnorm(1,mean=0,sd=factor3)+par3}
     b<-par4*rnorm(1,mean=0,sd=factor4)+par4
     while(b<=0){
           b<-par4*rnorm(1,mean=0,sd=factor4)+par4}
     smacro.two.out(PKtime,A,a,B,b,defun,par1,par2,par3,par4,Dose,i)        
     }       
  }
  else if (pick == 5){
     cat("\n\nEnter error factor for A\n")
     factor1<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for a\n")
     factor2<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for B\n")
     factor3<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for b\n")
     factor4<-scan(nlines=1,quiet=TRUE)
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" )
     A<-par1*runif(1,min=-factor1,max=factor1)+par1
     while(A<=0){
           A<-par1*runif(1,min=-factor1,max=factor1)+par1}
     a<-par2*runif(1,min=-factor2,max=factor2)+par2
     while(a<=0){
           a<-par2*runif(1,min=-factor2,max=factor2)+par2}
     B<-par3*runif(1,min=-factor3,max=factor3)+par3
     while(B<=0){
           B<-par3*runif(1,min=-factor3,max=factor3)+par3}
     b<-par4*runif(1,min=-factor4,max=factor4)+par4
     while(b<=0){
           b<-par4*runif(1,min=-factor4,max=factor4)+par4}
     smacro.two.out(PKtime,A,a,B,b,defun,par1,par2,par3,par4,Dose,i)          
     }       
  }
}

#Three exponential term
smacro.three <- function()
{
  cat("How many subject do you want?\n")
  Subject<-scan(nlines=1,quiet=TRUE)
  PKtime<-data.frame(time=c(0))
  PKtime<-edit(PKtime) 
  cat("\n")
  show(PKtime)

  par<-data.frame(Parameter=c("A","a","B","b","C","c"),Initial=c(0))
  par<-edit(par)
  cat("\n")
  
  par1<-par[1,2]
  par2<-par[2,2]
  par3<-par[3,2]
  par4<-par[4,2]
  par5<-par[5,2]
  par6<-par[6,2]
 
  file.menu <- c("no error",
                 "error=normal error", 
                 "error=uniform error",
                 "error=normal error*true value",
                 "error=uniform error*true value")
  pick <- menu(file.menu, title = "<< Error type >>")
  if (pick==1){
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" )
     A<-par1
     a<-par2   
     B<-par3
     b<-par4 
     C<-par5
     c<-par6 
     smacro.three.out(PKtime,A,a,B,b,C,c,defun,par1,par2,par3,par4,par5,par6,Dose,i)  
     }  
  }   
  else if (pick == 2){
     cat("\n\nEnter error factor for A\n")
     factor1<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for a\n")
     factor2<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for B\n")
     factor3<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for b\n")
     factor4<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for C\n")
     factor5<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for c\n")
     factor6<-scan(nlines=1,quiet=TRUE)
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" )
     A<-par1+rnorm(1,mean=0,sd=factor1)
     while(A<=0){
           A<-par1+rnorm(1,mean=0,sd=factor1)}
     a<-par2+rnorm(1,mean=0,sd=factor2)
     while(a<=0){
           a<-par2+rnorm(1,mean=0,sd=factor2)}
     B<-par3+rnorm(1,mean=0,sd=factor3)
     while(B<=0){
           B<-par3+rnorm(1,mean=0,sd=factor3)}
     b<-par4+rnorm(1,mean=0,sd=factor4)
     while(b<=0){
           b<-par4+rnorm(1,mean=0,sd=factor4)}
     C<-par5+rnorm(1,mean=0,sd=factor5)
     while(C<=0){
           C<-par5+rnorm(1,mean=0,sd=factor5)}
     c<-par6+rnorm(1,mean=0,sd=factor6)
     while(c<=0){
           c<-par6+rnorm(1,mean=0,sd=factor6)}
     smacro.three.out(PKtime,A,a,B,b,C,c,defun,par1,par2,par3,par4,par5,par6,Dose,i)      
     }       
  }
  else if (pick == 3){
     cat("\n\nEnter error factor for A\n")
     factor1<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for a\n")
     factor2<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for B\n")
     factor3<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for b\n")
     factor4<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for C\n")
     factor5<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for c\n")
     factor6<-scan(nlines=1,quiet=TRUE)
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" )
     A<-par1+runif(1,min=-factor1,max=factor1)
     while(A<=0){
           A<-par1+runif(1,min=-factor1,max=factor1)}
     a<-par2+runif(1,min=-factor2,max=factor2)
     while(a<=0){
           a<-par2+runif(1,min=-factor2,max=factor2)}
     B<-par3+runif(1,min=-factor3,max=factor3)
     while(B<=0){
           B<-par1+runif(1,min=-factor3,max=factor3)}
     b<-par4+runif(1,min=-factor4,max=factor4)
     while(b<=0){
           b<-par4+runif(1,min=-factor4,max=factor4)}
     C<-par5+runif(1,min=-factor5,max=factor5)
     while(C<=0){
           C<-par5+runif(1,min=-factor5,max=factor5)}
     c<-par6+runif(1,min=-factor6,max=factor6)
     while(c<=0){
           c<-par6+runif(1,min=-factor6,max=factor6)}
     smacro.three.out(PKtime,A,a,B,b,C,c,defun,par1,par2,par3,par4,par5,par6,Dose,i)      
     }       
  }
  else if (pick == 4){
     cat("\n\nEnter error factor for A\n")
     factor1<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for a\n")
     factor2<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for B\n")
     factor3<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for b\n")
     factor4<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for C\n")
     factor5<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for c\n")
     factor6<-scan(nlines=1,quiet=TRUE)
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" )
     A<-par1*rnorm(1,mean=0,sd=factor1)+par1
     while(A<=0){
           A<-par1*rnorm(1,mean=0,sd=factor1)+par1}
     a<-par2*rnorm(1,mean=0,sd=factor2)+par2
     while(a<=0){
           a<-par2*rnorm(1,mean=0,sd=factor2)+par2}
     B<-par3*rnorm(1,mean=0,sd=factor3)+par3
     while(B<=0){
           B<-par3*rnorm(1,mean=0,sd=factor3)+par3}
     b<-par4*rnorm(1,mean=0,sd=factor4)+par4
     while(b<=0){
           b<-par4*rnorm(1,mean=0,sd=factor4)+par4}
     C<-par5*rnorm(1,mean=0,sd=factor5)+par5
     while(C<=0){
           C<-par5*rnorm(1,mean=0,sd=factor5)+par5}
     c<-par6*rnorm(1,mean=0,sd=factor6)+par6
     while(c<=0){
           c<-par6*rnorm(1,mean=0,sd=factor6)+par6}
     smacro.three.out(PKtime,A,a,B,b,C,c,defun,par1,par2,par3,par4,par5,par6,Dose,i)      
     }       
  }
  else if (pick == 5){
     cat("\n\nEnter error factor for A\n")
     factor1<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for a\n")
     factor2<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for B\n")
     factor3<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for b\n")
     factor4<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for C\n")
     factor5<-scan(nlines=1,quiet=TRUE)
     cat("\nEnter error factor for c\n")
     factor6<-scan(nlines=1,quiet=TRUE)
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" )
     A<-par1*runif(1,min=-factor1,max=factor1)+par1
     while(A<=0){
           A<-par1*runif(1,min=-factor1,max=factor1)+par1}
     a<-par2*runif(1,min=-factor2,max=factor2)+par2
     while(a<=0){
           a<-par2*runif(1,min=-factor2,max=factor2)+par2}
     B<-par3*runif(1,min=-factor3,max=factor3)+par3
     while(B<=0){
           B<-par3*runif(1,min=-factor3,max=factor3)+par3}
     b<-par4*runif(1,min=-factor4,max=factor4)+par4
     while(b<=0){
           b<-par4*runif(1,min=-factor4,max=factor4)+par4}
     C<-par5*runif(1,min=-factor5,max=factor5)+par5
     while(C<=0){
           C<-par5*runif(1,min=-factor5,max=factor5)+par5}
     c<-par6*runif(1,min=-factor6,max=factor6)+par6
     while(c<=0){
           c<-par6*runif(1,min=-factor6,max=factor6)+par6}
     smacro.three.out(PKtime,A,a,B,b,C,c,defun,par1,par2,par3,par4,par5,par6,Dose,i)      
     }       
  }
}
