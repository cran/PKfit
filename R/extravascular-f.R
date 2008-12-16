### PKindex is usually the Dataset.


### Normal fitting
### One compartment PK model extravascualr single dose first-order absorption
### With lag time 
ffirst.lag<- function(PKindex,
                      Dose=NULL, 
                      ka=NULL,
                      Vm=NULL,Km=NULL, 
                      Vd=NULL,
                      kel=NULL,        
                      Tlag=TRUE,
                      MMe=FALSE) 
{
   #options(warn=-1)
        
   ## Input dose and Tlag and initial value for ka, kel and Vd

   if (is.null(Dose)) {
     cat("Enter Dose value\n")
     Dose <- scan(nlines=1,quiet=TRUE)
   } 
   else {
     cat("Dose from arguments is = ",Dose,"\n")
   }

   if ( Tlag ){
       cat("\nEnter lag time value\n")
       Tlag<-scan(nlines=1,quiet=TRUE)
       cat("\n")
   } 
   
   
   if (MMe){
      if (is.null(ka) || is.null(Vm) || is.null(Km) || is.null(Vd) ) {
        par<-data.frame(Parameter=c("ka","Vm","Km","Vd"),Initial=c(0))
        par<-edit(par)
        repeat{
           if ( par[1,2] == 0 || par[2,2] ==0 || par[3,2]==0 ||par[4,2]==0){
             cat("\n")
             cat("**********************************\n")
             cat(" Parameter value can not be zero. \n")
             cat(" Press enter to continue.         \n")
             cat("**********************************\n\n")
             readline()
             cat("\n")
             par<-edit(par)}   
           else{
             break
             return(edit(par))}
        } 
        cat("\n")       
        show(par)

      }
   } 
   else{
      if (is.null(ka) || is.null(kel) || is.null(Vd) ) {
        par<-data.frame(Parameter=c("ka","kel","Vd"),Initial=c(0))
        par<-edit(par)
        repeat{
           if ( par[1,2] == 0 || par[2,2] ==0 || par[3,2]==0){
             cat("\n")
             cat("**********************************\n")
             cat(" Parameter value can not be zero. \n")
             cat(" Press enter to continue.         \n")
             cat("**********************************\n\n")
             readline()
             cat("\n")
             par<-edit(par)}   
           else{
             break
             return(edit(par))}
        } 
        cat("\n")       
        show(par)

      }
   }
   
   
   cat("\n")
   
   if ( MMe ) { 
      ## User-supplied function w Michaelis-Mention elimination & w lag time
      defun<- function(time, y, parms) { 
      if(time <= Tlag) {
        dy1dt<-0
        dy2dt<-0
      }
      else  {
        dy1dt <- -parms["ka"] * y[1]
        dy2dt <-  parms["ka"] * y[1]/parms["Vd"] - (parms["Vm"]/parms["Vd"])*y[2]/(parms["Km"]/parms["Vd"]+y[2])
      }
      list(c(dy1dt,dy2dt)) 
      } 
    
      modfun2 <- function(time,ka,Vm,Km,Vd) { 
      out <- lsoda(c(Dose,0),c(0,time),defun,parms=c(ka=ka,Vm=Vm,Km=Km,Vd=Vd),
                   rtol=1e-8,atol=1e-8) 
      out[-1,3] 
      } 
   } 
   else {
      ## User-supplied function w/o MM elimination w lag time
      defun<- function(time, y, parms) { 
      if(time<=Tlag) {
       dy1dt<-0
       dy2dt<-0
      }
      else {
       dy1dt <- -parms["ka"]*y[1]
       dy2dt <-  parms["ka"]*y[1]/parms["Vd"]-parms["kel"]*y[2]
      }
      list(c(dy1dt,dy2dt)) 
      } 
    
      modfun1 <- function(time,ka,kel,Vd) { 
      out <- lsoda(c(Dose,0),c(0,time),defun,parms=c(ka=ka,kel=kel,Vd=Vd),
                   rtol=1e-5,atol=1e-8) 
      out[-1,3] 
      }  
   }
   ## Select weighting schemes
   file.menu <- c("equal weight", 
                  "1/Cp",
                  "1/Cp^2")           
   pick <- menu(file.menu, title = "<< Weighting Schemes >>")

   with(entertitle(),{

   for( i in 1:length(unique(PKindex$Subject)))  {
      cat("\n\n               << Subject",i,">>\n\n" )  
      objfun <- function(par) {
        if (MMe) {
           out <- modfun2(PKindex$time[PKindex$Subject==i], par[1], par[2],par[3],par[4]) 
        } 
        else  {
           out <- modfun1(PKindex$time[PKindex$Subject==i], par[1], par[2], par[3])
        }
     gift <- which( PKindex$conc[PKindex$Subject==i] != 0 )
     switch(pick,
            sum((PKindex$conc[PKindex$Subject==i][gift]-out[gift])^2),
            sum((PKindex$conc[PKindex$Subject==i][gift]-out[gift])^2/PKindex$conc[gift]),
            sum(((PKindex$conc[PKindex$Subject==i][gift] - out[gift])/PKindex$conc[gift])^2))
     }
     if (MMe) {
          gen<-genoud(objfun,nvars=4,max=FALSE,pop.size=30,max.generations=20,
               wait.generations=10,starting.value=c(par[1,2],par[2,2],par[3,2],par[4,2]),
               BFGS=FALSE,print.level=0,boundary.enforcement=2,
               Domains=matrix(c(0.01,1,1,1,10,100,100,100),4,2),
               MemoryMatrix=TRUE) 
      } 
      else {
          gen<-genoud(objfun,nvars=3,max=FALSE,pop.size=20,max.generations=15,
               wait.generations=10,starting.value=c(par[1,2],par[2,2],par[3,2]),
               BFGS=FALSE,print.level=0,boundary.enforcement=2,
               Domains=matrix(c(0.01,0.01,1,10,1,100),3,2),
               MemoryMatrix=TRUE)
      }
      cat("<< The value of parameter obtained from genetic algorithm >>\n\n")
      if (MMe) {
        namegen<-c("ka","Vm","Km","Vd")
        outgen<-c(gen$par[1],gen$par[2],gen$par[3],gen$par[4])
      } 
      else {
        ## No MM elimination
        namegen<-c("ka","kel","Vd")
        outgen<-c(gen$par[1],gen$par[2],gen$par[3])
      }
      print(data.frame(Parameter=namegen,Value=outgen))  
      F<-objfun(gen$par)
     
      if (MMe) {
        opt<-optim(c(gen$par[1],gen$par[2],gen$par[3],gen$par[4]),objfun, method="Nelder-Mead")
        nameopt<-c("ka","Vm","Km","Vd")
        outopt<-c(opt$par[1],opt$par[2],opt$par[3],opt$par[4])
      }
      else {
        opt<-optim(c(gen$par[1],gen$par[2],gen$par[3]),objfun, method="Nelder-Mead") 
        nameopt<-c("ka","kel","Vd")
        outopt<-c(opt$par[1],opt$par[2],opt$par[3])
      }
      
      cat("\n<< The value of parameter fitted by Nelder-Mead Simplex slgorithm >>\n\n")
      print(data.frame(Parameter=nameopt,Value=outopt))
      cat("\n<< Residual sum-of-squares and parameter values fitted by nls >>\n\n")
      
      if (MMe) {
        fm<-nls(conc ~ modfun2(time, ka, Vm, Km, Vd), data=PKindex,subset=Subject==i,
            start=list(ka=opt$par[1],Vm=opt$par[2],Km=opt$par[3],Vd=opt$par[4]),trace=TRUE,
            nls.control(maxiter=500,tol=1.0))
        cat("\n")
        plotting.non(PKindex, fm, i, pick, xaxis, yaxis)
      } 
      else {
        ## No MM elimination
        fm <-nls(conc ~ modfun1(time, ka, kel, Vd), data=PKindex,subset=Subject==i,
             start=list(ka=opt$par[1],kel=opt$par[2],Vd=opt$par[3]),trace=TRUE,
             nls.control(tol=1))
        cat("\n")
        coef<-data.frame(coef(fm)["kel"])
        plotting.lin(PKindex, fm, i, pick, coef, xaxis, yaxis)
      }
   }
   })
   cat("\n")   
} 

## Legacy function
ffirst.lagm <- function(PKindex,...) {
  ffirst.lag(PKindex,...,Tlag=TRUE,MMe=TRUE)
}

## Legacy function
ffirst.nolagm <- function(PKindex,...) {
  ffirst.lag(PKindex,...,Tlag=FALSE,MMe=TRUE)
}

### Normal fitting
### One compartment PK model extravascualr single dose first-order absorption 
### Without lag time
ffirst.nolag<- function(PKindex,
                        Dose=NULL, 
                        ka=NULL,
                        Vm=NULL,Km=NULL, 
                        Vd=NULL,
                        kel=NULL,        
                        Tlag=FALSE,
                        MMe=TRUE) 
{
   #options(warn=-1)
        
   ## Input dose and Tlag and initial value for ka, kel and Vd

   if (is.null(Dose)) {
     cat("Enter Dose value\n")
     Dose <- scan(nlines=1,quiet=TRUE)
   } 
   else {
     cat("Dose from arguments is = ",Dose,"\n")
   }
   
   
   if (MMe){
      if (is.null(ka) || is.null(Vm) || is.null(Km) || is.null(Vd) ) {
        par<-data.frame(Parameter=c("ka","Vm","Km","Vd"),Initial=c(0))
        par<-edit(par)
        repeat{
           if ( par[1,2] == 0 || par[2,2] ==0 || par[3,2]==0 ||par[4,2]==0){
             cat("\n")
             cat("**********************************\n")
             cat(" Parameter value can not be zero. \n")
             cat(" Press enter to continue.         \n")
             cat("**********************************\n\n")
             readline()
             cat("\n")
             par<-edit(par)}   
           else{
             break
             return(edit(par))}
        } 
        cat("\n")       
        show(par)

      }
   } 
   else{
      if (is.null(ka) || is.null(kel) || is.null(Vd) ) {
        par<-data.frame(Parameter=c("ka","kel","Vd"),Initial=c(0))
        par<-edit(par)
        repeat{
           if ( par[1,2] == 0 || par[2,2] ==0 || par[3,2]==0){
             cat("\n")
             cat("**********************************\n")
             cat(" Parameter value can not be zero. \n")
             cat(" Press enter to continue.         \n")
             cat("**********************************\n\n")
             readline()
             cat("\n")
             par<-edit(par)}   
           else{
             break
             return(edit(par))}
        } 
        cat("\n")       
        show(par)

      }
   }
   
   
   cat("\n")
   
   if ( MMe ) {

      ## User-supplied function with MM elimination w/o lag time
      defun<- function(time, y, parms) { 
      dy1dt <- -parms["ka"] * y[1]
      dy2dt <-  parms["ka"] * y[1]/parms["Vd"] - (parms["Vm"]/parms["Vd"])*y[2]/(parms["Km"]/parms["Vd"]+y[2])
      list(c(dy1dt,dy2dt)) 
      } 
    
      modfun4 <- function(time,ka,Vm,Km,Vd) { 
      out <- lsoda(c(Dose,0),c(0,time),defun,parms=c(ka=ka,Vm=Vm,Km=Km,Vd=Vd),
                   rtol=1e-5,atol=1e-7) 
      out[-1,3] 
      } 
   } 
   else {
      ## User-supplied function w/o MM elimination w/o lag time
      defun <- function(time, y, parms) { 
      dy1dt <- -parms["ka"] * y[1]
      dy2dt <-  parms["ka"] * y[1]/parms["Vd"] - parms["kel"] * y[2]
      list(c(dy1dt,dy2dt)) 
      } 
    
      modfun3 <- function(time,ka,kel,Vd) { 
      out <- lsoda(c(Dose, 0),c(0,time),defun,parms=c(ka=ka,kel=kel,Vd=Vd),
                   rtol=1e-5,atol=1e-7) 
      out[-1,3] 
      }
   }
   ## Select weighting schemes
   file.menu <- c("equal weight", 
                  "1/Cp",
                  "1/Cp^2")           
   pick <- menu(file.menu, title = "<< Weighting Schemes >>")

   with(entertitle(),{

   for( i in 1:length(unique(PKindex$Subject)))  {
      cat("\n\n               << Subject",i,">>\n\n" )  
      objfun <- function(par) {
        if (MMe) {
           out <- modfun4(PKindex$time[PKindex$Subject==i], par[1], par[2],par[3],par[4]) 
        } 
        else  {
           out <- modfun3(PKindex$time[PKindex$Subject==i], par[1], par[2], par[3])
        }
     gift <- which( PKindex$conc[PKindex$Subject==i] != 0 )
     switch(pick,
            sum((PKindex$conc[PKindex$Subject==i][gift]-out[gift])^2),
            sum((PKindex$conc[PKindex$Subject==i][gift]-out[gift])^2/PKindex$conc[gift]),
            sum(((PKindex$conc[PKindex$Subject==i][gift] - out[gift])/PKindex$conc[gift])^2))
     }
     if (MMe) {
          gen<-genoud(objfun,nvars=4,max=FALSE,pop.size=30,max.generations=20,
               wait.generations=10,starting.value=c(par[1,2],par[2,2],par[3,2],par[4,2]),
               BFGS=FALSE,print.level=0,boundary.enforcement=2,
               Domains=matrix(c(0.01,1,1,1,10,100,100,100),4,2),
               MemoryMatrix=TRUE) 
      } 
      else {
          gen<-genoud(objfun,nvars=3,max=FALSE,pop.size=20,max.generations=15,
               wait.generations=10,starting.value=c(par[1,2],par[2,2],par[3,2]),
               BFGS=FALSE,print.level=0,boundary.enforcement=2,
               Domains=matrix(c(0.01,0.01,1,10,1,100),3,2),
               MemoryMatrix=TRUE)
      }
      cat("<< The value of parameter obtained from genetic algorithm >>\n\n")
      if (MMe) {
        namegen<-c("ka","Vm","Km","Vd")
        outgen<-c(gen$par[1],gen$par[2],gen$par[3],gen$par[4])
      } 
      else {
        ## No MM elimination
        namegen<-c("ka","kel","Vd")
        outgen<-c(gen$par[1],gen$par[2],gen$par[3])
      }
      print(data.frame(Parameter=namegen,Value=outgen))  
      F<-objfun(gen$par)
     
      if (MMe) {
        opt<-optim(c(gen$par[1],gen$par[2],gen$par[3],gen$par[4]),objfun, method="Nelder-Mead")
        nameopt<-c("ka","Vm","Km","Vd")
        outopt<-c(opt$par[1],opt$par[2],opt$par[3],opt$par[4])
      }
      else {
        opt<-optim(c(gen$par[1],gen$par[2],gen$par[3]),objfun, method="Nelder-Mead") 
        nameopt<-c("ka","kel","Vd")
        outopt<-c(opt$par[1],opt$par[2],opt$par[3])
      }
      
      cat("\n<< The value of parameter fitted by Nelder-Mead Simplex slgorithm >>\n\n")
      print(data.frame(Parameter=nameopt,Value=outopt))
      cat("\n<< Residual sum-of-squares and parameter values fitted by nls >>\n\n")
      
      if (MMe) {
        fm<-nls(conc ~ modfun4(time, ka, Vm, Km, Vd), data=PKindex,subset=Subject==i,
            start=list(ka=opt$par[1],Vm=opt$par[2],Km=opt$par[3],Vd=opt$par[4]),trace=TRUE,
            nls.control(maxiter=500,tol=1.0))
        cat("\n")
        plotting.non(PKindex, fm, i, pick, xaxis, yaxis)
      } 
      else {
        ## No MM elimination
        fm <-nls(conc ~ modfun3(time, ka, kel, Vd), data=PKindex,subset=Subject==i,
             start=list(ka=opt$par[1],kel=opt$par[2],Vd=opt$par[3]),trace=TRUE,
             nls.control(tol=1))
        cat("\n")
        coef<-data.frame(coef(fm)["kel"])
        plotting.lin(PKindex, fm, i, pick, coef, xaxis, yaxis)
      }
   }
   })
   cat("\n")   
} 





###Simulation
###One compartment PK model extravascualr single dose first-order absorption 
###optional Michaelis-Menten Elimination
###optional lag time
sfirst.nolag <- function(Subject=NULL,  # N Subj's 
                         PKtime=NULL,   # times for sampling
                         Dose=NULL,     # single dose
                         ka=NULL,     
                         Vd=NULL,
                         kel=NULL, ## If not MM elimination
                         MMe=FALSE, ## michaelis-menten elimination?
                         Vm=NULL,Km=NULL,
                         Tlag=FALSE)
{
   #options(warn=-1)
   
   if (is.null(Subject) || !is.integer(Subject)) {
     cat("How many subject do you want?\n")
     Subject<-scan(nlines=1,quiet=TRUE)
   }

   if (is.null(PKtime)) { 
     ## need to verify is a numeric ordered vector.
     PKtime<-data.frame(time=c(0))
     PKtime<-edit(PKtime)
   }
  
   cat("\n")
   show(PKtime)

   if (is.null(Dose)) {
     cat("\nEnter dose value\n")
     Dose<-scan(nlines=1,quiet=TRUE)
   }
  
   if ( Tlag ){
     cat("\nEnter value for lag time\n")
     Tlag<-scan(nlines=1,quiet=TRUE)
   }
  
   if ( MMe ) { 
      if (is.null(ka) || is.null(Vm) || is.null(Km)|| is.null(Vd)) {
         par<-data.frame(Parameter=c("ka","Vm","Km","Vd"),Initial=c(0))
         par<-edit(par)
         repeat{
           if ( par[1,2] == 0 || par[2,2] ==0 || par[3,2]==0 || par[4,2]==0){
             cat("\n")
             cat("**********************************\n")
             cat(" Parameter value can not be zero. \n")
             cat(" Press enter to continue.         \n")
             cat("**********************************\n\n")
             readline()
             cat("\n")
             par<-edit(par)}   
           else{
             break
             return(edit(par))}
        } 
        cat("\n")       
        show(par)
         
         
         cat("\n")
         par1<-par[1,2]
         par2<-par[2,2]
         par3<-par[3,2]
         par4<-par[4,2]
      } 
      else {
         par1 <- ka
         par2 <- Vm
         par3 <- Km
         par4 <- Vd 
      } 
   }
   else{
      if ( is.null(ka) || is.null(kel) || is.null(Vd)){
         par<-data.frame(Parameter=c("ka","kel","Vd"),Initial=c(0))
         par<-edit(par)
         repeat{
           if ( par[1,2] == 0 || par[2,2] ==0 || par[3,2]==0){
             cat("\n")
             cat("**********************************\n")
             cat(" Parameter value can not be zero. \n")
             cat(" Press enter to continue.         \n")
             cat("**********************************\n\n")
             readline()
             cat("\n")
             par<-edit(par)}   
           else{
             break
             return(edit(par))}
        } 
        cat("\n")       
        show(par)
         
         
         cat("\n")
         par1<-par[1,2]
         par2<-par[2,2]
         par3<-par[3,2]
      } 
      else {
         par1 <- ka
         par2 <- kel
         par3 <- Vd
      }
   }   
  
   if (MMe ) {
      if (Tlag){
        ## User-supplied function w Michaelis-Mention elimination & w lag time
        defun<- function(time, y, parms) { 
           if(time <= Tlag) {
              dy1dt<-0
              dy2dt<-0
           }
           else {
              dy1dt <- -parms["ka"] * y[1]
              dy2dt <-  parms["ka"] * y[1]/parms["Vd"] - (parms["Vm"]/parms["Vd"])*y[2]/(parms["Km"]/parms["Vd"]+y[2])
           }
           list(c(dy1dt,dy2dt)) 
         } 
      }
      else{
        ## User-supplied function with MM elimination w/o lag time
        defun<- function(time, y, parms) { 
          dy1dt <- -parms["ka"] * y[1]
          dy2dt <-  parms["ka"] * y[1]/parms["Vd"] - (parms["Vm"]/parms["Vd"])*y[2]/(parms["Km"]/parms["Vd"]+y[2])
          list(c(dy1dt,dy2dt)) 
        } 
      }  
   } 
   else {
      if (Tlag){
        ## User-supplied function w/o MM elimination w lag time
        defun<- function(time, y, parms) { 
           if(time<=Tlag) {
              dy1dt<-0
              dy2dt<-0
           }
           else {
              dy1dt <- -parms["ka"]*y[1]
              dy2dt <-  parms["ka"]*y[1]/parms["Vd"]-parms["kel"]*y[2]
           }
           list(c(dy1dt,dy2dt)) 
         } 
      }
      else{
        ## User-supplied function w/o MM elimination w/o lag time
        defun <- function(time, y, parms) { 
           dy1dt <- -parms["ka"] * y[1]
           dy2dt <-  parms["ka"] * y[1]/parms["Vd"] - parms["kel"] * y[2]
           list(c(dy1dt,dy2dt)) 
        } 
      }  
   }
   
   file.menu <- c("Simulation with Error",
                  "Monte Carlo Simulation")
   pick <- menu(file.menu, title = "<< Simulation Type >>")
   if (pick ==1){
      cat("\n\n")
      file.menu <- c("No Error",
                     "Error=Normal Error", 
                     "Error=Uniform Error",
                     "Error=Normal Error*True Value",
                     "Error=Uniform Error*True Value")
      pick <- menu(file.menu, title = "<< Error Type >>")   
      
      type<-switch(pick, 
                  "No Error",
                  "Error=Normal Error", 
                  "Error=Uniform Error",
                  "Error=Normal Error*True Value",
                  "Error=Uniform Error*True Value")
      
      if (pick ==1){
          PKindex<-vector(Subject,mode="list")
          for( i in 1:Subject)  {
             cat("\n\n             << Subject",i,">>\n\n" )   
           
             if ( MMe ){
                ka<-par1
                Vm<-par2  
                Km<-par3
                Vd<-par4 
                PKindex[[i]]<-sfirst.mm.out(PKtime,ka,Vm,Km,Vd,defun,par1,par2,par3,par4,Dose,i,type,Tlag)      
             } 
             else{
                ka<-par1
                kel<-par2
                Vd<-par3       
                PKindex[[i]]<-sfirst1.out(PKtime,ka,kel,Vd,defun,par1,par2,par3,Dose,i,type,Tlag)   
             }
          }       
          PKindex<- as.data.frame(do.call("rbind",PKindex))
          rownames(PKindex) <- seq(nrow(PKindex)) 
          savefile(PKindex)
       }   
      else {
          if ( !MMe ){
             cat("\n\nEnter error factor for ka\n")
             factor1<-scan(nlines=1,quiet=TRUE)
             repeat{
              if ( factor1 == 0 ){
                cat("\n")
                cat("**********************************\n")
                cat(" Parameter value can not be zero. \n")
                cat(" Press enter to continue.         \n")
                cat("**********************************\n\n")
                readline()
                cat("\n")
                factor1<-scan(nlines=1,quiet=TRUE)}   
           else{
                break
                return(factor1)}
           } 
             
             cat("\nEnter error factor for kel\n")
             factor2<-scan(nlines=1,quiet=TRUE)
             repeat{
              if ( factor2 == 0 ){
                cat("\n")
                cat("**********************************\n")
                cat(" Parameter value can not be zero. \n")
                cat(" Press enter to continue.         \n")
                cat("**********************************\n\n")
                readline()
                cat("\n")
                factor2<-scan(nlines=1,quiet=TRUE)}   
           else{
                break
                return(factor2)}
           } 
             
             cat("\nEnter error factor for Vd\n")
             factor3<-scan(nlines=1,quiet=TRUE)
             repeat{
              if ( factor3 == 0 ){
                cat("\n")
                cat("**********************************\n")
                cat(" Parameter value can not be zero. \n")
                cat(" Press enter to continue.         \n")
                cat("**********************************\n\n")
                readline()
                cat("\n")
                factor3<-scan(nlines=1,quiet=TRUE)}   
           else{
                break
                return(factor3)}
           } 
           
             PKindex<-vector(Subject,mode="list")
             for( i in 1:Subject)  {
                cat("\n\n             << Subject",i,">>\n\n" )  
                switch(pick-1,
                      {ka<-par1+rnorm(1,mean=0,sd=factor1)
                       while(ka<=0){
                          ka<-par1+rnorm(1,mean=0,sd=factor1)}
                       kel<-par2+rnorm(1,mean=0,sd=factor2)
                       while(kel<=0){
                          kel<-par2+rnorm(1,mean=0,sd=factor2)}
                       Vd<-par3+rnorm(1,mean=0,sd=factor3)
                       while(Vd<=0){
                          Vd<-par3+rnorm(1,mean=0,sd=factor3)}},
                  
                      {ka<-par1+runif(1,min=-factor1,max=factor1)
                       while(ka<=0){
                          ka<-par1+runif(1,min=-factor1,max=factor1)}
                       kel<-par2+runif(1,min=-factor2,max=factor2)
                       while(kel<=0){
                          kel<-par2+runif(1,min=-factor2,max=factor2)}
                       Vd<-par3+runif(1,min=-factor3,max=factor3)
                       while(Vd<=0){
                          Vd<-par3+runif(1,min=-factor3,max=factor3)}}, 
                    
                      {ka<-par1*rnorm(1,mean=0,sd=factor1)+par1
                       while(ka<=0){
                          ka<-par1*rnorm(1,mean=0,sd=factor1)+par1}
                       kel<-par2*rnorm(1,mean=0,sd=factor2)+par2
                       while(kel<=0){
                          kel<-par2*rnorm(1,mean=0,sd=factor2)+par2}
                       Vd<-par3*rnorm(1,mean=0,sd=factor3)+par3
                       while(Vd<=0){
                          Vd<-par3*rnorm(1,mean=0,sd=factor3)+par3}},
                  
                      {ka<-par1*runif(1,min=-factor1,max=factor1)+par1
                       while(ka<=0){
                          ka<-par1*runif(1,min=-factor1,max=factor1)+par1}
                       kel<-par2*runif(1,min=-factor2,max=factor2)+par2
                       while(kel<=0){
                          kel<-par2*runif(1,min=-factor2,max=factor2)+par2}
                       Vd<-par3*runif(1,min=-factor3,max=factor3)+par3
                       while(Vd<=0){
                          Vd<-par3*runif(1,min=-factor3,max=factor3)+par3}}   
                )      
             PKindex[[i]]<-sfirst1.out(PKtime,ka,kel,Vd,defun,par1,par2,par3,Dose,i,type,Tlag)      
             }
          }
          else{
             cat("\n\nEnter error factor for ka\n")
             factor1<-scan(nlines=1,quiet=TRUE)
             repeat{
              if ( factor1 == 0 ){
                cat("\n")
                cat("**********************************\n")
                cat(" Parameter value can not be zero. \n")
                cat(" Press enter to continue.         \n")
                cat("**********************************\n\n")
                readline()
                cat("\n")
                factor1<-scan(nlines=1,quiet=TRUE)}   
           else{
                break
                return(factor1)}
           } 
             
             cat("\nEnter error factor for Vm\n")
             factor2<-scan(nlines=1,quiet=TRUE)
             repeat{
              if ( factor2 == 0 ){
                cat("\n")
                cat("**********************************\n")
                cat(" Parameter value can not be zero. \n")
                cat(" Press enter to continue.         \n")
                cat("**********************************\n\n")
                readline()
                cat("\n")
                factor2<-scan(nlines=1,quiet=TRUE)}   
           else{
                break
                return(factor2)}
           } 
             
             cat("\nEnter error factor for Km\n")
             factor3<-scan(nlines=1,quiet=TRUE)
             repeat{
              if ( factor3 == 0 ){
                cat("\n")
                cat("**********************************\n")
                cat(" Parameter value can not be zero. \n")
                cat(" Press enter to continue.         \n")
                cat("**********************************\n\n")
                readline()
                cat("\n")
                factor3<-scan(nlines=1,quiet=TRUE)}   
           else{
                break
                return(factor3)}
           } 
             
             cat("\nEnter error factor for Vd\n")
             factor4<-scan(nlines=1,quiet=TRUE)
             repeat{
              if ( factor4 == 0 ){
                cat("\n")
                cat("**********************************\n")
                cat(" Parameter value can not be zero. \n")
                cat(" Press enter to continue.         \n")
                cat("**********************************\n\n")
                readline()
                cat("\n")
                factor4<-scan(nlines=1,quiet=TRUE)}   
           else{
                break
                return(factor4)}
           } 
             
             PKindex<-vector(Subject,mode="list")
             for( i in 1:Subject)  {
               cat("\n\n             << Subject",i,">>\n\n" )  
               switch(pick-1,
                     {ka<-par1+rnorm(1,mean=0,sd=factor1)
                      while(ka<=0){
                         ka<-par1+rnorm(1,mean=0,sd=factor1)}
                      Vm<-par2+rnorm(1,mean=0,sd=factor2)
                      while(Vm<=0){
                         Vm<-par2+rnorm(1,mean=0,sd=factor2)}
                      Km<-par3+rnorm(1,mean=0,sd=factor3)
                      while(Km<=0){
                         Km<-par3+rnorm(1,mean=0,sd=factor3)}
                      Vd<-par4+rnorm(1,mean=0,sd=factor4)
                      while(Vd<=0){
                         Vd<-par4+rnorm(1,mean=0,sd=factor4)}},
                  
                     {ka<-par1+runif(1,min=-factor1,max=factor1)
                      while(ka<=0){
                         ka<-par1+runif(1,min=-factor1,max=factor1)}
                      Vm<-par2+runif(1,min=-factor2,max=factor2)
                      while(Vm<=0){
                         Vm<-par2+runif(1,min=-factor2,max=factor2)}
                      Km<-par3+runif(1,min=-factor3,max=factor3)
                      while(Km<=0){
                         Km<-par3+runif(1,min=-factor3,max=factor3)}
                      Vd<-par4+runif(1,min=-factor4,max=factor4)
                      while(Vd<=0){
                         Vd<-par4+runif(1,min=-factor4,max=factor4)}}, 
                    
                     {ka<-par1*rnorm(1,mean=0,sd=factor1)+par1
                      while(ka<=0){
                         ka<-par1*rnorm(1,mean=0,sd=factor1)+par1}
                      Vm<-par2*rnorm(1,mean=0,sd=factor2)+par2
                      while(Vm<=0){
                         Vm<-par2*rnorm(1,mean=0,sd=factor2)+par2}
                      Km<-par3*rnorm(1,mean=0,sd=factor3)+par3
                      while(Km<=0){
                         Km<-par3*rnorm(1,mean=0,sd=factor3)+par3}
                      Vd<-par4*rnorm(1,mean=0,sd=factor4)+par4
                      while(Vd<=0){
                         Vd<-par4*rnorm(1,mean=0,sd=factor4)+par4}},
                  
                     {ka<-par1*runif(1,min=-factor1,max=factor1)+par1
                      while(ka<=0){
                         ka<-par1*runif(1,min=-factor1,max=factor1)+par1}
                      Vm<-par2*runif(1,min=-factor2,max=factor2)+par2
                      while(Vm<=0){
                         Vm<-par2*runif(1,min=-factor2,max=factor2)+par2}
                      Km<-par3*runif(1,min=-factor3,max=factor3)+par3
                      while(Km<=0){
                         Km<-par3*runif(1,min=-factor3,max=factor3)+par3}
                      Vd<-par4*runif(1,min=-factor4,max=factor4)+par4
                      while(Vd<=0){
                         Vd<-par4*runif(1,min=-factor4,max=factor4)+par4}}   
                )      
                PKindex[[i]]<-sfirst.mm.out(PKtime,ka,Vm,Km,Vd,defun,par1,par2,par3,par4,Dose,i,type,Tlag)       
             } 
         }    
         PKindex<- as.data.frame(do.call("rbind",PKindex))
         rownames(PKindex) <- seq(nrow(PKindex)) 
         savefile(PKindex)
     } 
  } 
  else if (pick ==2){ 
     cat("\n\n")
     file.menu <- c("Error=Normal Error", 
                    "Error=Uniform Error",
                    "Error=Normal Error*True Value",
                    "Error=Uniform Error*True Value")
     pick <- menu(file.menu, title = "<< Error Type >>")
     
     type<-switch(pick, 
                  "Error=Normal Error", 
                  "Error=Uniform Error",
                  "Error=Normal Error*True Value",
                  "Error=Uniform Error*True Value")
     
     cat("\n\nHow many times do you want to iteration ?\n")
     re<-scan(nlines=1,quiet=TRUE)
     
       if (MMe ){
         if (Tlag){
          cat("\nEnter error factor for ka\n")
          factor1<-scan(nlines=1,quiet=TRUE)
          repeat{
              if ( factor1 == 0 ){
                cat("\n")
                cat("**********************************\n")
                cat(" Parameter value can not be zero. \n")
                cat(" Press enter to continue.         \n")
                cat("**********************************\n\n")
                readline()
                cat("\n")
                factor1<-scan(nlines=1,quiet=TRUE)}   
           else{
                break
                return(factor1)}
           } 
          
          cat("\nEnter error factor for Vm\n")
          factor2<-scan(nlines=1,quiet=TRUE)
          repeat{
              if ( factor2 == 0 ){
                cat("\n")
                cat("**********************************\n")
                cat(" Parameter value can not be zero. \n")
                cat(" Press enter to continue.         \n")
                cat("**********************************\n\n")
                readline()
                cat("\n")
                factor2<-scan(nlines=1,quiet=TRUE)}   
           else{
                break
                return(factor2)}
           } 
          
          cat("\nEnter error factor for Km\n")
          factor3<-scan(nlines=1,quiet=TRUE)
          repeat{
              if ( factor3 == 0 ){
                cat("\n")
                cat("**********************************\n")
                cat(" Parameter value can not be zero. \n")
                cat(" Press enter to continue.         \n")
                cat("**********************************\n\n")
                readline()
                cat("\n")
                factor3<-scan(nlines=1,quiet=TRUE)}   
           else{
                break
                return(factor3)}
           } 
          
          cat("\nEnter error factor for Vd\n")
          factor4<-scan(nlines=1,quiet=TRUE)
          repeat{
              if ( factor4 == 0 ){
                cat("\n")
                cat("**********************************\n")
                cat(" Parameter value can not be zero. \n")
                cat(" Press enter to continue.         \n")
                cat("**********************************\n\n")
                readline()
                cat("\n")
                factor4<-scan(nlines=1,quiet=TRUE)}   
           else{
                break
                return(factor4)}
           } 
          
          cat("\n")
          cat("************************************************************************\n")
          cat("Summary Table                                                           \n")
          cat("Model: 1-Compartment, Extravascular, Single-Dose, 1-Ordered Absorption, \n")
          cat("       and Michaelis-Menten Elimination with Lag Time Model             \n")
          cat("Subject #:", Subject,"                                                  \n")
          cat("Error Type:", type,"                                                    \n")
          cat("Simulation #:", re,"                                                    \n\n")
          sim<-matrix(c(par1,par2,par3,par4,factor1,factor2,factor3,factor4),4,2)
          dimnames(sim)<-list(c("ka","Vm","Km","Vd"),c("Original","Error factor"))
          show(sim)   
          cat("************************************************************************\n\n")
         }
         else{
          cat("\nEnter error factor for ka\n")
          factor1<-scan(nlines=1,quiet=TRUE)
          cat("\nEnter error factor for Vm\n")
          factor2<-scan(nlines=1,quiet=TRUE)
          cat("\nEnter error factor for Km\n")
          factor3<-scan(nlines=1,quiet=TRUE)
          cat("\nEnter error factor for Vd\n")
          factor4<-scan(nlines=1,quiet=TRUE)
          cat("\n")
          cat("*************************************************************************\n")
          cat("Summary Table                                                            \n")
          cat("Model: 1-Compartment, Extravascular, Single-Dose, 1t-Ordered Absorption, \n")
          cat("       and Michaelis-Menten Elimination without Lag Time Model           \n")
          cat("Subject #:", Subject,"                                                   \n")
          cat("Error Type:", type,"                                                     \n")
          cat("Simulation #:", re,"                                                     \n\n")
          sim<-matrix(c(par1,par2,par3,par4,factor1,factor2,factor3,factor4),4,2)
          dimnames(sim)<-list(c("ka","Vm","Km","Vd"),c("Original","Error factor"))
          show(sim)   
          cat("*************************************************************************\n\n")
         }
       }
       else{
         if (Tlag){
         cat("\nEnter error factor for ka\n")
         factor1<-scan(nlines=1,quiet=TRUE)
         cat("\nEnter error factor for kel\n")
         factor2<-scan(nlines=1,quiet=TRUE)
         cat("\nEnter error factor for Vd\n")
         factor3<-scan(nlines=1,quiet=TRUE)
         cat("\n")
         cat("****************************************************\n")
         cat("Summary Table                                       \n")
         cat("Model: 1-Compartment, Extravascular, Single-Dose,   \n")
         cat("       and 1-Ordered Absorption with Lag Time Model \n")
         cat("Subject #:", Subject,"                              \n")
         cat("Error Type:", type,"                                \n")
         cat("Simulation #:", re,"                                \n\n")
         sim<-matrix(c(par1,par2,par3,factor1,factor2,factor3),3,2)
         dimnames(sim)<-list(c("ka","kel","Vd"),c("Original","Error factor"))
         show(sim)   
         cat("****************************************************\n\n")
         }
         else{
         cat("\nEnter error factor for ka\n")
         factor1<-scan(nlines=1,quiet=TRUE)
         cat("\nEnter error factor for kel\n")
         factor2<-scan(nlines=1,quiet=TRUE)
         cat("\nEnter error factor for Vd\n")
         factor3<-scan(nlines=1,quiet=TRUE)
         cat("\n")
         cat("*******************************************************\n")
         cat("Summary Table                                          \n")
         cat("Model: 1-Compartment, Extravascular, Single-Dose,      \n")
         cat("       and 1-Ordered Absorption without Lag Time Model \n")
         cat("Subject #:", Subject,"                                 \n")
         cat("Error Type:", type,"                                   \n")
         cat("Simulation #:", re,"                                   \n\n")
         sim<-matrix(c(par1,par2,par3,factor1,factor2,factor3),3,2)
         dimnames(sim)<-list(c("ka","kel","Vd"),c("Original","Error factor"))
         show(sim)   
         cat("*******************************************************\n\n")
         }
       }

     
     PKindex<-vector(Subject,mode="list")
     for( i in 1:Subject)  {
       cat("\n\n             << Subject",i,">>\n\n" ) 
       C1.lsoda<-list()

         for (j in 1:re){
     
           if ( !MMe ){
              switch(pick, 
                    {ka<-par1+rnorm(1,mean=0,sd=factor1)
                    while(ka<=0){
                       ka<-par1+rnorm(1,mean=0,sd=factor1)}
                    kel<-par2+rnorm(1,mean=0,sd=factor2)
                    while(kel<=0){
                       kel<-par2+rnorm(1,mean=0,sd=factor2)}
                    Vd<-par3+rnorm(1,mean=0,sd=factor3)
                    while(Vd<=0){
                       Vd<-par3+rnorm(1,mean=0,sd=factor3)}},
                  
                   {ka<-par1+runif(1,min=-factor1,max=factor1)
                    while(ka<=0){
                       ka<-par1+runif(1,min=-factor1,max=factor1)}
                    kel<-par2+runif(1,min=-factor2,max=factor2)
                    while(kel<=0){
                       kel<-par2+runif(1,min=-factor2,max=factor2)}
                    Vd<-par3+runif(1,min=-factor3,max=factor3)
                    while(Vd<=0){
                       Vd<-par3+runif(1,min=-factor3,max=factor3)}}, 
                    
                   {ka<-par1*rnorm(1,mean=0,sd=factor1)+par1
                    while(ka<=0){
                       ka<-par1*rnorm(1,mean=0,sd=factor1)+par1}
                    kel<-par2*rnorm(1,mean=0,sd=factor2)+par2
                    while(kel<=0){
                       kel<-par2*rnorm(1,mean=0,sd=factor2)+par2}
                    Vd<-par3*rnorm(1,mean=0,sd=factor3)+par3
                    while(Vd<=0){
                       Vd<-par3*rnorm(1,mean=0,sd=factor3)+par3}},
                  
                   {ka<-par1*runif(1,min=-factor1,max=factor1)+par1
                    while(ka<=0){
                       ka<-par1*runif(1,min=-factor1,max=factor1)+par1}
                    kel<-par2*runif(1,min=-factor2,max=factor2)+par2
                    while(kel<=0){
                      kel<-par2*runif(1,min=-factor2,max=factor2)+par2}
                    Vd<-par3*runif(1,min=-factor3,max=factor3)+par3
                    while(Vd<=0){
                       Vd<-par3*runif(1,min=-factor3,max=factor3)+par3}}   
                    )    
              time1<-PKtime$time
              parms<-c(ka=ka,kel=kel,Vd=Vd)  
              XX<-data.frame(lsoda(c(Dose,0),c(0,time1), defun, parms))  
              good<-ifelse(XX[2:(length(time1)+1),3]<=1e-5,
                           0,
                           XX[2:(length(time1)+1),3])
              C1.lsoda[[j]]<-data.frame(XX[2:(length(time1)+1),1],good)
              colnames(C1.lsoda[[j]])<-list("time","concentration") 
           }
           else{
              switch(pick, 
                    {ka<-par1+rnorm(1,mean=0,sd=factor1)
                    while(ka<=0){
                       ka<-par1+rnorm(1,mean=0,sd=factor1)}
                    Vm<-par2+rnorm(1,mean=0,sd=factor2)
                    while(Vm<=0){
                       Vm<-par2+rnorm(1,mean=0,sd=factor2)}
                    Km<-par3+rnorm(1,mean=0,sd=factor3)
                    while(Km<=0){
                       Km<-par3+rnorm(1,mean=0,sd=factor3)}
                    Vd<-par4+rnorm(1,mean=0,sd=factor4)
                    while(Vd<=0){
                       Vd<-par4+rnorm(1,mean=0,sd=factor4)}},
                  
                   {ka<-par1+runif(1,min=-factor1,max=factor1)
                    while(ka<=0){
                       ka<-par1+runif(1,min=-factor1,max=factor1)}
                    Vm<-par2+runif(1,min=-factor2,max=factor2)
                    while(Vm<=0){
                       Vm<-par2+runif(1,min=-factor2,max=factor2)}
                    Km<-par3+runif(1,min=-factor3,max=factor3)
                    while(Km<=0){
                       Km<-par3+runif(1,min=-factor3,max=factor3)}
                    Vd<-par4+runif(1,min=-factor4,max=factor4)
                    while(Vd<=0){
                      Vd<-par4+runif(1,min=-factor4,max=factor4)}}, 
                    
                   {ka<-par1*rnorm(1,mean=0,sd=factor1)+par1
                    while(ka<=0){
                       ka<-par1*rnorm(1,mean=0,sd=factor1)+par1}
                    Vm<-par2*rnorm(1,mean=0,sd=factor2)+par2
                    while(Vm<=0){
                       Vm<-par2*rnorm(1,mean=0,sd=factor2)+par2}
                    Km<-par3*rnorm(1,mean=0,sd=factor3)+par3
                    while(Km<=0){
                       Km<-par3*rnorm(1,mean=0,sd=factor3)+par3}
                    Vd<-par4*rnorm(1,mean=0,sd=factor4)+par4
                    while(Vd<=0){
                       Vd<-par4*rnorm(1,mean=0,sd=factor4)+par4}},
                  
                   {ka<-par1*runif(1,min=-factor1,max=factor1)+par1
                    while(ka<=0){
                       ka<-par1*runif(1,min=-factor1,max=factor1)+par1}
                    Vm<-par2*runif(1,min=-factor2,max=factor2)+par2
                    while(Vm<=0){
                       Vm<-par2*runif(1,min=-factor2,max=factor2)+par2}
                    Km<-par3*runif(1,min=-factor3,max=factor3)+par3
                    while(Km<=0){
                       Km<-par3*runif(1,min=-factor3,max=factor3)+par3}
                    Vd<-par4*runif(1,min=-factor4,max=factor4)+par4
                    while(Vd<=0){
                       Vd<-par4*runif(1,min=-factor4,max=factor4)+par4}}  
                    )
              time1<-PKtime$time
              parms<-c(ka=ka,Vm=Vm,Km=Km,Vd=Vd) 
              XX<-data.frame(lsoda(c(Dose,0),c(0,time1), defun, parms))  
              good<-ifelse(XX[2:(length(time1)+1),3]<=1e-5,
                           0,
                           XX[2:(length(time1)+1),3])
              C1.lsoda[[j]]<-data.frame(XX[2:(length(time1)+1),1],good)
              colnames(C1.lsoda[[j]])<-list("time","concentration") 
           } 
         }        
         PKindex[[i]]<-montecarlo(C1.lsoda,time1,i,re) 
     }  
     PKindex<-as.data.frame(do.call("rbind",PKindex))
     rownames(PKindex)<-seq(nrow(PKindex)) 
     savefile(PKindex) 
  }
} 



## Legacy function
sfirst.lag <- function(...) {
  sfirst.nolag(...,Tlag=TRUE,MMe=FALSE)
}


## Legacy function
sfirst.lagm <- function(...) {
  sfirst.nolag(...,Tlag=TRUE,MMe=TRUE)
}


## Legacy function
sfirst.nolagm <- function(...) {
  sfirst.nolag(...,Tlag=FALSE,MMe=TRUE)
}
