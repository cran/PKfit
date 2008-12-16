### PKindex is usually the Dataset.


### Normal fitting
### Two compartment PK model iv bolus single dose
fbolus2<- function(PKindex,
                   Dose=NULL, 
                   kel=NULL,
                   k12=NULL,  
                   k21=NULL,      
                   Vd=NULL) 
{
   #options(warn=-1)
   
   ## Input dose and initial value for kel, k12, k21 and Vd

   if (is.null(Dose)) {
     cat("Enter Dose value\n")
     Dose <- scan(nlines=1,quiet=TRUE)
   } 
   else {
     cat("Dose from arguments is = ",Dose,"\n")
   }
   
   if (is.null(kel) || is.null(k12) || is.null(k21) || is.null(Vd) ) {
        par<-data.frame(Parameter=c("kel","k12","k21","Vd"),Initial=c(0))
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
   }
   
   cat("\n")
   
   defun<- function(time, y, parms) { 
     dCp1dt <- -parms["kel"]*y[1]-parms["k12"]*y[1]+parms["k21"]*y[2] 
     dCp2dt <-  parms["k12"]*y[1]-parms["k21"]*y[2]
     list(c(dCp1dt,dCp2dt)) 
   } 
    
   modfun1 <- function(time,kel,k12,k21, Vd) { 
     out <- lsoda(c(Dose/Vd,0),c(0,time),defun,parms=c(kel=kel,k12=k12,k21=k21,Vd=Vd),
                  rtol=1e-5,atol=1e-5) 
     out[-1,2] 
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
        out <- modfun1(PKindex$time[PKindex$Subject==i], par[1], par[2], par[3], par[4])
        gift <- which( PKindex$conc[PKindex$Subject==i] != 0 )
        switch(pick,
               sum((PKindex$conc[PKindex$Subject==i][gift]-out[gift])^2),
               sum((PKindex$conc[PKindex$Subject==i][gift]-out[gift])^2/PKindex$conc[gift]),
               sum(((PKindex$conc[PKindex$Subject==i][gift] - out[gift])/PKindex$conc[gift])^2))
     }
     
     gen<-genoud(objfun,nvars=4,max=FALSE,pop.size=30,max.generations=20,
          wait.generations=10,
          starting.value=c(par[1,2],par[2,2],par[3,2],par[4,2]),
          BFGS=FALSE,print.level=0,boundary.enforcement=2,
          Domains=matrix(c(0.01,0.01,0.01,1,1,10,1,100),4,2),
          MemoryMatrix=TRUE)    
          
     cat("<< The value of parameter obtained from genetic algorithm >>\n\n")
     
     namegen<-c("kel","k12","k21","Vd")
     outgen<-c(gen$par[1],gen$par[2],gen$par[3],gen$par[4])
     print(data.frame(Parameter=namegen,Value=outgen))  
     F<-objfun(gen$par)
     
     opt <- optim(c(gen$par[1],gen$par[2],gen$par[3],gen$par[4]),objfun,method="Nelder-Mead")   
     nameopt<-c("kel","k12","k21","Vd")
     outopt<-c(opt$par[1],opt$par[2],opt$par[3],opt$par[4])
     
     cat("\n<< The value of parameter fitted by Nelder-Mead Simplex slgorithm >>\n\n")
     print(data.frame(Parameter=nameopt,Value=outopt))
     cat("\n<< Residual sum-of-squares and parameter values fitted by nls >>\n\n")
     
     fm<-nls(conc ~ modfun1(time,kel,k12,k21,Vd), data=PKindex,subset=Subject==i,
         start=list(kel=opt$par[1],k12=opt$par[2],k21=opt$par[3],Vd=opt$par[4]),trace=TRUE,
         nls.control(maxiter=50,tol=1,minFactor=1/1024))
     cat("\n")
     coef<-data.frame(coef(fm)["kel"])     
     plotting.lin(PKindex, fm, i, pick, coef, xaxis, yaxis)
     }
  })
  cat("\n")   
}


### Two compartment PK model iv infusion single dose    
finfu2<- function(PKindex,
                  Tinf=NULL,
                  Dose=NULL, 
                  kel=NULL,
                  k12=NULL,  
                  k21=NULL,      
                  Vd=NULL) 
{
   #options(warn=-1)
   
   ## Input dose and Tinf and initial value for kel, k12, k21 and Vd

   if (is.null(Dose)) {
     cat("Enter Dose value\n")
     Dose <- scan(nlines=1,quiet=TRUE)
   } 
   else {
     cat("Dose from arguments is = ",Dose,"\n")
   }
   
   if (is.null(Tinf)) {
     cat("\nEnter infusion duration\n")
     Tinf<-scan(nlines=1,quiet=TRUE)
     cat("\n")
   } 
   else {
     cat("Tinf from arguments is = ",Tinf,"\n")
   }
   
   if (is.null(kel) || is.null(k12) || is.null(k21) || is.null(Vd) ) {
        par<-data.frame(Parameter=c("kel","k12","k21","Vd"),Initial=c(0))
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
   }
   
   cat("\n")
   
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
      
   modfun2 <- function(time,kel,k12,k21,Vd) { 
      out <- lsoda(c(0,0),c(0,time),defun,parms=c(kel=kel,k12=k12,k21=k21,Vd=Vd),
                   rtol=1e-8,atol=1e-8) 
      out[-1,2] 
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
        out <- modfun2(PKindex$time[PKindex$Subject==i], par[1], par[2], par[3], par[4])
        gift <- which( PKindex$conc[PKindex$Subject==i] != 0 )
        switch(pick,
               sum((PKindex$conc[PKindex$Subject==i][gift]-out[gift])^2),
               sum((PKindex$conc[PKindex$Subject==i][gift]-out[gift])^2/PKindex$conc[gift]),
               sum(((PKindex$conc[PKindex$Subject==i][gift] - out[gift])/PKindex$conc[gift])^2))
     }
     
     gen<-genoud(objfun,nvars=4,max=FALSE,pop.size=30,max.generations=20,
          wait.generations=10,
          starting.value=c(par[1,2],par[2,2],par[3,2],par[4,2]),
          BFGS=FALSE,print.level=0,boundary.enforcement=2,
          Domains=matrix(c(0.01,0.01,0.01,1,1,10,1,100),4,2),
          MemoryMatrix=TRUE)    
          
     cat("<< The value of parameter obtained from genetic algorithm >>\n\n")
     
     namegen<-c("kel","k12","k21","Vd")
     outgen<-c(gen$par[1],gen$par[2],gen$par[3],gen$par[4])
     print(data.frame(Parameter=namegen,Value=outgen))  
     F<-objfun(gen$par)
     
     opt <- optim(c(gen$par[1],gen$par[2],gen$par[3],gen$par[4]),objfun,method="Nelder-Mead")   
     nameopt<-c("kel","k12","k21","Vd")
     outopt<-c(opt$par[1],opt$par[2],opt$par[3],opt$par[4])
     
     cat("\n<< The value of parameter fitted by Nelder-Mead Simplex slgorithm >>\n\n")
     print(data.frame(Parameter=nameopt,Value=outopt))
     cat("\n<< Residual sum-of-squares and parameter values fitted by nls >>\n\n")
     
     fm<-nls(conc ~ modfun2(time,kel,k12,k21,Vd), data=PKindex,subset=Subject==i,
         start=list(kel=opt$par[1],k12=opt$par[2],k21=opt$par[3],Vd=opt$par[4]),trace=TRUE,
         nls.control(maxiter=50,tol=1,minFactor=1/1024))
     cat("\n")
     coef<-data.frame(coef(fm)["kel"])     
     plotting.lin(PKindex, fm, i, pick, coef, xaxis, yaxis)
     }
  })
  cat("\n")   
}


### Two compartment PK model extravascular single dose first order absorption
ffirst2<- function(PKindex,
                  ka=NULL,
                  Dose=NULL, 
                  kel=NULL,
                  k12=NULL,  
                  k21=NULL,      
                  Vd=NULL) 
{
   #options(warn=-1)
   
   ## Input dose and initial value for ka, kel, k12, k21 and Vd

   if (is.null(Dose)) {
     cat("Enter Dose value\n")
     Dose <- scan(nlines=1,quiet=TRUE)
   } 
   else {
     cat("Dose from arguments is = ",Dose,"\n")
   }
   
   if (is.null(ka) || is.null(kel) || is.null(k12) || is.null(k21) || is.null(Vd) ) {
        par<-data.frame(Parameter=c("ka","kel","k12","k21","Vd"),Initial=c(0))
        par<-edit(par)
        repeat{
           if ( par[1,2] == 0 || par[2,2] ==0 || par[3,2]==0 || par[4,2]==0 || par[5,2]==0){
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
   
   cat("\n")
   
   defun<- function(time, y, parms) { 
     dCp1dt <- -parms["ka"]*y[1]
     dCp2dt <- parms["ka"]*y[1]/parms["Vd"]-parms["kel"]*y[2]+parms["k21"]*y[3]-parms["k12"]*y[2]
     dCp3dt <- parms["k12"]*y[2]-parms["k21"]*y[3]
     list(c(dCp1dt,dCp2dt,dCp3dt)) 
   } 
    
   modfun3 <- function(time,ka,kel,k12,k21,Vd) { 
     out <- lsoda(c(Dose,0,0),c(0,time),defun,parms=c(ka=ka,kel=kel,k12=k12,k21=k21,Vd=Vd),
                  rtol=1e-6,atol=1e-8) 
     out[-1,3] 
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
        out <- modfun3(PKindex$time[PKindex$Subject==i],par[1],par[2],par[3],par[4],par[5])
        gift <- which( PKindex$conc[PKindex$Subject==i] != 0 )
        switch(pick,
               sum((PKindex$conc[PKindex$Subject==i][gift]-out[gift])^2),
               sum((PKindex$conc[PKindex$Subject==i][gift]-out[gift])^2/PKindex$conc[gift]),
               sum(((PKindex$conc[PKindex$Subject==i][gift] - out[gift])/PKindex$conc[gift])^2))
     }
     
     gen<-genoud(objfun,nvars=5,max=FALSE,pop.size=15,max.generations=10,
          wait.generations=5,
          starting.value=c(par[1,2],par[2,2],par[3,2],par[4,2],par[5,2]),
          BFGS=FALSE,print.level=0,boundary.enforcement=2,
          Domains=matrix(c(0.01,0.01,0.01,0.01,1,10,1,10,1,100),5,2),
          MemoryMatrix=TRUE)     
          
     cat("<< The value of parameter obtained from genetic algorithm >>\n\n")
     
     namegen<-c("ka","kel","k12","k21","Vd")
     outgen<-c(gen$par[1],gen$par[2],gen$par[3],gen$par[4],gen$par[5])
     print(data.frame(Parameter=namegen,Value=outgen)) 
     F<-objfun(gen$par)
     
     opt<-optim(c(gen$par[1],gen$par[2],gen$par[3],gen$par[4],gen$par[5]),objfun, method="Nelder-Mead")
     nameopt<-c("ka","kel","k12","k21","Vd")
     outopt<-c(opt$par[1],opt$par[2],opt$par[3],opt$par[4],opt$par[5])
     
     cat("\n<< The value of parameter fitted by Nelder-Mead Simplex slgorithm >>\n\n")
     print(data.frame(Parameter=nameopt,Value=outopt))
     cat("\n<< Residual sum-of-squares and parameter values fitted by nls >>\n\n")
     
     fm<-nls(conc ~ modfun3(time,ka,kel,k12,k21,Vd), data=PKindex,subset=Subject==i,
         start=list(ka=opt$par[1],kel=opt$par[2],k12=opt$par[3],k21=opt$par[4],Vd=opt$par[5]),trace=TRUE,
         nls.control(maxiter=50,tol=1e-2,minFactor=1/1024))
     cat("\n")
     coef<-data.frame(coef(fm)["kel"])     
     plotting.lin(PKindex, fm, i, pick, coef, xaxis, yaxis)
     }
  })
  cat("\n")   
}

###Simulation
###Two compartment PK model iv bolus single dose
sbolus2 <- function(Subject=NULL,  # N Subj's 
                    PKtime=NULL,   # times for sampling
                    Dose=NULL,     # single dose
                    Vd=NULL,
                    kel=NULL, 
                    k12=NULL,
                    k21=NULL)
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
   
   if (is.null(kel) || is.null(k12) || is.null(k21) || is.null(Vd)) {
       par<-data.frame(Parameter=c("kel","k12","k21","Vd"),Initial=c(0))
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
       par1 <- kel
       par2 <- k12
       par3 <- k21
       par4 <- Vd
    }
    
    defun<- function(time, y, parms) { 
      dCp1dt <- -parms["kel"]*y[1]-parms["k12"]*y[1]+parms["k21"]*y[2] 
      dCp2dt <-  parms["k12"]*y[1]-parms["k21"]*y[2]
      list(c(dCp1dt,dCp2dt)) 
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
           kel<-par1
           k12<-par2
           k21<-par3
           Vd<-par4   
           PKindex[[i]]<-sbolus2.out(PKtime,kel,k12,k21,Vd,defun,par1,par2,par3,par4,Dose,i,type)
         }   
         PKindex<- as.data.frame(do.call("rbind",PKindex))
         rownames(PKindex) <- seq(nrow(PKindex)) 
         savefile(PKindex)    
       }   
       else {
         cat("\n\nEnter error factor for kel\n")
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
         
         cat("\nEnter error factor for k12\n")
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
         
         cat("\nEnter error factor for k21\n")
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
                    {kel<-par1+rnorm(1,mean=0,sd=factor1)
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
                        Vd<-par4+rnorm(1,mean=0,sd=factor4)}},
                        
                    {kel<-par1+runif(1,min=-factor1,max=factor1)
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
                        Vd<-par4+runif(1,min=-factor4,max=factor4)}},
                        
                     {kel<-par1*rnorm(1,mean=0,sd=factor1)+par1
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
                         Vd<-par4*rnorm(1,mean=0,sd=factor4)+par4}},
                         
                      {kel<-par1*runif(1,min=-factor1,max=factor1)+par1
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
                          Vd<-par4*runif(1,min=-factor4,max=factor4)+par4}}
               )
               PKindex[[i]]<-sbolus2.out(PKtime,kel,k12,k21,Vd,defun,par1,par2,par3,par4,Dose,i,type)         
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
     
     cat("\nEnter error factor for kel\n")
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
     
     cat("\nEnter error factor for k12\n")
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
     
     cat("\nEnter error factor for k21\n")
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
     cat("****************************************************\n")
     cat("Summary Table                                       \n")
     cat("Model: 2-Compartment, IV-Bolus, & Single-Dose Model \n") 
     cat("Subject #:", Subject,"                              \n")
     cat("Error Type:", type,"                                \n")
     cat("Simulation #:", re,"                                \n\n")
     sim<-matrix(c(par1,par2,par3,par4,factor1,factor2,factor3,factor4),4,2)
     dimnames(sim)<-list(c("kel","k12","k21","Vd"),c("Original","Error factor"))
     show(sim)   
     cat("****************************************************\n\n")
     
     PKindex<-vector(Subject,mode="list")
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" ) 
     C1.lsoda<-list()
        for (j in 1:re){
            switch(pick,
                   {kel<-par1+rnorm(1,mean=0,sd=factor1)
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
                        Vd<-par4+rnorm(1,mean=0,sd=factor4)}},
                        
                    {kel<-par1+runif(1,min=-factor1,max=factor1)
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
                        Vd<-par4+runif(1,min=-factor4,max=factor4)}},
                        
                     {kel<-par1*rnorm(1,mean=0,sd=factor1)+par1
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
                         Vd<-par4*rnorm(1,mean=0,sd=factor4)+par4}},
                         
                      {kel<-par1*runif(1,min=-factor1,max=factor1)+par1
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
                          Vd<-par4*runif(1,min=-factor4,max=factor4)+par4}}
               )
               time1<-PKtime$time
               parms<-c(kel=kel,k12=k12,k21=k21,Vd=Vd)  
               XX<-data.frame(lsoda(c(Dose/Vd,0),c(0,time1),defun,parms)) 
               C1.lsoda[[j]]<-data.frame(XX[2:(length(time1)+1),1],XX[2:(length(time1)+1),2])
               colnames(C1.lsoda[[j]])<-list("time","concentration") 
          }   
       PKindex[[i]]<-montecarlo(C1.lsoda,time1,i,re) 
      }  
    PKindex<-as.data.frame(do.call("rbind",PKindex))
    rownames(PKindex)<-seq(nrow(PKindex)) 
    savefile(PKindex)
   }  
}      
                   

###Two compartment PK model iv infusion single dose 
sinfu2<- function(Subject=NULL,  # N Subj's 
                  PKtime=NULL,   # times for sampling
                  Dose=NULL,     # single dose
                  Tinf=NULL,
                  Vd=NULL,
                  kel=NULL, 
                  k12=NULL,
                  k21=NULL)
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
   
   if (is.null(Tinf)) {
    cat("\nEnter value for infusion time\n")
    Tinf<-scan(nlines=1,quiet=TRUE)
    cat("\n")
  }
   
   if (is.null(kel) || is.null(k12) || is.null(k21) || is.null(Vd)) {
       par<-data.frame(Parameter=c("kel","k12","k21","Vd"),Initial=c(0))
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
       par1 <- kel
       par2 <- k12
       par3 <- k21
       par4 <- Vd
    }
    
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
           kel<-par1
           k12<-par2
           k21<-par3
           Vd<-par4   
           PKindex[[i]]<-sinfu2.out(PKtime,kel,k12,k21,Vd,defun,par1,par2,par3,par4,Dose,i,type)  
         }   
         PKindex<- as.data.frame(do.call("rbind",PKindex))
         rownames(PKindex) <- seq(nrow(PKindex)) 
         savefile(PKindex)    
       }   
       else {
         cat("\n\nEnter error factor for kel\n")
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
         
         cat("\nEnter error factor for k12\n")
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
         
         cat("\nEnter error factor for k21\n")
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
                    {kel<-par1+rnorm(1,mean=0,sd=factor1)
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
                        Vd<-par4+rnorm(1,mean=0,sd=factor4)}},
                        
                    {kel<-par1+runif(1,min=-factor1,max=factor1)
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
                        Vd<-par4+runif(1,min=-factor4,max=factor4)}},
                        
                     {kel<-par1*rnorm(1,mean=0,sd=factor1)+par1
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
                         Vd<-par4*rnorm(1,mean=0,sd=factor4)+par4}},
                         
                      {kel<-par1*runif(1,min=-factor1,max=factor1)+par1
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
                          Vd<-par4*runif(1,min=-factor4,max=factor4)+par4}}
               )
               PKindex[[i]]<-sinfu2.out(PKtime,kel,k12,k21,Vd,defun,par1,par2,par3,par4,Dose,i,type)           
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
     
     cat("\nEnter error factor for kel\n")
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
     
     cat("\nEnter error factor for k12\n")
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
     
     cat("\nEnter error factor for k21\n")
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
     cat("*******************************************************\n")
     cat("Summary Table                                          \n")
     cat("Model: 2-Compartment, IV-Infusion, & Single-Dose Model \n") 
     cat("Subject #:", Subject,"                                 \n")
     cat("Error Type:", type,"                                   \n")
     cat("Simulation #:", re,"                                   \n\n")
     sim<-matrix(c(par1,par2,par3,par4,factor1,factor2,factor3,factor4),4,2)
     dimnames(sim)<-list(c("kel","k12","k21","Vd"),c("Original","Error factor"))
     show(sim)   
     cat("*******************************************************\n\n")
     
     PKindex<-vector(Subject,mode="list")
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" ) 
     C1.lsoda<-list()
        for (j in 1:re){
            switch(pick,
                   {kel<-par1+rnorm(1,mean=0,sd=factor1)
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
                        Vd<-par4+rnorm(1,mean=0,sd=factor4)}},
                        
                    {kel<-par1+runif(1,min=-factor1,max=factor1)
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
                        Vd<-par4+runif(1,min=-factor4,max=factor4)}},
                        
                     {kel<-par1*rnorm(1,mean=0,sd=factor1)+par1
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
                         Vd<-par4*rnorm(1,mean=0,sd=factor4)+par4}},
                         
                      {kel<-par1*runif(1,min=-factor1,max=factor1)+par1
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
                          Vd<-par4*runif(1,min=-factor4,max=factor4)+par4}}
               )
               time1<-PKtime$time
               parms<-c(kel=kel,k12=k12,k21=k21,Vd=Vd)  
               XX<-data.frame(lsoda(c(0,0),c(0,time1),defun,parms)) 
               C1.lsoda[[j]]<-data.frame(XX[2:(length(time1)+1),1],XX[2:(length(time1)+1),2])
               colnames(C1.lsoda[[j]])<-list("time","concentration") 
          }   
      PKindex[[i]]<-montecarlo(C1.lsoda,time1,i,re) 
     }  
    PKindex<-as.data.frame(do.call("rbind",PKindex))
    rownames(PKindex)<-seq(nrow(PKindex)) 
    savefile(PKindex)
  }  
}   


###Two compartment PK model extravascular single dose first order absorption
sfirst2<- function(Subject=NULL,  # N Subj's 
                  PKtime=NULL,   # times for sampling
                  Dose=NULL,     # single dose
                  ka=NULL,
                  Vd=NULL,
                  kel=NULL, 
                  k12=NULL,
                  k21=NULL)
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
   
   
   if (is.null(ka) || is.null(kel) || is.null(k12) || is.null(k21) || is.null(Vd)) {
       par<-data.frame(Parameter=c("ka","kel","k12","k21","Vd"),Initial=c(0))
       par<-edit(par)
       repeat{
           if ( par[1,2] == 0 || par[2,2] ==0 || par[3,2]==0 || par[4,2]==0 || par[5,2]==0){
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
       par5<-par[5,2]
    } 
    else {
       par1 <- ka
       par2 <- kel
       par3 <- k12
       par4 <- k21
       par5 <- Vd
    }
    
    defun<- function(time, y, parms) { 
       dCp1dt <- -parms["ka"]*y[1]
       dCp2dt <- parms["ka"]*y[1]/parms["Vd"]-parms["kel"]*y[2]+parms["k21"]*y[3]-parms["k12"]*y[2]
       dCp3dt <- parms["k12"]*y[2]-parms["k21"]*y[3]
       list(c(dCp1dt,dCp2dt,dCp3dt)) 
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
           ka<-par1
           kel<-par2
           k12<-par3
           k21<-par4
           Vd<-par5   
           PKindex[[i]]<-sfirst2.out(PKtime,ka,kel,k12,k21,Vd,defun,par1,par2,par3,par4,par5,Dose,i,type)  
         }   
         PKindex<- as.data.frame(do.call("rbind",PKindex))
         rownames(PKindex) <- seq(nrow(PKindex)) 
         savefile(PKindex)    
       }   
       else {
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
         
         cat("\nEnter error factor for k12\n")
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
         
         cat("\nEnter error factor for k21\n")
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
         
         cat("\nEnter error factor for Vd\n")
         factor5<-scan(nlines=1,quiet=TRUE)
         repeat{
              if ( factor5 == 0 ){
                cat("\n")
                cat("**********************************\n")
                cat(" Parameter value can not be zero. \n")
                cat(" Press enter to continue.         \n")
                cat("**********************************\n\n")
                readline()
                cat("\n")
                factor5<-scan(nlines=1,quiet=TRUE)}   
           else{
                break
                return(factor5)}
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
                     k12<-par3+rnorm(1,mean=0,sd=factor3)
                     while(k12<=0){
                        k12<-par3+rnorm(1,mean=0,sd=factor3)}
                     k21<-par4+rnorm(1,mean=0,sd=factor4)
                     while(k21<=0){
                        k21<-par4+rnorm(1,mean=0,sd=factor4)}
                     Vd<-par5+rnorm(1,mean=0,sd=factor5)
                     while(Vd<=0){
                        Vd<-par5+rnorm(1,mean=0,sd=factor5)}},
                        
                    {ka<-par1+runif(1,min=-factor1,max=factor1)
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
                        Vd<-par5+runif(1,min=-factor5,max=factor5)}},
                        
                     {ka<-par1*rnorm(1,mean=0,sd=factor1)+par1
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
                         Vd<-par5*rnorm(1,mean=0,sd=factor5)+par5}},
                         
                      {ka<-par1*runif(1,min=-factor1,max=factor1)+par1
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
                          Vd<-par5*runif(1,min=-factor5,max=factor5)+par5}}
               )
               PKindex[[i]]<-sfirst2.out(PKtime,ka,kel,k12,k21,Vd,defun,par1,par2,par3,par4,par5,Dose,i,type)               
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
     
     cat("\nEnter error factor for k12\n")
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
     
     cat("\nEnter error factor for k21\n")
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
     
     cat("\nEnter error factor for Vd\n")
     factor5<-scan(nlines=1,quiet=TRUE)
     repeat{
              if ( factor5 == 0 ){
                cat("\n")
                cat("**********************************\n")
                cat(" Parameter value can not be zero. \n")
                cat(" Press enter to continue.         \n")
                cat("**********************************\n\n")
                readline()
                cat("\n")
                factor5<-scan(nlines=1,quiet=TRUE)}   
           else{
                break
                return(factor5)}
           }
     
     cat("\n")
     cat("******************************************************\n")
     cat("Summary Table                                         \n")
     cat("Model: 2-Compartment, Extravascular,                  \n") 
     cat("       Single-Dose & 1-Ordered without Lag Time Model \n") 
     cat("Subject #:", Subject,"                                \n")
     cat("Error Type:", type,"                                  \n")
     cat("Simulation #:", re,"                                  \n\n")
     sim<-matrix(c(par1,par2,par3,par4,par5,factor1,factor2,factor3,factor4,factor5),5,2)
     dimnames(sim)<-list(c("ka","kel","k12","k21","Vd"),c("Original","Error factor"))
     show(sim)   
     cat("******************************************************\n\n")
     
     PKindex<-vector(Subject,mode="list")
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" ) 
     C1.lsoda<-list()
        for (j in 1:re){
            switch(pick,
                   {ka<-par1+rnorm(1,mean=0,sd=factor1)
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
                        Vd<-par5+rnorm(1,mean=0,sd=factor5)}},
                        
                    {ka<-par1+runif(1,min=-factor1,max=factor1)
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
                        Vd<-par5+runif(1,min=-factor5,max=factor5)}},
                        
                     {ka<-par1*rnorm(1,mean=0,sd=factor1)+par1
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
                         Vd<-par5*rnorm(1,mean=0,sd=factor5)+par5}},
                         
                      {ka<-par1*runif(1,min=-factor1,max=factor1)+par1
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
                          Vd<-par5*runif(1,min=-factor5,max=factor5)+par5}}
               )
               time1<-PKtime$time
               parms<-c(ka=ka,kel=kel,k12=k12,k21=k21,Vd=Vd)
               XX<-data.frame(lsoda(c(Dose,0,0),c(0,time1),defun,parms)) 
               C1.lsoda[[j]]<-data.frame(XX[2:(length(time1)+1),1],XX[2:(length(time1)+1),3])
               colnames(C1.lsoda[[j]])<-list("time","concentration") 
          }      
    PKindex[[i]]<-montecarlo(C1.lsoda,time1,i,re) 
    }  
    PKindex<-as.data.frame(do.call("rbind",PKindex))
    rownames(PKindex)<-seq(nrow(PKindex)) 
    savefile(PKindex)
  }  
}   
