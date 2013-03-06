### PKindex is the target Dataset.


### Normal fitting
### One compartment PK model extravascular single 
### dose zero-order absorption without lag time
### optional Michaelis-Menten Elimination
fzero.nolag <- function(PKindex,
                        Dose=NULL, 
                        Tabs=NULL,
                        Vm=NULL,Km=NULL, ## MMe=TRUE
                        Vd=NULL,
                        kel=NULL,        ## MMe=FALSE
                        MMe=FALSE) 
{ 
   #options(warn=-1)
        
   ## Input dose and initial value for Tabs, kel and Vd
   
   if (is.null(Dose)) {
     cat("Enter Dose\n")
     Dose <- scan(nlines=1,quiet=TRUE)
   } 
   else {
     cat("Dose from arguments is = ",Dose,"\n")
   }
   
   if (MMe){
      if ( is.null(Tabs) || is.null(Vm) || is.null(Km) || is.null(Vd) ) {
        par<-data.frame(Parameter=c("Tabs","Vm","Km","Vd"),Initial=c(0))
        par<-edit(par)
        repeat{
           if ( par[1,2] == 0 || par[2,2] ==0 || par[3,2]==0 || par[4,2]==0){
             cat("\n")
             cat("**********************************\n")
             cat(" Parameter value can not be zero. \n")
             cat(" Press Enter to continue.         \n")
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
   else {
      ## No MM elimination
      if ( is.null(Tabs) || is.null(kel) || is.null(Vd)) {
        par<-data.frame(Parameter=c("Tabs","kel","Vd"),Initial=c(0))
        par<-edit(par)
        repeat{
           if ( par[1,2] == 0 || par[2,2] ==0 || par[3,2]==0){
             cat("\n")
             cat("**********************************\n")
             cat(" Parameter value can not be zero. \n")
             cat(" Press Enter to continue.         \n")
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
   
   if (!MMe) {
      ## User-supplied function w/o Michaelis-Mention elimination
      defun<- function(time, y, parms) { 
      if(time<=parms["Tabs"]) 
        dCpdt <- (Dose/parms["Tabs"])/parms["Vd"] - parms["kel"] * y[1]
      else
        dCpdt <- - parms["kel"] * y[1]
      list(dCpdt) 
      } 
   
      modfun1 <- function(time,Tabs,kel,Vd) { 
         out <- lsoda(0,c(0,time),defun,parms=c(Tabs=Tabs,kel=kel,Vd=Vd),
                      rtol=1e-8,atol=1e-8) 
         out[-1,2]
      } 
   } 
   else {
     ## User-supplied function with MM elimination
      defun<- function(time, y, parms) { 
      if(time<=parms["Tabs"]) 
        dCpdt <- (Dose/parms["Tabs"])/parms["Vd"]-(parms["Vm"]/parms["Vd"])*y[1]/(parms["Km"]/parms["Vd"]+y[1])            
      else
        dCpdt <- -(parms["Vm"]/parms["Vd"])*y[1]/(parms["Km"]/parms["Vd"]+y[1])            
      list(dCpdt) 
      }   
   
      modfun2 <- function(time,Tabs,Vm,Km,Vd) { 
      out <- lsoda(0,c(0,time),defun,parms=c(Tabs=Tabs,Vm=Vm,Km=Km,Vd=Vd),
                   rtol=1e-6,atol=1e-8) 
      out[-1,2]
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
           out<-modfun2(PKindex$time[PKindex$Subject==i],par[1],par[2],par[3],par[4])
        } 
        else {
           ## No MM elimination
           out<-modfun1(PKindex$time[PKindex$Subject==i],par[1],par[2],par[3])
        }
      gift <- which( PKindex$conc[PKindex$Subject==i] != 0 )
      switch(pick,
             sum((PKindex$conc[PKindex$Subject==i][gift]-out[gift])^2),
             sum((PKindex$conc[PKindex$Subject==i][gift]-out[gift])^2/PKindex$conc[gift]),
             sum(((PKindex$conc[PKindex$Subject==i][gift] - out[gift])/PKindex$conc[gift])^2))
      }
      
      if (MMe) {
         gen<-genoud(objfun,nvars=4,max=FALSE,pop.size=30,max.generations=20,
              wait.generations=10,
              starting.value=c(par[1,2],par[2,2],par[3,2],par[4,2]),
              BFGS=FALSE,print.level=0,boundary.enforcement=2,
              Domains=matrix(c(1,1,1,1,50,100,100,100),4,2),
              MemoryMatrix=TRUE)
      } 
      else {
         gen<-genoud(objfun,nvars=3,max=FALSE,pop.size=30,max.generations=20,
              wait.generations=10,
              starting.value=c(par[1,2],par[2,2],par[3,2]),
              BFGS=FALSE,print.level=0,boundary.enforcement=2,
              Domains=matrix(c(1,0.01,1,50,1,100),3,2),
              MemoryMatrix=TRUE)
      }
      
     cat("<< PK parameters obtained from genetic algorithm >>\n\n")
     if (MMe) {
        namegen<-c("Tabs","Vm","Km","Vd")
        outgen<-c(gen$par[1],gen$par[2],gen$par[3],gen$par[4])
     } 
     else {
        ## No MM elimination
        namegen<-c("Tabs","kel","Vd")
        outgen<-c(gen$par[1],gen$par[2],gen$par[3])
     }
     print(data.frame(Parameter=namegen,Value=outgen))  
     F<-objfun(gen$par)  
     
     if (MMe) {
        opt<-optim(c(gen$par[1],gen$par[2],gen$par[3],gen$par[4]),objfun, method="Nelder-Mead")  
        nameopt<-c("Tabs","Vm","Km","Vd")
        outopt<-c(opt$par[1],opt$par[2],opt$par[3],opt$par[4])
     }
     else {
        opt<-optim(c(gen$par[1],gen$par[2],gen$par[3]),objfun, method="Nelder-Mead")  
        nameopt<-c("Tabs","kel","Vd")
        outopt<-c(opt$par[1],opt$par[2],opt$par[3])
     }
     
     cat("\n<< PK parameters obtained from Nelder-Mead Simplex algorithm >>\n\n")
     print(data.frame(Parameter=nameopt,Value=outopt))
     cat("\n<< Residual sum-of-square (RSS) and final PK parameters with nls >>\n\n")

     if (MMe) {
        fm<-nls(conc~modfun2(time,Tabs,Vm,Km,Vd),data=PKindex, algorithm="default",subset=Subject==i,
            start=list(Tabs=opt$par[1],Vm=opt$par[2],Km=opt$par[3],Vd=opt$par[4]),trace=TRUE,
            nls.control(maxiter=50,tol=10,minFactor=1/2048))
        cat("\n")
        plotting.non(PKindex, fm, i, pick,xaxis,yaxis)
     } 
     else {
        ## No MM elimination
        fm<-nls(conc~modfun1(time,Tabs,kel,Vd),data=PKindex, algorithm="default",subset=Subject==i,
            start=list(Tabs=opt$par[1],kel=opt$par[2],Vd=opt$par[3]),trace=TRUE,
            nls.control(tol=1.0,minFactor=1/2048))
        cat("\n")
        coef<-data.frame(coef(fm)["kel"])
        plotting.lin(PKindex, fm, i, pick, coef, xaxis, yaxis)
     }
  }
  })
  cat("\n")   
}


## Legacy function
fzero.nolagm <- function(PKindex,...) {
  fzero.nolag(PKindex,...,MMe=TRUE)
}


###Simulation
###One compartment PK model iv infusion single dose
###optional Michaelis-Menten Elimination
szero.nolag<- function(Subject=NULL,  # N Subj's 
                   PKtime=NULL,   ## times for sampling
                   Dose=NULL,     ## single dose
                   Tabs=NULL,    
                   Vd=NULL,
                   kel=NULL,      ## If not MM elimination
                   MMe=FALSE,     ## michaelis-menten elimination?
                   Vm=NULL,Km=NULL)
{
   #options(warn=-1)
   
   if (is.null(Subject) || !is.integer(Subject)) {
     cat("How many subjects do you want?\n")
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
    cat("\nEnter dose\n")
    Dose<-scan(nlines=1,quiet=TRUE)
   }
   
   if ( !MMe) { 
    if ( is.null(Tabs) || is.null(kel) || is.null(Vd)) {
       par<-data.frame(Parameter=c("Tabs","kel","Vd"),Initial=c(0))
       par<-edit(par)
       repeat{
           if ( par[1,2] == 0 || par[2,2] ==0 || par[3,2]==0){
             cat("\n")
             cat("**********************************\n")
             cat(" Parameter value can not be zero. \n")
             cat(" Press Enter to continue.         \n")
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
       par1 <- Tabs
       par2 <- kel
       par3 <- Vd
     } 
  }
  else{
    if ( is.null(Tabs) || is.null(Vm) || is.null(Km) || is.null(Vd)){
       par<-data.frame(Parameter=c("Tabs","Vm","Km","Vd"),Initial=c(0))
       par<-edit(par)
       repeat{
           if ( par[1,2] == 0 || par[2,2] ==0 || par[3,2]==0 || par[4,2]==0){
             cat("\n")
             cat("**********************************\n")
             cat(" Parameter value can not be zero. \n")
             cat(" Press Enter to continue.         \n")
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
       par1 <- Tabs
       par2 <- Vm
       par3 <- Km
       par4 <- Vd
    }
  } 
   
  if ( ! MMe){
    defun<- function(time, y, parms) { 
     if(time<=parms["Tabs"]) 
        dCpdt <- (Dose/parms["Tabs"])/parms["Vd"] - parms["kel"] * y[1]
     else
        dCpdt <- - parms["kel"] * y[1]
     list(dCpdt) 
    } 
  }
  else{
    defun<- function(time, y, parms) { 
     if(time<=parms["Tabs"]) 
       dCpdt <- (Dose/parms["Tabs"])/parms["Vd"]-(parms["Vm"]/parms["Vd"])*y[1]/(parms["Km"]/parms["Vd"]+y[1])            
     else
       dCpdt <- -(parms["Vm"]/parms["Vd"])*y[1]/(parms["Km"]/parms["Vd"]+y[1])            
     list(dCpdt) 
    }
  }    
  
  file.menu <- c("Simulation with Error",
                 "Monte Carlo Simulation")
  pick <- menu(file.menu, title = "<< Simulation Types >>")
  if (pick ==1){
     cat("\n\n")
     file.menu <- c("No Error",
                    "Error = Normal Error", 
                    "Error = Uniform Error",
                    "Error = Normal Error*True Value",
                    "Error = Uniform Error*True Value")
     pick <- menu(file.menu, title = "<< Error Types >>")      
     
     type<-switch(pick, 
                  "No Error",
                  "Error = Normal Error", 
                  "Error = Uniform Error",
                  "Error = Normal Error*True Value",
                  "Error = Uniform Error*True Value")
     
     
       if (pick ==1){
         PKindex<-vector(Subject,mode="list")
         for( i in 1:Subject)  {
           cat("\n\n             << Subject",i,">>\n\n" )   
           
           if (! MMe ){
             Tabs<-par1
             kel<-par2
             Vd<-par3  
             PKindex[[i]]<-szero.out(PKtime,Tabs,kel,Vd,defun,par1,par2,par3,Dose,i,type)     
           } 
           else{
             Tabs<-par1
             Vm<-par2
             Km<-par3
             Vd<-par4     
             PKindex[[i]]<-szero.mm.out(PKtime,Tabs,Vm,Km,Vd,defun,par1,par2,par3,par4,Dose,i,type) 
           }
         }       
         PKindex<- as.data.frame(do.call("rbind",PKindex))
         rownames(PKindex) <- seq(nrow(PKindex)) 
         savefile(PKindex)
       }   
       
       else {
         if (! MMe){
           cat("\n\nEnter error factor for Tabs\n")
           factor1<-scan(nlines=1,quiet=TRUE)
           repeat{
              if ( factor1 == 0 ){
                cat("\n")
                cat("**********************************\n")
                cat(" Parameter value can not be zero. \n")
                cat(" Press Enter to continue.         \n")
                cat("**********************************\n\n")
                readline()
                cat("\n")
                factor1<-scan(nlines=1,quiet=TRUE)}   
           else{
                break
                return(factor1)}
           } 
           
           cat("\n\nEnter error factor for kel\n")
           factor2<-scan(nlines=1,quiet=TRUE)
           repeat{
              if ( factor2 == 0 ){
                cat("\n")
                cat("**********************************\n")
                cat(" Parameter value can not be zero. \n")
                cat(" Press Enter to continue.         \n")
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
                cat(" Press Enter to continue.         \n")
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
                   {Tabs<-par1+rnorm(1,mean=0,sd=factor1)
                    while(Tabs<=0){
                       Tabs<-par1+rnorm(1,mean=0,sd=factor1)}
                    kel<-par2+rnorm(1,mean=0,sd=factor2)
                    while(kel<=0){
                       kel<-par2+rnorm(1,mean=0,sd=factor2)}
                    Vd<-par3+rnorm(1,mean=0,sd=factor3)
                    while(Vd<=0){
                       Vd<-par3+rnorm(1,mean=0,sd=factor3)}},
                  
                   {Tabs<-par1+runif(1,min=-factor1,max=factor1)
                    while(Tabs<=0){
                       Tabs<-par1+runif(1,min=-factor1,max=factor1)}
                    kel<-par2+runif(1,min=-factor2,max=factor2)
                    while(kel<=0){
                       kel<-par2+runif(1,min=-factor2,max=factor2)}
                    Vd<-par3+runif(1,min=-factor3,max=factor3)
                    while(Vd<=0){
                       Vd<-par3+runif(1,min=-factor3,max=factor3)}}, 
                    
                   {Tabs<-par1*rnorm(1,mean=0,sd=factor1)+par1
                    while(Tabs<=0){
                       Tabs<-par1*rnorm(1,mean=0,sd=factor1)+par1}
                    kel<-par2*rnorm(1,mean=0,sd=factor2)+par2
                    while(kel<=0){
                       kel<-par2*rnorm(1,mean=0,sd=factor2)+par2}
                    Vd<-par3*rnorm(1,mean=0,sd=factor3)+par3
                    while(Vd<=0){
                       Vd<-par3*rnorm(1,mean=0,sd=factor3)+par3}},
                  
                   {Tabs<-par1*runif(1,min=-factor1,max=factor1)+par1
                    while(Tabs<=0){
                       Tabs<-par1*runif(1,min=-factor1,max=factor1)+par1}
                    kel<-par2*runif(1,min=-factor2,max=factor2)+par2
                    while(kel<=0){
                       kel<-par2*runif(1,min=-factor2,max=factor2)+par2}
                    Vd<-par3*runif(1,min=-factor3,max=factor3)+par3
                    while(Vd<=0){
                      Vd<-par3*runif(1,min=-factor3,max=factor3)+par3}}   
             )      
             PKindex[[i]]<-szero.out(PKtime,Tabs,kel,Vd,defun,par1,par2,par3,Dose,i,type)    
           }
         }
         else{
           cat("\n\nEnter error factor for Tabs\n")
           factor1<-scan(nlines=1,quiet=TRUE)
           repeat{
              if ( factor1 == 0 ){
                cat("\n")
                cat("**********************************\n")
                cat(" Parameter value can not be zero. \n")
                cat(" Press Enter to continue.         \n")
                cat("**********************************\n\n")
                readline()
                cat("\n")
                factor1<-scan(nlines=1,quiet=TRUE)}   
           else{
                break
                return(factor1)}
           } 
           
           cat("\n\nEnter error factor for Vm\n")
           factor2<-scan(nlines=1,quiet=TRUE)
           repeat{
              if ( factor2 == 0 ){
                cat("\n")
                cat("**********************************\n")
                cat(" Parameter value can not be zero. \n")
                cat(" Press Enter to continue.         \n")
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
                cat(" Press Enter to continue.         \n")
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
                cat(" Press Enter to continue.         \n")
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
                   {Tabs<-par1+rnorm(1,mean=0,sd=factor1)
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
                       Vd<-par4+rnorm(1,mean=0,sd=factor4)}},
                  
                   {Tabs<-par1+runif(1,min=-factor1,max=factor1)
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
                       Vd<-par4+runif(1,min=-factor4,max=factor4)}}, 
                    
                   {Tabs<-par1*rnorm(1,mean=0,sd=factor1)+par1
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
                       Vd<-par4*rnorm(1,mean=0,sd=factor4)+par4}},
                  
                   {Tabs<-par1*runif(1,min=-factor1,max=factor1)+par1
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
                       Vd<-par4*runif(1,min=-factor4,max=factor4)+par4}}   
             )      
             PKindex[[i]]<-szero.mm.out(PKtime,Tabs,Vm,Km,Vd,defun,par1,par2,par3,par4,Dose,i,type) 
           } 
         }    
         PKindex<- as.data.frame(do.call("rbind",PKindex))
         rownames(PKindex) <- seq(nrow(PKindex)) 
         savefile(PKindex)
   } 
  } 
  else if (pick ==2){ 
     cat("\n\n")
     file.menu <- c("Error = Normal Error", 
                    "Error = Uniform Error",
                    "Error = Normal Error*True Value",
                    "Error = Uniform Error*True Value")
                    
     type<-switch(pick, 
                  "Error = Normal Error", 
                  "Error = Uniform Error",
                  "Error = Normal Error*True Value",
                  "Error = Uniform Error*True Value")
                    
     pick <- menu(file.menu, title = "<< Error Types >>")
     cat("\n\nHow many interations do you want to run?\n")
     re<-scan(nlines=1,quiet=TRUE)
     
       if (! MMe){
         cat("\nEnter error factor for Tabs\n")
         factor1<-scan(nlines=1,quiet=TRUE)
         repeat{
              if ( factor1 == 0 ){
                cat("\n")
                cat("**********************************\n")
                cat(" Parameter value can not be zero. \n")
                cat(" Press Enter to continue.         \n")
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
                cat(" Press Enter to continue.         \n")
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
                cat(" Press Enter to continue.         \n")
                cat("**********************************\n\n")
                readline()
                cat("\n")
                factor3<-scan(nlines=1,quiet=TRUE)}   
           else{
                break
                return(factor3)}
           } 
         
         cat("\n")
         cat("********************************************************\n")
         cat("Summary Table                                           \n")
         cat("Model: 1-Compartment, Extravascular, Single-Dose,       \n")
         cat("       & Zero-Ordered Absorption without Lag Time Model \n")
         cat("Subject #:", Subject,"                                  \n")
         cat("Error Type:", type,"                                    \n")
         cat("Simulation #:", re,"                                    \n\n")
         sim<-matrix(c(par1,par2,par3,factor1,factor2,factor3),3,2)
         dimnames(sim)<-list(c("Tabs","kel","Vd"),c("Selected","Error factor"))
         show(sim)   
         cat("********************************************************\n\n")
       }
       else{
         cat("\nEnter error factor for Tabs\n")
         factor1<-scan(nlines=1,quiet=TRUE)
         cat("\nEnter error factor for Vm\n")
         factor2<-scan(nlines=1,quiet=TRUE)
         cat("\nEnter error factor for Km\n")
         factor3<-scan(nlines=1,quiet=TRUE)
         cat("\nEnter error factor for Vd\n")
         factor4<-scan(nlines=1,quiet=TRUE)
         cat("\n")
         cat("***************************************************************************\n")
         cat("Summary Table                                                              \n")
         cat("Model: 1-Compartment, Extravascular, Single-Dose, Zero-Ordered Absorption, \n") 
         cat("       & Michaelis-Menten Elimination without Lag Time Model               \n")
         cat("Subject #:", Subject,"                                                     \n")
         cat("Error Type:", type,"                                                       \n")
         cat("Simulation #:", re,"                                                       \n\n")
         sim<-matrix(c(par1,par2,par3,par4,factor1,factor2,factor3,factor4),4,2)
         dimnames(sim)<-list(c("Tabs","Vm","Km","Vd"),c("Selected","Error factor"))
         show(sim)   
         cat("**************************************************************************\n\n")
       }
       
       PKindex<-vector(Subject,mode="list")
       for( i in 1:Subject)  {
       cat("\n\n             << Subject",i,">>\n\n" ) 
       C1.lsoda<-list()

         for (j in 1:re){
     
           if (! MMe){
              switch(pick, 
                    {Tabs<-par1+rnorm(1,mean=0,sd=factor1)
                    while(Tabs<=0){
                       Tabs<-par1+rnorm(1,mean=0,sd=factor1)}
                    kel<-par2+rnorm(1,mean=0,sd=factor2)
                    while(kel<=0){
                       kel<-par2+rnorm(1,mean=0,sd=factor2)}
                    Vd<-par3+rnorm(1,mean=0,sd=factor3)
                    while(Vd<=0){
                       Vd<-par3+rnorm(1,mean=0,sd=factor3)}},
                  
                   {Tabs<-par1+runif(1,min=-factor1,max=factor1)
                    while(Tabs<=0){
                       Tabs<-par1+runif(1,min=-factor1,max=factor1)}
                    kel<-par2+runif(1,min=-factor2,max=factor2)
                    while(kel<=0){
                       kel<-par2+runif(1,min=-factor2,max=factor2)}
                    Vd<-par3+runif(1,min=-factor3,max=factor3)
                    while(Vd<=0){
                       Vd<-par3+runif(1,min=-factor3,max=factor3)}}, 
                    
                   {Tabs<-par1*rnorm(1,mean=0,sd=factor1)+par1
                    while(Tabs<=0){
                       Tabs<-par1*rnorm(1,mean=0,sd=factor1)+par1}
                    kel<-par2*rnorm(1,mean=0,sd=factor2)+par2
                    while(kel<=0){
                       kel<-par2*rnorm(1,mean=0,sd=factor2)+par2}
                    Vd<-par3*rnorm(1,mean=0,sd=factor3)+par3
                    while(Vd<=0){
                       Vd<-par3*rnorm(1,mean=0,sd=factor3)+par3}},
                  
                   {Tabs<-par1*runif(1,min=-factor1,max=factor1)+par1
                    while(Tabs<=0){
                       Tabs<-par1*runif(1,min=-factor1,max=factor1)+par1}
                    kel<-par2*runif(1,min=-factor2,max=factor2)+par2
                    while(kel<=0){
                       kel<-par2*runif(1,min=-factor2,max=factor2)+par2}
                    Vd<-par3*runif(1,min=-factor3,max=factor3)+par3
                    while(Vd<=0){
                      Vd<-par3*runif(1,min=-factor3,max=factor3)+par3}}     
              )    
              time1<-PKtime$time
              parms<-c(Tabs=Tabs,kel=kel,Vd=Vd) 
              XX<-data.frame(lsoda(0, c(0,time1), defun, parms))   
              C1.lsoda[[j]]<-data.frame(XX[2:(length(time1)+1),1],XX[2:(length(time1)+1),2])
              colnames(C1.lsoda[[j]])<-list("time","concentration")
           }
           else{
              switch(pick, 
                    {Tabs<-par1+rnorm(1,mean=0,sd=factor1)
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
                       Vd<-par4+rnorm(1,mean=0,sd=factor4)}},
                  
                   {Tabs<-par1+runif(1,min=-factor1,max=factor1)
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
                       Vd<-par4+runif(1,min=-factor4,max=factor4)}}, 
                    
                   {Tabs<-par1*rnorm(1,mean=0,sd=factor1)+par1
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
                       Vd<-par4*rnorm(1,mean=0,sd=factor4)+par4}},
                  
                   {Tabs<-par1*runif(1,min=-factor1,max=factor1)+par1
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
                       Vd<-par4*runif(1,min=-factor4,max=factor4)+par4}}         
              )
              time1<-PKtime$time
              parms<-c(Tabs=Tabs,Vm=Vm,Km=Km,Vd=Vd) 
              XX<-data.frame(lsoda(0,c(0,time1), defun, parms)) 
              C1.lsoda[[j]]<-data.frame(XX[2:(length(time1)+1),1],XX[2:(length(time1)+1),2])
              colnames(C1.lsoda[[j]])<-list("time","concentration")     
           } 
         }        
         PKindex[[i]]<-montecarlo(C1.lsoda,time,i,re) 
     }  
     PKindex<-as.data.frame(do.call("rbind",PKindex))
     rownames(PKindex)<-seq(nrow(PKindex)) 
     savefile(PKindex) 
  }
} 



## Legacy function
szero.nolagm <- function(...) {
  szero.nolag(...,MMe=TRUE)
}



