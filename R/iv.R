### PKindex is usually the Dataset.


### Normal fitting
### One compartment PK model iv bolus single dose
### optional Michaelis-Menten Elimination
fbolus1 <- function(PKindex,
                    Dose=NULL,
                    Vm=NULL,Km=NULL, ## MMe=TRUE
                    Vd=NULL,
                    kel=NULL,        ## MMe=FALSE
                    MMe=FALSE)
{
   #options(warn=-1)
   
   ## Input dose and initial value for kel and Vd
   
   if (is.null(Dose)) {
     cat("Enter Dose value\n")
     Dose <- scan(nlines=1,quiet=TRUE)
   } 
   else {
     cat("Dose from arguments is = ",Dose,"\n")
   }
   
   if (MMe){
      if (is.null(Vm) || is.null(Km) || is.null(Vd) ) {
        par<-data.frame(Parameter=c("Vm","Km","Vd"),Initial=c(0))
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
   else {
      ## No MM elimination
      if (is.null(kel) || is.null(Vd)) {
        par<-data.frame(Parameter=c("kel","Vd"),Initial=c(0))
        par<-edit(par)
        repeat{
           if ( par[1,2] == 0 || par[2,2] ==0 ){
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
   
   if (!MMe) {
      ## User-supplied function w/o Michaelis-Mention elimination
      defun <- function(time, y, parms) { 
      dCpdt <- -parms["kel"] * y[1] 
      list(dCpdt) 
      } 
    
      modfun1 <- function(time,kel, Vd) {  
      out <- lsoda(Dose/Vd,c(0,time),defun,parms=c(kel=kel,Vd=Vd),
                   rtol=1e-3,atol=1e-5) 
      out[-1,2] 
      }
   } 
   else {
      ## User-supplied function with MM elimination
      defun<- function(time, y, parms) { 
      dCpdt <- -(parms["Vm"]/parms["Vd"])*y[1]/(parms["Km"]/parms["Vd"]+y[1]) 
      list(dCpdt)
      }

      modfun2 <- function(time,Vm,Km,Vd) { 
      out <- lsoda(Dose/Vd,c(0,time),defun,parms=c(Vm=Vm,Km=Km,Vd=Vd),
                   rtol=1e-3,atol=1e-5)
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
           out <- modfun2(PKindex$time[PKindex$Subject==i], par[1], par[2],par[3])
        } 
        else {
           ## No MM elimination
           out <- modfun1(PKindex$time[PKindex$Subject==i], par[1], par[2])
        }
        gift <- which( PKindex$conc[PKindex$Subject==i] != 0 )
        switch(pick,
             sum((PKindex$conc[PKindex$Subject==i][gift]-out[gift])^2),
             sum((PKindex$conc[PKindex$Subject==i][gift]-out[gift])^2/PKindex$conc[gift]),
             sum(((PKindex$conc[PKindex$Subject==i][gift]-out[gift])/PKindex$conc[gift])^2))
        }
        if (MMe) {
          gen<-genoud(objfun,nvars=3,max=FALSE,pop.size=30,
               max.generations=20,wait.generations=10,
               starting.value=c(par[1,2],par[2,2],par[3,2]),
               BFGS=FALSE,print.level=0,boundary.enforcement=2,
               Domains=matrix(c(1,1,1,100,100,100),3,2),
               MemoryMatrix=TRUE) 
        }
        else {
          gen<-genoud(objfun,nvars=2,max=FALSE,pop.size=30,
               max.generations=20,wait.generations=10,
               starting.value=c(par[1,2],par[2,2]),
               BFGS=FALSE,print.level=0,boundary.enforcement=2,
               Domains=matrix(c(0.01,0.01,100,100),2,2),
               MemoryMatrix=TRUE)  
         }
       cat("<< The value of parameter obtained from genetic algorithm >>\n\n")
       if (MMe) {
          namegen<-c("Vm","Km","Vd")
          outgen<-c(gen$par[1],gen$par[2],gen$par[3])
       } 
       else {
          ## No MM elimination
          namegen<-c("kel","Vd")
          outgen<-c(gen$par[1],gen$par[2])
       }
       print(data.frame(Parameter=namegen,Value=outgen))  
       F<-objfun(gen$par)
       
       if (MMe) {
        opt<-optim(c(gen$par[1],gen$par[2],gen$par[3]),objfun, method="Nelder-Mead")
        nameopt<-c("Vm","Km","Vd")
        outopt<-c(opt$par[1],opt$par[2],opt$par[3])
       }
       else {
        opt<-optim(c(gen$par[1],gen$par[2]),objfun,method="Nelder-Mead")  
        nameopt<-c("kel","Vd")
        outopt<-c(opt$par[1],opt$par[2])
       }
       cat("\n<< The value of parameter fitted by Nelder-Mead Simplex slgorithm >>\n\n")
       print(data.frame(Parameter=nameopt,Value=outopt))
       cat("\n<< Residual sum-of-squares and parameter values fitted by nls >>\n\n")
       if (MMe) {
         fm<-nls(conc~modfun2(time,Vm,Km,Vd),data=PKindex,subset=Subject==i,
                 start=list(Vm=opt$par[1],Km=opt$par[2],Vd=opt$par[3]),trace=TRUE,
                 nls.control(tol=1))
         cat("\n")
         plotting.non(PKindex, fm, i, pick, xaxis, yaxis)
       } 
       else {
        ## No MM elimination
         fm<-nls(conc ~ modfun1(time, kel, Vd), data=PKindex,
                 start=list(kel=opt$par[1],Vd=opt$par[2]),trace=TRUE,subset=Subject==i,
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
fbolus.mm <- function(PKindex,...) {
  fbolus1(PKindex,...,MMe=TRUE)
}

###Simulation
###One compartment PK model iv infusion single dose
###optional Michaelis-Menten Elimination
sbolus1 <- function(Subject=NULL,  # N Subj's 
                   PKtime=NULL,   # times for sampling
                   Dose=NULL,     # single dose
                   Vd=NULL,
                   kel=NULL, ## If not MM elimination
                   MMe=FALSE, ## michaelis-menten elimination?
                   Vm=NULL,Km=NULL)
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
   if ( !MMe) { 
     if (is.null(kel) || is.null(Vd)) {
       par<-data.frame(Parameter=c("kel","Vd"),Initial=c(0))
       par<-edit(par)
       repeat{
  if ( par[1,2] == 0 || par[2,2] ==0){
       cat("\n")
       cat("**********************************\n")
       cat(" Parameter value can not be zero. \n")
       cat(" Press enter to continue.         \n")
       cat("**********************************\n\n")
       readline()
       cat("\n")
       par<-edit(par)
       }   
  else{
       break
       return(edit(par))}
  } 
  cat("\n")       
  show(par)     
       cat("\n")
       par1<-par[1,2]
       par2<-par[2,2]
     } 
     else {
       par1 <- kel
       par2 <- Vd
     } 
   }
   else{
     if (is.null(Vm) || is.null(Km) || is.null(Vd)){
       par<-data.frame(Parameter=c("Vm","Km","Vd"),Initial=c(0))
       par<-edit(par)
       repeat{
  if ( par[1,2] == 0 || par[2,2]==0 || par[3,2]==0 ){
       cat("\n")
       cat("**********************************\n")
       cat(" Parameter value can not be zero. \n")
       cat(" Press enter to continue.         \n")
       cat("**********************************\n\n")
       readline()
       cat("\n")
       par<-edit(par)
       }   
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
       par1 <- Vm
       par2 <- Km
       par3 <- Vd
     }
   } 
   
   if ( ! MMe){
     defun<- function(time, y, parms) { 
     dCpdt <- -parms["kel"]* y[1]
     list(c(dCpdt)) 
     }
   }
   else{
     defun<- function(time, y, parms) { 
     dCpdt <- -(parms["Vm"]/parms["Vd"])*y[1]/(parms["Km"]/parms["Vd"]+y[1])
     list(c(dCpdt)) 
     }   
   }
   cat("\n")
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
       
           if (! MMe ){
             kel<-par1
             Vd<-par2   
             PKindex[[i]]<-sbolus1.out(PKtime,kel,Vd,defun,par1,par2,Dose,i,type) 
           } 
           else{
             Vm<-par1
             Km<-par2
             Vd<-par3       
             PKindex[[i]]<-sbolus.mm.out(PKtime,Vm,Km,Vd,defun,par1,par2,par3,Dose,i,type)
           }
         }       
         PKindex<- as.data.frame(do.call("rbind",PKindex))
         rownames(PKindex) <- seq(nrow(PKindex)) 
         savefile(PKindex)
     }
     else {
         if (! MMe){
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
          
           
           cat("\nEnter error factor for Vd\n")
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
           
           PKindex<-vector(Subject,mode="list")
           for( i in 1:Subject)  {
             cat("\n\n             << Subject",i,">>\n\n" )  
             switch(pick-1,
                   {kel<-par1+rnorm(1,mean=0,sd=factor1)
                    while(kel<=0){
                       kel<-par1+rnorm(1,mean=0,sd=factor1)}
                    Vd<-par2+rnorm(1,mean=0,sd=factor2)
                    while(Vd<=0){
                       Vd<-par2+rnorm(1,mean=0,sd=factor2)}},
                  
                   {kel<-par1+runif(1,min=-factor1,max=factor1)
                    while(kel<=0){
                       kel<-par1+runif(1,min=-factor1,max=factor1)}
                    Vd<-par2+runif(1,min=-factor2,max=factor2)
                    while(Vd<=0){
                       Vd<-par2+runif(1,min=-factor2,max=factor2)}}, 
                    
                   {kel<-par1*rnorm(1,mean=0,sd=factor1)+par1
                    while(kel<=0){
                       kel<-par1*rnorm(1,mean=0,sd=factor1)+par1}
                    Vd<-par2*rnorm(1,mean=0,sd=factor2)+par2
                    while(Vd<=0){
                       Vd<-par2*rnorm(1,mean=0,sd=factor2)+par2}},
                  
                   {kel<-par1*runif(1,min=-factor1,max=factor1)+par1
                    while(kel<=0){
                       kel<-par1*runif(1,min=-factor1,max=factor1)+par1}
                    Vd<-par2*runif(1,min=-factor2,max=factor2)+par2
                    while(Vd<=0){
                      Vd<-par2*runif(1,min=-factor2,max=factor2)+par2}}   
             )      
             PKindex[[i]]<-sbolus1.out(PKtime,kel,Vd,defun,par1,par2,Dose,i,type)
           }
         }
         else{
           cat("\n\nEnter error factor for Vm\n")
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
           
           cat("\nEnter error factor for Km\n")
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
                   {Vm<-par1+rnorm(1,mean=0,sd=factor1)
                    while(Vm<=0){
                       Vm<-par1+rnorm(1,mean=0,sd=factor1)}
                    Km<-par2+rnorm(1,mean=0,sd=factor2)
                    while(Km<=0){
                       Km<-par2+rnorm(1,mean=0,sd=factor2)}
                    Vd<-par3+rnorm(1,mean=0,sd=factor3)
                    while(Vd<=0){
                       Vd<-par3+rnorm(1,mean=0,sd=factor3)}},
                  
                   {Vm<-par1+runif(1,min=-factor1,max=factor1)
                    while(Vm<=0){
                       Vm<-par1+runif(1,min=-factor1,max=factor1)}
                    Km<-par2+runif(1,min=-factor2,max=factor2)
                    while(Km<=0){
                       Km<-par2+runif(1,min=-factor2,max=factor2)}
                    Vd<-par3+runif(1,min=-factor3,max=factor3)
                    while(Vd<=0){
                       Vd<-par3+runif(1,min=-factor3,max=factor3)}}, 
                    
                   {Vm<-par1*rnorm(1,mean=0,sd=factor1)+par1
                    while(Vm<=0){
                       Vm<-par1*rnorm(1,mean=0,sd=factor1)+par1}
                    Km<-par2*rnorm(1,mean=0,sd=factor2)+par2
                    while(Km<=0){
                       Km<-par2*rnorm(1,mean=0,sd=factor2)+par2}
                    Vd<-par3*rnorm(1,mean=0,sd=factor3)+par3
                    while(Vd<=0){
                       Vd<-par3*rnorm(1,mean=0,sd=factor3)+par3}},
                  
                   {Vm<-par1*runif(1,min=-factor1,max=factor1)+par1
                    while(Vm<=0){
                       Vm<-par1*runif(1,min=-factor1,max=factor1)+par1}
                    Km<-par2*runif(1,min=-factor2,max=factor2)+par2
                    while(Km<=0){
                       Km<-par2*runif(1,min=-factor2,max=factor2)+par2}
                    Vd<-par3*runif(1,min=-factor3,max=factor3)+par3
                    while(Vd<=0){
                       Vd<-par3*runif(1,min=-factor3,max=factor3)+par3}}   
             )      
             PKindex[[i]]<-sbolus.mm.out(PKtime,Vm,Km,Vd,defun,par1,par2,par3,Dose,i,type) 
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
    
     
     if (! MMe){
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
    
         cat("\nEnter error factor for Vd\n")
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
         
         
         cat("\n")
         cat("****************************************************\n")
         cat("Summary Table                                       \n")
         cat("Model: 1-Compartment, IV-Bolus, & Single-Dose Model \n")
         cat("Subject #:", Subject,"                              \n")
         cat("Error Type:", type,"                                \n") 
         cat("Simulation #:", re,"                                \n\n")
         sim<-matrix(c(par1,par2,factor1,factor2),2,2)
         dimnames(sim)<-list(c("kel","Vd"),c("Original","Error factor"))
         show(sim)   
         cat("****************************************************\n\n")
       }
       else{
         cat("\nEnter error factor for Vm\n")
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
         
         cat("\nEnter error factor for Km\n")
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
         
         cat("\n")
         cat("**********************************************\n")
         cat("Summary Table                                 \n")
         cat("Model: 1-Compartment, IV-Bolus, Single-Dose,  \n") 
         cat("       & Michaelis-Menten Elimination Model   \n") 
         cat("Subject #:", Subject,"                        \n")
         cat("Error Type:", type,"                          \n") 
         cat("Simulation #:", re,"                          \n\n")
         sim<-matrix(c(par1,par2,par3,factor1,factor2,factor3),3,2)
         dimnames(sim)<-list(c("Vm","Km","Vd"),c("Original","Error factor"))
         show(sim)   
         cat("**********************************************\n\n")
       }  
       PKindex<-vector(Subject,mode="list")
       for( i in 1:Subject)  {
       cat("\n\n             << Subject",i,">>\n\n" ) 
       C1.lsoda<-list()

         for (j in 1:re){
            if (! MMe){
              switch(pick, 
                    {kel<-par1+rnorm(1,mean=0,sd=factor1)
                     while(kel<=0){
                        kel<-par1+rnorm(1,mean=0,sd=factor1)}
                     Vd<-par2+rnorm(1,mean=0,sd=factor2)
                     while(Vd<=0){
                        Vd<-par2+rnorm(1,mean=0,sd=factor2)}},
               
                    {kel<-par1+runif(1,min=-factor1,max=factor1)
                     while(kel<=0){
                        kel<-par1+runif(1,min=-factor1,max=factor1)}
                     Vd<-par2+runif(1,min=-factor2,max=factor2)
                     while(Vd<=0){
                        Vd<-par2+runif(1,min=-factor2,max=factor2)}},
                        
                    {kel<-par1*rnorm(1,mean=0,sd=factor1)+par1
                    while(kel<=0){
                       kel<-par1*rnorm(1,mean=0,sd=factor1)+par1}
                    Vd<-par2*rnorm(1,mean=0,sd=factor2)+par2
                    while(Vd<=0){
                       Vd<-par2*rnorm(1,mean=0,sd=factor2)+par2}},
                  
                   {kel<-par1*runif(1,min=-factor1,max=factor1)+par1
                    while(kel<=0){
                       kel<-par1*runif(1,min=-factor1,max=factor1)+par1}
                    Vd<-par2*runif(1,min=-factor2,max=factor2)+par2
                    while(Vd<=0){
                      Vd<-par2*runif(1,min=-factor2,max=factor2)+par2}}     
              )    
              time1<-PKtime$time
              parms<-c(kel=kel,Vd=Vd) 
              XX<-data.frame(lsoda(Dose/Vd,c(0,time1),defun,parms))
              C1.lsoda[[j]]<-data.frame(XX[2:(length(time1)+1),1],XX[2:(length(time1)+1),2])
              colnames(C1.lsoda[[j]])<-list("time","concentration") 
           }
           else{
              switch(pick, 
                    {Vm<-par1+rnorm(1,mean=0,sd=factor1)
                     while(Vm<=0){
                        Vm<-par1+rnorm(1,mean=0,sd=factor1)}
                     Km<-par2+rnorm(1,mean=0,sd=factor2)
                     while(Km<=0){
                        Km<-par2+rnorm(1,mean=0,sd=factor2)}
                     Vd<-par3+rnorm(1,mean=0,sd=factor3)
                     while(Vd<=0){
                        Vd<-par3+rnorm(1,mean=0,sd=factor3)}},
               
                    {Vm<-par1+runif(1,min=-factor1,max=factor1)
                     while(Vm<=0){
                        Vm<-par1+runif(1,min=-factor1,max=factor1)}
                     Km<-par2+runif(1,min=-factor2,max=factor2)
                     while(Km<=0){
                        Km<-par2+runif(1,min=-factor2,max=factor2)}
                     Vd<-par3+runif(1,min=-factor3,max=factor3)
                     while(Vd<=0){
                        Vd<-par3+runif(1,min=-factor3,max=factor3)}},
                        
                    {Vm<-par1*rnorm(1,mean=0,sd=factor1)+par1
                    while(Vm<=0){
                       Vm<-par1*rnorm(1,mean=0,sd=factor1)+par1}
                    Km<-par2*rnorm(1,mean=0,sd=factor2)+par2
                    while(Km<=0){
                       Km<-par2*rnorm(1,mean=0,sd=factor2)+par2}
                    Vd<-par3*rnorm(1,mean=0,sd=factor3)+par3
                    while(Vd<=0){
                       Vd<-par3*rnorm(1,mean=0,sd=factor3)+par3}},
                  
                   {Vm<-par1*runif(1,min=-factor1,max=factor1)+par1
                    while(Vm<=0){
                       Vm<-par1*runif(1,min=-factor1,max=factor1)+par1}
                    Km<-par2*runif(1,min=-factor2,max=factor2)+par2
                    while(Km<=0){
                       Km<-par2*runif(1,min=-factor2,max=factor2)+par2}
                    Vd<-par3*runif(1,min=-factor3,max=factor3)+par3
                    while(Vd<=0){
                       Vd<-par3*runif(1,min=-factor3,max=factor3)+par3}}       
              )
              time1<-PKtime$time
              parms<-c(Vm=Vm,Km=Km,Vd=Vd)  
              XX<-data.frame(lsoda(Dose/Vd,c(0,time1),defun,parms))
              C1.lsoda[[j]]<-data.frame(XX[2:(length(time1)+1),1],XX[2:(length(time1)+1),2])
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
sbolus.mm <- function(...) {
  sbolus1(...,MMe=TRUE)
}

  