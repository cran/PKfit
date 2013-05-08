###Simulation
###One compartment PK model iv infusion single dose
###optional Michaelis-Menten Elimination
sinfu1 <- function(Subject=NULL,  # N Subj's 
                   PKtime=NULL,   # times for sampling
                   Dose=NULL,     # single dose
                   Tinf=NULL,     # infusion time (length)
                   Vd=NULL,
                   kel=NULL,  ## If not MM elimination
                   MMe=FALSE, ## michaelis-menten elimination?
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

  if (is.null(Tinf)) {
    cat("\nEnter value for infusion time\n")
    Tinf<-scan(nlines=1,quiet=TRUE)
    cat("\n")
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
       par1 <- Vm
       par2 <- Km
       par3 <- Vd
    }
  } 

  if ( ! MMe){
    defun<- function(time,y,parms) { 
      if(time<=Tinf)         
        dCpdt <- (Dose/Tinf)/parms["Vd"]-parms["kel"]*y[1]
      else
        dCpdt <- -parms["kel"]*y[1]
      list(dCpdt) 
    } 
  }
  else{
    defun<- function(time, y, parms) { 
      if(time<=Tinf)  
        dCpdt <- (Dose/Tinf)/parms["Vd"]-(parms["Vm"]/parms["Vd"])*y[1]/(parms["Km"]/parms["Vd"]+y[1]) 
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
             kel<-par1
             Vd<-par2   
             PKindex[[i]]<-sinfu1.out(PKtime,kel,Vd,defun,par1,par2,Dose,i,type) 
           } 
           else{
             Vm<-par1
             Km<-par2
             Vd<-par3       
             PKindex[[i]]<-sinfu.mm.out(PKtime,Vm,Km,Vd,defun,par1,par2,par3,Dose,i,type)
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
                cat(" Press Enter to continue.         \n")
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
                cat(" Press Enter to continue.         \n")
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
             PKindex[[i]]<-sinfu1.out(PKtime,kel,Vd,defun,par1,par2,Dose,i,type)   
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
                cat(" Press Enter to continue.         \n")
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
             PKindex[[i]]<-sinfu.mm.out(PKtime,Vm,Km,Vd,defun,par1,par2,par3,Dose,i,type) 
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
     pick <- menu(file.menu, title = "<< Error Types >>")
     
     type<-switch(pick, 
                  "Error = Normal Error", 
                  "Error = Uniform Error",
                  "Error = Normal Error*True Value",
                  "Error = Uniform Error*True Value")
     
     cat("\n\nHow many interations do you want to run?\n")
     re<-scan(nlines=1,quiet=TRUE)
     
       if (! MMe){
         cat("\nEnter error factor for kel\n")
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
         
         cat("\nEnter error factor for Vd\n")
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
         
         cat("\n")
         cat("*******************************************************\n")
         cat("Summary Table                                          \n")
         cat("Model: 1-Compartment, IV-Infusion, & Single-Dose Model \n")
         cat("Subject #:", Subject,"                                 \n")
         cat("Error Type:", type,"                                   \n")
         cat("Simulation #:", re,"                                   \n\n")
         sim<-matrix(c(par1,par2,factor1,factor2),2,2)
         dimnames(sim)<-list(c("kel","Vd"),c("Selected","Error factor"))
         show(sim)   
         cat("*******************************************************\n\n")
       }
       else{
         cat("\nEnter error factor for Vm\n")
         factor1<-scan(nlines=1,quiet=TRUE)
         cat("\nEnter error factor for Km\n")
         factor2<-scan(nlines=1,quiet=TRUE)
         cat("\nEnter error factor for Vd\n")
         factor3<-scan(nlines=1,quiet=TRUE)
         cat("\n")
         cat("************************************************\n")
         cat("Summary Table                                   \n")
         cat("Model: 1-Compartment, IV-Infusion, Single-Dose, \n") 
         cat("       & Michaelis-Menten Elimination Model     \n")
         cat("Subject #:", Subject,"                          \n")
         cat("Error Type:", type,"                            \n")
         cat("Simulation #:", re,"                            \n\n")
         sim<-matrix(c(par1,par2,par3,factor1,factor2,factor3),3,2)
         dimnames(sim)<-list(c("Vm","Km","Vd"),c("Selected","Error factor"))
         show(sim)   
         cat("************************************************\n\n")
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
              XX<-data.frame(lsoda(0,c(0,time1),defun,parms,rtol=1e-6,atol=1e-6))
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
              XX<-data.frame(lsoda(0,c(0,time1),defun,parms,rtol=1e-6,atol=1e-6))
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
