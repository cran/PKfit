###Simulation
###One compartment PK model iv bolus single dose
###optional Michaelis-Menten Elimination

sbolus1 <- function(Subject=NULL,  # N Subj's 
                    PKtime=NULL,   # times for sampling
                    Dose=NULL,     # single dose
                    Vd=NULL,
                    kel=NULL,      ## If not MM elimination
                    MMe=FALSE,     ## michaelis-menten elimination?
                    Vm=NULL,Km=NULL)
{
   options(warn=-1)
   cat("\n First enter all paramaters used for simulation profile.\n");readline(" Press Enter to continue...\n\n")

   if (!MMe) {
      par<-data.frame(Parameter=c("Total_subject#","Dose","kel","Vd"),Initial=c(24,300,0.31,11.7))
      par<-edit(par)                                                       
      repeat{
          if (par[1,2] <= 0 || par[2,2] <=0 || par[3,2]<=0 || par[4,2]<=0){
            cat("\n")
            cat("**********************************\n")
            cat(" Any parameter value can not be equal to or less than zero.\n")
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
      Subject<-par[1,2]
      Dose<-par[2,2]
      par1<-par[3,2]
      par2<-par[4,2]
   }
   else{
      par<-data.frame(Parameter=c("Total_subject#","Dose","Vm","Km","Vd"),Initial=c(24,300,0.17,4.84,38.2))
      par<-edit(par)                                                       
      repeat{
          if (par[1,2] <= 0 || par[2,2] <=0 || par[3,2]<=0 || par[4,2]<=0 || par[5,2]<=0){
            cat("\n")
            cat("**********************************\n")
            cat(" Any parameter value can not be equal to or less than zero.\n")
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
      Subject<-par[1,2]
      Dose<-par[2,2]
      par1<-par[3,2]
      par2<-par[4,2]
      par3<-par[5,2]
}
   readline("\n Next enter time points. Press Enter to continue.\n\n")
   ### PKtime<-data.frame(time=c(0))
   PKtime<-data.frame(time=c(0,.1,.2,.3,.4,.5,.6,.8,1,2,4,6,8,12,14,16,18,24))
   PKtime<-edit(PKtime)
   cat("\n")
   show(PKtime);cat("\n\n")   
   
   if (!MMe){
     defun<- function(time, y, parms) { 
     dCpdt <- -parms["kel"]* y[1]
     list(c(dCpdt)) 
     }
   }
   else{
     defun<- function(time, y, parms) { 
     dCpdt <- -parms["Vm"]*y[1]/(parms["Km"]+y[1])
     list(c(dCpdt)) 
     }   
   }
   cat("\n")
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
         PKindex<-vector(Subject,mode="list");cat("\n")
         for( i in 1:Subject)  {
           cat("\n     << Subject:- #",i,">>\n\n" )
       
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
         if (!MMe){
           par.err<-data.frame(Parameter=c("kel","Vd"),error_factor=c(.15,.15))
           par.err<-edit(par.err)
           repeat{
           if (par.err[1,2] <=0 || par.err[2,2] <=0){
               cat("\n")
               cat("**********************************\n")
               cat(" Parameter values can be equal or less than zero. \n")
               cat(" Press Enter to continue.         \n")
               cat("**********************************\n\n")
               readline()
               cat("\n")
               par.err<-edit(par.err)}   
             else{
               break
               return(edit(par.err))}
           }
           
           factor1<-par.err[1,2]
           factor2<-par.err[2,2]
           
           PKindex<-vector(Subject,mode="list");cat("\n")
           for( i in 1:Subject)  {
             cat("\n     << Subject:- #",i,">>\n\n" )
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
           par.err<-data.frame(Parameter=c("Vm","Km","Vd"),error_factor=c(.15,.15,.15))
           par.err<-edit(par.err)
           repeat{
           if (par.err[1,2] <=0 || par.err[2,2] <=0|| par.err[3,2] <=0){
               cat("\n")
               cat("**********************************\n")
               cat(" Parameter values can be equal or less than zero. \n")
               cat(" Press Enter to continue.         \n")
               cat("**********************************\n\n")
               readline()
               cat("\n")
               par.err<-edit(par.err)}   
             else{
               break
               return(edit(par.err))}
           }
           
           factor1<-par.err[1,2]
           factor2<-par.err[2,2]
           factor3<-par.err[3,2]
           
           PKindex<-vector(Subject,mode="list");cat("\n")
           for( i in 1:Subject)  {
             cat("\n     << Subject:- #",i,">>\n\n" )
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
     
     cat("\n\nHow many iterations to run for each subject?\n")
     re<-scan(nlines=1,quiet=TRUE)
     
     if (!MMe){
         par.err<-data.frame(Parameter=c("kel","Vd"),error_factor=c(.15,.15))
         par.err<-edit(par.err)
         repeat{
         if (par.err[1,2] <=0 || par.err[2,2] <=0){
             cat("\n")
             cat("**********************************\n")
             cat(" Parameter values can be equal or less than zero. \n")
             cat(" Press Enter to continue.         \n")
             cat("**********************************\n\n")
             readline()
             cat("\n")
             par.err<-edit(par.err)}   
           else{
             break
             return(edit(par.err))}
         }
         
         factor1<-par.err[1,2]
         factor2<-par.err[2,2]
         
         cat("\n")
         cat("****************************************************\n")
         cat(" Summary Table                                       \n")
         cat(" Model: 1-compartment, iv bolus, single-dose model \n")
         cat(" Subject #:", Subject,"                              \n")
         cat(" Error Type:", type,"                                \n") 
         cat(" Simulation #:", re,"                                \n\n")
         sim<-matrix(c(par1,par2,factor1,factor2),2,2)
         dimnames(sim)<-list(c("kel","Vd"),c("Selected","Error factor"))
         show(sim)   
         cat("****************************************************\n\n")
       }
       else{
         par.err<-data.frame(Parameter=c("Vm","Km","Vd"),error_factor=c(.15,.15,.15))
         par.err<-edit(par.err)
         repeat{
         if (par.err[1,2] <=0 || par.err[2,2] <=0|| par.err[3,2] <=0){
             cat("\n")
             cat("**********************************\n")
             cat(" Parameter values can be equal or less than zero. \n")
             cat(" Press Enter to continue.         \n")
             cat("**********************************\n\n")
             readline()
             cat("\n")
             par.err<-edit(par.err)}   
           else{
             break
             return(edit(par.err))}
         }
         
         factor1<-par.err[1,2]
         factor2<-par.err[2,2]
         factor3<-par.err[3,2]
         
         cat("\n")
         cat("**********************************************\n")
         cat(" Summary Table                                 \n")
         cat(" Model: 1-compartment, iv bolus, single-dose,  \n") 
         cat("        & Michaelis-Menten elimination model   \n") 
         cat(" Subject #:", Subject,"                        \n")
         cat(" Error Type:", type,"                          \n") 
         cat(" Simulation #:", re,"                          \n\n")
         sim<-matrix(c(par1,par2,par3,factor1,factor2,factor3),3,2)
         dimnames(sim)<-list(c("Vm","Km","Vd"),c("Selected","Error factor"))
         show(sim)   
         cat("**********************************************\n\n")
       }  
       PKindex<-vector(Subject,mode="list");cat("\n")
       for( i in 1:Subject)  {
       cat("\n     << Subject:- #",i,">>\n\n" )
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
              XX<-data.frame(lsoda(Dose/Vd,c(0,time1),defun,parms,rtol=1e-6,atol=1e-6))
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
              XX<-data.frame(lsoda(Dose/Vd,c(0,time1),defun,parms,rtol=1e-6,atol=1e-6))
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
