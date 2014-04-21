###Simulation
###One compartment PK model extravascualr single dose first-order absorption 
###optional Michaelis-Menten Elimination
###optional lag time
sfirst.nolag <- function(Subject=NULL,  # N Subj's 
                         PKtime=NULL,   # times for sampling
                         Dose=NULL,     # single dose
                         Tlag_time=NULL,
                         ka=NULL,     
                         Vd=NULL,
                         kel=NULL, ## If not MM elimination
                         MMe=FALSE, ## michaelis-menten elimination?
                         Vm=NULL,Km=NULL,
                         Tlag=FALSE)
{
   options(warn=-1)
    
   cat("\n First enter all paramaters used for simulation profile.\n");readline(" Press Enter to continue...\n\n")
   if (MMe) {
      if(Tlag){
        par<-data.frame(Parameter=c("Total_subject#","Dose","Tlag","ka","Vm","Km","Vd"),
                          Initial=c(24,300,0.25,0.36,0.13,0.31,11.7))
        par<-edit(par)                                                       
        repeat{
            if (par[1,2] <= 0 || par[2,2] <=0 || par[3,2]<=0 || par[4,2]<=0 || par[5,2]<=0|| par[6,2]<=0|| par[7,2]<=0){
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
        Tlag_time<-par[3,2]
        par1<-par[4,2]
        par2<-par[5,2]
        par3<-par[6,2]
        par4<-par[7,2]
        }
       else{
        par<-data.frame(Parameter=c("Total_subject#","Dose","ka","Vm","Km","Vd"),Initial=c(24,300,0.36,0.13,0.31,11.7))
        par<-edit(par)                                                       
        repeat{
            if (par[1,2] <= 0 || par[2,2] <=0 || par[3,2]<=0 || par[4,2]<=0 || par[5,2]<=0|| par[6,2]<=0){
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
        par4<-par[6,2]       
       }
   }

   else{
      if(Tlag){
        par<-data.frame(Parameter=c("Total_subject#","Dose","Tlag","ka","kel","Vd"),Initial=c(24,300,0.25,0.36,0.13,11.7))
        par<-edit(par)                                                       
        repeat{
            if (par[1,2] <= 0 || par[2,2] <=0 || par[3,2]<=0 || par[4,2]<=0 || par[5,2]<=0|| par[6,2]<=0){
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
        Tlag_time<-par[3,2]
        par1<-par[4,2]
        par2<-par[5,2]
        par3<-par[6,2]
      }
      else{
        par<-data.frame(Parameter=c("Total_subject#","Dose","ka","kel","Vd"),Initial=c(24,300,0.36,0.13,11.7))
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
     }
     
   readline("\n Next enter time points. Press Enter to continue.\n\n")
   PKtime<-data.frame(time=c(0))
   ### PKtime<-data.frame(time=c(0,.1,.2,.3,.4,.6,.8,1,2,4,6,8,12,14,16,18,24))
   PKtime<-edit(PKtime)
   cat("\n")
   show(PKtime);cat("\n\n")
  
   if (MMe) {
      if (Tlag){
        ## User-supplied function w Michaelis-Mention elimination & w lag time
        defun<- function(time, y, parms) { 
           if(time <= Tlag_time) {
              dy1dt<-0
              dy2dt<-0
           }
           else {
              dy1dt <- -parms["ka"] * y[1]
              dy2dt <-  parms["ka"] * y[1]/parms["Vd"] - parms["Vm"]*y[2]/(parms["Km"]+y[2])
           }
           list(c(dy1dt,dy2dt)) 
         } 
      }
      else{
        ## User-supplied function with MM elimination w/o lag time
        defun<- function(time, y, parms) { 
          dy1dt <- -parms["ka"] * y[1]
          dy2dt <-  parms["ka"] * y[1]/parms["Vd"] - parms["Vm"]*y[2]/(parms["Km"]+y[2])
          list(c(dy1dt,dy2dt)) 
        } 
      }  
   } 
   else {
      if (Tlag){
        ## User-supplied function w/o MM elimination w lag time
        defun<- function(time, y, parms) { 
           if(time<=Tlag_time) {
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
   pick <- menu(file.menu, title = "<< Simulation Types >>")
   if (pick ==1){    ### normal simulation starting here
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
          for(i in 1:Subject)  {
             cat("\n     << Subject:- #",i,">>\n\n" )
           
             if ( MMe ){
                ka<-par1
                Vm<-par2  
                Km<-par3
                Vd<-par4 
                PKindex[[i]]<-sfirst.mm.out(PKtime,ka,Vm,Km,Vd,defun,par1,par2,par3,par4,Dose,i,type,Tlag)      
             } 
             else{
                ka <-par1
                kel<-par2
                Vd <-par3       
                PKindex[[i]]<-sfirst1.out(PKtime,ka,kel,Vd,defun,par1,par2,par3,Dose,i,type,Tlag)   
             }
          }       
          PKindex<- as.data.frame(do.call("rbind",PKindex))
          rownames(PKindex) <- seq(nrow(PKindex)) 
          savefile(PKindex)
       }   
      else {
         if (!MMe){
         par.err<-data.frame(Parameter=c("ka","kel","Vd"),error_factor=c(.15,.15,.15))
         par.err<-edit(par.err)
         repeat{
         if (par.err[1,2] <= 0 || par.err[2,2] <=0 || par.err[3,2]<=0){
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
             for(i in 1:Subject)  {
                cat("\n     << Subject:- #",i,">>\n\n" )
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
      else{   ### if(MMe) here
         par.err<-data.frame(Parameter=c("ka","Vm","Km","Vd"),error_factor=c(.15,.15,.15,.15))
         par.err<-edit(par.err)
         repeat{
         if (par.err[1,2] <= 0 || par.err[2,2] <=0 || par.err[3,2]<=0 || par.err[4,2]<=0){
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
         factor4<-par.err[4,2]
             
             PKindex<-vector(Subject,mode="list");cat("\n")
             for(i in 1:Subject)  {
               cat("\n     << Subject:- #",i,">>\n\n" )
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
  else if (pick ==2){    ### starting monte-carlo simulation here 
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
     
       if (MMe){
         par.err<-data.frame(Parameter=c("ka","Vm","Km","Vd"),error_factor=c(.15,.15,.15,.15))
         par.err<-edit(par.err)
         repeat{
         if (par.err[1,2] <= 0 || par.err[2,2] <=0 || par.err[3,2]<=0 || par.err[4,2]<=0){
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
         factor4<-par.err[4,2]       
          
          cat("\n")
          cat("************************************************************************\n")
          cat("Summary Table                                                           \n")
          cat("Model: 1-Compartment, Extravascular, Single-Dose, 1-Ordered Absorption, \n")
          cat("       and Michaelis-Menten Elimination with Lag Time Model             \n")
          cat("Subject #:", Subject,"                                                  \n")
          cat("Error Type:", type,"                                                    \n")
          cat("Simulation #:", re,"                                                    \n\n")
          sim<-matrix(c(par1,par2,par3,par4,factor1,factor2,factor3,factor4),4,2)
          dimnames(sim)<-list(c("ka","Vm","Km","Vd"),c("Selected","Error factor"))
          show(sim)   
          cat("************************************************************************\n\n")
       }
       else{
         par.err<-data.frame(Parameter=c("ka","kel","Vd"),error_factor=c(.15,.15,.15))
         par.err<-edit(par.err)
         repeat{
         if (par.err[1,2] <= 0 || par.err[2,2] <=0 || par.err[3,2]<=0){
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
         cat("****************************************************\n")
         cat("Summary Table                                       \n")
         cat("Model: 1-Compartment, Extravascular, Single-Dose,   \n")
         cat("       and 1-Ordered Absorption with Lag Time Model \n")
         cat("Subject #:", Subject,"                              \n")
         cat("Error Type:", type,"                                \n")
         cat("Simulation #:", re,"                                \n\n")
         sim<-matrix(c(par1,par2,par3,factor1,factor2,factor3),3,2)
         dimnames(sim)<-list(c("ka","kel","Vd"),c("Selected","Error factor"))
         show(sim)   
         cat("****************************************************\n\n")
         }
     cat("\n")
     PKindex<-vector(Subject,mode="list");cat("\n")
     for(i in 1:Subject)  {
       cat("\n     << Subject:- #",i,">>\n\n" )
       C1.lsoda<-list()

         for (j in 1:re){
     
           if (!MMe ){
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
              XX<-data.frame(lsoda(c(Dose,0),c(0,time1), defun, parms,rtol=1e-6,atol=1e-6))
              good<-ifelse(XX[2:(length(time1)+1),3]<=1e-5,
                           0,
                           XX[2:(length(time1)+1),3])
              C1.lsoda[[j]]<-data.frame(XX[2:(length(time1)+1),1],good)
              colnames(C1.lsoda[[j]])<-list("time","concentration") 
           }
           else{   ### if(MMe) here
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
              XX<-data.frame(lsoda(c(Dose,0),c(0,time1), defun, parms,rtol=1e-6,atol=1e-6))
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