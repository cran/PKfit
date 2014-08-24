###
### copy this file to your working path and modify it as what you need;
###
### One compartment PK model extravascualr single dose first-order absorption 
### for multiple-dosed simulation: different doses, 1-compt PK model
### 
require(PKfit)
require(deSolve)

Subject = NULL   # N Subj's 
PKtime  = NULL   # times for sampling
ka      = NULL     
Vd      = NULL
kel     = NULL
F       = NULL
MD      = TRUE

options(warn=-1)
graphics.off()
    
par<-data.frame(Parameter=c("Total_subject#","ka","kel","Vd","F"),
                Initial=c(12,0.361,0.132,11.7,0.85))  ### (1) set default PK parameters for sim here
        cat("\n")       
        show(par)
        cat("\n")
        Subject<-par[1,2]
        par1   <-par[2,2]
        par2   <-par[3,2]
        par3   <-par[4,2]
        par4   <-par[5,2]
     
### PKtime for multiple-dose simulation
PKtime<-data.frame(time=seq(0,72,0.1))  ### (2) set time line for sim

### (3) define ODE for 1-compt, 1st-ordered abs/1st-ordered elimin. PK model
defun <- function(time, y, parms) { 
        dy1dt <- -parms["ka"] * y[1]
        dy2dt <-  parms["ka"] * y[1]/parms["Vd"] -parms["kel"]*y[2]
        list(c(dy1dt,dy2dt)) 
} 
      
### Dose   time of Dose
###  10     0 (the 1st day)
###  20    24 (the 2nd day)
###  30    48 (the 3rd day)

time<-PKtime$time
Dose1<- 10              ### (4) change doses if necessary
Dose2<- 20
Dose3<- 30
dosing.time<-c(24,48)   ### (5) change dosing time if necessary; also 'events' at below
###
### move the following 6 lines right before lsoda().
###
### yini<-c(dy1dt=Dose1,dy2dt=0)   ### set values for time = zero. for ODEs
### events <- data.frame(var="dy1dt",time=dosing.time,value=c(Dose2,Dose3),method="add")
           
sim.type<- 1   ### (6) select general simulation or monte-carlo sim here

###
par(mfrow=c(2,1),las=1)
   
if (sim.type ==1){    ### normal simulation starting here

err.type<- 3   ### (7) select error type for general simulation
pick<-err.type
      
      type<-switch(err.type, 
                  "No Error",
                  "Error = Normal Error", 
                  "Error = Uniform Error",
                  "Error = Normal Error*True Value",
                  "Error = Uniform Error*True Value")
      
      if (err.type ==1){
          Subject<-1   ### set subject=1 with no error sim.
          PKindex<-vector(Subject,mode="list")
          for(i in 1:Subject)  {
           cat("\n\n     << Subject:- #",i,">>" )
           cat("\n\n")                                                               
           cat("*******************************************************\n")        
           cat(" Summary Table                                       \n\n")   
           cat(" Model: 1-compartment, extravascular, multiple-dose,  \n")
           cat("     Error Type:", type,"\n\n")     
           sim<-matrix(c(ka,kel,Vd,F,par1,par2,par3,par4),4,2)
           dimnames(sim)<-list(c("ka","kel","Vd","F"),c("Simulated Values","Input Values"))}
           show(sim)                                                               
           cat("*******************************************************\n\n")      

             ka <-par1
             kel<-par2
             Vd <-par3
             F  <-par4
             parms<-c(ka=ka,kel=kel,Vd=Vd,F=F)
             
             yini<-c(dy1dt=Dose1*F,dy2dt=0)   ### set values for time = zero. for ODEs
             events <- data.frame(var="dy1dt",time=dosing.time,value=c(Dose2*F,Dose3*F),method="add")
             
             C1.lsoda<-data.frame(lsode(yini,c(0,time),defun,parms,rtol=1e-08,atol=1e-08,
                                  events=list(data=events)))
                
             good<-ifelse(C1.lsoda[2:(length(time)+1),3]<=1e-5,
                          0,
                          C1.lsoda[2:(length(time)+1),3])
             PKindex<-data.frame(i,
                                 C1.lsoda[2:(length(time)+1),1],
                                 good)
             colnames(PKindex)<-list("Subject","time","conc")
             show(PKindex)
             x<-C1.lsoda[2:(length(time)+1),1]
             y<-good
             dev.new();plotting.sim(i,x,y,MD)
          
       }   
  else {
  par.err<-data.frame(Parameter=c("ka","kel","Vd","F"),error_factor=c(.15,.15,.15,.15))  ### (8) change error factors for general sim.
         factor1<-par.err[1,2]
         factor2<-par.err[2,2]
         factor3<-par.err[3,2]
         factor4<-par.err[4,2]
         
             PKindex<-vector(Subject,mode="list")
             pick<- pick-1  # reset pick = 1, 2, 3 or 4 ('0' is zero error & excluded)
             for(i in 1:Subject){

                ka  <- err_types(pick,par1,factor1)
                kel <- err_types(pick,par2,factor2)
                Vd  <- err_types(pick,par3,factor3)
                F   <- err_types(pick,par4,factor4)
                
           cat("\n\n     << Subject:- #",i,">>" )
           cat("\n\n")                                                               
           cat("*******************************************************\n")        
           cat(" Summary Table                                       \n\n")   
           cat(" Model: 1-compartment, extravascular, multiple-dose,  \n")
           cat("     Error Type:", type,"\n\n")     
 
           sim<-matrix(c(ka,kel,Vd,F,par1,par2,par3,par4),4,2)
           dimnames(sim)<-list(c("ka","kel","Vd","F"),c("Simulated Values","Input Values"))
           show(sim)                                                               
           cat("*******************************************************\n\n")      
           parms<-c(ka=ka,kel=kel,Vd=Vd,F=F)
           
           yini<-c(dy1dt=Dose1*F,dy2dt=0)   ### set values for time = zero. for ODEs
           events <- data.frame(var="dy1dt",time=dosing.time,value=c(Dose2*F,Dose3*F),method="add")
             
           C1.lsoda<-data.frame(lsode(yini,c(0,time),defun,parms,rtol=1e-08,atol=1e-08,
                                events=list(data=events)))
           
           good<-ifelse(C1.lsoda[2:(length(time)+1),3]<=1e-5,
                        0,
                        C1.lsoda[2:(length(time)+1),3])
           PKindex<-data.frame(i,
                               C1.lsoda[2:(length(time)+1),1],
                               good)
           colnames(PKindex)<-list("Subject","time","conc")
           show(PKindex)
           x<-C1.lsoda[2:(length(time)+1),1]
           y<-good
           dev.new();plotting.sim(i,x,y,MD)
     } 
   }
}
if (sim.type ==2){    ### doing monte-carlo simulation here 

err.type<- 2   ### (7) select error type for monte-carlo simulation here
pcik<-err.type
     
     type<-switch(err.type, 
                  "Error = Normal Error", 
                  "Error = Uniform Error",
                  "Error = Normal Error*True Value",
                  "Error = Uniform Error*True Value")
     

     re<- 20   ### (8) set iteration # for monte-carlo sim.
    
         par.err<-data.frame(Parameter=c("ka","kel","Vd","F"),error_factor=c(.15,.15,.15,.15))  ### (8) change error factors for MC sim
         factor1<-par.err[1,2]
         factor2<-par.err[2,2]
         factor3<-par.err[3,2]
         factor4<-par.err[4,2]

         cat("\n\n")                                                               
         cat("*******************************************************\n")        
         cat("Summary Table - Monte-Carlo simulation runs          \n\n")   
         cat(" Model: 1-compartment, extravascular, multiple-dose,  \n")
         sim<-matrix(c(par1,par2,par3,par4,factor1,factor2,factor3,factor4),4,2)
         dimnames(sim)<-list(c("ka","kel","Vd","F"),c("Input Values","Error factor"))
         show(sim)                                                               
         cat("*******************************************************")
                  
     PKindex<-vector(Subject,mode="list")
     for(i in 1:Subject)  {
       cat("\n\n     << Subject:- #",i,">>\n" )
       C1.lsoda<-list()

       for (j in 1:re){
       
             ka  <- err_types(pick,par1,factor1)
             kel <- err_types(pick,par2,factor2)
             Vd  <- err_types(pick,par3,factor3)
             F   <- err_types(pick,par4,factor4)
                             
             parms<-c(ka=ka,kel=kel,Vd=Vd,F=F)
             
             yini<-c(dy1dt=Dose1*F,dy2dt=0)   ### set values for time = zero. for ODEs
             events <- data.frame(var="dy1dt",time=dosing.time,value=c(Dose2*F,Dose3*F),method="add")
             
             XX<-data.frame(lsode(yini,c(0,time),defun,parms,rtol=1e-08,atol=1e-08,
                                  events=list(data=events)))
             good<-ifelse(XX[2:(length(time)+1),3]<=1e-5,
                          0,
                          XX[2:(length(time)+1),3])
             C1.lsoda[[j]]<-data.frame(XX[2:(length(time)+1),1],good)
             colnames(C1.lsoda[[j]])<-list("time","concentration") 
         }   
             dev.new()     
             PKindex[[i]]<-montecarlo(C1.lsoda,time,i,re) 
     }  
  }
