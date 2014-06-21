###Simulation
###One compartment PK model iv infusion single dose
###optional Michaelis-Menten Elimination

sinfu1 <- function(Subject=NULL,  # N Subj's 
                   PKtime=NULL,   # times for sampling
                   Dose=NULL,     # single dose
                   Tinf=NULL,     # infusion time (length)
                   Vd=NULL,
                   kel=NULL,      ## If not MM elimination
                   MMe=FALSE,     ## michaelis-menten elimination?
                   Vm=NULL,Km=NULL,
                   nDose=NULL,    # for multiple-dose
                   Tau=NULL,
                   MD=FALSE)
{
   options(warn=-1)
   sim.outputs_to_txt<-sim.outputs_to_txt
   sim.plots_to_pdf<-sim.plots_to_pdf
     
cat("\n First enter all paramaters used for simulation profile.\n");readline(" Press Enter to continue...\n\n")
if(MD) {
if (!MMe) { 
  if(file.exists("sinfu1_md.csv")){
      par<-read.csv(file="sinfu1_md.csv",row.names=NULL,header=TRUE)
      par<-edit(par)}
  else{   
      par<-data.frame(Parameter=c("Total_subject#","Dose","Tau","#Dose","Tinf","kel","Vd"),Initial=c(24,300,12,10,1.0,0.12,11.7))
      par<-edit(par)                                                       
      repeat{
          if (par[1,2]<=0 || par[2,2]<=0 || par[3,2]<=0 || par[4,2]<=0|| par[5,2]<=0|| par[6,2]<=0|| par[7,2]<=0){
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
     }
      write.csv(par,file="sinfu1_md.csv",row.names=FALSE,col.names=TRUE)
      cat("\n")       
      show(par)
      
      cat("\n")
      Subject<-par[1,2]
      Dose   <-par[2,2]
      Tau    <-par[3,2]
      nDose  <-par[4,2]
      Tinf<-par[5,2]
      par1<-par[6,2]
      par2<-par[7,2]
     
  }
else{
  if(file.exists("sinfu1_mm_md.csv")){
    par<-read.csv(file="sinfu1_mm_md.csv",row.names=NULL,header=TRUE)
    par<-edit(par)}
  else{ 
    par<-data.frame(Parameter=c("Total_subject#","Dose","Tau","#Dose","Tinf","Vm","Km","Vd"),Initial=c(24,300,12,10,0.5,2.12,4.74,11.7))
    par<-edit(par)                                                       
    repeat{
        if (par[1,2]<=0 || par[2,2]<=0 || par[3,2]<=0 || par[4,2]<=0|| par[5,2]<=0|| par[6,2]<=0|| par[7,2]<=0|| par[8,2]<=0){
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
   }
    write.csv(par,file="sinfu1_mm_md.csv",row.names=FALSE,col.names=TRUE)
    cat("\n")       
    show(par)
    
    cat("\n")
    Subject<-par[1,2]
    Dose <-par[2,2]
    Tau  <-par[3,2]
    nDose<-par[4,2]
    Tinf<-par[5,2]
    par1<-par[6,2]
    par2<-par[7,2]
    par3<-par[8,2]
}

### no PKtime require for multiple-dose simulation
PKtime<-data.frame(time=seq(0,Tau*nDose-1,0.1))  ### to see more output details, smaller time steps. -YJ
   
  if (!MMe){
    defun<- function(time,y,parms) {
     if(time%%Tau<=Tinf)
           dCpdt <- (Dose/Tinf)/parms["Vd"]-parms["kel"]*y[1]
         else
           dCpdt <- -parms["kel"]*y[1]
      list(dCpdt)
    }
  }
  else{
    defun<- function(time, y, parms) {
      if(time%%Tau<=Tinf)
        dCpdt <- (Dose/Tinf)/parms["Vd"]-parms["Vm"]*y[1]/(parms["Km"]+y[1]) 
      else
        dCpdt <- -parms["Vm"]*y[1]/(parms["Km"]+y[1])
      list(dCpdt)
    }
  }
}
else{  
if (!MMe) { 
  if(file.exists("sinfu1.csv")){
      par<-read.csv(file="sinfu1.csv",row.names=NULL,header=TRUE)
      par<-edit(par)}
  else{   
      par<-data.frame(Parameter=c("Total_subject#","Dose","Tinf","kel","Vd"),Initial=c(24,300,1.0,0.12,11.7))
      par<-edit(par)                                                       
      repeat{
          if (par[1,2]<=0 || par[2,2]<=0 || par[3,2]<=0 || par[4,2]<=0|| par[5,2]<=0){
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
     }
      write.csv(par,file="sinfu1.csv",row.names=FALSE,col.names=TRUE)
      cat("\n")       
      show(par)
      
      cat("\n")
      Subject<-par[1,2]
      Dose<-par[2,2]
      Tinf<-par[3,2]
      par1<-par[4,2]
      par2<-par[5,2]
     
  }
else{
  if(file.exists("sinfu1_mm.csv")){
    par<-read.csv(file="sinfu1_mm.csv",row.names=NULL,header=TRUE)
    par<-edit(par)}
  else{ 
    par<-data.frame(Parameter=c("Total_subject#","Dose","Tinf","Vm","Km","Vd"),Initial=c(24,300,0.5,2.31,4.74,11.7))
    par<-edit(par)                                                       
    repeat{
        if (par[1,2]<=0 || par[2,2]<=0 || par[3,2]<=0 || par[4,2]<=0|| par[5,2]<=0|| par[6,2]<=0){
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
   }
    write.csv(par,file="sinfu1_mm.csv",row.names=FALSE,col.names=TRUE)
    cat("\n")       
    show(par)
    
    cat("\n")
    Subject<-par[1,2]
    Dose<-par[2,2]
    Tinf<-par[3,2]
    par1<-par[4,2]
    par2<-par[5,2]
    par3<-par[6,2]
}

   readline("\n Next enter or load time points. Press Enter to continue.\n\n")
   if(file.exists("sim_times.csv")){
      PKtime<-read.csv(file="sim_times.csv",row.names=NULL,header=TRUE)
      PKtime<-edit(PKtime)}
   else{
      ### PKtime<-data.frame(time=c(0))
      PKtime<-data.frame(time=c(0,.1,.2,.3,.4,.6,.8,1,2,4,6,8,12,14,16,18,24,48,72))
      PKtime<-edit(PKtime)}
   
   write.csv(PKtime,file="sim_times.csv",row.names=FALSE,col.names=TRUE)
   cat("\n")
   show(PKtime);cat("\n\n")
   
  if (!MMe){
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
        dCpdt <- (Dose/Tinf)/parms["Vd"]-parms["Vm"]*y[1]/(parms["Km"]+y[1]) 
      else
        dCpdt <- -parms["Vm"]*y[1]/(parms["Km"]+y[1])
      list(dCpdt)
    }
  }   
}
  file.menu <- c("Simulation with Error",
                 "Monte Carlo Simulation")
  pick <- menu(file.menu, title = "<< Simulation Types >>")

###
dev.new()
par(mfrow=c(2,1),las=1)
pdf_activate=FALSE  ### set pdf device activate? as FALSE at beginning
### pdf(sim.plots_to_pdf,paper="a4")
###
###
### log to outputs.txt here
###
zz <- file(sim.outputs_to_txt, open="wt")
sink(zz,split=TRUE)   ### use sink(zz.split=TURE) will output to the txt file, as well as the screen at the same time. YJ
description_version()
sink()
### 
  
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
     
     sink(zz,split=TRUE)
     
       if (pick ==1){
         Subject<-1          ### does not make any sense to simulate more than 1 subject with this option.
         PKindex<-vector(Subject,mode="list")
         for( i in 1:Subject)  {
           cat("\n\n     << Subject:- #",i,">>" )
           cat("\n\n")                                                               
           cat("*******************************************************\n")        
           cat(" Summary Table                                          \n")   
           if(MD)
           cat(" Model: 1-compartment, iv infusion, multiple-dose model\n\n")
           else    
           cat(" Model: 1-compartment, iv infusion, single-dose model \n")
           cat("      Subject #:", Subject,"                       \n")
           cat("           Dose:", Dose,"                          \n") 
           if(MD) {
           cat("Dosing Interval:", Tau,"                           \n")
           cat("      # of Dose:", nDose,"                         \n")}
           cat("  Infusion time:", Tinf,"                          \n")
           cat("     Error Type:", type,"\n\n")     
           if(MMe){
             sim<-matrix(c(Vm,Km,Vd,par1,par2,par3),3,2)
             dimnames(sim)<-list(c("Vm","Km","Vd"),c("Simulated Values","Input Values"))}
           else{
             sim<-matrix(c(kel,Vd,par1,par2),2,2)
             dimnames(sim)<-list(c("kel","Vd"),c("Simulated Values","Input Values"))}
           show(sim)                                                               
           cat("*******************************************************\n\n")      

       
           if (!MMe){
             kel<-par1
             Vd <-par2
             ### PKindex[[i]]<-sinfu1.out(PKtime,kel,Vd,defun,par1,par2,Dose,i,type,MD) 
             time<-PKtime$time
             parms<-c(kel=kel,Vd=Vd)
             if(MD){
                dosing.time<-seq(0,Tau*nDose,Tau)
                yini<-c(dCpdt=0)
                events <- data.frame(var="dCpdt",time=dosing.time,value=Dose/Tinf/Vd,method="add")
                C1.lsoda<-data.frame(lsode(yini,c(0,time),defun,parms,rtol=1e-08,atol=1e-08,
                                  events=list(data=events)))
             }
             else{
               C1.lsoda<-data.frame(lsoda(0,c(0,time),defun,parms,rtol=1e-08,atol=1e-08))}
             
             good<-ifelse(C1.lsoda[2:(length(time)+1),2]<=0,
                          0,
                          C1.lsoda[2:(length(time)+1),2])
             PKindex.this.subj<-data.frame(i,
                                 C1.lsoda[2:(length(time)+1),1],
                                 good)
             colnames(PKindex.this.subj)<-list("Subject","time","conc")
             show(PKindex.this.subj)
             x<-C1.lsoda[2:(length(time)+1),1]
             y<-good
             plotting.sim(i,x,y,MD);PKindex[[i]]<-PKindex.this.subj   ### dump this subj's data to subj i.
           } 
           else{
             Vm<-par1
             Km<-par2
             Vd<-par3       
             ### PKindex[[i]]<-sinfu.mm.out(PKtime,Vm,Km,Vd,defun,par1,par2,par3,Dose,i,type,MD)
             time<-PKtime$time
             parms<-c(Vm=Vm,Km=Km,Vd=Vd)
             if(MD){
                dosing.time<-seq(0,Tau*nDose,Tau)
                yini<-c(dCpdt=0)
                events <- data.frame(var="dCpdt",time=dosing.time,value=Dose/Tinf/Vd,method="add")
                C1.lsoda<-data.frame(lsode(yini,c(0,time),defun,parms,rtol=1e-08,atol=1e-08,
                                  events=list(data=events)))
             }
             else{
             C1.lsoda<-data.frame(lsoda(0,c(0,time),defun,parms,atol=1e-10))}
              
             good<-ifelse(C1.lsoda[2:(length(time)+1),2]<=0,
                          0,
                          C1.lsoda[2:(length(time)+1),2])
             PKindex.this.subj<-data.frame(i,
                                 C1.lsoda[2:(length(time)+1),1],
                                 good)
             colnames(PKindex.this.subj)<-list("Subject","time","conc")
             show(PKindex.this.subj)
             x<-C1.lsoda[2:(length(time)+1),1]
             y<-good
             plotting.sim(i,x,y,MD);PKindex[[i]]<-PKindex.this.subj
           }
            ###         
            ### here revert between pdf() and graphic device                          ### added by YJ
            ### 
                      if(pdf_activate){
                         dev.copy()                          ## copy to pdf file 2nd plots to end
                         dev.set(which=x11c)                 ## back from graphic device now to continue...
                                      }
                      else{
                         x11c<-dev.cur()                     ## the current graphics device
                         pdf(sim.plots_to_pdf,paper="a4")
                         pdf_activate=TRUE                   ## set pdf_activate=TRUE from now on
                         dev.set(which=x11c)                 ## go to graphics device...
                         dev.copy()                          ## copy the first plot here
                         dev.set(which=x11c)                 ## back from graphics device
                          }
            ###
            ###  end plotting here...           
         }       
         PKindex<- as.data.frame(do.call("rbind",PKindex))
         rownames(PKindex) <- seq(nrow(PKindex)) 
         ### savefile(PKindex)
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
           
           PKindex<-vector(Subject,mode="list")
           for(i in 1:Subject)  {
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
             cat("\n\n     << Subject:- #",i,">>" )
             cat("\n\n")                                                               
             cat("*******************************************************\n")        
             cat(" Summary Table                                          \n")   
             if(MD)
             cat(" Model: 1-compartment, iv infusion, multiple-dose model\n\n")
             else    
             cat(" Model: 1-compartment, iv infusion, single-dose model \n")
             cat("      Subject #:", Subject,"                       \n")
             cat("           Dose:", Dose,"                          \n") 
             if(MD) {
             cat("Dosing Interval:", Tau,"                           \n")
             cat("      # of Dose:", nDose,"                         \n")}  
             cat("  Infusion time:", Tinf,"                          \n")     
             cat("     Error Type:", type,"\n\n")     
             if(MMe){
             sim<-matrix(c(Vm,Km,Vd,par1,par2,par3),3,2)
             dimnames(sim)<-list(c("Vm","Km","Vd"),c("Simulated Values","Input Values"))
             }
             else{
             sim<-matrix(c(kel,Vd,par1,par2),2,2)
             dimnames(sim)<-list(c("kel","Vd"),c("Simulated Values","Input Values"))
             }
             show(sim)                                                               
             cat("*******************************************************\n\n")                            

             ### PKindex[[i]]<-sinfu1.out(PKtime,kel,Vd,defun,par1,par2,Dose,i,type,MD)
             time<-PKtime$time
             parms<-c(kel=kel,Vd=Vd)
             if(MD){
                dosing.time<-seq(0,Tau*nDose,Tau)
                yini<-c(dCpdt=0)
                events <- data.frame(var="dCpdt",time=dosing.time,value=Dose/Tinf/Vd,method="add")
                C1.lsoda<-data.frame(lsode(yini,c(0,time),defun,parms,rtol=1e-08,atol=1e-08,
                                  events=list(data=events)))
             }
             else{
               C1.lsoda<-data.frame(lsoda(0,c(0,time),defun,parms,rtol=1e-08,atol=1e-08))}
             
             good<-ifelse(C1.lsoda[2:(length(time)+1),2]<=0,
                          0,
                          C1.lsoda[2:(length(time)+1),2])
             PKindex.this.subj<-data.frame(i,
                                 C1.lsoda[2:(length(time)+1),1],
                                 good)
             colnames(PKindex.this.subj)<-list("Subject","time","conc")
             show(PKindex.this.subj)
             x<-C1.lsoda[2:(length(time)+1),1]
             y<-good
             plotting.sim(i,x,y,MD);PKindex[[i]]<-PKindex.this.subj
             ###         
             ### here revert between pdf() and graphic device                          ### added by YJ
             ### 
                       if(pdf_activate){
                          dev.copy()                          ## copy to pdf file 2nd plots to end
                          dev.set(which=x11c)                 ## back from graphic device now to continue...
                                       }
                       else{
                          x11c<-dev.cur()                     ## the current graphics device
                          pdf(sim.plots_to_pdf,paper="a4")
                          pdf_activate=TRUE                   ## set pdf_activate=TRUE from now on
                          dev.set(which=x11c)                 ## go to graphics device...
                          dev.copy()                          ## copy the first plot here
                          dev.set(which=x11c)                 ## back from graphics device
                           }
             ###
             ###  end plotting here...                
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
           
           PKindex<-vector(Subject,mode="list")
           for( i in 1:Subject)  {
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
             cat("\n\n     << Subject:- #",i,">>" )
             cat("\n\n")                                                               
             cat("*******************************************************\n")        
             cat(" Summary Table                                          \n")   
             if(MD)
             cat(" Model: 1-compartment, iv infusion, multiple-dose model\n\n")
             else    
             cat(" Model: 1-compartment, iv infusion, single-dose model \n")
             cat("      Subject #:", Subject,"                       \n")
             cat("           Dose:", Dose,"                          \n") 
             if(MD) {
             cat("Dosing Interval:", Tau,"                           \n")
             cat("      # of Dose:", nDose,"                         \n")}  
             cat("  Infusion time:", Tinf,"                          \n")     
             cat("     Error Type:", type,"\n\n")     
             if(MMe){
             sim<-matrix(c(Vm,Km,Vd,par1,par2,par3),3,2)
             dimnames(sim)<-list(c("Vm","Km","Vd"),c("Simulated Values","Input Values"))}
             else{
             sim<-matrix(c(kel,Vd,par1,par2),2,2)
             dimnames(sim)<-list(c("kel","Vd"),c("Simulated Values","Input Values"))}
             show(sim)                                                               
             cat("*******************************************************\n\n")

             ### PKindex[[i]]<-sinfu.mm.out(PKtime,Vm,Km,Vd,defun,par1,par2,par3,Dose,i,type,MD)
             time<-PKtime$time
             parms<-c(Vm=Vm,Km=Km,Vd=Vd)
             if(MD){
                dosing.time<-seq(0,Tau*nDose,Tau)
                yini<-c(dCpdt=0)
                events <- data.frame(var="dCpdt",time=dosing.time,value=Dose/Tinf/Vd,method="add")
                C1.lsoda<-data.frame(lsode(yini,c(0,time),defun,parms,rtol=1e-08,atol=1e-08,
                                  events=list(data=events)))
             }
             else{
             C1.lsoda<-data.frame(lsoda(0,c(0,time),defun,parms,atol=1e-10))}
              
             good<-ifelse(C1.lsoda[2:(length(time)+1),2]<=0,
                          0,
                          C1.lsoda[2:(length(time)+1),2])
             PKindex.this.subj<-data.frame(i,
                                 C1.lsoda[2:(length(time)+1),1],
                                 good)
             colnames(PKindex.this.subj)<-list("Subject","time","conc")
             show(PKindex.this.subj)
             x<-C1.lsoda[2:(length(time)+1),1]
             y<-good
             plotting.sim(i,x,y,MD);PKindex[[i]]<-PKindex.this.subj
             ###         
             ### here revert between pdf() and graphic device                          ### added by YJ
             ### 
                       if(pdf_activate){
                          dev.copy()                          ## copy to pdf file 2nd plots to end
                          dev.set(which=x11c)                 ## back from graphic device now to continue...
                                       }
                       else{
                          x11c<-dev.cur()                     ## the current graphics device
                          pdf(sim.plots_to_pdf,paper="a4")
                          pdf_activate=TRUE                   ## set pdf_activate=TRUE from now on
                          dev.set(which=x11c)                 ## go to graphics device...
                          dev.copy()                          ## copy the first plot here
                          dev.set(which=x11c)                 ## back from graphics device
                           }
             ###
             ###  end plotting here...             
           } 
         }    
         PKindex<- as.data.frame(do.call("rbind",PKindex))
         rownames(PKindex) <- seq(nrow(PKindex)) 
         ### savefile(PKindex)
   } 
  } 
  else if (pick ==2){   ### start mcsim from here
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
       }
       
         sink(zz,split=TRUE)
                
         cat("\n\n")
         cat("****************************************************\n")
         cat("Summary Table - Monte-Carlo simulation runs         \n")
         if(MD)
         cat("Model: 1-compartment, iv infusion, multiple-dose model \n\n") 
         else
         cat("Model: 1-compartment, iv infusion, single-dose model   \n\n") 
         cat("      Subject #:", Subject,"                      \n")
         if(MD) {
         cat("Dosing Interval:", Tau,"                          \n")
         cat("      # of Dose:", nDose,"                        \n")}
         cat("  Infusion time:", Tinf,"                          \n")
         cat("   Simulation #:", re,"                            \n\n")
         cat("     Error Type:", type,"                          \n")
         if(MMe){
             sim<-matrix(c(Vm,Km,Vd,factor1,factor2,factor3),3,2)
             dimnames(sim)<-list(c("Vm","Km","Vd"),c("Input Values","Error factor"))
             }
             else{
             sim<-matrix(c(kel,Vd,factor1,factor2),2,2)
             dimnames(sim)<-list(c("kel","Vd"),c("Input Values","Error factor"))
             }
         show(sim)   
         cat("****************************************************")          
         
     PKindex<-vector(Subject,mode="list")
     for( i in 1:Subject)  {
       cat("\n\n     << Subject:- #",i,">>\n" )
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
              if(MD){
                 dosing.time<-seq(0,Tau*nDose,Tau)
                 yini<-c(dCpdt=0)
                 events <- data.frame(var="dCpdt",time=dosing.time,value=Dose/Tinf/Vd,method="add")
                 XX<-data.frame(lsode(yini,c(0,time1),defun,parms,rtol=1e-08,atol=1e-08,
                                      events=list(data=events)))
              }
              else{ 
              XX<-data.frame(lsoda(0,c(0,time1),defun,parms,rtol=1e-08,atol=1e-08))}
              C1.lsoda[[j]]<-data.frame(XX[2:(length(time1)+1),1],XX[2:(length(time1)+1),2])
              colnames(C1.lsoda[[j]])<-list("time","concentration")        
           } 
         }        
            PKindex[[i]]<-montecarlo(C1.lsoda,time1,i,re) 
            ###         
            ### here revert between pdf() and graphic device                          ### added by YJ
            ### 
                      if(pdf_activate){
                         dev.copy()                          ## copy to pdf file 2nd plots to end
                         dev.set(which=x11c)                 ## back from graphic device now to continue...
                                      }
                      else{
                         x11c<-dev.cur()                     ## the current graphics device
                         pdf(sim.plots_to_pdf,paper="a4")
                         pdf_activate=TRUE                   ## set pdf_activate=TRUE from now on
                         dev.set(which=x11c)                 ## go to graphics device...
                         dev.copy()                          ## copy the first plot here
                         dev.set(which=x11c)                 ## back from graphics device
                          }
            ###
            ###  end plotting here...         
     }  
     PKindex<-as.data.frame(do.call("rbind",PKindex))
     rownames(PKindex)<-seq(nrow(PKindex)) 
     ### savefile(PKindex) 
  }
   sink()
   close(zz)
   cat(paste("\n\n Two outputs,",sim.outputs_to_txt,"&",sim.plots_to_pdf,",\n have been generated at",getwd(),"\n\n"))
   readline(" Press any key to continue...")   
   savefile(PKindex)
   dev.off()
   graphics.off()
   if(MD) sone.iv.route.MD()
   else   sone.iv.route.SD()
} 
