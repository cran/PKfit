### Simulation
### One compartment PK model iv bolus single dose
### optional Michaelis-Menten Elimination

sbolus1 <- function(Subject=NULL,  # N Subj's 
                    PKtime=NULL,   # times for sampling
                    Dose=NULL,     # single dose
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
if(MD){
   if (!MMe) {
      if(file.exists("sbolus1_md.csv")){
           par<-read.csv(file="sbolus1_md.csv",row.names=NULL,header=TRUE)
           par<-edit(par)}
      else{   
          par<-data.frame(Parameter=c("Total_subject#","Dose","Tau","#Dose","kel","Vd"),Initial=c(24,300,12,10,0.21,11.7))
          par<-edit(par)                                                       
          repeat{
              if (par[1,2]<= 0 || par[2,2]<=0 || par[3,2]<=0 || par[4,2]<=0|| par[5,2]<=0 || par[6,2]<=0){
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
      cat("\n")       
      show(par)
   
      cat("\n")
      Subject<-par[1,2]
      Dose   <-par[2,2]
      Tau    <-par[3,2]
      nDose  <-par[4,2]
      par1   <-par[5,2]
      par2   <-par[6,2]
 
      write.csv(par,file="sbolus1_md.csv",row.names=FALSE,col.names=TRUE)
  }
   else{
      if(file.exists("sbolus1_mm_md.csv")){
           par<-read.csv(file="sbolus1_mm_md.csv",row.names=NULL,header=TRUE)
           par<-edit(par)}
      else{     
      par<-data.frame(Parameter=c("Total_subject#","Dose","Tau","#Dose","Vm","Km","Vd"),Initial=c(24,300,12,10,1.84,4.51,11.2))
      par<-edit(par)                                                       
      repeat{
          if (par[1,2]<= 0 || par[2,2]<=0 || par[3,2]<=0 || par[4,2]<=0 || par[5,2]<=0|| par[6,2]<=0 || par[7,2]<=0){
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
    }                                    #### move one '}' here from the end of write.csv(...)
      cat("\n")       
      show(par)
      write.csv(par,file="sbolus1_mm_md.csv",row.names=FALSE,col.names=TRUE)
      cat("\n")
      Subject<-par[1,2]
      Dose   <-par[2,2]
      Tau    <-par[3,2]
      nDose  <-par[4,2]
      par1<-par[5,2]
      par2<-par[6,2]
      par3<-par[7,2]
   }
      
### no PKtime require for multiple-dose simulation
### PKtime<-c(0:(Tau*nDose-1))
### show(Tau);cat("\n");show(nDose);readline()    ### for debugging; move one "}" to upper... -YJ
   PKtime<-data.frame(time=seq(0.1,Tau*nDose-1,0.1))
   
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
}
else{
   if (!MMe) {
      if(file.exists("sbolus1.csv")){
           par<-read.csv(file="sbolus1.csv",row.names=NULL,header=TRUE)
           par<-edit(par)}
      else{   
          par<-data.frame(Parameter=c("Total_subject#","Dose","kel","Vd"),Initial=c(24,300,0.21,11.7))
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
      write.csv(par,file="sbolus1.csv",row.names=FALSE,col.names=TRUE)
  }
   else{
      if(file.exists("sbolus1_mm.csv")){
           par<-read.csv(file="sbolus1_mm.csv",row.names=NULL,header=TRUE)
           par<-edit(par)}
      else{     
      par<-data.frame(Parameter=c("Total_subject#","Dose","Vm","Km","Vd"),Initial=c(24,300,2.17,4.84,11.7))
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
      write.csv(par,file="sbolus1_mm.csv",row.names=FALSE,col.names=TRUE)
}
   readline("\n Next enter or load time points. Press Enter to continue.\n\n")
   if(file.exists("sim_times.csv")){
      PKtime<-read.csv(file="sim_times.csv",row.names=NULL,header=TRUE)
      PKtime<-edit(PKtime)}
   else{
      ## PKtime<-data.frame(time=c(0))
      PKtime<-data.frame(time=c(0,.1,.2,.3,.4,.6,.8,1,2,4,6,8,12,14,16,18,24,48,72))
      PKtime<-edit(PKtime)}
   
   write.csv(PKtime,file="sim_times.csv",row.names=FALSE,col.names=TRUE)
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
}
   cat("\n")
   file.menu <- c("Simulation with Error",
                  "Monte Carlo Simulation")
   pick <- menu(file.menu, title = "<< Simulation Types >>")

###
dev.new()
par(mfrow=c(2,1),las=1)
pdf_activate=FALSE  ### set pdf device activate? as FALSE at beginning
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
     
     if (pick ==1){    ### we will use this module for tdm too with mutiple-dose scheme.
         Subject<-1          ### does not make any sense to simulate more than 1 subject with this option.
         PKindex<-vector(Subject,mode="list")
         for( i in 1:Subject)  {
             cat("\n\n     << Subject:- #",i,">>" )
             cat("\n\n")
             cat("****************************************************\n")
             cat(" Summary Table                                     \n\n")
             if(MD)
             cat(" Model: 1-compartment, iv bolus, multiple-dose model\n\n")
             else
             cat(" Model: 1-compartment, iv bolus, single-dose model \n\n")             
             cat("      Subject #:", Subject,"                       \n")
             cat("           Dose:", Dose,"                          \n")
             if(MD) {
             cat("Dosing Interval:", Tau,"                           \n")
             cat("      # of Dose:", nDose,"                         \n")}
             cat("     Error Type:", type,"                        \n\n")
             if(MMe){
             sim<-matrix(c(Vm,Km,Vd,par1,par2,par3),3,2)
             dimnames(sim)<-list(c("Vm","Km","Vd"),c("Simulated Values","Input Values"))
             }
             else{
             sim<-matrix(c(kel,Vd,par1,par2),2,2)
             dimnames(sim)<-list(c("kel","Vd"),c("Simulated Values","Input Values"))
             }
             show(sim)
             cat("****************************************************\n\n")
       
           if (!MMe ){
             kel<-par1
             Vd <-par2   
             ### PKindex[[i]]<-sbolus1.out(PKtime,kel,Vd,defun,par1,par2,Dose,i,type,MD)
             time<-PKtime$time
             parms<-c(kel=kel,Vd=Vd)
             if(MD){
             dosing.time<-seq(Tau,Tau*nDose,Tau)
             yini<-c(dCpdt=Dose/Vd)
             events <- data.frame(var="dCpdt",time=dosing.time,value=Dose/Vd,method="add")
             C1.lsoda<-data.frame(lsode(yini,c(0,time),defun,parms,rtol=1e-08,atol=1e-08,
                                  events=list(data=events)))
             }
             else{
             C1.lsoda<-data.frame(lsoda(Dose/Vd,c(0,time),defun,parms,rtol=1e-08,atol=1e-08))}
             
             good <- ifelse(C1.lsoda[2:(length(time)+1),2]<=0,
                            0,
                            C1.lsoda[2:(length(time)+1),2])
             PKindex.this.subj <- data.frame(i,
                                   C1.lsoda[2:(length(time)+1),1],
                                   good)
             colnames(PKindex.this.subj)<-list("Subject","time","conc")
             show(PKindex.this.subj)
             x <- C1.lsoda[2:(length(time)+1),1]
             y <- good
             plotting.sim(i,x,y,MD);PKindex[[i]]<-PKindex.this.subj   ### dump this subj's data to subj i.
           } 
           else{
             Vm<-par1
             Km<-par2
             Vd<-par3       
             ### PKindex[[i]]<-sbolus.mm.out(PKtime,Vm,Km,Vd,defun,par1,par2,par3,Dose,i,type,MD)
             time<-PKtime$time
             parms<-c(Vm=Vm,Km=Km,Vd=Vd)
             if(MD){
             dosing.time<-seq(Tau,Tau*nDose,Tau)
             yini<-c(dCpdt=Dose/Vd)
             events <- data.frame(var="dCpdt",time=dosing.time,value=Dose/Vd,method="add")
             C1.lsoda<-data.frame(lsode(yini,c(0,time),defun,parms,rtol=1e-08,atol=1e-08,
                                  events=list(data=events)))             
             }
             else{  
             C1.lsoda<-data.frame(lsoda(Dose/Vd,c(0,time),defun,parms,atol=1e-08,rtol=1e-08))}
             
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
           
           PKindex<-vector(Subject,mode="list");cat("\n")
           pick<- pick-1  # reset pick = 1, 2, 3 or 4 ('0' is zero error & excluded)
           
           for( i in 1:Subject)  {

             ka  <- err_types(pick,par1,factor1)  # reset pick = 1, 2, 3 or 4 ('0' is zero error & excluded)
             kel <- err_types(pick,par2,factor2)
             
             cat("\n\n     << Subject:- #",i,">>" )
             cat("\n\n")
             cat("****************************************************\n")
             cat(" Summary Table                                      \n")
             if(MD)
             cat(" Model: 1-compartment, iv bolus, multiple-dose model \n")
             else
             cat(" Model: 1-compartment, iv bolus, single-dose model \n")             
             cat("      Subject #:", Subject,"                      \n")
             cat("           Dose:", Dose,"                         \n")
             if(MD) {
             cat("Dosing Interval:", Tau,"                          \n")
             cat("      # of Dose:", nDose,"                        \n")}
             cat("     Error Type:", type,"\n\n")
             sim<-matrix(c(kel,Vd,par1,par2),2,2)
             dimnames(sim)<-list(c("kel","Vd"),c("Simulated Values","Input Values"))
             show(sim)
             cat("****************************************************\n\n")             
            ### PKindex[[i]]<-sbolus1.out(PKtime,kel,Vd,defun,par1,par2,Dose,i,type,MD)
            
             time<-PKtime$time
             parms<-c(kel=kel,Vd=Vd)
             if(MD){
             dosing.time<-seq(Tau,Tau*nDose,Tau)
             yini<-c(dCpdt=Dose/Vd)
             events <- data.frame(var="dCpdt",time=dosing.time,value=Dose/Vd,method="add")
             C1.lsoda<-data.frame(lsode(yini,c(0,time),defun,parms,rtol=1e-08,atol=1e-08,
                                  events=list(data=events)))
             }
             else{
             C1.lsoda<-data.frame(lsoda(Dose/Vd,c(0,time),defun,parms,rtol=1e-08,atol=1e-08))}
             
             good <- ifelse(C1.lsoda[2:(length(time)+1),2]<=0,
                            0,
                            C1.lsoda[2:(length(time)+1),2])
             PKindex.this.subj <- data.frame(i,
                                   C1.lsoda[2:(length(time)+1),1],
                                   good)
             colnames(PKindex.this.subj)<-list("Subject","time","conc")
             show(PKindex.this.subj)
             x <- C1.lsoda[2:(length(time)+1),1]
             y <- good
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
           
           PKindex<-vector(Subject,mode="list");cat("\n")
           pick<- pick-1  # reset pick = 1, 2, 3 or 4 ('0' is zero error & excluded)
           
           for( i in 1:Subject)  {

             Vm <- err_types(pick,par1,factor1)
             Km <- err_types(pick,par2,factor2)
             Vd <- err_types(pick,par3,factor3)
             
             cat("\n\n     << Subject:- #",i,">>" )
             cat("\n\n")
             cat("****************************************************\n")
             cat(" Summary Table                                      \n")
             if(MD)
             cat(" Model: 1-compartment, iv bolus, multiple-dose model\n\n")
             else
             cat(" Model: 1-compartment, iv bolus, single-dose model \n\n")             
             cat("      Subject #:", Subject,"                       \n")
             cat("           Dose:", Dose,"                          \n")
             if(MD) {
             cat("Dosing Interval:", Tau,"                           \n")
             cat("      # of Dose:", nDose,"                         \n")}
             cat("     Error Type:", type,"                        \n\n")
             if(MMe){
             sim<-matrix(c(Vm,Km,Vd,par1,par2,par3),3,2)
             dimnames(sim)<-list(c("Vm","Km","Vd"),c("Simulated Values","Input Values"))
             }
             else{
             sim<-matrix(c(kel,Vd,par1,par2),2,2)
             dimnames(sim)<-list(c("kel","Vd"),c("Simulated Values","Input Values"))
             }
             show(sim)
             cat("****************************************************\n\n")             
             ### PKindex[[i]]<-sbolus.mm.out(PKtime,Vm,Km,Vd,defun,par1,par2,par3,Dose,i,type,MD)
             time<-PKtime$time
             parms<-c(Vm=Vm,Km=Km,Vd=Vd)
             if(MD){
             dosing.time<-seq(Tau,Tau*nDose,Tau)
             yini<-c(dCpdt=Dose/Vd)
             events <- data.frame(var="dCpdt",time=dosing.time,value=Dose/Vd,method="add")
             C1.lsoda<-data.frame(lsode(yini,c(0,time),defun,parms,rtol=1e-08,atol=1e-08,
                                  events=list(data=events)))             
             }
             else{  
             C1.lsoda<-data.frame(lsoda(Dose/Vd,c(0,time),defun,parms,atol=1e-08,rtol=1e-08))}
             
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
         cat("Model: 1-Compartment, IV-Bolus, Multiple-Dose Model \n\n") 
         else
         cat("Model: 1-Compartment, IV-Bolus, Single-Dose Model   \n\n") 
         cat("      Subject #:", Subject,"                      \n")
         if(MD) {
         cat("Dosing Interval:", Tau,"                          \n")
         cat("      # of Dose:", nDose,"                        \n")}
         cat("   Simulation #:", re,"                            \n\n")
         cat("     Error Type:", type,"\n")
         if(MMe){
           sim<-matrix(c(par1,par2,par3,factor1,factor2,factor3),3,2)
           dimnames(sim)<-list(c("Vm","Km","Vd"),c("Input Values","Error factor"))}
         else{
           sim<-matrix(c(par1,par2,factor1,factor2),2,2)
           dimnames(sim)<-list(c("kel","Vd"),c("Input Values","Error factor"))}
         show(sim)   
         cat("****************************************************") 
               
       PKindex<-vector(Subject,mode="list")
       for( i in 1:Subject)  {
       cat("\n\n     << Subject:- #",i,">>\n" )
       C1.lsoda<-list()

         for (j in 1:re){
            if (! MMe){
            
             ka  <- err_types(pick,par1,factor1)
             kel <- err_types(pick,par2,factor2)
             
              time1<-PKtime$time
              parms<-c(kel=kel,Vd=Vd)
              if(MD){
                 dosing.time<-seq(Tau,Tau*nDose,Tau)
                 yini<-c(dCpdt=Dose/Vd)
                 events <- data.frame(var="dCpdt",time=dosing.time,value=Dose/Vd,method="add")
                 XX<-data.frame(lsode(yini,c(0,time1),defun,parms,rtol=1e-08,atol=1e-08,
                                      events=list(data=events)))
               }
              else{
              XX<-data.frame(lsoda(Dose/Vd,c(0,time1),defun,parms,rtol=1e-08,atol=1e-08))}
              C1.lsoda[[j]]<-data.frame(XX[2:(length(time1)+1),1],XX[2:(length(time1)+1),2])
              colnames(C1.lsoda[[j]])<-list("time","concentration") 
           }
           else{

             Vm <- err_types(pick,par1,factor1)
             Km <- err_types(pick,par2,factor2)
             Vd <- err_types(pick,par3,factor3)
             
              time1<-PKtime$time
              parms<-c(Vm=Vm,Km=Km,Vd=Vd)
              if(MD){
                 dosing.time<-seq(Tau,Tau*nDose,Tau)
                 yini<-c(dCpdt=Dose/Vd)
                 events <- data.frame(var="dCpdt",time=dosing.time,value=Dose/Vd,method="add")
                 XX<-data.frame(lsode(yini,c(0,time1),defun,parms,rtol=1e-08,atol=1e-08,
                                      events=list(data=events)))
               }
              else{
              XX<-data.frame(lsoda(Dose/Vd,c(0,time1),defun,parms,rtol=1e-6,atol=1e-10))}
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
     ###savefile(PKindex) 
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
