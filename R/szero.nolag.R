### Simulation
### One compartment PK model zero-ordered abs., single dose
### optional Michaelis-Menten Elimination
### 
szero.nolag<- function(Subject=NULL,  ## N Subj's 
                       PKtime=NULL,   ## times for sampling
                       Dose=NULL,     ## single dose
                       Tabs=NULL,    
                       Vd=NULL,
                       kel=NULL,      ## If not MM elimination
                       MMe=FALSE,     ## michaelis-menten elimination?
                       Vm=NULL,Km=NULL,
                       nDose=NULL,    ## for multiple-dose
                       Tau=NULL,
                       MD=FALSE)
{
   options(warn=-1)
   sim.outputs_to_txt<-sim.outputs_to_txt
   sim.plots_to_pdf<-sim.plots_to_pdf
      
cat("\n First enter all paramaters used for simulation profile.\n");readline(" Press Enter to continue...\n\n")
if(MD){
if (!MMe) { 
  if(file.exists("szero_nolag_md.csv")){
      par<-read.csv(file="szero_nolag_md.csv",row.names=NULL,header=TRUE)
      par<-edit(par)}
  else{ 
      par<-data.frame(Parameter=c("Total_subject#","Dose","Tau","#Dose","Tabs","kel","Vd"),
                        Initial=c(24,300,12,10,1.0,0.21,11.7))
      par<-edit(par)                                                       
      repeat{
          if (par[1,2]<=0||par[2,2]<=0||par[3,2]<=0||par[4,2]<=0||par[5,2]<=0||par[6,2]<=0||par[7,2]<=0){
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
      write.csv(par,file="szero_nolag_md.csv",row.names=FALSE,col.names=TRUE)
      cat("\n")       
      show(par)
      
      cat("\n")
      Subject<-par[1,2]
      Dose   <-par[2,2]
      Tau    <-par[3,2]
      nDose  <-par[4,2]
      par1<-par[5,2]
      par2<-par[6,2]
      par3<-par[7,2]
  }
else{
  if(file.exists("szero_nolag_mm_md.csv")){
    par<-read.csv(file="szero_nolag_mm_md.csv",row.names=NULL,header=TRUE)
    par<-edit(par)}
  else{ 
    par<-data.frame(Parameter=c("Total_subject#","Dose","Tau","#Dose","Tabs","Vm","Km","Vd"),
                      Initial=c(24,300,12,10,0.5,2.31,4.74,11.7))
    par<-edit(par)                                                       
    repeat{
        if (par[1,2]<=0||par[2,2]<=0||par[3,2]<=0||par[4,2]<=0||par[5,2]<=0||par[6,2]<=0||par[7,2]<=0||par[8,2]<=0){
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
    write.csv(par,file="szero_nolag_mm_md.csv",row.names=FALSE,col.names=TRUE)
    cat("\n")       
    show(par)
    
    cat("\n")
    Subject<-par[1,2]
    Dose   <-par[2,2]
    Tau    <-par[3,2]
    nDose  <-par[4,2]
    par1<-par[5,2]
    par2<-par[6,2]
    par3<-par[7,2]
    par4<-par[8,2]
} 

### no PKtime require for multiple-dose simulation
### PKtime<-c(0:(Tau*nDose-1))
### show(Tau);cat("\n");show(nDose);readline()    ### for debugging; move one "}" to upper... -YJ
PKtime<-data.frame(time=seq(0,Tau*nDose-1,0.1))
   
  if (!MMe){
    defun<- function(time, y, parms) { 
     if(time<Tau){                  ### for this 1st dose. --YJ
     if(time<=parms["Tabs"]) 
        dCpdt <- (Dose/parms["Tabs"])/parms["Vd"] - parms["kel"] * y[1]
     else
        dCpdt <- - parms["kel"] * y[1]
     }
     else{                         ### for 2nd, 3rd,.... dose. -YJ
     if(time%%Tau<=1 && time%%Tau<=parms["Tabs"]) 
        dCpdt <- (Dose/parms["Tabs"])/parms["Vd"] - parms["kel"] * y[1]
     else dCpdt <- - parms["kel"] * y[1]
     }
     list(dCpdt) 
    } 
  }
  else{
    defun<- function(time, y, parms) { 
     if(time<Tau){
     if(time<=parms["Tabs"]) 
       dCpdt <- (Dose/parms["Tabs"])/parms["Vd"]-parms["Vm"]*y[1]/(parms["Km"]+y[1])            
     else
       dCpdt <- -parms["Vm"]*y[1]/(parms["Km"]+y[1])
     }
     else{
     if(time%%Tau<=1 && time%%Tau<=parms["Tabs"])
       dCpdt <- (Dose/parms["Tabs"])/parms["Vd"]-parms["Vm"]*y[1]/(parms["Km"]+y[1])            
     else
       dCpdt <- -parms["Vm"]*y[1]/(parms["Km"]+y[1]) 
     }          
     list(dCpdt) 
    }
  }
}
else{     ### single-dose sim from here.
if (!MMe) { 
  if(file.exists("szero_nolag.csv")){
      par<-read.csv(file="szero_nolag.csv",row.names=NULL,header=TRUE)
      par<-edit(par)}
  else{ 
      par<-data.frame(Parameter=c("Total_subject#","Dose","Tabs","kel","Vd"),Initial=c(24,300,1.0,0.21,11.7))
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
      write.csv(par,file="szero_nolag.csv",row.names=FALSE,col.names=TRUE)
      cat("\n")       
      show(par)
      
      cat("\n")
      Subject<-par[1,2]
      Dose<-par[2,2]
      par1<-par[3,2]
      par2<-par[4,2]
      par3<-par[5,2]
  }
else{
  if(file.exists("szero_nolag_mm.csv")){
    par<-read.csv(file="szero_nolag_mm.csv",row.names=NULL,header=TRUE)
    par<-edit(par)}
  else{ 
    par<-data.frame(Parameter=c("Total_subject#","Dose","Tabs","Vm","Km","Vd"),Initial=c(24,300,0.5,2.31,4.74,11.7))
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
    write.csv(par,file="szero_nolag_mm.csv",row.names=FALSE,col.names=TRUE)
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
       dCpdt <- (Dose/parms["Tabs"])/parms["Vd"]-parms["Vm"]*y[1]/(parms["Km"]+y[1])            
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
         Subject<-1       ### no error does not req. more than one subj.
         PKindex<-vector(Subject,mode="list")
         for( i in 1:Subject)  {
             cat("\n\n     << Subject:- #",i,">>" )
             cat("\n\n")
             cat("****************************************************\n")
             cat(" Summary Table                                     \n\n")
             if(MD)
             cat(" Model: 1-compartment, extravascular, multiple-dose, \n")
             else
             cat(" Model: 1-compartment, extravascular, single-dose, \n")
             if(MMe)
             cat("   & zero-ordered Abs. with MM elim., no lag time\n\n")
             else
             cat("   & zero-ordered Abs. with 1st-ordered elim, no lag time\n\n")        
             cat("      Subject #:", Subject,"                       \n")
             cat("           Dose:", Dose,"                          \n")
             if(MD) {
             cat("Dosing Interval:", Tau,"                           \n")
             cat("      # of Dose:", nDose,"                         \n")}
             cat("     Error Type:", type,"                        \n\n")
             if(MMe){
             sim<-matrix(c(Tabs,Vm,Km,Vd,par1,par2,par3,par4),4,2)
             dimnames(sim)<-list(c("Tabs","Vm","Km","Vd"),c("Simulated Values","Input Values"))
             }
             else{
             sim<-matrix(c(Tabs,kel,Vd,par1,par2,par3),3,2)
             dimnames(sim)<-list(c("Tabs","kel","Vd"),c("Simulated Values","Input Values"))
             }
             show(sim)
             cat("****************************************************\n\n")
           
           if (! MMe ){
             Tabs<-par1
             kel <-par2
             Vd  <-par3  
             ### PKindex[[i]]<-szero.out(PKtime,Tabs,kel,Vd,defun,par1,par2,par3,Dose,i,type,MD)     
             time<-PKtime$time
             parms<-c(Tabs=Tabs,kel=kel,Vd=Vd)  
             C1.lsoda<-data.frame(lsoda(0, c(0,time), defun, parms,atol=1e-08,rtol=1e-08))
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
             Tabs<-par1
             Vm  <-par2
             Km  <-par3
             Vd  <-par4
             ### PKindex[[i]]<-szero.mm.out(PKtime,Tabs,Vm,Km,Vd,defun,par1,par2,par3,par4,Dose,i,type,MD)
             time<-PKtime$time
             parms<-c(Tabs=Tabs,Vm=Vm,Km=Km,Vd=Vd) 
             C1.lsoda<-data.frame(lsoda(0,c(0,time), defun, parms,atol=1e-08,rtol=1e-08))
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
         if (! MMe){
           par.err<-data.frame(Parameter=c("Tabs","kel","Vd"),error_factor=c(.15,.15,.15))
           par.err<-edit(par.err)
           repeat{
           if (par.err[1,2]<=0 || par.err[2,2]<=0|| par.err[3,2]<=0){
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
           
              Tabs<- err_types(pick,par1,factor1)
              kel <- err_types(pick,par2,factor2)
              Vd  <- err_types(pick,par3,factor3)
              
             cat("\n\n     << Subject:- #",i,">>" )
             cat("\n\n")
             cat("****************************************************\n")
             cat(" Summary Table                                     \n\n")
             if(MD)
             cat(" Model: 1-compartment, extravascular, multiple-dose, \n")
             else
             cat(" Model: 1-compartment, extravascular, single-dose, \n")
             if(MMe)
             cat("   & zero-ordered Abs. with MM elim., no lag time\n\n")
             else
             cat("   & zero-ordered Abs. with 1st-ordered elim, no lag time\n\n") 
             cat("      Subject #:", Subject,"                       \n")
             cat("           Dose:", Dose,"                          \n")
             if(MD) {
             cat("Dosing Interval:", Tau,"                           \n")
             cat("      # of Dose:", nDose,"                         \n")}
             cat("     Error Type:", type,"                        \n\n")
             if(MMe){
             sim<-matrix(c(Tabs,Vm,Km,Vd,par1,par2,par3,par4),4,2)
             dimnames(sim)<-list(c("Tabs","Vm","Km","Vd"),c("Simulated Values","Input Values"))
             }
             else{
             sim<-matrix(c(Tabs,kel,Vd,par1,par2,par3),3,2)
             dimnames(sim)<-list(c("Tabs","kel","Vd"),c("Simulated Values","Input Values"))
             }
             show(sim)
             cat("****************************************************\n\n")
                           
             ### PKindex[[i]]<-szero.out(PKtime,Tabs,kel,Vd,defun,par1,par2,par3,Dose,i,type,MD)
             
             time<-PKtime$time
             parms<-c(Tabs=Tabs,kel=kel,Vd=Vd)  
             C1.lsoda<-data.frame(lsoda(0, c(0,time), defun, parms,atol=1e-08,rtol=1e-08))
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
           par.err<-data.frame(Parameter=c("Tabs","Vm","Km","Vd"),error_factor=c(.15,.15,.15,.15))
           par.err<-edit(par.err)
           repeat{
           if (par.err[1,2] <=0 || par.err[2,2] <=0|| par.err[3,2] <=0|| par.err[4,2] <=0){
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
           pick<- pick-1  # reset pick = 1, 2, 3 or 4 ('0' is zero error & excluded)
           
           for( i in 1:Subject)  {

              Tabs<- err_types(pick,par1,factor1)
              Vm  <- err_types(pick,par2,factor2)
              Km  <- err_types(pick,par3,factor3)
              Vd  <- err_types(pick,par4,factor4)

             cat("\n\n     << Subject:- #",i,">>" )
             cat("\n\n")
             cat("****************************************************\n")
             cat(" Summary Table                                     \n\n")
             if(MD)
             cat(" Model: 1-compartment, extravascular, multiple-dose, \n")
             else
             cat(" Model: 1-compartment, extravascular, single-dose, \n")
             if(MMe)
             cat("   & zero-ordered Abs. with MM elim., no lag time\n\n")
             else
             cat("   & zero-ordered Abs. with 1st-ordered elim, no lag time\n\n") 
             cat("      Subject #:", Subject,"                       \n")
             cat("           Dose:", Dose,"                          \n")
             if(MD) {
             cat("Dosing Interval:", Tau,"                           \n")
             cat("      # of Dose:", nDose,"                         \n")}
             cat("     Error Type:", type,"                        \n\n")
             if(MMe){
             sim<-matrix(c(Tabs,Vm,Km,Vd,par1,par2,par3,par4),4,2)
             dimnames(sim)<-list(c("Tabs","Vm","Km","Vd"),c("Simulated Values","Input Values"))
             }
             else{
             sim<-matrix(c(Tabs,kel,Vd,par1,par2,par3),3,2)
             dimnames(sim)<-list(c("Tabs","kel","Vd"),c("Simulated Values","Input Values"))
             }
             show(sim)
             cat("****************************************************\n\n")
                          
             ### PKindex[[i]]<-szero.mm.out(PKtime,Tabs,Vm,Km,Vd,defun,par1,par2,par3,par4,Dose,i,type,MD)
             time<-PKtime$time
             parms<-c(Tabs=Tabs,Vm=Vm,Km=Km,Vd=Vd) 
             C1.lsoda<-data.frame(lsoda(0,c(0,time), defun, parms,atol=1e-08,rtol=1e-08))
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
             ### here revert between pdf() and graphic device     ### added by YJ
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
     cat("\n\nHow many iterations to run for each subject?\n")
     re<-scan(nlines=1,quiet=TRUE)
     
       if (!MMe){
           par.err<-data.frame(Parameter=c("Tabs","kel","Vd"),error_factor=c(.15,.15,.15))
           par.err<-edit(par.err)
           repeat{
           if (par.err[1,2]<=0 || par.err[2,2]<=0|| par.err[3,2]<=0){
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
       else{
         par.err<-data.frame(Parameter=c("Tabs","Vm","Km","Vd"),error_factor=c(.15,.15,.15,.15))
         par.err<-edit(par.err)
         repeat{
         if (par.err[1,2] <=0 || par.err[2,2] <=0|| par.err[3,2] <=0|| par.err[4,2] <=0){
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
       }
         sink(zz,split=TRUE)
         
         cat("\n")
         cat("********************************************************\n")
         cat("Summary Table - Monte-Carlo simulation runs             \n")
         if(MD)
         cat(" Model: 1-compartment, extravascular, multiple-dose,    \n")
         else
         cat(" Model: 1-compartment, extravascular, single-dose,      \n")
         if(MMe)
         cat("   & zero-ordered Abs. with MM elim., no lag time\n\n")
         else
         cat("   & zero-ordered Abs. with 1st-ordered elim, no lag time\n\n") 
         cat("      Subject #:", Subject,"                       \n")
         cat("           Dose:", Dose,"                          \n")
         if(MD) {
         cat("Dosing Interval:", Tau,"                           \n")
         cat("      # of Dose:", nDose,"                         \n")}
         cat("     Error Type:", type,"                        \n\n")
         cat("   Simulation #:", re,"                          \n\n")
         if(MMe){
         sim<-matrix(c(par1,par2,par3,par4,factor1,factor2,factor3,factor4),4,2)
         dimnames(sim)<-list(c("Tabs","Vm","Km","Vd"),c("Input Values","Error factor"))}
         else{
         sim<-matrix(c(par1,par2,par3,factor1,factor2,factor3),3,2)
         dimnames(sim)<-list(c("Tabs","kel","Vd"),c("Input Values","Error factor"))}
         show(sim)   
         cat("********************************************************")
       
       PKindex<-vector(Subject,mode="list")
       for( i in 1:Subject)  {
       cat("\n\n             << Subject",i,">>\n\n" ) 
       C1.lsoda<-list()

         for (j in 1:re){
     
           if (!MMe){

              Tabs<- err_types(pick,par1,factor1) # pick = 1,2,3 or 4
              kel <- err_types(pick,par2,factor2)
              Vd  <- err_types(pick,par3,factor3)
                         
              time1<-PKtime$time
              parms<-c(Tabs=Tabs,kel=kel,Vd=Vd) 
              XX<-data.frame(lsoda(0, c(0,time1), defun, parms))   
              C1.lsoda[[j]]<-data.frame(XX[2:(length(time1)+1),1],XX[2:(length(time1)+1),2])
              colnames(C1.lsoda[[j]])<-list("time","concentration")
           }
           else{
           
              Tabs<- err_types(pick,par1,factor1)
              Vm  <- err_types(pick,par2,factor2)
              Km  <- err_types(pick,par3,factor3)
              Vd  <- err_types(pick,par4,factor4)

              time1<-PKtime$time
              parms<-c(Tabs=Tabs,Vm=Vm,Km=Km,Vd=Vd) 
              XX<-data.frame(lsoda(0,c(0,time1), defun, parms)) 
              C1.lsoda[[j]]<-data.frame(XX[2:(length(time1)+1),1],XX[2:(length(time1)+1),2])
              colnames(C1.lsoda[[j]])<-list("time","concentration")     
           } 
         }        
            PKindex[[i]]<-montecarlo(C1.lsoda,time1,i,re) ### changed from 'time' to 'time1' since v1.2.1
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
   if(MD) sone.noniv.route.MD()
     else sone.noniv.route.SD()
} 
