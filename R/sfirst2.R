### Two compartment PK model extravascular single dose first order absorption
sfirst2<- function(Subject=NULL,   # N Subj's 
                   PKtime=NULL,    # times for sampling
                   Dose=NULL,      # single dose
                   ka=NULL,
                   Vd=NULL,
                   kel=NULL, 
                   k12=NULL,
                   k21=NULL,
                   nDose=NULL,    # for multiple-dose
                   Tau=NULL,
                   MD=FALSE)
{
   options(warn=-1)
   sim.outputs_to_txt<-sim.outputs_to_txt
   sim.plots_to_pdf<-sim.plots_to_pdf   
   
cat("\n First enter all paramaters used for simulation profile.\n");readline(" Press Enter to continue...\n\n")
if(MD){
if(file.exists("sfirst2_md.csv")){
   par<-read.csv(file="sfirst2_md.csv",row.names=NULL,header=TRUE)
   par<-edit(par)}
else{  
   par<-data.frame(Parameter=c("Total_subject#","Dose","Tau","#Dose","ka","kel","k12","k21","Vd"),
                     Initial=c(24,300,8,10,0.51,0.36,0.13,0.31,11.7))
   par<-edit(par)                                                       
   repeat{
       if (par[1,2]<=0||par[2,2]<=0||par[3,2]<=0||par[4,2]<=0||par[5,2]<=0||par[6,2]<=0||par[7,2]<=0||
           par[8,2]<=0||par[9,2]<=0){
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
   write.csv(par,file="sfirst2_md.csv",row.names=FALSE,col.names=TRUE) 
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
   par5<-par[9,2]

### no PKtime require for multiple-dose simulation
PKtime<-data.frame(time=seq(0,Tau*nDose-1,0.1))  ### set '0.1' to see more output details.
    
    defun<- function(time, y, parms) { 
       dCp1dt <- -parms["ka"]*y[1]
       dCp2dt <- parms["ka"]*y[1]/parms["Vd"]-parms["kel"]*y[2]+parms["k21"]*y[3]-parms["k12"]*y[2]
       dCp3dt <- parms["k12"]*y[2]-parms["k21"]*y[3]
       list(c(dCp1dt,dCp2dt,dCp3dt)) 
    } 
}
else{
if(file.exists("sfirst2.csv")){
   par<-read.csv(file="sfirst2.csv",row.names=NULL,header=TRUE)
   par<-edit(par)}
else{  
   par<-data.frame(Parameter=c("Total_subject#","Dose","ka","kel","k12","k21","Vd"),Initial=c(24,300,0.51,0.36,0.13,0.31,11.7))
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
 }
   write.csv(par,file="sfirst2.csv",row.names=FALSE,col.names=TRUE) 
   cat("\n")       
   show(par)
   
   cat("\n")
   Subject<-par[1,2]
   Dose<-par[2,2]
   par1<-par[3,2]
   par2<-par[4,2]
   par3<-par[5,2]
   par4<-par[6,2]
   par5<-par[7,2]

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
    
    defun<- function(time, y, parms) { 
       dCp1dt <- -parms["ka"]*y[1]    ### y[1] here is amount of drug, not drug plasma conc. here -YJ
       dCp2dt <- parms["ka"]*y[1]/parms["Vd"]-parms["kel"]*y[2]+parms["k21"]*y[3]-parms["k12"]*y[2]
       dCp3dt <- parms["k12"]*y[2]-parms["k21"]*y[3]
       list(c(dCp1dt,dCp2dt,dCp3dt)) 
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
         for(i in 1:Subject)  {
           cat("\n\n     << Subject:- #",i,">>" )
           cat("\n\n")                                                               
           cat("*******************************************************\n")        
           cat(" Summary Table                                       \n\n")   
           if(MD){
           cat(" Model: 2-compartment, extravascular,                  \n") 
           cat("        multiple-dose, & 1st-ordered without lag time\n\n")}
           else{    
           cat(" Model: 2-compartment, extravascular,                  \n") 
           cat("        single-dose, & 1st-ordered without lag time\n\n")}
           cat("      Subject #:", Subject,"                       \n")
           cat("           Dose:", Dose,"                          \n") 
           if(MD) {
           cat("Dosing Interval:", Tau,"                           \n")
           cat("      # of Dose:", nDose,"                         \n")}
           cat("     Error Type:", type,"\n\n")     
           sim<-matrix(c(ka,kel,k12,k21,Vd,par1,par2,par3,par4,par5),5,2)
           dimnames(sim)<-list(c("ka","kel","k12","k21","Vd"),c("Simulated Values","Input Values"))
           show(sim)                                                               
           cat("*******************************************************\n\n")
           ka <-par1
           kel<-par2
           k12<-par3
           k21<-par4
           Vd <-par5   
           ### PKindex[[i]]<-sfirst2.out(PKtime,ka,kel,k12,k21,Vd,defun,par1,par2,par3,par4,par5,Dose,i,type,MD)
           time<-PKtime$time
           parms<-c(ka=ka,kel=kel,k12=k12,k21=k21,Vd=Vd)
           if(MD){
              dosing.time<-seq(Tau,Tau*nDose,Tau)
              yini<-c(dCp1dt=Dose,dCp2dt=0,dCp3dt=0)
              events <- data.frame(var="dCp1dt",time=dosing.time,value=Dose,method="add")
              C1.lsoda<-data.frame(lsode(yini,c(0,time),defun,parms,rtol=1e-08,atol=1e-08,
                                   events=list(data=events)))
           }
           else{
             parms<-c(ka=ka,kel=kel,k12=k12,k21=k21,Vd=Vd) 
             C1.lsoda<-data.frame(lsoda(c(Dose,0,0),c(0,time),defun,parms,rtol=1e-08,atol=1e-08))}
           good<-ifelse(C1.lsoda[2:(length(time)+1),3]<=0,
                        0,
                        C1.lsoda[2:(length(time)+1),3])
           PKindex.this.subj<-data.frame(i,
                               C1.lsoda[2:(length(time)+1),1],
                               good)
           colnames(PKindex.this.subj)<-list("Subject","time","conc")
           show(PKindex.this.subj)
           x<-C1.lsoda[2:(length(time)+1),1]
           y<-good
           plotting.sim(i,x,y,MD);PKindex[[i]]<-PKindex.this.subj   ### dump this subj's data to subj i.
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
         par.err<-data.frame(Parameter=c("ka","kel","k12","k21","Vd"),error_factor=c(.15,.15,.15,.15,.15))
         par.err<-edit(par.err)
         repeat{
         if (par.err[1,2] <= 0 || par.err[2,2] <=0 || par.err[3,2]<=0 || par.err[4,2]<=0 || par.err[5,2]<=0){
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
         factor5<-par.err[5,2]
         
         PKindex<-vector(Subject,mode="list")
         pick<- pick-1  # reset pick = 1, 2, 3 or 4 ('0' is zero error & excluded)
         
         for(i in 1:Subject)  {
         
             ka  <- err_types(pick,par1,factor1)
             kel <- err_types(pick,par2,factor2)
             k12 <- err_types(pick,par3,factor3)
             k21 <- err_types(pick,par4,factor4)
             Vd  <- err_types(pick,par5,factor5)
             

               cat("\n\n     << Subject:- #",i,">>" )
               cat("\n\n")                                                               
               cat("*******************************************************\n")        
               cat(" Summary Table                                       \n\n")   
               if(MD){
               cat(" Model: 2-compartment, extravascular,                  \n") 
               cat("        multiple-dose, & 1st-ordered without lag time\n\n")}
               else{    
               cat(" Model: 2-compartment, extravascular,                  \n") 
               cat("        single-dose, & 1st-ordered without lag time\n\n")}
               cat("      Subject #:", Subject,"                       \n")
               cat("           Dose:", Dose,"                          \n") 
               if(MD) {
               cat("Dosing Interval:", Tau,"                           \n")
               cat("      # of Dose:", nDose,"                         \n")}
               cat("     Error Type:", type,"\n\n")     
               sim<-matrix(c(ka,kel,k12,k21,Vd,par1,par2,par3,par4,par5),5,2)
               dimnames(sim)<-list(c("ka","kel","k12","k21","Vd"),c("Simulated Values","Input Values"))
               show(sim)                                                               
               cat("*******************************************************\n\n")
                          
               ### PKindex[[i]]<-sfirst2.out(PKtime,ka,kel,k12,k21,Vd,defun,par1,par2,par3,par4,par5,Dose,i,type,MD)
           time<-PKtime$time
           parms<-c(ka=ka,kel=kel,k12=k12,k21=k21,Vd=Vd)
           if(MD){
              dosing.time<-seq(Tau,Tau*nDose,Tau)
              yini<-c(dCp1dt=Dose,dCp2dt=0,dCp3dt=0)
              events<- data.frame(var="dCp1dt",time=dosing.time,value=Dose,method="add")
              C1.lsoda<-data.frame(lsode(yini,c(0,time),defun,parms,rtol=1e-08,atol=1e-08,
                                   events=list(data=events)))
           }
           else{
             parms<-c(ka=ka,kel=kel,k12=k12,k21=k21,Vd=Vd) 
             C1.lsoda<-data.frame(lsoda(c(Dose,0,0),c(0,time),defun,parms,rtol=1e-08,atol=1e-08))}
           good<-ifelse(C1.lsoda[2:(length(time)+1),3]<=0,
                        0,
                        C1.lsoda[2:(length(time)+1),3])
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
     pick <- menu(file.menu, title = "<< Error Types >>")
     
     type<-switch(pick, 
                  "Error = Normal Error", 
                  "Error = Uniform Error",
                  "Error = Normal Error*True Value",
                  "Error = Uniform Error*True Value")
     
     cat("\n\nHow many iterations to run for each subject?\n")
     re<-scan(nlines=1,quiet=TRUE)
     
     par.err<-data.frame(Parameter=c("ka","kel","k12","k21","Vd"),error_factor=c(.15,.15,.15,.15,.15))
     par.err<-edit(par.err)
     repeat{
     if (par.err[1,2] <= 0 || par.err[2,2] <=0 || par.err[3,2]<=0 || par.err[4,2]<=0 || par.err[5,2]<=0){
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
     factor5<-par.err[5,2]
     
     sink(zz,split=TRUE)
     
     cat("\n\n")
     cat("******************************************************\n")
     cat(" Summary Table - Monte-Carlo simulation runs          \n")
     if(MD){
     cat(" Model: 2-compartment, extravascular,                  \n") 
     cat("        multiple-dose, & 1st-ordered without lag time\n\n")}
     else{    
     cat(" Model: 2-compartment, extravascular,                  \n") 
     cat("        single-dose, & 1st-ordered without lag time\n\n")}
     cat("    Subject #:", Subject,"                             \n")
     cat("   Error Type:", type,"                                \n")
     cat(" Simulation #:", re,"                                 \n\n")
     sim<-matrix(c(par1,par2,par3,par4,par5,factor1,factor2,factor3,factor4,factor5),5,2)
     dimnames(sim)<-list(c("ka","kel","k12","k21","Vd"),c("Input Values","Error factor"))
     show(sim)   
     cat("******************************************************")
     
     PKindex<-vector(Subject,mode="list")
     for(i in 1:Subject)  {
     cat("\n\n     << Subject:- #",i,">>\n" )
     C1.lsoda<-list()
        for (j in 1:re){

             ka  <- err_types(pick,par1,factor1)
             kel <- err_types(pick,par2,factor2)
             k12 <- err_types(pick,par3,factor3)
             k21 <- err_types(pick,par4,factor4)
             Vd  <- err_types(pick,par5,factor5)

               time1<-PKtime$time
               parms<-c(ka=ka,kel=kel,k12=k12,k21=k21,Vd=Vd)
               if(MD){
                  dosing.time<-seq(Tau,Tau*nDose,Tau)
                  yini<-c(dCp1dt=Dose,dCp2dt=0,dCp3dt=0)
                  events<- data.frame(var="dCp1dt",time=dosing.time,value=Dose,method="add")
                  XX<-data.frame(lsode(yini,c(0,time1),defun,parms,rtol=1e-08,atol=1e-08,
                                   events=list(data=events)))
               }
               else{
                  XX<-data.frame(lsoda(c(Dose,0,0),c(0,time1),defun,parms,rtol=1e-08,atol=1e-08))}
               C1.lsoda[[j]]<-data.frame(XX[2:(length(time1)+1),1],XX[2:(length(time1)+1),3])
               colnames(C1.lsoda[[j]])<-list("time","concentration") 
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
   if(MD) stwo.MD.all()
   else   stwo.SD.all()
}   
