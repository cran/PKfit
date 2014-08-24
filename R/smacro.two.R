### Two exponential term
### 'B' term can be negative in this case; such as one-compartment, PO, 1st-ordered pk model. Cp(t)=A(exp(-kel*t)-exp(-ka*t))
### will need to create another code for this.  --YJ

smacro.two <-function(Subject=NULL,  # N Subj's 
                      PKtime=NULL,   # times for sampling
                      A=NULL,
                      alpha=NULL,
                      B=NULL,
                      beta=NULL,
                      MD=FALSE)
{
   options(warn=-1)
   sim.outputs_to_txt<-sim.outputs_to_txt
   sim.plots_to_pdf<-sim.plots_to_pdf
      
cat("\n First enter all paramaters used for simulation profile.\n");readline(" Press Enter to continue...\n\n")
if(file.exists("smacro_two.csv")){
   par<-read.csv(file="smacro_two.csv",row.names=NULL,header=TRUE)
   par<-edit(par)}
else{ 
   par<-data.frame(Parameter=c("Total_subject#","A","alpha","B","beta"),Initial=c(24,11.6,0.12,16.7,0.51))
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
  }
   write.csv(par,file="smacro_two.csv",row.names=FALSE,col.names=TRUE) 
   cat("\n")       
   show(par)
   
   cat("\n")
   Subject<-par[1,2]
   par1<-par[2,2]
   par2<-par[3,2]
   par3<-par[4,2]
   par4<-par[5,2]

   readline("\n Next enter or load time points. Press Enter to continue.\n\n")
   if(file.exists("sim_times.csv")){
      PKtime<-read.csv(file="sim_times.csv",row.names=NULL,header=TRUE)
      PKtime<-edit(PKtime)}
   else{
      PKtime<-data.frame(time=c(0))
      ### PKtime<-data.frame(time=c(0,.1,.2,.3,.4,.6,.8,1,2,4,6,8,12,14,16,18,24,48,72))
      PKtime<-edit(PKtime)}
   
   write.csv(PKtime,file="sim_times.csv",row.names=FALSE,col.names=TRUE)
   cat("\n")
   show(PKtime);cat("\n\n")
    
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
            A    <-par1
            alpha<-par2  
            B    <-par3
            beta <-par4
            PKindex[[i]]<-smacro.two.out(PKtime,A,alpha,B,beta,defun,par1,par2,par3,par4,i,type,MD)
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
         par.err<-data.frame(Parameter=c("A","alpha","B","beta"),error_factor=c(.15,.15,.15,.15))
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
         
         PKindex<-vector(Subject,mode="list")
         pick<- pick-1  # reset pick = 1, 2, 3 or 4 ('0' is zero error & excluded)
         
         for( i in 1:Subject)  { 

             A     <- err_types(pick,par1,factor1)
             alpha <- err_types(pick,par2,factor2)
             B     <- err_types(pick,par3,factor3)
             beta  <- err_types(pick,par4,factor4)
                      
            cat("\n\n       << Subject:- #",i,">>" )    ### used to implement. -YJ
            PKindex[[i]]<-smacro.two.out(PKtime,A,alpha,B,beta,defun,par1,par2,par3,par4,i,type,MD)
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
      else if (pick ==2){ ### starting monte-carlo sim here
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
         par.err<-data.frame(Parameter=c("A","alpha","B","beta"),error_factor=c(.15,.15,.15,.15))
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
      
      sink(zz,split=TRUE)
      
      cat("\n\n")
      cat("**********************************\n")
      cat(" Summary Table                     \n")
      cat(" Model: two-exponential term model \n") 
      cat(" Subject #:", Subject,"            \n")
      cat(" Error Type:", type,"              \n")
      cat(" Simulation #:", re,"              \n\n")
      sim<-matrix(c(par1,par2,par3,par4,factor1,factor2,factor3,factor4),4,2)
      dimnames(sim)<-list(c("A","alpha","B","beta"),c("Input Values","Error factor"))
      show(sim)   
      cat("**********************************")
      PKindex<-vector(Subject,mode="list")
      for( i in 1:Subject)  {
        cat("\n\n   << Subject:- #",i," >>\n\n" )    ### used to implement. -YJ
        C1.lsoda<-list()
           for (j in 1:re){
           
             A     <- err_types(pick,par1,factor1)
             alpha <- err_types(pick,par2,factor2)
             B     <- err_types(pick,par3,factor3)
             beta  <- err_types(pick,par4,factor4)
                        
              time1<-PKtime$time
              defun<- A*exp(-alpha*time1)+B*exp(-beta*time1) 
              XX<-data.frame(time1,defun) 
              C1.lsoda[[j]]<-data.frame(XX$time1,XX$defun) 
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
   PK.sim()
}   
