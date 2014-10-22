### PKindex is the target Dataset.
### Normal fitting
### One compartment PK model extravascualr single dose first-order absorption
### With lag time

ffirst.lag<- function(PKindex,
                      Dose=NULL, 
                      ka=NULL,
                      Vm=NULL,Km=NULL, 
                      Vd=NULL,
                      kel=NULL,        
                      Tlag=TRUE,
                      MMe=FALSE) 
{
   options(warn=-1)
   modfun1<-NULL
   modfun2<-NULL
   
   fit.outputs_to_txt<-fit.outputs_to_txt
   fit.plots_to_pdf<-fit.plots_to_pdf
      
   ## Input dose and Tlag and initial value for ka, kel and Vd

   if (MMe){
      if(file.exists("ffirst_lag_mm.csv")){
        par.init<-read.csv(file="ffirst_lag_mm.csv",row.names=NULL,header=TRUE)
        par.init<-edit(par.init)}
      else{
        par.init<-data.frame(Parameter=c("Dose","Tlag","ka","Vm","Km","Vd"),Initial=c(0,0,0,0,0,0))
        par.init<-edit(par.init)
        repeat{
           if (par.init[1,2]<=0||par.init[2,2]<=0||par.init[3,2]<=0||par.init[4,2]<=0||par.init[5,2]<=0||par.init[6,2]<=0){
             cat("\n")
             cat("**********************************\n")
             cat(" Parameter initial values can not be zero. \n")
             cat(" Press Enter to continue.         \n")
             cat("**********************************\n\n")
             readline()
             cat("\n")
             par.init<-edit(par.init)}   
           else{
             break
             return(edit(par.init))}
        } 
      }
           write.csv(par.init,file="ffirst_lag_mm.csv",row.names=FALSE,col.names=TRUE)
   } 
   else{
      if(file.exists("ffirst_lag.csv")){
        par.init<-read.csv(file="ffirst_lag.csv",row.names=NULL,header=TRUE)
        par.init<-edit(par.init)}
      else{
        par.init<-data.frame(Parameter=c("Dose","Tlag","ka","kel","Vd"),Initial=c(0,0,0,0,0))
        par.init<-edit(par.init)
        repeat{
           if (par.init[1,2]<=0 || par.init[2,2]<=0 || par.init[3,2]<=0|| par.init[4,2]<=0|| par.init[5,2]<=0){
             cat("\n")
             cat("**********************************\n")
             cat(" Parameter initial values can not be zero. \n")
             cat(" Press Enter to continue.         \n")
             cat("**********************************\n\n")
             readline()
             cat("\n")
             par.init<-edit(par.init)}   
           else{
             break
             return(edit(par.init))}
        } 
      }
           write.csv(par.init,file="ffirst_lag.csv",row.names=FALSE,col.names=TRUE)
   }
        cat("\n")
        Dose<-par.init[1,2]
        Tlag<-par.init[2,2]       
        show(par.init);cat("\n")
   
   cat("\n")
   
   if (MMe) { 
      ## User-supplied function w Michaelis-Mention elimination & w lag time
      defun<- function(time, y, parms) { 
      if(time <= Tlag) {
        dy1dt<-0
        dy2dt<-0
      }
      else  {
        dy1dt <- -parms["ka"] * y[1]
        dy2dt <-  parms["ka"] * y[1]/parms["Vd"] - parms["Vm"]*y[2]/(parms["Km"]+y[2])
      }
      list(c(dy1dt,dy2dt)) 
      } 
    
      modfun2 <<- function(time,ka,Vm,Km,Vd) { 
      out <- lsoda(c(Dose,0),c(0,time),defun,parms=c(ka=ka,Vm=Vm,Km=Km,Vd=Vd),rtol=1e-6,atol=1e-10) 
      out[-1,3] 
      } 
   } 
   else {
      ## User-supplied function w/o MM elimination w lag time
      defun<- function(time, y, parms) { 
      if(time<=Tlag) {
       dy1dt<-0
       dy2dt<-0
      }
      else {
       dy1dt <- -parms["ka"]*y[1]
       dy2dt <-  parms["ka"]*y[1]/parms["Vd"]-parms["kel"]*y[2]
      }
      list(c(dy1dt,dy2dt)) 
      } 
    
      modfun1 <<- function(time,ka,kel,Vd) { 
      out <- lsoda(c(Dose,0),c(0,time),defun,parms=c(ka=ka,kel=kel,Vd=Vd),rtol=1e-08,atol=1e-08) 
      out[-1,3] 
      }  
   }
   ## Select weighting schemes
   file.menu <- c("equal weight", 
                  "1/Cp",
                  "1/Cp^2")           
   pick <- menu(file.menu, title = "<< Weighting Schemes >>")

   with(entertitle(),{
### give warning below
###
cat("\n The following steps may go wrong. If so, please check\n")
cat(" your data, your model, initial values and/or weightings.\n\n")
readline(" Press Enter to continue...");cat("\n\n")
cat(" Please wait...\n\n")
###
dev.new()
par(mfrow=c(2,2),las=1)
pdf_activate=FALSE  ### set pdf device activate? as FALSE at beginning
###
###
### log to outputs.txt here
###
zz <- file(fit.outputs_to_txt, open="wt")
sink(zz,split=TRUE)   ### use sink(zz.split=TURE) will output to the txt file, as well as the screen at the same time. YJ
description_version();cat("\n\n")
sink()  ### turn off temporarily to avoid logging too many warnings... -YJ

   for( i in 1:length(unique(PKindex$Subject)))  {
      objfun <- function(par) {
        if (MMe) {
           out <- modfun2(PKindex$time[PKindex$Subject==i], par[1], par[2],par[3],par[4]) 
        } 
        else  {
           out <- modfun1(PKindex$time[PKindex$Subject==i], par[1], par[2], par[3])
        }
        gift <- which(PKindex$conc[PKindex$Subject==i] != 0)
        ### sum(((PKindex$conc[PKindex$Subject==i][gift]-out[gift])/PKindex$conc[gift])^2)
        switch(pick,
               sum((PKindex$conc[PKindex$Subject==i][gift]-out[gift])^2),
               sum((PKindex$conc[PKindex$Subject==i][gift]-out[gift])^2/PKindex$conc[gift]),
               sum(((PKindex$conc[PKindex$Subject==i][gift]-out[gift])/PKindex$conc[gift])^2))
     }
     
      if (MMe) {
        opt<-optimx(c(par.init[3,2],par.init[4,2],par.init[5,2],par.init[6,2]),objfun,
                   method="Nelder-Mead",control=list(maxit=5000))
        nameopt<-c("ka","Vm","Km","Vd")
        outopt<-c(opt$p1,opt$p2,opt$p3,opt$p4)
      }
      else {
        opt<-optimx(c(par.init[3,2],par.init[4,2],par.init[5,2]),objfun,
                  method="Nelder-Mead",control=list(maxit=5000)) 
        nameopt<-c("ka","kel","Vd")
        outopt<-c(opt$p1,opt$p2,opt$p3)
      }
      
      if (MMe){
              if(opt$p1<0) {opt$p1<-0.0001}
              if(opt$p2<0) {opt$p2<-0.0001}
              if(opt$p3<0) {opt$p3<-0.0001}
              if(opt$p4<0) {opt$p4<-0.0001}
       }
       else {
              if(opt$p1<0) {opt$p1<-0.0001}
              if(opt$p2<0) {opt$p2<-0.0001}
              if(opt$p3<0) {opt$p3<-0.0001}
       }      
     
     conc<-PKindex$conc[PKindex$Subject==i]
     
     if(pick==1) weights<- ifelse(conc==0.,1,1/conc^0)  ### equal weight
     if(pick==2) weights<- ifelse(conc==0.,1,1/conc^1)  ### 1/Cp
     if(pick==3) weights<- ifelse(conc==0.,1,1/conc^2)  ### 1/Cp^2
           
     if (MMe) {
        fm<-nlsLM(conc ~ modfun2(time, ka, Vm, Km, Vd),data=subset(PKindex,Subject==i),
                  start=list(ka=opt$p1,Vm=opt$p2,Km=opt$p3,Vd=opt$p4),weights=weights,
                  control=nls.lm.control(maxiter=500,maxfev=5000,factor=100),lower=c(1e-06,1e-06,1e-06,1e-06))  ### set 'lower=c(...)' may cause crashed.  --YJ
     sink(zz,split=TRUE)
     cat(" ********************************\n\n")
     cat("      --- Subject:- #",i,"---    \n\n")
     cat(" ********************************\n\n")
     cat("--- input data ---\n")
     conc<-PKindex$conc[PKindex$Subject==i]
     time<-PKindex$time[PKindex$Subject==i]
     this_subj<-data.frame(time, conc)
     show(this_subj);cat("\n")     # show input data    
     cat("--- initial values for parameters ---\n")
     show(par.init);cat("\n")    # show initial values here
     cat("--- weighting scheme: ")
     switch(pick,                  ## show weighting scheme
       cat("equal weight"),
       cat("1/Cp"),
       cat("1/Cp^2"));cat("\n\n")
     cat("--- model selection: one-compartment, extravascular PK model\n   with 1st-ordered lag-time abs. and M-M elim.\n\n")
     cat("<< PK parameter obtained from Nelder-Mead Simplex algorithm >>\n\n")
     print(data.frame(Parameter=nameopt,Value=outopt));cat("\n")     
     plotting.non(PKindex, fm, i, pick, xaxis, yaxis)
     sink()  ### turn off temporarily to avoid logging too many warnings... -YJ     
###         
### here revert between pdf() and graphic device                          ### added by YJ
### 
          if(pdf_activate){
             dev.copy()                      ## copy to pdf file 2nd plots to end
             dev.set(which=x11c)             ## back from graphic device now to continue...
                          }
          else{
             x11c<-dev.cur()                     ## the current graphics device
             pdf(fit.plots_to_pdf,paper="a4")
             pdf_activate=TRUE               ## set pdf_activate=TRUE from now on
             dev.set(which=x11c)             ## go to graphics device...
             dev.copy()                      ## copy the first plot here
             dev.set(which=x11c)             ## back from graphics device
              }
###
###  end plotting here...
###        
        
      } 
      else {
        ## No MM elimination
        fm <-nlsLM(conc ~ modfun1(time, ka, kel, Vd), data=PKindex,subset=Subject==i,
             start=list(ka=opt$p1,kel=opt$p2,Vd=opt$p3),weights=weights,
             control=nls.lm.control(maxiter=500,maxfev=5000,factor=100),lower=c(1e-06,1e-06,1e-06))    ### set 'lower=c(...)' may cause crashed.  --YJ
        coef<-data.frame(coef(fm)["kel"])
     sink(zz,split=TRUE)
     cat(" ********************************\n\n")
     cat("      --- Subject:- #",i,"---    \n\n")
     cat(" ********************************\n\n")
     cat("--- input data ---\n")
     conc<-PKindex$conc[PKindex$Subject==i]
     time<-PKindex$time[PKindex$Subject==i]
     this_subj<-data.frame(time, conc)
     show(this_subj);cat("\n")     # show input data    
     cat("--- initial values for parameters ---\n")
     show(par.init);cat("\n")    # show initial values here
     cat("--- weighting scheme: ")
     switch(pick,                  ## show weighting scheme
       cat("equal weight"),
       cat("1/Cp"),
       cat("1/Cp^2"));cat("\n\n")
     cat("--- model selection: one-compartment, extravascular PK model\n   with 1st-ordered lag-time abs./elim.\n\n")
     cat("<< PK parameter obtained from Nelder-Mead Simplex algorithm >>\n\n")
     print(data.frame(Parameter=nameopt,Value=outopt));cat("\n")     
     plotting.lin(PKindex, fm, i, pick, coef, xaxis, yaxis)
     sink()  ### turn off temporarily to avoid logging too many warnings... -YJ     
###         
### here revert between pdf() and graphic device                          ### added by YJ
### 
          if(pdf_activate){
             dev.copy()                      ## copy to pdf file 2nd plots to end
             dev.set(which=x11c)             ## back from graphic device now to continue...
                          }
          else{
             x11c<-dev.cur()                     ## the current graphics device
             pdf(fit.plots_to_pdf,paper="a4")
             pdf_activate=TRUE               ## set pdf_activate=TRUE from now on
             dev.set(which=x11c)             ## go to graphics device...
             dev.copy()                      ## copy the first plot here
             dev.set(which=x11c)             ## back from graphics device
              }
###
###  end plotting here...
###
        
      }
   }
  sink()           # reset sink()
  close(zz)        # close outputs.txt
  cat(paste(" Two outputs,",fit.outputs_to_txt,"&",fit.plots_to_pdf,",\n have been generated at",getwd(),"\n\n"))
  readline(" Press any key to continue...")
  dev.off()        # close pdf()
  graphics.off()   # close plot windows   
   })
  cat("\n")
  ### run()   
  PK.fit(PKindex)
} 