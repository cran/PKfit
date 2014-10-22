### PKindex is the target Dataset.
### Normal fitting
### Three compartment PK model iv bolus single dose
fbolus3<- function(PKindex,
                   Dose=NULL, 
                   kel=NULL,
                   k12=NULL,  
                   k21=NULL,
                   k13=NULL,
                   k31=NULL,      
                   Vd=NULL) 
{
   options(warn=-1)
   modfun<-NULL
   
   fit.outputs_to_txt<-fit.outputs_to_txt
   fit.plots_to_pdf<-fit.plots_to_pdf
      
   ## Input dose and initial value for kel, k12, k21 and Vd
  
   if(file.exists("fbolus3.csv")){
        par.init<-read.csv(file="fbolus3.csv",row.names=NULL,header=TRUE)
        par.init<-edit(par.init)}
   else{
        par.init<-data.frame(Parameter=c("Dose","kel","k12","k21","k13","k31","Vd"),Initial=c(0,0,0,0,0,0,0))
        par.init<-edit(par.init)
        repeat{
           if (par.init[1,2]<=0||par.init[2,2]<=0||par.init[3,2]<=0||par.init[4,2]<=0||par.init[5,2]<=0||par.init[6,2]<=0||par.init[7,2]<=0){
             cat("\n")
             cat("**********************************\n")
             cat(" Parameter initial value can not be zero. \n")
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
        write.csv(par.init,file="fbolus3.csv",row.names=FALSE,col.names=TRUE)
        cat("\n")
        Dose<-par.init[1,2]        
        show(par.init)
   
   cat("\n")
   
   defun<- function(time, y, parms) { 
     dCp1dt <- -parms["kel"]*y[1]-parms["k12"]*y[1]+parms["k21"]*y[2]-parms["k13"]*y[1]+parms["k31"]*y[3] 
     dCp2dt <-  parms["k12"]*y[1]-parms["k21"]*y[2]
     dCp3dt <-  parms["k13"]*y[1]-parms["k31"]*y[3]
     list(c(dCp1dt,dCp2dt,dCp3dt)) 
   } 

   modfun <<- function(time,kel,k12,k21,k13,k31,Vd) { 
      out <- lsoda(y=c(Dose/Vd,0,0),c(0,time),defun,parms=c(kel=kel,k12=k12,k21=k21,k13=k13,k31=k31,Vd=Vd),
                  rtol=1e-06,atol=1e-10) 
     ### plot(out)
     out[-1,2] 
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
description_version()
cat("\n\n")
sink()  ### turn off temporarily to avoid logging too many warnings... -YJ
   
   for( i in 1:length(unique(PKindex$Subject)))  {
     objfun <- function(par) {
        out <- modfun(PKindex$time[PKindex$Subject==i], par[1], par[2], par[3], par[4], par[5], par[6])
        gift <- which(PKindex$conc[PKindex$Subject==i] != 0)
        ### sum(((PKindex$conc[PKindex$Subject==i][gift]-out[gift])/PKindex$conc[gift])^2)
        switch(pick,
               sum((PKindex$conc[PKindex$Subject==i][gift]-out[gift])^2),
               sum((PKindex$conc[PKindex$Subject==i][gift]-out[gift])^2/PKindex$conc[gift]),
               sum(((PKindex$conc[PKindex$Subject==i][gift]-out[gift])/PKindex$conc[gift])^2))
     }
     
###      
     opt <- optimx(c(par.init[2,2],par.init[3,2],par.init[4,2],par.init[5,2],par.init[6,2],par.init[7,2]),
                 objfun,method="Nelder-Mead",control=list(maxit=5000))
     nameopt<-c("kel","k12","k21","k13","k31","Vd")
     outopt<-c(opt$p1,opt$p2,opt$p3,opt$p4,opt$p5,opt$p6)
     
              if(opt$p1<0) {opt$p1<-0.0001}
              if(opt$p2<0) {opt$p2<-0.0001}
              if(opt$p3<0) {opt$p3<-0.0001}
              if(opt$p4<0) {opt$p4<-0.0001}
              if(opt$p5<0) {opt$p5<-0.0001}
              if(opt$p6<0) {opt$p6<-0.0001}
     
     conc<-PKindex$conc[PKindex$Subject==i]
     
     if(pick==1) weights<- ifelse(conc==0.,1,1/conc^0)  ### equal weight
     if(pick==2) weights<- ifelse(conc==0.,1,1/conc^1)  ### 1/Cp
     if(pick==3) weights<- ifelse(conc==0.,1,1/conc^2)  ### 1/Cp^2
          
     fm<-nlsLM(conc ~ modfun(time,kel,k12,k21,k13,k31,Vd),data=subset(PKindex,Subject==i),weights=weights,
         start=list(kel=opt$p1,k12=opt$p2,k21=opt$p3,k13=opt$p4,k31=opt$p5,Vd=opt$p6),
         control=nls.lm.control(maxiter=500,maxfev=5000,factor=100),lower=c(1e-06,1e-06,1e-06,1e-06,1e-06,1e-06))  ### set 'lower=c(...)' may cause crashed.  --YJ
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
     cat("--- model selection: a three-compartment, iv bolus pk model with\n    1st-ordered elim.\n\n")
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
             x11c<-dev.cur()                 ## the current graphics device
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
