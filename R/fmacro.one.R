### PKindex is the target Dataset.


### Normal fitting
### One exponential term

fmacro.one <- function(PKindex,
                       A=NULL,
                       a=NULL) 
{
   options(warn=-1)
   defun<-NULL    ### only exponential terms use 'defun'

   fit.outputs_to_txt<-fit.outputs_to_txt
   fit.plots_to_pdf<-fit.plots_to_pdf
      
   ## Input and initial value for A and a
   
   if(file.exists("fmacro_one.csv")){
        par.init<-read.csv(file="fmacro_one.csv",row.names=NULL,header=TRUE)
        par.init<-edit(par.init)}
   else{
       par.init<-data.frame(Parameter=c("A","alpha"),Initial=c(0,0))
       ### par.init<-data.frame(Parameter=c("A","alpha"),Initial=c(10,0.1))
       par.init<-edit(par.init)
       repeat{
           if ( par.init[1,2]<= 0 || par.init[2,2]<=0){
             cat("\n")
             cat("**********************************\n")
             cat(" Parameter value can not be zero. \n")
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
        write.csv(par.init,file="fmacro_one.csv",row.names=FALSE,col.names=TRUE)
        cat("\n")       
        show(par.init)
   
   cat("\n")
   
   defun<<-function(time,A,alpha){
     pred<- A*exp(-alpha*time)
   }
   
   ## Select weighting schemes
   file.menu <-c("equal weight",  ## may cause error with equal weight; wait to be solved. --YJ
                 "1/Cp",
                 "1/Cp^2")
   pick <- menu(file.menu, title="<< weighting schemes >>")

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
         out<-defun(PKindex$time[PKindex$Subject==i],par[1],par[2])
         gift<- which( PKindex$conc[PKindex$Subject==i] != 0 )
         sum(((PKindex$conc[PKindex$Subject==i][gift]-out[gift])/PKindex$conc[gift])^2)
         ### switch(pick,
         ###        sum((PKindex$conc[PKindex$Subject==i][gift]-out[gift])^2),
         ###        sum((PKindex$conc[PKindex$Subject==i][gift]-out[gift])^2/PKindex$conc[gift]),
         ###        sum(((PKindex$conc[PKindex$Subject==i][gift]-out[gift])/PKindex$conc[gift])^2)
         ###        )
      }
      
      opt<-optim(c(par.init[1,2],par.init[2,2]),objfun,method="Nelder-Mead",control=list(maxit=5000))       
      nameopt<-c("A","alpha")
      outopt<-c(opt$par[1],opt$par[2])
      
            if(opt$par[1]<0) {opt$par[1]<-0.01}
            if(opt$par[2]<0) {opt$par[2]<-0.01}

     conc<-PKindex$conc[PKindex$Subject==i]
     
     if(pick==1) weights=(1/conc^0)  ### equal weight
     if(pick==2) weights=(1/conc^1)  ### 1/Cp
     if(pick==3) weights=(1/conc^2)  ### 1/Cp^2
            
     fm<-nlsLM(conc~defun(time,A,alpha),data=subset(PKindex,Subject==i),start=list(A=opt$par[1],alpha=opt$par[2]),
              weights=weights, control=nls.lm.control(maxiter=500))  ### set 'lower=c(...)' may cause crashed.  --YJ
     coef<-data.frame(coef(fm)["alpha"])     
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
     cat("--- model selection: a one-exponential macroconstant\n\n")       
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
