### PKindex is the target Dataset.
### Normal fitting
### One compartment PK model extravascular single 
### dose zero-order absorption without lag time (absorption rate, Ro = Dose/Tabs)
### optional Michaelis-Menten Elimination

fzero.nolag <- function(PKindex,
                        Dose=NULL, 
                        Tabs=NULL,
                        Vm=NULL,Km=NULL, ## MMe=TRUE
                        Vd=NULL,
                        kel=NULL,        ## MMe=FALSE
                        MMe=FALSE) 
{ 
   options(warn=-1)
   modfun1<-NULL
   modfun2<-NULL
        
   ## Input dose and initial value for Tabs, kel and Vd
   
   if (MMe){
      if (is.null(Dose)||is.null(Tabs)||is.null(Vm)||is.null(Km)||is.null(Vd)) {
        par.init<-data.frame(Parameter=c("Dose","Tabs","Vm","Km","Vd"),Initial=c(0,0,0,0,0))
        par.init<-edit(par.init)
        repeat{
           if (par.init[1,2] <= 0 || par.init[2,2] <= 0 || par.init[3,2]<= 0 || par.init[4,2]<= 0|| par.init[5,2]<= 0){
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
   } 
   else {
      ## No MM elimination
      if (is.null(Dose)||is.null(Tabs) || is.null(kel) || is.null(Vd)) {
        par.init<-data.frame(Parameter=c("Dose","Tabs","kel","Vd"),Initial=c(0,0,0,0))
        par.init<-edit(par.init)
        repeat{
           if (par.init[1,2] <= 0 || par.init[2,2] <= 0 || par.init[3,2]<= 0|| par.init[4,2]<= 0){
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
   }
        cat("\n")
        Dose<-par.init[1,2]       
        show(par.init);cat("\n")
   
   if (!MMe) {
      ## User-supplied function w/o Michaelis-Mention elimination
      defun<- function(time, y, parms) { 
      if(time<=parms["Tabs"]) 
        dCpdt <- (Dose/parms["Tabs"])/parms["Vd"] - parms["kel"] * y[1]
      else
        dCpdt <- - parms["kel"] * y[1]
      list(dCpdt) 
      } 
   
      modfun1 <<- function(time,Tabs,kel,Vd) { 
         out <- lsoda(0,c(0,time),defun,parms=c(Tabs=Tabs,kel=kel,Vd=Vd),
                      rtol=1e-6,atol=1e-6) 
         out[-1,2]
      } 
   } 
   else {
     ## User-supplied function with MM elimination
      defun<- function(time, y, parms) { 
      if(time<=parms["Tabs"]) 
        dCpdt <- (Dose/parms["Tabs"])/parms["Vd"]-(parms["Vm"]/parms["Vd"])*y[1]/(parms["Km"]/parms["Vd"]+y[1])            
      else
        dCpdt <- -(parms["Vm"]/parms["Vd"])*y[1]/(parms["Km"]/parms["Vd"]+y[1])            
      list(dCpdt) 
      }   
   
      modfun2 <<- function(time,Tabs,Vm,Km,Vd) { 
      out <- lsoda(0,c(0,time),defun,parms=c(Tabs=Tabs,Vm=Vm,Km=Km,Vd=Vd),
                   rtol=1e-6,atol=1e-6) 
      out[-1,2]
      } 
    }
    
   ## Select weighting schemes
   file.menu <- c("equal weight", 
                  "1/Cp",
                  "1/Cp^2")           
   pick <- menu(file.menu, title = "<< Weighting Schemes >>")

   with(entertitle(),{
###
### windows(record=TRUE)
dev.new()
par(mfrow=c(2,2),las=1)
pdf_activate=FALSE  ### set pdf device activate? as FALSE at beginning
###
### give warning below
###
cat("\n The following steps may go wrong. If so, please check\n")
cat("  your data, check your model and check initial values.\n\n")
readline(" Press Enter to continue...")
cat("\n\n")
###
### log to outputs.txt here
###
zz <- file("pkfit_fitting_outputs.txt", open="wt")
sink(zz,split=TRUE)   ### use sink(zz.split=TURE) will output to the txt file, as well as the screen at the same time. YJ
cat("\n\n");cat("--- input data ---\n")
show(PKindex);cat("\n\n")     # show input data    
cat("--- initial values for parameters ---\n")
show(par.init);cat("\n")    # show initial values here
cat("--- weighting scheme: ")
switch(pick,                  ## show weighting scheme
  cat("equal weight\n"),
  cat("1/Cp\n"),
  cat("1/Cp^2\n"));cat("\n")
if(MMe){
cat("--- model selection: a one-compartment pk model with\n    zero-ordered abs., M-M elim.")}
else{
cat("--- model selection: a one-compartment pk model with\n    zero-ordered abs., 1st-ordered elim.")}
   
   for( i in 1:length(unique(PKindex$Subject)))  {
     cat("\n\n               << Subject",i,">>\n\n" )  
     objfun <- function(par) {
        if (MMe) {
           out<-modfun2(PKindex$time[PKindex$Subject==i],par[1],par[2],par[3],par[4])
        } 
        else {
           ## No MM elimination
           out<-modfun1(PKindex$time[PKindex$Subject==i],par[1],par[2],par[3])
        }
      gift <- which( PKindex$conc[PKindex$Subject==i] != 0 )
      switch(pick,
             sum((PKindex$conc[PKindex$Subject==i][gift]-out[gift])^2),
             sum((PKindex$conc[PKindex$Subject==i][gift]-out[gift])^2/PKindex$conc[gift]),
             sum(((PKindex$conc[PKindex$Subject==i][gift]-out[gift])/PKindex$conc[gift])^2))
      }
      
###      
     if (MMe) {
        opt<-optim(c(par.init[2,2],par.init[3,2],par.init[4,2],par.init[5,2]),objfun, method="Nelder-Mead")  
        nameopt<-c("Tabs","Vm","Km","Vd")
        outopt<-c(opt$par[1],opt$par[2],opt$par[3],opt$par[4])
     }
     else {
        opt<-optim(c(par.init[2,2],par.init[3,2],par.init[4,2]),objfun, method="Nelder-Mead")  
        nameopt<-c("Tabs","kel","Vd")
        outopt<-c(opt$par[1],opt$par[2],opt$par[3])
     }
     
     cat("\n<< PK parameters obtained from Nelder-Mead Simplex algorithm >>\n\n")
     print(data.frame(Parameter=nameopt,Value=outopt))
      if (MMe){
              if(opt$par[1]<0) {opt$par[1]<-0.01}
              if(opt$par[2]<0) {opt$par[2]<-0.01}
              if(opt$par[3]<0) {opt$par[3]<-0.01}
              if(opt$par[4]<0) {opt$par[4]<-0.01}
       }
       else {
              if(opt$par[1]<0) {opt$par[1]<-0.01}
              if(opt$par[2]<0) {opt$par[2]<-0.01}
              if(opt$par[3]<0) {opt$par[3]<-0.01}
       }      
     
     cat("\n<< Residual sum-of-square (RSS) and final PK parameters with nlsLM >>\n\n")

     if (MMe) {
        fm<-nlsLM(conc~modfun2(time,Tabs,Vm,Km,Vd),data=subset(PKindex,Subject==i),
            start=list(Tabs=opt$par[1],Vm=opt$par[2],Km=opt$par[3],Vd=opt$par[4]),
            control=nls.lm.control(maxiter=500),lower=c(0,0,0,1e-06)) ### lower of Vd should not be zero due to Dose/Vd. --YJ
        plotting.non(PKindex, fm, i, pick,xaxis,yaxis)
###
### copied from the original plotting.lin()
###
     main<-paste(c("Subject# ", i),collapse=" ")
     j<-1:length(PKindex$time[PKindex$Subject==i])
     xx<-PKindex$time[PKindex$Subject==i]
     yy<-PKindex$conc[PKindex$Subject==i]
     cal<-predict(fm,list(time=xx))
     wei <- switch(pick,
               ifelse(yy[j]==0.0, 0, yy[j]-cal[j]),
               ifelse(yy[j]==0.0, 0, sqrt(1/(yy[j]))*(yy[j]-cal[j])),
               ifelse(yy[j]==0.0, 0, sqrt(1/((yy[j])^2))*(yy[j]-cal[j])))
     
    #Linear plot
     plot(yy~xx,data=PKindex,type='p',main=main, 
          xlab=xaxis, ylab=yaxis,pch=15,col="black",bty="l",
          font.lab=2,cex.lab=1,cex.axis=1,cex.main=1) 
     lines(xx,predict(fm,list(time=xx)),type="l",lty=1,
           col="firebrick3",lwd="2")
     mtext("Linear",side=3,cex=0.88)
       
    #Semi-log plot
     plot(xx,yy,log="y",type='p',main=main,
          xlab=xaxis, ylab=yaxis,pch=15,col="black",bty="l",
          font.lab=2,cex.lab=1,cex.axis=1,cex.main=1) 
     lines(xx,predict(fm,list(time=xx)),type="l",lty=1,
           col="firebrick3",lwd="2")
     mtext("Semi-log",side=3,cex=0.88)
        
    #Residual plot, time vs weighted residual
     plot(xx,wei,pch=15,col="blue",bty="l",xlab=xaxis,
          ylab="Weighted Residual",main="Residual Plots",cex.lab=1,
          cex.axis=1,cex.main=1,font.lab=2)
     abline(h=0,lwd=2,col="black",lty=2)
       
    #Residual plot, calculated concentration vs weigthed residual
     plot(cal,wei,pch=15,col="blue",bty="l",xlab="Calc Cp(i)",
          ylab="Weighted Residual",main="Weighted Residual Plots",cex.lab=1,
          cex.axis=1,cex.main=1,font.lab=2)
     abline(h=0,lwd=2,col="black",lty=2) 
###         
### here revert between pdf() and graphic device                          ### added by YJ
### 
          if(pdf_activate){
             dev.copy()                      ## copy to pdf file 2nd plots to end
             dev.set(which=x11c)             ## back from graphic device now to continue...
                          }
          else{
             x11c<-dev.cur()                 ## the current graphics device
             pdf(file="pkfit_plots.pdf",     ## activate pdf log file from now on... starting with ref. product
                  paper="a4")
###             description_plot()              ## bear output logo
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
        fm<-nlsLM(conc~modfun1(time,Tabs,kel,Vd),data=subset(PKindex,Subject==i),
            start=list(Tabs=opt$par[1],kel=opt$par[2],Vd=opt$par[3]),
         control=nls.lm.control(maxiter=500),lower=c(0,0,1e-06)) ### lower of Vd should not be zero due to Dose/Vd. --YJ
        coef<-data.frame(coef(fm)["kel"])
        plotting.lin(PKindex, fm, i, pick, coef, xaxis, yaxis)
###
### copied from the original plotting.lin()
###
     main<-paste(c("Subject# ", i),collapse=" ")
     j<-1:length(PKindex$time[PKindex$Subject==i])
     xx<-PKindex$time[PKindex$Subject==i]
     yy<-PKindex$conc[PKindex$Subject==i]
     cal<-predict(fm,list(time=xx))
     wei <- switch(pick,
               ifelse(yy[j]==0.0, 0, yy[j]-cal[j]),
               ifelse(yy[j]==0.0, 0, sqrt(1/(yy[j]))*(yy[j]-cal[j])),
               ifelse(yy[j]==0.0, 0, sqrt(1/((yy[j])^2))*(yy[j]-cal[j])))
     
    #Linear plot
     plot(yy~xx,data=PKindex,type='p',main=main, 
          xlab=xaxis, ylab=yaxis,pch=15,col="black",bty="l",
          font.lab=2,cex.lab=1,cex.axis=1,cex.main=1) 
     lines(xx,predict(fm,list(time=xx)),type="l",lty=1,
           col="firebrick3",lwd="2")
     mtext("Linear",side=3,cex=0.88)
       
    #Semi-log plot
     plot(xx,yy,log="y",type='p',main=main,
          xlab=xaxis, ylab=yaxis,pch=15,col="black",bty="l",
          font.lab=2,cex.lab=1,cex.axis=1,cex.main=1) 
     lines(xx,predict(fm,list(time=xx)),type="l",lty=1,
           col="firebrick3",lwd="2")
     mtext("Semi-log",side=3,cex=0.88)
        
    #Residual plot, time vs weighted residual
     plot(xx,wei,pch=15,col="blue",bty="l",xlab=xaxis,
          ylab="Weighted Residual",main="Residual Plots",cex.lab=1,
          cex.axis=1,cex.main=1,font.lab=2)
     abline(h=0,lwd=2,col="black",lty=2)
       
    #Residual plot, calculated concentration vs weigthed residual
     plot(cal,wei,pch=15,col="blue",bty="l",xlab="Calc Cp(i)",
          ylab="Weighted Residual",main="Weighted Residual Plots",cex.lab=1,
          cex.axis=1,cex.main=1,font.lab=2)
     abline(h=0,lwd=2,col="black",lty=2) 
###         
### here revert between pdf() and graphic device                          ### added by YJ
### 
          if(pdf_activate){
             dev.copy()                      ## copy to pdf file 2nd plots to end
             dev.set(which=x11c)             ## back from graphic device now to continue...
                          }
          else{
             x11c<-dev.cur()                 ## the current graphics device
             pdf(file="pkfit_plots.pdf",     ## activate pdf log file from now on... starting with ref. product
                  paper="a4")
###             description_plot()              ## bear output logo
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
  cat(" All outputs (pkfit_fitting_outputs.txt & pkfit_plots.pdf)\n can be found at",getwd(),"\n")
  readline(" Press any key to continue...")
  dev.off()        # close pdf()
  graphics.off()   # close plot windows
  })
  cat("\n")
  run()
}