### PKindex is the target Dataset.
### no finish yet...  YJ (04/07/2013)

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
   #options(warn=-1)
        
   ## Input dose and Tlag and initial value for ka, kel and Vd

   if (is.null(Dose)) {
     cat("Enter the Dose:\n")
     Dose <- scan(nlines=1,quiet=TRUE)
   } 
   else {
     cat("Dose is = ",Dose,"\n")
   }

   if ( Tlag ){
       cat("\nEnter the lag time value:\n")
       Tlag<-scan(nlines=1,quiet=TRUE)
       cat("\n")
   } 
   
   
   if (MMe){
      if (is.null(ka) || is.null(Vm) || is.null(Km) || is.null(Vd) ) {
        par<-data.frame(Parameter=c("ka","Vm","Km","Vd"),Initial=c(0))
        par<-edit(par)
        repeat{
           if ( par[1,2] == 0 || par[2,2] ==0 || par[3,2]==0 ||par[4,2]==0){
             cat("\n")
             cat("**********************************\n")
             cat(" Parameter value can not be zero. \n")
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

      }
   } 
   else{
      if (is.null(ka) || is.null(kel) || is.null(Vd) ) {
        par<-data.frame(Parameter=c("ka","kel","Vd"),Initial=c(0))
        par<-edit(par)
        repeat{
           if ( par[1,2] == 0 || par[2,2] ==0 || par[3,2]==0){
             cat("\n")
             cat("**********************************\n")
             cat(" Parameter value can not be zero. \n")
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

      }
   }
   
   
   cat("\n")
   
   if ( MMe ) { 
      ## User-supplied function w Michaelis-Mention elimination & w lag time
      defun<- function(time, y, parms) { 
      if(time <= Tlag) {
        dy1dt<-0
        dy2dt<-0
      }
      else  {
        dy1dt <- -parms["ka"] * y[1]
        dy2dt <-  parms["ka"] * y[1]/parms["Vd"] - (parms["Vm"]/parms["Vd"])*y[2]/(parms["Km"]/parms["Vd"]+y[2])
      }
      list(c(dy1dt,dy2dt)) 
      } 
    
      modfun2 <- function(time,ka,Vm,Km,Vd) { 
      out <- lsoda(c(Dose,0),c(0,time),defun,parms=c(ka=ka,Vm=Vm,Km=Km,Vd=Vd),rtol=1e-6,atol=1e-6) 
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
    
      modfun1 <- function(time,ka,kel,Vd) { 
      out <- lsoda(c(Dose,0),c(0,time),defun,parms=c(ka=ka,kel=kel,Vd=Vd),rtol=1e-6,atol=1e-6) 
      out[-1,3] 
      }  
   }
   ## Select weighting schemes
   file.menu <- c("equal weight", 
                  "1/Cp",
                  "1/Cp^2")           
   pick <- menu(file.menu, title = "<< Weighting Schemes >>")

   with(entertitle(),{
###
windows(record=TRUE)
par(mfrow=c(2,2),las=1)
pdf_activate=FALSE  ### set pdf device activate? as FALSE at beginning

###
### log to outputs.txt here
###
zz <- file("pkfit_fitting_outputs.txt", open="wt")
sink(zz,split=TRUE)   ### use sink(zz.split=TURE) will output to the txt file, as well as the screen at the same time. YJ
###
### give warning below
###
cat("\n The following steps may go wrong, if so please check your model,\n")
cat(" check your data and check your initial values next time.\n\n")
readline(" Press Enter to continue...")
cat("\n\n")
   

   for( i in 1:length(unique(PKindex$Subject)))  {
      cat("\n\n               << Subject",i,">>\n\n" )  
      objfun <- function(par) {
        if (MMe) {
           out <- modfun2(PKindex$time[PKindex$Subject==i], par[1], par[2],par[3],par[4]) 
        } 
        else  {
           out <- modfun1(PKindex$time[PKindex$Subject==i], par[1], par[2], par[3])
        }
     gift <- which( PKindex$conc[PKindex$Subject==i] != 0 )
     switch(pick,
            sum((PKindex$conc[PKindex$Subject==i][gift]-out[gift])^2),
            sum((PKindex$conc[PKindex$Subject==i][gift]-out[gift])^2/PKindex$conc[gift]),
            sum(((PKindex$conc[PKindex$Subject==i][gift] - out[gift])/PKindex$conc[gift])^2))
     }
###      if (MMe) {
###           gen<-genoud(objfun,nvars=4,max=FALSE,pop.size=30,max.generations=20,
###                wait.generations=10,starting.value=c(par[1,2],par[2,2],par[3,2],par[4,2]),
###                BFGS=FALSE,print.level=0,boundary.enforcement=2,
###                Domains=matrix(c(0.01,1,1,1,10,100,100,100),4,2),
###                MemoryMatrix=TRUE) 
###       } 
###       else {
###           gen<-genoud(objfun,nvars=3,max=FALSE,pop.size=20,max.generations=15,
###                wait.generations=10,starting.value=c(par[1,2],par[2,2],par[3,2]),
###                BFGS=FALSE,print.level=0,boundary.enforcement=2,
###                Domains=matrix(c(0.01,0.01,1,10,1,100),3,2),
###                MemoryMatrix=TRUE)
###       }
###       cat("<< PK parameters obtained from genetic algorithm >>\n\n")
###       if (MMe) {
###         namegen<-c("ka","Vm","Km","Vd")
###         outgen<-c(gen$par[1],gen$par[2],gen$par[3],gen$par[4])
###       } 
###       else {
###         ## No MM elimination
###         namegen<-c("ka","kel","Vd")
###         outgen<-c(gen$par[1],gen$par[2],gen$par[3])
###       }
###       print(data.frame(Parameter=namegen,Value=outgen))  
###       F<-objfun(gen$par)
     
      if (MMe) {
        opt<-optim(c(par[1,2],par[2,2],par[3,2],par[4,2]),objfun, method="Nelder-Mead")
        nameopt<-c("ka","Vm","Km","Vd")
        outopt<-c(opt$par[1],opt$par[2],opt$par[3],opt$par[4])
      }
      else {
        opt<-optim(c(par[1,2],par[2,2],par[3,2]),objfun, method="Nelder-Mead") 
        nameopt<-c("ka","kel","Vd")
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
      cat("\n<< Residual sum-of-square (RSS) and final PK parameters with nls >>\n\n")
      
      if (MMe) {
        fm<-nls(conc ~ modfun2(time, ka, Vm, Km, Vd),data=subset(PKindex,Subject==i),
            start=list(ka=opt$par[1],Vm=opt$par[2],Km=opt$par[3],Vd=opt$par[4]),trace=TRUE,
            nls.control(tol=1)) ### it seems MM should use 'nls.control(tol=1)'; otherwise if use 'port' it can crash...  --YJ
        cat("\n")
        plotting.non(PKindex, fm, i, pick, xaxis, yaxis)
###
### copied from the original plotting.lin()
###
         main<-paste(c("Plots for Subject# ", i),collapse=" ")
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
        fm <-nls(conc ~ modfun1(time, ka, kel, Vd), data=PKindex,subset=Subject==i,
             start=list(ka=opt$par[1],kel=opt$par[2],Vd=opt$par[3]),trace=TRUE,
             nls.control(maxiter=5000,tol=1e-06,minFactor=1/1024/1024),algorithm = "port",lower=c(0,0,0))
        cat("\n")
        coef<-data.frame(coef(fm)["kel"])
        plotting.lin(PKindex, fm, i, pick, coef, xaxis, yaxis)
###
### copied from the original plotting.lin()
###
         main<-paste(c("Plots for Subject# ", i),collapse=" ")
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
  dev.off()        # close pdf()
  graphics.off()   # close plot windows
  sink()           # reset sink()
  close(zz)        # close outputs.txt
   
   })
   cat("\n")
   run()   
} 