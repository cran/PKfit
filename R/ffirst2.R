### Two compartment PK model extravascular single dose first order absorption
ffirst2<- function(PKindex,
                  ka=NULL,
                  Dose=NULL, 
                  kel=NULL,
                  k12=NULL,  
                  k21=NULL,      
                  Vd=NULL) 
{
   #options(warn=-1)
   
   ## Input dose and initial value for ka, kel, k12, k21 and Vd

   if (is.null(Dose)) {
     cat("Enter Dose\n")
     Dose <- scan(nlines=1,quiet=TRUE)
   } 
   else {
     cat("Dose from arguments is = ",Dose,"\n")
   }
   
   if (is.null(ka) || is.null(kel) || is.null(k12) || is.null(k21) || is.null(Vd) ) {
        par<-data.frame(Parameter=c("ka","kel","k12","k21","Vd"),Initial=c(0))
        par<-edit(par)
        repeat{
           if ( par[1,2] == 0 || par[2,2] ==0 || par[3,2]==0 || par[4,2]==0 || par[5,2]==0){
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
   
   cat("\n")
   
   defun<- function(time, y, parms) { 
     dCp1dt <- -parms["ka"]*y[1]
     dCp2dt <- parms["ka"]*y[1]/parms["Vd"]-parms["kel"]*y[2]+parms["k21"]*y[3]-parms["k12"]*y[2]
     dCp3dt <- parms["k12"]*y[2]-parms["k21"]*y[3]
     list(c(dCp1dt,dCp2dt,dCp3dt)) 
   } 
    
   modfun3 <- function(time,ka,kel,k12,k21,Vd) { 
     out <- lsoda(c(Dose,0,0),c(0,time),defun,parms=c(ka=ka,kel=kel,k12=k12,k21=k21,Vd=Vd),
                  rtol=1e-6,atol=1e-6) 
     out[-1,3] 
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
        out <- modfun3(PKindex$time[PKindex$Subject==i],par[1],par[2],par[3],par[4],par[5])
        gift <- which( PKindex$conc[PKindex$Subject==i] != 0 )
        switch(pick,
               sum((PKindex$conc[PKindex$Subject==i][gift]-out[gift])^2),
               sum((PKindex$conc[PKindex$Subject==i][gift]-out[gift])^2/PKindex$conc[gift]),
               sum(((PKindex$conc[PKindex$Subject==i][gift] - out[gift])/PKindex$conc[gift])^2))
     }
     
###      gen<-genoud(objfun,nvars=5,max=FALSE,pop.size=15,max.generations=10,    ### seems OK to abandon genetic algorithm here (waste a lot of time to do it) YJ
###           wait.generations=5,
###           starting.value=c(par[1,2],par[2,2],par[3,2],par[4,2],par[5,2]),
###           BFGS=FALSE,print.level=0,boundary.enforcement=2,
###           Domains=matrix(c(0.01,0.01,0.01,0.01,1,10,1,10,1,100),5,2),
###           MemoryMatrix=TRUE)     
###           
###      cat("<< PK parameters obtained from genetic algorithm >>\n\n")
###      
###      namegen<-c("ka","kel","k12","k21","Vd")
###      outgen<-c(gen$par[1],gen$par[2],gen$par[3],gen$par[4],gen$par[5])
###      print(data.frame(Parameter=namegen,Value=outgen)) 
###      F<-objfun(gen$par)
###      if(gen$par[1]<0) {gen$par[1]<-0.01}
###      if(gen$par[2]<0) {gen$par[2]<-0.01}
###      if(gen$par[3]<0) {gen$par[3]<-0.01}
###      if(gen$par[4]<0) {gen$par[4]<-0.01}
###      if(gen$par[5]<0) {gen$par[5]<-0.01}
###      
###      opt<-optim(c(gen$par[1],gen$par[2],gen$par[3],gen$par[4],gen$par[5]),objfun, method="Nelder-Mead")
     opt<-optim(c(par[1,2],par[2,2],par[3,2],par[4,2],par[5,2]),objfun, method="Nelder-Mead")
     nameopt<-c("ka","kel","k12","k21","Vd")
     outopt<-c(opt$par[1],opt$par[2],opt$par[3],opt$par[4],opt$par[5])
     
     cat("\n<< PK parameter obtained from Nelder-Mead Simplex algorithm >>\n\n")
     print(data.frame(Parameter=nameopt,Value=outopt))
     ### to avoid negative values from Nelder-Mead Simplex; otherwise, it will crash before nls().  -YJ
     if(opt$par[1]<0) {opt$par[1]<-0.01}
     if(opt$par[2]<0) {opt$par[2]<-0.01}
     if(opt$par[3]<0) {opt$par[3]<-0.01}
     if(opt$par[4]<0) {opt$par[4]<-0.01}
     if(opt$par[5]<0) {opt$par[5]<-0.01}
     
     cat("\n<< Residual sum-of-squares and final PK parameters values with nls >>\n\n")
     fm<-nls(conc ~ modfun3(time,ka,kel,k12,k21,Vd),data=subset(PKindex,Subject==i),  ### for nls(); nlsLM() cannot work out. -YJ
         start=list(ka=opt$par[1],kel=opt$par[2],k12=opt$par[3],k21=opt$par[4],Vd=opt$par[5]),trace=TRUE,
         nls.control(maxiter=5000,tol=1e-06,minFactor=1/1024/1024),algorithm = "port",lower=c(0,0,0,0,0))  ### I added nls(...,algorithm = "port", lower=c(0,0,0,0,0)) seems stable now.
                                                                                                           ### only 'port' can set up lower limits in nls().
     cat("\n")
     coef<-data.frame(coef(fm)["kel"])     ### for nls() only? is the same for nlsLM()? try it... YJ
     
     
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
  dev.off()        # close pdf()
  graphics.off()   # close plot windows
  sink()           # reset sink()
  close(zz)        # close outputs.txt
  })
  cat("\n")
  run()
}

