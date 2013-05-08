### ---------------plot for nonlinear----------------
plotting.non <- function (PKindex, fm, i, pick, xaxis, yaxis,
                          separateWindows=TRUE) 
{  
  #options(warn=-1)
  j <- 1:length(PKindex$time[PKindex$Subject==i])
  x <- PKindex$time[PKindex$Subject==i]
  y <- PKindex$conc[PKindex$Subject==i]
   
  ## Calculated concentration
  cal<-predict(fm,list(time=x))

  if (!(pick %in% 1:3)) {
    stop("Your pick is invalid.")
  }
  wei <- switch(pick,
                ifelse(y[j]==0.0, 0, y[j]-cal[j]),
                ifelse(y[j]==0.0, 0, sqrt(1/(y[j]))*(y[j]-cal[j])),
                ifelse(y[j]==0.0, 0, sqrt(1/((y[j])^2))*(y[j]-cal[j])))
  
 #calculate AUC and AUMC       
  add <- function(time,conc) {
    auc<-0 ; aumc<-0
    for(i in 2:length(time)) {
      auc[i]<-1/2*(time[i]-time[i-1])*(conc[i]+conc[i-1])
      auc[i]<-auc[i]+auc[i-1]
      aumc[i]<-1/2*(time[i]-time[i-1])*(conc[i]*time[i]+conc[i-1]*time[i-1])
      aumc[i]<-aumc[i]+aumc[i-1]
    }
    return(list(auc=auc,aumc=aumc))
  }
  
  add1 <- add(x,y)
  AUC<-c(NA,add1$auc[-1])
  AUMC<-c(NA,add1$aumc[-1])
              
  cat("<< Output >>\n\n")      
  output <- data.frame(x,y,cal,wei,AUC,AUMC)
  colnames(output) <- list("Time","Observed","Calculated","Wt. Residuals","AUC","AUMC")
  show(output)
  cat("\n") 
        
  aicllsbc(fm)
  cat("\n\n")
    
###   ## Divide plot console into four parts
###   if (separateWindows) {
###     windows(record=TRUE)
###   }
###   par(mfrow=c(2,2))
###   main<-paste(c("Subject", i ,"plot"),collapse=" ")
###      
###  #Linear plot
###   plot(y~x,data=PKindex,type='p',main=main, 
###        xlab=xaxis, ylab=yaxis,pch=15,col="black",bty="l",
###        font.lab=2,cex.lab=1,cex.axis=1,cex.main=1) 
###   lines(x,predict(fm,list(time=x)),type="l",lty=1,
###         col="firebrick3",lwd="2")
###   mtext("Linear",side=3,cex=0.88)
###   
###  #Semi-log plot
###   plot(x,y,log="y",type='p',main=main,
###        xlab=xaxis, ylab=yaxis,pch=15,col="black",bty="l",
###        font.lab=2,cex.lab=1,cex.axis=1,cex.main=1) 
###   lines(x,predict(fm,list(time=x)),type="l",lty=1,
###         col="firebrick3",lwd="2")
###   mtext("Semi-log",side=3,cex=0.88)
###   
###   ## Residual plot, time vs weighted residual----- 
###   plot(x,wei,pch=15,col="blue",bty="l",xlab=xaxis,
###        ylab="Weighted Residual",main="Residual Plots",cex.lab=1,
###        cex.axis=1,cex.main=1,font.lab=2)
###   abline(h=0,lwd=2,col="black",lty=2)
###   
###   ## Residual plot, calcukated concentration vs weigthed residual-----
###   plot(cal,wei,pch=15,col="blue",bty="l",xlab="Calc Cp(i)",
###        ylab="Weighted Residual",main="Residual Plots",cex.lab=1,
###        cex.axis=1,cex.main=1,font.lab=2)
###   abline(h=0,lwd=2,col="black",lty=2)
}    
