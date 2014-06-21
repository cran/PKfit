#---------------plot for linear----------------
plotting.lin <- function (PKindex, fm, i, pick, coef, xaxis, yaxis,
                          separateWindows=TRUE)
{               
  #options(warn=-1)
  j<-1:length(PKindex$time[PKindex$Subject==i])
  x<-PKindex$time[PKindex$Subject==i]
  y<-PKindex$conc[PKindex$Subject==i]
  
  ### if(x[1]!=0){    ### when (0,0) is not included
  ###   tmpXY<-data.frame(x,y)
  ###   tmpZZ<-data.frame(x=0,y=0)
  ###   tmpXY<-rbind(tmpZZ,tmpXY)
  ###   time<-tmpXY$x;x<-tmpXY$x
  ###   conc<-tmpXY$y;y<-tmpXY$y
  ###   j<-j+1
  ### }
     
 #Calculated concentration
  cal<-predict(fm,list(time=x))
    
 #Weighted residuals   
 if (!(pick %in% 1:3)) {
    stop("Pick is illegal")
  }
  
 wei <- switch(pick,
               ifelse(y[j]==0.0, 0, y[j]-cal[j]),
               ifelse(y[j]==0.0, 0, sqrt(1/(y[j]))*(y[j]-cal[j])),
               ifelse(y[j]==0.0, 0, sqrt(1/((y[j])^2))*(y[j]-cal[j])))
  
 #calculate AUC and AUMC   
  add<-function(time,conc){
     auc<-0 ; aumc<-0
     for(i in 2:length(time)) {
     auc[i]<-1/2*(time[i]-time[i-1])*(conc[i]+conc[i-1])
     auc[i]<-auc[i]+auc[i-1]
     aumc[i]<-1/2*(time[i]-time[i-1])*(conc[i]*time[i]+conc[i-1]*time[i-1])
     aumc[i]<-aumc[i]+aumc[i-1]
     }
     return(list(auc=auc,aumc=aumc))
  }
  add1<-add(x,y)
  if (x[1]==0){
    AUC<-add1$auc
    AUMC<-add1$aumc
    }
  else {
  AUC<-c(NA,add1$auc[-1])
  AUMC<-c(NA,add1$aumc[-1])
  }
        
 #Output   
  cat("<< Residual sum-of-squares and final PK parameters values with nlsLM >>\n\n")
  output<-data.frame(x,y,cal,wei,AUC,AUMC)
  colnames(output)<-list("Time","Observed","Calculated","Wt. Residuals","AUC","AUMC")
  show(fm);cat("\n");show(output)  
  
 #AUC (0 to infinity)              
  ### cat(paste(c("\n<< AUC (0 to",max(y),"computed by trapezoidal rule >>\n\n"))
  auc.infinity<-y[length(y)]/coef[1,1]
  auc<-AUC[length(y)]+auc.infinity
  ### show(auc) 
  
 #AUMC (0 to infinity) 
  ### cat(paste(c("\n<< AUMC (0 to",max(y),"computed by trapezoidal rule >>\n\n")
  aumc.infinity<-(x[length(x)]*y[length(y)])/coef[1,1]+y[length(y)]/((coef[1,1])^2)
  aumc<-AUMC[length(y)]+aumc.infinity
  ### show(aumc)
  sumNCA<-data.frame(Parameters=c("Cmax","Cmin","AUC0-inf","AUMC0-inf"),
                     values=c(max(y),min(y),auc,aumc))
  cat("\n");show(sumNCA)
         
  aicllsbc(fm)
  cat("\n")

par(mfrow=c(2,2))
main<-paste(c("Subject# ", i),collapse=" ")
j<-1:length(PKindex$time[PKindex$Subject==i])
xxstep<-seq(from=min(PKindex$time[PKindex$Subject==i]),to=max(PKindex$time[PKindex$Subject==i]),by=0.01)
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
lines(xxstep,predict(fm,list(time=xxstep)),type="l",lty=1,
      col="firebrick3",lwd="2")
mtext("Linear",side=3,cex=0.88)
  
#Semi-log plot
plot(xx,yy,log="y",type='p',main=main,
     xlab=xaxis, ylab=yaxis,pch=15,col="black",bty="l",
     font.lab=2,cex.lab=1,cex.axis=1,cex.main=1) 
lines(xxstep,predict(fm,list(time=xxstep)),type="l",lty=1,
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

#Divide plot console into four parts
###   if (separateWindows) {
###     windows(record=TRUE)
###   }
###   par(mfrow=c(2,2))
###   main<-paste(c("Subject# ", i ,"plots"),collapse=" ")
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
###  #Residual plot, time vs weighted residual
###   plot(x,wei,pch=15,col="blue",bty="l",xlab=xaxis,
###        ylab="Weighted Residual",main="Residual Plots",cex.lab=1,
###        cex.axis=1,cex.main=1,font.lab=2)
###   abline(h=0,lwd=2,col="black",lty=2)
###     
###  #Residual plot, calculated concentration vs weigthed residual
###   plot(cal,wei,pch=15,col="blue",bty="l",xlab="Calc Cp(i)",
###        ylab="Weighted Residual",main="Residual Plots",cex.lab=1,
###        cex.axis=1,cex.main=1,font.lab=2)
###   abline(h=0,lwd=2,col="black",lty=2)          
}    
