### ----------------plot for simulation----------------
### Draw 2 windows, and consider that 2 more windows are coming
plotting.sim <- function(i,x,y,separateWindows=TRUE) 
{
  #options(warn=-1)
  ### par(Ask=TRUE,las=1)   ### YJ
windows(record=TRUE)      ### not to log into a pdf file; keep as it was here.  -YJ

  par(mfrow=c(2,1))
  main<-paste(c("Plot of Subject#_ ", i),collapse=" ")
   
  ## linear plot
  plot(y~x,type='b',main=main, 
       xlab="Time",ylab="Simulated drug plasma Conc.",pch=15,col="black",bty="l",
       font.lab=2,cex.lab=1,cex.axis=1,cex.main=1) 
  text("Linear",side=3,cex=0.88)

  ## semi-log plot
  plot(x,y,log="y",type='b',main=main,
       xlab="Time",ylab="Simulated drug plasma Conc.",pch=15,col="black",bty="l",
       font.lab=2,cex.lab=1,cex.axis=1,cex.main=1) 
  mtext("Semi-log",side=3,cex=0.88) 
  
  cat("\n\n")
}
