### ----------------plot for simulation----------------
### Draw 2 windows, and consider that 2 more windows are coming
plotting.sim <- function(i,x,y,MD=FALSE) 
{

### dev.new()
par(mfrow=c(2,1), ask = FALSE)

main<-paste(c("Subject:-", i),collapse="")
   
if(MD){
  ## linear plot
  plot(y~x,type='l',main=main,     ### plot lines only('l'ine only)
       xlab="Time after dosing",ylab="Simulated drug plasma Conc.",pch=15,col="black",bty="l",
       font.lab=2,cex.lab=1,cex.axis=1,cex.main=1) 
  text("Linear",side=3,cex=0.88)

  ## semi-log plot
  plot(x,y,log="y",type='l',main=main,
       xlab="Time after dosing",ylab="Simulated drug plasma Conc.",pch=15,col="black",bty="l",
       font.lab=2,cex.lab=1,cex.axis=1,cex.main=1) 
  mtext("Semi-log",side=3,cex=0.88) 
}
else{
  ## linear plot
  plot(y~x,type='b',main=main,    ### plot lines and symbols ('b'oth)
       xlab="Time after dosing",ylab="Simulated drug plasma Conc.",pch=15,col="black",bty="l",
       font.lab=2,cex.lab=1,cex.axis=1,cex.main=1) 
  text("Linear",side=3,cex=0.88)

  ## semi-log plot
  plot(x,y,log="y",type='b',main=main,
  ### plot(x,y,log="y",type='l',main=main,
       xlab="Time after dosing",ylab="Simulated drug plasma Conc.",pch=15,col="black",bty="l",
       font.lab=2,cex.lab=1,cex.axis=1,cex.main=1) 
  mtext("Semi-log",side=3,cex=0.88) 
}
}
