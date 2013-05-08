#for extravascular zero order absorption without lag time
szero.out<-function(PKtime,Tabs,kel,Vd,defun,par1,par2,par3,Dose,i,type) 
{
  time<-PKtime$time
  parms<-c(Tabs=Tabs,kel=kel,Vd=Vd)  
  C1.lsoda<-data.frame(lsoda(0, c(0,time), defun, parms))   
  cat("\n")
  cat("********************************************************\n")
  cat("Summary Table                                           \n")
  cat("Model: 1-Compartment, Extravascular, Single-Dose,       \n")
  cat("       & Zero-Ordered Absorption without Lag Time Model \n")
  cat("Error Type:", type,"                                    \n\n")
  sim<-matrix(c(Tabs,kel,Vd,par1,par2,par3),3,2)
  dimnames(sim)<-list(c("Tabs","kel","Vd"),c("Value","Selected"))
  show(sim)  
  cat("********************************************************\n\n")
  
  good<-ifelse(C1.lsoda[2:(length(time)+1),2]<=0,
               0,
               C1.lsoda[2:(length(time)+1),2])
  PKindex<-data.frame(i,
                      C1.lsoda[2:(length(time)+1),1],
                      good)
  colnames(PKindex)<-list("Subject","time","conc")
  show(PKindex)
  x<-C1.lsoda[2:(length(time)+1),1]
  y<-good
  plotting.sim(i,x,y)
  return(PKindex) 
}
