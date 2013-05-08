#for extravascular first order absorption with/without lag time
sfirst1.out<-function(PKtime,ka,kel,Vd,defun,par1,par2,par3,Dose,i,type,Tlag)  
{
  time<-PKtime$time
  parms<-c(ka=ka,kel=kel,Vd=Vd)  
  C1.lsoda<-data.frame(lsoda(c(Dose,0),c(0,time), defun, parms))  
  
  if (missing(Tlag)){
  
  cat("\n\n")
  cat("*****************************************************\n")
  cat("Summary Table                                        \n")
  cat("Model: 1-Compartment, Extravascular, Single-Dose,    \n")
  cat("       & 1-Ordered Absorption without Lag Time Model \n")
  cat("Error Type:", type,"                                 \n\n")
  sim<-matrix(c(ka,kel,Vd,par1,par2,par3),3,2)
  dimnames(sim)<-list(c("ka","kel","Vd"),c("Value","Selected"))
  show(sim)
  cat("*****************************************************\n\n")
  }
  
  else{
  
  cat("\n\n")
  cat("**************************************************\n")
  cat("Summary Table                                     \n")
  cat("Model: 1-Compartment, Extravascular, Single-Dose, \n")
  cat("       & 1-Ordered Absorption with Lag Time Model \n")
  cat("Error Type:", type,"                              \n\n")
  sim<-matrix(c(ka,kel,Vd,par1,par2,par3),3,2)
  dimnames(sim)<-list(c("ka","kel","Vd"),c("Value","Selected"))
  show(sim)
  cat("**************************************************\n\n")
  }  
  
  good<-ifelse(C1.lsoda[2:(length(time)+1),3]<=1e-5,
               0,
               C1.lsoda[2:(length(time)+1),3])
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
