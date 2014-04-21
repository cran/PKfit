#for extravascular zero order absorption without lag time nonlinear elimination
szero.mm.out<-function(PKtime,Tabs,Vm,Km,Vd,defun,par1,par2,par3,par4,Dose,i,type) 
{
  time<-PKtime$time
  parms<-c(Tabs=Tabs,Vm=Vm,Km=Km,Vd=Vd) 
  C1.lsoda<-data.frame(lsoda(0,c(0,time), defun, parms))   
  cat("\n\n")
  cat("***************************************************************************\n")
  cat(" Summary Table                                                              \n")
  cat(" Model: 1-compartment, extravascular, single-dose, zero-ordered absorption, \n") 
  cat("        & Michaelis-Menten elimination without lag time model               \n")
  cat(" Error Type:", type,"                                                       \n\n")
  sim<-matrix(c(Tabs,Vm,Km,Vd,par1,par2,par3,par4),4,2)
  dimnames(sim)<-list(c("Tabs","Vm","Km","Vd"),c("Value","Selected"))
  show(sim)
  cat("***************************************************************************\n\n")
  
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
