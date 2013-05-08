#for two compartment iv bolus
sbolus2.out<-function(PKtime,kel,k12,k21,Vd,defun,par1,par2,par3,par4,Dose,i,type)
{ 
  time<-PKtime$time
  parms<-c(kel=kel,k12=k12,k21=k21,Vd=Vd) 
  C1.lsoda<-data.frame(lsoda(c(Dose/Vd,0),c(0,time),defun,parms)) 
  cat("\n\n")
  cat("****************************************************\n")
  cat("Summary Table                                       \n")
  cat("Model: 2-Compartment, IV-Bolus, & Single-Dose Model \n") 
  cat("Error Type:", type,"                                \n\n")
  sim<-matrix(c(kel,k12,k21,Vd,par1,par2,par3,par4),4,2)
  dimnames(sim)<-list(c("kel","k12","k21","Vd"),c("Value","Selected"))
  show(sim)
  cat("****************************************************\n\n")
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
