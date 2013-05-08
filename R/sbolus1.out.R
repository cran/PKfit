#for iv bolus 
sbolus1.out<-function(PKtime,kel,Vd,defun,par1,par2,Dose,i,type)
{
  time<-PKtime$time
  parms<-c(kel=kel,Vd=Vd) 
  C1.lsoda<-data.frame(lsoda(Dose/Vd,c(0,time),defun,parms)) 
  cat("\n")
  cat("****************************************************\n")
  cat("Summary Table                                       \n")
  cat("Model: 1-Compartment, IV-Bolus, & Single-Dose Model \n")
  cat("Error Type:", type,"                                \n\n") 
  sim<-matrix(c(kel,Vd,par1,par2),2,2)
  dimnames(sim)<-list(c("kel","Vd"),c("Value","Selected"))
  show(sim)
  cat("****************************************************\n\n")
  
  good <- ifelse(C1.lsoda[2:(length(time)+1),2]<=0,
                 0,
                 C1.lsoda[2:(length(time)+1),2])
  PKindex <- data.frame(i,
                        C1.lsoda[2:(length(time)+1),1],
                        good)
  colnames(PKindex)<-list("Subject","time","conc")
  show(PKindex)
  x <- C1.lsoda[2:(length(time)+1),1]
  y <- good
  plotting.sim(i,x,y)
  return(PKindex) 
}