#for one exponential term
smacro.one.out<-function(PKtime,A,a,defun,par1,par2,Dose,i,type) 
{
  time<-PKtime$time
  defun<- A*exp(-a*time)
  output<-data.frame(time,defun) 
  cat("\n\n")
  cat("*******************************************\n")
  cat("Summary Table                              \n")
  cat("Model: One-exponential Term Model          \n") 
  cat("Error Type:", type,"                       \n\n")
  sim<-matrix(c(A,a,par1,par2),2,2)
  dimnames(sim)<-list(c("A","a"),c("Value","Selected"))
  show(sim)
  cat("******************************************\n\n")
  
  PKindex<-data.frame(i,output)
  colnames(PKindex)<-list("Subject","time","conc")
  show(PKindex)
  x<-PKindex[,2]
  y<-PKindex[,3]
  plotting.sim(i,x,y)
  return(PKindex) 
}
