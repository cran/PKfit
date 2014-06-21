#for one exponential term
smacro.one.out<-function(PKtime,A,alpha,defun,par1,par2,i,type,MD=FALSE) 
{
  time<-PKtime$time
  defun<- A*exp(-alpha*time)
  output<-data.frame(time,defun) 
  cat("\n\n")
  cat("*******************************************\n")
  cat(" Summary Table                              \n")
  cat(" Model: one-exponential term model          \n") 
  cat(" Error Type:", type,"                       \n\n")
  sim<-matrix(c(A,alpha,par1,par2),2,2)
  dimnames(sim)<-list(c("A","alpha"),c("Simulated Values","Input Values"))
  show(sim)
  cat("******************************************\n\n")
  
  PKindex<-data.frame(i,output)
  colnames(PKindex)<-list("Subject","time","conc")
  show(PKindex)
  x<-PKindex[,2]
  y<-PKindex[,3]
  plotting.sim(i,x,y,MD)
  return(PKindex) 
}
