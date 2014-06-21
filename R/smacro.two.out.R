#for two exponential term
smacro.two.out<-function(PKtime,A,alpha,B,beta,defun,par1,par2,par3,par4,i,type,MD=FALSE)
{
  time<-PKtime$time
  defun<- A*exp(-alpha*time)+B*exp(-beta*time)
  output<-data.frame(time,defun) 
  cat("\n\n")
  cat("*******************************************\n")
  cat(" Summary Table                              \n")
  cat(" Model: two-exponential term model          \n") 
  cat(" Error Type:", type,"                       \n\n")
  sim<-matrix(c(A,alpha,B,beta,par1,par2,par3,par4),4,2)
  dimnames(sim)<-list(c("A","alpha","B","beta"),c("Simulated Values","Input Values"))
  show(sim)
  cat("*******************************************\n\n")
  PKindex<-data.frame(i,output)
  colnames(PKindex)<-list("Subject","time","conc")
  show(PKindex)
  x<-PKindex[,2]
  y<-PKindex[,3]
  plotting.sim(i,x,y,MD)
  return(PKindex) 
}
