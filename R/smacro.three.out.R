#for three exponential term
smacro.three.out<-function(PKtime,A,alpha,B,beta,C,gamma,defun,par1,par2,par3,par4,par5,par6,i,type)      
{         
  time<-PKtime$time
  defun<- A*exp(-alpha*time)+B*exp(-beta*time)+C*exp(-gamma*time)
  output<-data.frame(time,defun) 
  cat("\n\n")
  cat("*******************************************\n")
  cat(" Summary Table                              \n")
  cat(" Model: three-exponential term model        \n") 
  cat(" Error Type:", type,"                       \n\n")
  sim<-matrix(c(A,alpha,B,beta,C,gamma,par1,par2,par3,par4,par5,par6),6,2)
  dimnames(sim)<-list(c("A","alpha","B","beta","C","gamma"),c("Value","Selected"))
  show(sim)
  cat("*******************************************\n\n")
  PKindex<-data.frame(i,output)
  colnames(PKindex)<-list("Subject","time","conc")
  show(PKindex)
  x<-PKindex[,2]
  y<-PKindex[,3]
  plotting.sim(i,x,y)
  return(PKindex) 
}
