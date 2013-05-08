#for three exponential term
smacro.three.out<-function(PKtime,A,a,B,b,C,c,defun,par1,par2,par3,par4,par5,par6,Dose,i,type)      
{         
  time<-PKtime$time
  defun<- A*exp(-a*time)+B*exp(-b*time)+C*exp(-c*time)
  output<-data.frame(time,defun) 
  cat("\n\n")
  cat("*******************************************\n")
  cat("Summary Table                              \n")
  cat("Model: Three-exponential Term Model        \n") 
  cat("Error Type:", type,"                       \n\n")
  sim<-matrix(c(A,a,B,b,C,c,par1,par2,par3,par4,par5,par6),6,2)
  dimnames(sim)<-list(c("A","a","B","b","C","c"),c("Value","Selected"))
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
