#for two compartment extravascular first order absorption
sfirst2.out<-function(PKtime,ka,kel,k12,k21,Vd,defun,par1,par2,par3,par4,par5,Dose,i,type)
{                   
  time<-PKtime$time
  parms<-c(ka=ka,kel=kel,k12=k12,k21=k21,Vd=Vd) 
  C1.lsoda<-data.frame(lsoda(c(Dose,0,0),c(0,time),defun,parms,rtol=1e-6,atol=1e-6))
  cat("\n\n")
  cat("*******************************************************\n")
  cat(" Summary Table                                        \n\n")
  cat(" Model: 2-compartment, extravascular,                   \n") 
  cat("        single-dose, & 1st-ordered without lag time model \n") 
  cat(" Error Type:", type,"                                 \n\n")
  sim<-matrix(c(ka,kel,k12,k21,Vd,par1,par2,par3,par4,par5),5,2)
  dimnames(sim)<-list(c("ka","kel","k12","k21","Vd"),c("Value","Selected"))
  show(sim)
##  readline()  ### pasue here
  cat("*******************************************************\n\n")
  good<-ifelse(C1.lsoda[2:(length(time)+1),3]<=0,
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
