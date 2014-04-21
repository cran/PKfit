### Simulation
###One exponential term
smacro.one <-function(Subject=NULL,  # N Subj's 
                      PKtime=NULL,   # times for sampling
                      A=NULL,
                      alpha=NULL)
{
   options(warn=-1)
   cat("\n First enter all paramaters used for simulation profile.\n");readline(" Press Enter to continue...\n\n")
   par<-data.frame(Parameter=c("Total_subject#","A","alpha"),Initial=c(24,11.6,0.12))
   par<-edit(par)                                                       
   repeat{
       if (par[1,2] <= 0 || par[2,2] <=0 || par[3,2]<=0){
         cat("\n")
         cat("**********************************\n")
         cat(" Any parameter value can not be equal to or less than zero.\n")
         cat(" Press Enter to continue.         \n")
         cat("**********************************\n\n")
         readline()
         cat("\n")
         par<-edit(par)}   
       else{
         break
         return(edit(par))}
    } 
    cat("\n")       
    show(par)
   
   cat("\n")
   Subject<-par[1,2]
   par1<-par[2,2]
   par2<-par[3,2]

   readline("\n Next enter time points. Press Enter to continue.\n\n")
   PKtime<-data.frame(time=c(0))
   PKtime<-edit(PKtime)
   cat("\n")
    
    file.menu <- c("Simulation with Error",
                   "Monte Carlo Simulation")
    pick <- menu(file.menu, title = "<< Simulation Types >>") 
    
    if (pick ==1){
      cat("\n\n")
      file.menu <- c("No Error",
                     "Error = Normal Error", 
                     "Error = Uniform Error",
                     "Error = Normal Error*True Value",
                     "Error = Uniform Error*True Value")
      pick <- menu(file.menu, title = "<< Error types >>")      
      
      type<-switch(pick, 
                   "No Error",
                   "Error = Normal Error", 
                   "Error = Uniform Error",
                   "Error = Normal Error*True Value",
                   "Error = Uniform Error*True Value")
      
      if (pick ==1){
          PKindex<-vector(Subject,mode="list");cat("\n")
          for(i in 1:Subject)  {
            cat("\n     << Subject:- #",i,">>\n\n" )
            A<-par1
            alpha<-par2  
            PKindex[[i]]<-smacro.one.out(PKtime,A,alpha,defun,par1,par2,i,type)
            }  
      PKindex<- as.data.frame(do.call("rbind",PKindex))
      rownames(PKindex) <- seq(nrow(PKindex)) 
      savefile(PKindex)  
      }  
      else {
         par.err<-data.frame(Parameter=c("A","alpha"),error_factor=c(.15,.15))
         par.err<-edit(par.err)
         repeat{
         if (par.err[1,2] <= 0 || par.err[2,2] <=0){
             cat("\n")
             cat("**********************************\n")
             cat(" Parameter values can be equal or less than zero. \n")
             cat(" Press Enter to continue.         \n")
             cat("**********************************\n\n")
             readline()
             cat("\n")
             par.err<-edit(par.err)}   
           else{
             break
             return(edit(par.err))}
         }
         
         factor1<-par.err[1,2]
         factor2<-par.err[2,2]
         
         PKindex<-vector(Subject,mode="list");cat("\n")
         for(i in 1:Subject)  { 
            cat("\n     << Subject:- #",i,">>\n\n" )
            
            switch(pick-1,
                   {A<-par1+rnorm(1,mean=0,sd=factor1)
                    while(A<=0){
                       A<-par1+rnorm(1,mean=0,sd=factor1)}
                    alpha<-par2+rnorm(1,mean=0,sd=factor2)
                    while(alpha<=0){
                       alpha<-par2+rnorm(1,mean=0,sd=factor2)}},
                       
                    {A<-par1+runif(1,min=-factor1,max=factor1)
                     while(A<=0){
                        A<-par1+runif(1,min=-factor1,max=factor1)}
                     alpha<-par2+runif(1,min=-factor2,max=factor2)
                     while(alpha<=0){
                        alpha<-par2+runif(1,min=-factor2,max=factor2)}},
                        
                    {A<-par1*rnorm(1,mean=0,sd=factor1)+par1
                     while(A<=0){
                        A<-par1*rnorm(1,mean=0,sd=factor1)+par1}
                     alpha<-par2*rnorm(1,mean=0,sd=factor2)+par2
                     while(alpha<=0){
                        alpha<-par2*rnorm(1,mean=0,sd=factor2)+par2}},
                        
                    {A<-par1*runif(1,min=-factor1,max=factor1)+par1
                     while(A<=0){
                        A<-par1*runif(1,min=-factor1,max=factor1)+par1}
                     alpha<-par2*runif(1,min=-factor2,max=factor2)+par2
                     while(alpha<=0){
                        alpha<-par2*runif(1,min=-factor2,max=factor2)+par2}}   
                    )
            PKindex[[i]]<-smacro.one.out(PKtime,A,alpha,defun,par1,par2,i,type)    
            }  
         PKindex<- as.data.frame(do.call("rbind",PKindex))
         rownames(PKindex) <- seq(nrow(PKindex)) 
         savefile(PKindex)    
        }  
      }
      else if (pick ==2){ ### starting monte-carlo sim here
      cat("\n\n")
      file.menu <- c("Error = Normal Error", 
                     "error=uniform error",
                     "error=normal error*true value",
                     "error=uniform error*true value")
      pick <- menu(file.menu, title = "<< Error types >>")
      
      type<-switch(pick, 
                  "Error = Normal Error", 
                  "Error = Uniform Error",
                  "Error = Normal Error*True Value",
                  "Error = Uniform Error*True Value")
      
      cat("\n\nHow many iterations to run for each subject?\n")
      re<-scan(nlines=1,quiet=TRUE)
         par.err<-data.frame(Parameter=c("A","alpha"),error_factor=c(.15,.15))
         par.err<-edit(par.err)
         repeat{
         if (par.err[1,2] <= 0 || par.err[2,2] <=0){
             cat("\n")
             cat("**********************************\n")
             cat(" Parameter values can be equal or less than zero. \n")
             cat(" Press Enter to continue.         \n")
             cat("**********************************\n\n")
             readline()
             cat("\n")
             par.err<-edit(par.err)}   
           else{
             break
             return(edit(par.err))}
         }
         
         factor1<-par.err[1,2]
         factor2<-par.err[2,2]
      
      cat("\n")
      cat("**********************************\n")
      cat(" Summary Table                     \n")
      cat(" Model: one-exponential term model \n") 
      cat(" Subject #:", Subject,"            \n")
      cat(" Error Type:", type,"              \n")
      cat(" Simulation #:", re,"              \n\n")
      sim<-matrix(c(par1,par2,factor1,factor2),2,2)
      dimnames(sim)<-list(c("A","alpha"),c("Selected","Error factor"))
      show(sim)   
      cat("**********************************\n\n")
      PKindex<-vector(Subject,mode="list")
      for( i in 1:Subject)  {
        cat("\n\n             << Subject",i,">>\n\n" ) 
        C1.lsoda<-list()
           for (j in 1:re){
              switch(pick, 
                     {A<-par1+rnorm(1,mean=0,sd=factor1)
                      while(A<=0){
                         A<-par1+rnorm(1,mean=0,sd=factor1)}
                      alpha<-par2+rnorm(1,mean=0,sd=factor2)
                      while(alpha<=0){
                         alpha<-par2+rnorm(1,mean=0,sd=factor2)}},
                         
                     {A<-par1+runif(1,min=-factor1,max=factor1)
                      while(A<=0){
                         A<-par1+runif(1,min=-factor1,max=factor1)} 
                      alpha<-par2+runif(1,min=-factor2,max=factor2)
                      while(alpha<=0){
                      alpha<-par2+runif(1,min=-factor2,max=factor2)}},
                      
                      {A<-par1*rnorm(1,mean=0,sd=factor1)+par1
                       while(A<=0){
                          a<-par1*rnorm(1,mean=0,sd=factor1)+par1}
                       alpha<-par2*rnorm(1,mean=0,sd=factor2)+par2
                       while(alpha<=0){
                          alpha<-par2*rnorm(1,mean=0,sd=factor2)+par2}},
                          
                      {A<-par1*runif(1,min=-factor1,max=factor1)+par1
                       while(A<=0){
                          A<-par1*runif(1,min=-factor1,max=factor1)+par1}
                       alpha<-par2*runif(1,min=-factor2,max=factor2)+par2
                       while(alpha<=0){
                          alpha<-par2*runif(1,min=-factor2,max=factor2)+par2}}
              )
              time1<-PKtime$time
              defun<- A*exp(-alpha*time1) 
              XX<-data.frame(time1,defun) 
              C1.lsoda[[j]]<-data.frame(XX$time1,XX$defun) 
              colnames(C1.lsoda[[j]])<-list("time","concentration")
           }
         PKindex[[i]]<-montecarlo(C1.lsoda,time1,i,re) 
         }
    PKindex<-as.data.frame(do.call("rbind",PKindex))
    rownames(PKindex)<-seq(nrow(PKindex)) 
    savefile(PKindex)  
   } 
}      
