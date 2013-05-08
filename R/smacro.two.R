###Two exponential term
smacro.two <-function(Subject=NULL,  # N Subj's 
                      PKtime=NULL,   # times for sampling
                      A=NULL,
                      a=NULL,
                      B=NULL,
                      b=NULL)
{
   #options(warn=-1)
   
   if (is.null(Subject) || !is.integer(Subject)) {
     cat("How many subjects do you want?\n")
     Subject<-scan(nlines=1,quiet=TRUE)
   }

   if (is.null(PKtime)) { 
     ## need to verify is a numeric ordered vector.
     PKtime<-data.frame(time=c(0))
     PKtime<-edit(PKtime)
   }
  
   cat("\n")
   show(PKtime)
   
   if (is.null(A) || is.null(a) || is.null(B) || is.null(b)) {
       par<-data.frame(Parameter=c("A","a","B","b"),Initial=c(0))
       par<-edit(par)
       repeat{
           if ( par[1,2] == 0 || par[2,2] ==0 || par[3,2]==0 || par[4,2]==0){
             cat("\n")
             cat("**********************************\n")
             cat(" Parameter value can not be zero. \n")
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
       par1<-par[1,2]
       par2<-par[2,2]
       par3<-par[3,2]
       par4<-par[4,2]
    } 
    else {
       par1 <- A
       par2 <- a
       par3 <- B
       par4 <- b
    } 
    
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
      pick <- menu(file.menu, title = "<< Error Types >>")      
      
      type<-switch(pick, 
                  "No Error",
                  "Error = Normal Error", 
                  "Error = Uniform Error",
                  "Error = Normal Error*True Value",
                  "Error = Uniform Error*True Value")
      
      if (pick ==1){
          PKindex<-vector(Subject,mode="list")
          for( i in 1:Subject)  {
            cat("\n\n             << Subject",i,">>\n\n" ) 
            A<-par1
            a<-par2  
            B<-par3
            b<-par4
            PKindex[[i]]<-smacro.two.out(PKtime,A,a,B,b,defun,par1,par2,par3,par4,i,type)
            }  
      PKindex<- as.data.frame(do.call("rbind",PKindex))
      rownames(PKindex) <- seq(nrow(PKindex)) 
      savefile(PKindex)  
      }  
      else {
         cat("\n\nEnter error factor for A\n")
         factor1<-scan(nlines=1,quiet=TRUE)
         repeat{
              if ( factor1 == 0 ){
                cat("\n")
                cat("**********************************\n")
                cat(" Parameter value can not be zero. \n")
                cat(" Press Enter to continue.         \n")
                cat("**********************************\n\n")
                readline()
                cat("\n")
                factor1<-scan(nlines=1,quiet=TRUE)}   
           else{
                break
                return(factor1)}
           }
         
         cat("\nEnter error factor for a\n")
         factor2<-scan(nlines=1,quiet=TRUE)
         repeat{
              if ( factor2 == 0 ){
                cat("\n")
                cat("**********************************\n")
                cat(" Parameter value can not be zero. \n")
                cat(" Press Enter to continue.         \n")
                cat("**********************************\n\n")
                readline()
                cat("\n")
                factor2<-scan(nlines=1,quiet=TRUE)}   
           else{
                break
                return(factor2)}
           }
         
         cat("\nEnter error factor for B\n")
         factor3<-scan(nlines=1,quiet=TRUE)
         repeat{
              if ( factor3 == 0 ){
                cat("\n")
                cat("**********************************\n")
                cat(" Parameter value can not be zero. \n")
                cat(" Press Enter to continue.         \n")
                cat("**********************************\n\n")
                readline()
                cat("\n")
                factor3<-scan(nlines=1,quiet=TRUE)}   
           else{
                break
                return(factor3)}
           }
         
         cat("\nEnter error factor for b\n")
         factor4<-scan(nlines=1,quiet=TRUE)
         repeat{
              if ( factor4 == 0 ){
                cat("\n")
                cat("**********************************\n")
                cat(" Parameter value can not be zero. \n")
                cat(" Press Enter to continue.         \n")
                cat("**********************************\n\n")
                readline()
                cat("\n")
                factor4<-scan(nlines=1,quiet=TRUE)}   
           else{
                break
                return(factor4)}
           }
         
         PKindex<-vector(Subject,mode="list")
         for( i in 1:Subject)  { 
            cat("\n\n             << Subject",i,">>\n\n" )  

            switch(pick-1,
                   {A<-par1+rnorm(1,mean=0,sd=factor1)
                    while(A<=0){
                       A<-par1+rnorm(1,mean=0,sd=factor1)}
                    a<-par2+rnorm(1,mean=0,sd=factor2)
                    while(a<=0){
                       a<-par2+rnorm(1,mean=0,sd=factor2)}
                    B<-par3+rnorm(1,mean=0,sd=factor3)
                    while(B<=0){
                       B<-par3+rnorm(1,mean=0,sd=factor3)}
                    b<-par4+rnorm(1,mean=0,sd=factor4)
                    while(b<=0){
                       b<-par4+rnorm(1,mean=0,sd=factor4)}},
                       
                    {A<-par1+runif(1,min=-factor1,max=factor1)
                     while(A<=0){
                        A<-par1+runif(1,min=-factor1,max=factor1)}
                     a<-par2+runif(1,min=-factor2,max=factor2)
                     while(a<=0){
                        a<-par2+runif(1,min=-factor2,max=factor2)}
                     B<-par3+runif(1,min=-factor3,max=factor3)
                     while(B<=0){
                        B<-par1+runif(1,min=-factor3,max=factor3)}
                     b<-par4+runif(1,min=-factor4,max=factor4)
                     while(b<=0){
                        b<-par4+runif(1,min=-factor4,max=factor4)}},
                        
                    {A<-par1*rnorm(1,mean=0,sd=factor1)+par1
                     while(A<=0){
                        A<-par1*rnorm(1,mean=0,sd=factor1)+par1}
                     a<-par2*rnorm(1,mean=0,sd=factor2)+par2
                     while(a<=0){
                        a<-par2*rnorm(1,mean=0,sd=factor2)+par2}
                     B<-par3*rnorm(1,mean=0,sd=factor3)+par3
                     while(B<=0){
                        B<-par3*rnorm(1,mean=0,sd=factor3)+par3}
                     b<-par4*rnorm(1,mean=0,sd=factor4)+par4
                     while(b<=0){
                        b<-par4*rnorm(1,mean=0,sd=factor4)+par4}},
                        
                    {A<-par1*runif(1,min=-factor1,max=factor1)+par1
                     while(A<=0){
                        A<-par1*runif(1,min=-factor1,max=factor1)+par1}
                     a<-par2*runif(1,min=-factor2,max=factor2)+par2
                     while(a<=0){
                        a<-par2*runif(1,min=-factor2,max=factor2)+par2}
                     B<-par3*runif(1,min=-factor3,max=factor3)+par3
                     while(B<=0){
                        B<-par3*runif(1,min=-factor3,max=factor3)+par3}
                     b<-par4*runif(1,min=-factor4,max=factor4)+par4
                     while(b<=0){
                        b<-par4*runif(1,min=-factor4,max=factor4)+par4}}   
                    )
            PKindex[[i]]<-smacro.two.out(PKtime,A,a,B,b,defun,par1,par2,par3,par4,i,type)          
            }   
         PKindex<- as.data.frame(do.call("rbind",PKindex))
         rownames(PKindex) <- seq(nrow(PKindex)) 
         savefile(PKindex)     
        } 
      }
      else if (pick ==2){ 
      cat("\n\n")
      file.menu <- c("Error = Normal Error", 
                     "Error = Uniform Error",
                     "Error = Normal Error*True Value",
                     "Error = Uniform Error*True Value")
      pick <- menu(file.menu, title = "<< Error Types >>")
      
      type<-switch(pick, 
                  "Error = Normal Error", 
                  "Error = Uniform Error",
                  "Error = Normal Error*True Value",
                  "Error = Uniform Error*True Value")
      
      cat("\n\nHow many interations do you want to run?\n")
      re<-scan(nlines=1,quiet=TRUE)
      
      cat("\nEnter error factor for A\n")
      factor1<-scan(nlines=1,quiet=TRUE)
      repeat{
              if ( factor1 == 0 ){
                cat("\n")
                cat("**********************************\n")
                cat(" Parameter value can not be zero. \n")
                cat(" Press Enter to continue.         \n")
                cat("**********************************\n\n")
                readline()
                cat("\n")
                factor1<-scan(nlines=1,quiet=TRUE)}   
           else{
                break
                return(factor1)}
           }
      
      cat("\nEnter error factor for a\n")
      factor2<-scan(nlines=1,quiet=TRUE)
      repeat{
              if ( factor2 == 0 ){
                cat("\n")
                cat("**********************************\n")
                cat(" Parameter value can not be zero. \n")
                cat(" Press Enter to continue.         \n")
                cat("**********************************\n\n")
                readline()
                cat("\n")
                factor2<-scan(nlines=1,quiet=TRUE)}   
           else{
                break
                return(factor2)}
           }
      
      cat("\nEnter error factor for B\n")
      factor3<-scan(nlines=1,quiet=TRUE)
      repeat{
              if ( factor3 == 0 ){
                cat("\n")
                cat("**********************************\n")
                cat(" Parameter value can not be zero. \n")
                cat(" Press Enter to continue.         \n")
                cat("**********************************\n\n")
                readline()
                cat("\n")
                factor3<-scan(nlines=1,quiet=TRUE)}   
           else{
                break
                return(factor3)}
           }
      
      cat("\nEnter error factor for b\n")
      factor4<-scan(nlines=1,quiet=TRUE)
      repeat{
              if ( factor4 == 0 ){
                cat("\n")
                cat("**********************************\n")
                cat(" Parameter value can not be zero. \n")
                cat(" Press Enter to continue.         \n")
                cat("**********************************\n\n")
                readline()
                cat("\n")
                factor4<-scan(nlines=1,quiet=TRUE)}   
           else{
                break
                return(factor4)}
           }
      
      cat("\n")
      cat("**********************************\n")
      cat("Summary Table                     \n")
      cat("Model: Two-exponential Term Model \n") 
      cat("Subject #:", Subject,"            \n")
      cat("Error Type:", type,"              \n")
      cat("Simulation #:", re,"              \n\n")
      sim<-matrix(c(par1,par2,par3,par4,factor1,factor2,factor3,factor4),4,2)
      dimnames(sim)<-list(c("A","a","B","b"),c("Selected","Error factor"))
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
                    a<-par2+rnorm(1,mean=0,sd=factor2)
                    while(a<=0){
                       a<-par2+rnorm(1,mean=0,sd=factor2)}
                    B<-par3+rnorm(1,mean=0,sd=factor3)
                    while(B<=0){
                       B<-par3+rnorm(1,mean=0,sd=factor3)}
                    b<-par4+rnorm(1,mean=0,sd=factor4)
                    while(b<=0){
                       b<-par4+rnorm(1,mean=0,sd=factor4)}},
                       
                    {A<-par1+runif(1,min=-factor1,max=factor1)
                     while(A<=0){
                        A<-par1+runif(1,min=-factor1,max=factor1)}
                     a<-par2+runif(1,min=-factor2,max=factor2)
                     while(a<=0){
                        a<-par2+runif(1,min=-factor2,max=factor2)}
                     B<-par3+runif(1,min=-factor3,max=factor3)
                     while(B<=0){
                        B<-par1+runif(1,min=-factor3,max=factor3)}
                     b<-par4+runif(1,min=-factor4,max=factor4)
                     while(b<=0){
                        b<-par4+runif(1,min=-factor4,max=factor4)}},
                        
                    {A<-par1*rnorm(1,mean=0,sd=factor1)+par1
                     while(A<=0){
                        A<-par1*rnorm(1,mean=0,sd=factor1)+par1}
                     a<-par2*rnorm(1,mean=0,sd=factor2)+par2
                     while(a<=0){
                        a<-par2*rnorm(1,mean=0,sd=factor2)+par2}
                     B<-par3*rnorm(1,mean=0,sd=factor3)+par3
                     while(B<=0){
                        B<-par3*rnorm(1,mean=0,sd=factor3)+par3}
                     b<-par4*rnorm(1,mean=0,sd=factor4)+par4
                     while(b<=0){
                        b<-par4*rnorm(1,mean=0,sd=factor4)+par4}},
                        
                    {A<-par1*runif(1,min=-factor1,max=factor1)+par1
                     while(A<=0){
                        A<-par1*runif(1,min=-factor1,max=factor1)+par1}
                     a<-par2*runif(1,min=-factor2,max=factor2)+par2
                     while(a<=0){
                        a<-par2*runif(1,min=-factor2,max=factor2)+par2}
                     B<-par3*runif(1,min=-factor3,max=factor3)+par3
                     while(B<=0){
                        B<-par3*runif(1,min=-factor3,max=factor3)+par3}
                     b<-par4*runif(1,min=-factor4,max=factor4)+par4
                     while(b<=0){
                        b<-par4*runif(1,min=-factor4,max=factor4)+par4}}   
              )
              time1<-PKtime$time
              defun<- A*exp(-a*time1)+B*exp(-b*time1) 
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