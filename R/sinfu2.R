###Two compartment PK model iv infusion single dose 
sinfu2<- function(Subject=NULL,  # N Subj's 
                  PKtime=NULL,   # times for sampling
                  Dose=NULL,     # single dose
                  Tinf=NULL,
                  Vd=NULL,
                  kel=NULL, 
                  k12=NULL,
                  k21=NULL)
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
   
   if (is.null(Dose)) {
    cat("\nEnter dose\n")
    Dose<-scan(nlines=1,quiet=TRUE)
   }
   
   if (is.null(Tinf)) {
    cat("\nEnter value for infusion time\n")
    Tinf<-scan(nlines=1,quiet=TRUE)
    cat("\n")
  }
   
   if (is.null(kel) || is.null(k12) || is.null(k21) || is.null(Vd)) {
       par<-data.frame(Parameter=c("kel","k12","k21","Vd"),Initial=c(0))
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
       par1 <- kel
       par2 <- k12
       par3 <- k21
       par4 <- Vd
    }
    
    defun<- function(time, y, parms) { 
       if(time<=Tinf) {
           dCp1dt<- (Dose/Tinf)/parms["Vd"]+parms["k21"]*y[2]-parms["kel"]*y[1]-parms["k12"]*y[1]
           dCp2dt<- parms["k12"]*y[1]-parms["k21"]*y[2]
         }
       else {       
           dCp1dt <- -parms["kel"]*y[1]-parms["k12"]*y[1]+parms["k21"]*y[2] 
           dCp2dt <- parms["k12"]*y[1]-parms["k21"]*y[2]
         }
       list(c(dCp1dt,dCp2dt)) 
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
           kel<-par1
           k12<-par2
           k21<-par3
           Vd<-par4   
           PKindex[[i]]<-sinfu2.out(PKtime,kel,k12,k21,Vd,defun,par1,par2,par3,par4,Dose,i,type)  
         }   
         PKindex<- as.data.frame(do.call("rbind",PKindex))
         rownames(PKindex) <- seq(nrow(PKindex)) 
         savefile(PKindex)    
       }   
       else {
         cat("\n\nEnter error factor for kel\n")
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
         
         cat("\nEnter error factor for k12\n")
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
         
         cat("\nEnter error factor for k21\n")
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
         
         cat("\nEnter error factor for Vd\n")
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
                    {kel<-par1+rnorm(1,mean=0,sd=factor1)
                     while(kel<=0){
                        kel<-par1+rnorm(1,mean=0,sd=factor1)} 
                     k12<-par2+rnorm(1,mean=0,sd=factor2)
                     while(k12<=0){
                        k12<-par2+rnorm(1,mean=0,sd=factor2)}
                     k21<-par3+rnorm(1,mean=0,sd=factor3)
                     while(k21<=0){
                        k21<-par3+rnorm(1,mean=0,sd=factor3)}        
                     Vd<-par4+rnorm(1,mean=0,sd=factor4)
                     while(Vd<=0){
                        Vd<-par4+rnorm(1,mean=0,sd=factor4)}},
                        
                    {kel<-par1+runif(1,min=-factor1,max=factor1)
                     while(kel<=0){
                        kel<-par1+runif(1,min=-factor1,max=factor1)}
                     k12<-par2+runif(1,min=-factor2,max=factor2)
                     while(k12<=0){
                        k12<-par2+runif(1,min=-factor2,max=factor2)}
                     k21<-par3+runif(1,min=-factor3,max=factor3)
                     while(k21<=0){
                        k21<-par3+runif(1,min=-factor3,max=factor3)}
                     Vd<-par4+runif(1,min=-factor4,max=factor4)
                     while(Vd<=0){
                        Vd<-par4+runif(1,min=-factor4,max=factor4)}},
                        
                     {kel<-par1*rnorm(1,mean=0,sd=factor1)+par1
                      while(kel<=0){
                         kel<-par1*rnorm(1,mean=0,sd=factor1)+par1}
                      k12<-par2*rnorm(1,mean=0,sd=factor2)+par2
                      while(k12<=0){
                         k12<-par2*rnorm(1,mean=0,sd=factor2)+par2}
                      k21<-par3*rnorm(1,mean=0,sd=factor3)+par3
                      while(k21<=0){
                         k21<-par3*rnorm(1,mean=0,sd=factor3)+par3}
                      Vd<-par4*rnorm(1,mean=0,sd=factor4)+par4
                      while(Vd<=0){
                         Vd<-par4*rnorm(1,mean=0,sd=factor4)+par4}},
                         
                      {kel<-par1*runif(1,min=-factor1,max=factor1)+par1
                       while(kel<=0){
                          kel<-par1*runif(1,min=-factor1,max=factor1)+par1}
                       k12<-par2*runif(1,min=-factor2,max=factor2)+par2
                       while(k12<=0){
                          k12<-par2*runif(1,min=-factor2,max=factor2)+par2}
                       k21<-par3*runif(1,min=-factor3,max=factor3)+par3
                       while(k21<=0){
                       k21<-par3*runif(1,min=-factor3,max=factor3)+par3}
                       Vd<-par4*runif(1,min=-factor4,max=factor4)+par4
                       while(Vd<=0){
                          Vd<-par4*runif(1,min=-factor4,max=factor4)+par4}}
               )
               PKindex[[i]]<-sinfu2.out(PKtime,kel,k12,k21,Vd,defun,par1,par2,par3,par4,Dose,i,type)           
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
     
     cat("\nEnter error factor for kel\n")
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
     
     cat("\nEnter error factor for k12\n")
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
     
     cat("\nEnter error factor for k21\n")
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
     
     cat("\nEnter error factor for Vd\n")
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
     cat("*******************************************************\n")
     cat("Summary Table                                          \n")
     cat("Model: 2-Compartment, IV-Infusion, & Single-Dose Model \n") 
     cat("Subject #:", Subject,"                                 \n")
     cat("Error Type:", type,"                                   \n")
     cat("Simulation #:", re,"                                   \n\n")
     sim<-matrix(c(par1,par2,par3,par4,factor1,factor2,factor3,factor4),4,2)
     dimnames(sim)<-list(c("kel","k12","k21","Vd"),c("Selected","Error factor"))
     show(sim)   
     cat("*******************************************************\n\n")
     
     PKindex<-vector(Subject,mode="list")
     for( i in 1:Subject)  {
     cat("\n\n             << Subject",i,">>\n\n" ) 
     C1.lsoda<-list()
        for (j in 1:re){
            switch(pick,
                   {kel<-par1+rnorm(1,mean=0,sd=factor1)
                     while(kel<=0){
                        kel<-par1+rnorm(1,mean=0,sd=factor1)} 
                     k12<-par2+rnorm(1,mean=0,sd=factor2)
                     while(k12<=0){
                        k12<-par2+rnorm(1,mean=0,sd=factor2)}
                     k21<-par3+rnorm(1,mean=0,sd=factor3)
                     while(k21<=0){
                        k21<-par3+rnorm(1,mean=0,sd=factor3)}        
                     Vd<-par4+rnorm(1,mean=0,sd=factor4)
                     while(Vd<=0){
                        Vd<-par4+rnorm(1,mean=0,sd=factor4)}},
                        
                    {kel<-par1+runif(1,min=-factor1,max=factor1)
                     while(kel<=0){
                        kel<-par1+runif(1,min=-factor1,max=factor1)}
                     k12<-par2+runif(1,min=-factor2,max=factor2)
                     while(k12<=0){
                        k12<-par2+runif(1,min=-factor2,max=factor2)}
                     k21<-par3+runif(1,min=-factor3,max=factor3)
                     while(k21<=0){
                        k21<-par3+runif(1,min=-factor3,max=factor3)}
                     Vd<-par4+runif(1,min=-factor4,max=factor4)
                     while(Vd<=0){
                        Vd<-par4+runif(1,min=-factor4,max=factor4)}},
                        
                     {kel<-par1*rnorm(1,mean=0,sd=factor1)+par1
                      while(kel<=0){
                         kel<-par1*rnorm(1,mean=0,sd=factor1)+par1}
                      k12<-par2*rnorm(1,mean=0,sd=factor2)+par2
                      while(k12<=0){
                         k12<-par2*rnorm(1,mean=0,sd=factor2)+par2}
                      k21<-par3*rnorm(1,mean=0,sd=factor3)+par3
                      while(k21<=0){
                         k21<-par3*rnorm(1,mean=0,sd=factor3)+par3}
                      Vd<-par4*rnorm(1,mean=0,sd=factor4)+par4
                      while(Vd<=0){
                         Vd<-par4*rnorm(1,mean=0,sd=factor4)+par4}},
                         
                      {kel<-par1*runif(1,min=-factor1,max=factor1)+par1
                       while(kel<=0){
                          kel<-par1*runif(1,min=-factor1,max=factor1)+par1}
                       k12<-par2*runif(1,min=-factor2,max=factor2)+par2
                       while(k12<=0){
                          k12<-par2*runif(1,min=-factor2,max=factor2)+par2}
                       k21<-par3*runif(1,min=-factor3,max=factor3)+par3
                       while(k21<=0){
                       k21<-par3*runif(1,min=-factor3,max=factor3)+par3}
                       Vd<-par4*runif(1,min=-factor4,max=factor4)+par4
                       while(Vd<=0){
                          Vd<-par4*runif(1,min=-factor4,max=factor4)+par4}}
               )
               time1<-PKtime$time
               parms<-c(kel=kel,k12=k12,k21=k21,Vd=Vd)  
               XX<-data.frame(lsoda(c(0,0),c(0,time1),defun,parms,rtol=1e-6,atol=1e-6))
               C1.lsoda[[j]]<-data.frame(XX[2:(length(time1)+1),1],XX[2:(length(time1)+1),2])
               colnames(C1.lsoda[[j]])<-list("time","concentration") 
          }   
      PKindex[[i]]<-montecarlo(C1.lsoda,time1,i,re) 
     }  
    PKindex<-as.data.frame(do.call("rbind",PKindex))
    rownames(PKindex)<-seq(nrow(PKindex)) 
    savefile(PKindex)
  }  
}   
