###Two compartment PK model extravascular single dose first order absorption
sfirst2<- function(Subject=NULL,   # N Subj's 
                   PKtime=NULL,    # times for sampling
                   Dose=NULL,      # single dose
                   ka=NULL,
                   Vd=NULL,
                   kel=NULL, 
                   k12=NULL,
                   k21=NULL)
{
   options(warn=-1)
   
   cat("\n First enter all paramaters used for simulation profile.\n");readline(" Press Enter to continue...\n\n")
   par<-data.frame(Parameter=c("Total_subject#","Dose","ka","kel","k12","k21","Vd"),Initial=c(24,300,0.51,0.36,0.13,0.31,11.7))
   par<-edit(par)                                                       
   repeat{
       if (par[1,2] <= 0 || par[2,2] <=0 || par[3,2]<=0 || par[4,2]<=0 || par[5,2]<=0|| par[6,2]<=0|| par[7,2]<=0){
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
   Dose<-par[2,2]
   par1<-par[3,2]
   par2<-par[4,2]
   par3<-par[5,2]
   par4<-par[6,2]
   par5<-par[7,2]

   readline("\n Next enter time points. Press Enter to continue.\n\n")
   PKtime<-data.frame(time=c(0))
   PKtime<-edit(PKtime)
   cat("\n")
   show(PKtime);cat("\n\n")   
    
    defun<- function(time, y, parms) { 
       dCp1dt <- -parms["ka"]*y[1]    ### y[1] here is amount of drug, not drug plasma conc. here -YJ
       dCp2dt <- parms["ka"]*y[1]/parms["Vd"]-parms["kel"]*y[2]+parms["k21"]*y[3]-parms["k12"]*y[2]
       dCp3dt <- parms["k12"]*y[2]-parms["k21"]*y[3]
       list(c(dCp1dt,dCp2dt,dCp3dt)) 
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
         PKindex<-vector(Subject,mode="list");cat("\n")
         for(i in 1:Subject)  {
           cat("\n     << Subject:- #",i,">>\n\n" )
           ka <-par1
           kel<-par2
           k12<-par3
           k21<-par4
           Vd <-par5   
           PKindex[[i]]<-sfirst2.out(PKtime,ka,kel,k12,k21,Vd,defun,par1,par2,par3,par4,par5,Dose,i,type)  
         }   
         PKindex<- as.data.frame(do.call("rbind",PKindex))
         rownames(PKindex) <- seq(nrow(PKindex)) 
         savefile(PKindex)    
       }   
       else {
         par.err<-data.frame(Parameter=c("ka","kel","k12","k21","Vd"),error_factor=c(.15,.15,.15,.15,.15))
         par.err<-edit(par.err)
         repeat{
         if (par.err[1,2] <= 0 || par.err[2,2] <=0 || par.err[3,2]<=0 || par.err[4,2]<=0 || par.err[5,2]<=0){
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
         factor3<-par.err[3,2]
         factor4<-par.err[4,2]
         factor5<-par.err[5,2]
         
         PKindex<-vector(Subject,mode="list");cat("\n")
         for(i in 1:Subject)  {
             cat("\n     << Subject:- #",i,">>\n\n" )  
             switch(pick-1,
                    {ka<-par1+rnorm(1,mean=0,sd=factor1)
                     while(ka<=0){
                        ka<-par1+rnorm(1,mean=0,sd=factor1)}          
                     kel<-par2+rnorm(1,mean=0,sd=factor2)
                     while(kel<=0){
                        kel<-par2+rnorm(1,mean=0,sd=factor2)} 
                     k12<-par3+rnorm(1,mean=0,sd=factor3)
                     while(k12<=0){
                        k12<-par3+rnorm(1,mean=0,sd=factor3)}
                     k21<-par4+rnorm(1,mean=0,sd=factor4)
                     while(k21<=0){
                        k21<-par4+rnorm(1,mean=0,sd=factor4)}
                     Vd<-par5+rnorm(1,mean=0,sd=factor5)
                     while(Vd<=0){
                        Vd<-par5+rnorm(1,mean=0,sd=factor5)}},
                        
                    {ka<-par1+runif(1,min=-factor1,max=factor1)
                     while(ka<=0){
                        ka<-par1+runif(1,min=-factor1,max=factor1)}      
                     kel<-par2+runif(1,min=-factor2,max=factor2)
                     while(kel<=0){
                        kel<-par2+runif(1,min=-factor2,max=factor2)}
                     k12<-par3+runif(1,min=-factor3,max=factor3)
                     while(k12<=0){
                        k12<-par3+runif(1,min=-factor3,max=factor3)}
                     k21<-par4+runif(1,min=-factor4,max=factor4)
                     while(k21<=0){
                        k21<-par4+runif(1,min=-factor4,max=factor4)}
                     Vd<-par5+runif(1,min=-factor5,max=factor5)
                     while(Vd<=0){
                        Vd<-par5+runif(1,min=-factor5,max=factor5)}},
                        
                     {ka<-par1*rnorm(1,mean=0,sd=factor1)+par1
                      while(ka<=0){
                         ka<-par1*rnorm(1,mean=0,sd=factor1)+par1}         
                      kel<-par2*rnorm(1,mean=0,sd=factor2)+par2
                      while(kel<=0){
                         kel<-par2*rnorm(1,mean=0,sd=factor2)+par2}
                      k12<-par3*rnorm(1,mean=0,sd=factor3)+par3
                      while(k12<=0){
                         k12<-par3*rnorm(1,mean=0,sd=factor3)+par3}
                      k21<-par4*rnorm(1,mean=0,sd=factor4)+par4
                      while(k21<=0){
                         k21<-par4*rnorm(1,mean=0,sd=factor4)+par4}
                      Vd<-par5*rnorm(1,mean=0,sd=factor5)+par5
                      while(Vd<=0){
                         Vd<-par5*rnorm(1,mean=0,sd=factor5)+par5}},
                         
                      {ka<-par1*runif(1,min=-factor1,max=factor1)+par1
                       while(ka<=0){
                          ka<-par1*runif(1,min=-factor1,max=factor1)+par1}
                       kel<-par2*runif(1,min=-factor2,max=factor2)+par2
                       while(kel<=0){
                          kel<-par2*runif(1,min=-factor2,max=factor2)+par2}
                       k12<-par3*runif(1,min=-factor3,max=factor3)+par3
                       while(k12<=0){
                          k12<-par3*runif(1,min=-factor3,max=factor3)+par3}
                       k21<-par4*runif(1,min=-factor4,max=factor4)+par4
                       while(k21<=0){
                          k21<-par4*runif(1,min=-factor4,max=factor4)+par4}
                       Vd<-par5*runif(1,min=-factor5,max=factor5)+par5
                       while(Vd<=0){
                          Vd<-par5*runif(1,min=-factor5,max=factor5)+par5}}
               )
               PKindex[[i]]<-sfirst2.out(PKtime,ka,kel,k12,k21,Vd,defun,par1,par2,par3,par4,par5,Dose,i,type) 
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
     
     cat("\n\nHow many iterations to run for each subject?\n")
     re<-scan(nlines=1,quiet=TRUE)
     
     par.err<-data.frame(Parameter=c("ka","kel","k12","k21","Vd"),error_factor=c(.15,.15,.15,.15,.15))
     par.err<-edit(par.err)
     repeat{
     if (par.err[1,2] <= 0 || par.err[2,2] <=0 || par.err[3,2]<=0 || par.err[4,2]<=0 || par.err[5,2]<=0){
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
     factor3<-par.err[3,2]
     factor4<-par.err[4,2]
     factor5<-par.err[5,2]
     
     cat("\n")
     cat("******************************************************\n")
     cat(" Summary Table                                         \n")
     cat(" Model: a 2-compartment, extravascular,                  \n") 
     cat("        single-dose & 1st-ordered without lag time model \n") 
     cat(" Subject #:", Subject,"                                \n")
     cat(" Error Type:", type,"                                  \n")
     cat(" Simulation #:", re,"                                  \n\n")
     sim<-matrix(c(par1,par2,par3,par4,par5,factor1,factor2,factor3,factor4,factor5),5,2)
     dimnames(sim)<-list(c("ka","kel","k12","k21","Vd"),c("Selected","Error factor"))
     show(sim)   
     cat("******************************************************\n\n")
     
     PKindex<-vector(Subject,mode="list");cat("\n")
     for(i in 1:Subject)  {
     cat("\n     << Subject:- #",i,">>\n\n" )
     C1.lsoda<-list()
        for (j in 1:re){
            switch(pick,
                   {ka<-par1+rnorm(1,mean=0,sd=factor1)
                     while(ka<=0){
                        ka<-par1+rnorm(1,mean=0,sd=factor1)}          
                     kel<-par2+rnorm(1,mean=0,sd=factor2)
                     while(kel<=0){
                        kel<-par2+rnorm(1,mean=0,sd=factor2)} 
                     k12<-par3+rnorm(1,mean=0,sd=factor3)
                     while(k12<=0){
                        k12<-par3+rnorm(1,mean=0,sd=factor3)}
                     k21<-par4+rnorm(1,mean=0,sd=factor4)
                     while(k21<=0){
                        k21<-par4+rnorm(1,mean=0,sd=factor4)}
                     Vd<-par5+rnorm(1,mean=0,sd=factor5)
                     while(Vd<=0){
                        Vd<-par5+rnorm(1,mean=0,sd=factor5)}},
                        
                    {ka<-par1+runif(1,min=-factor1,max=factor1)
                     while(ka<=0){
                        ka<-par1+runif(1,min=-factor1,max=factor1)}      
                     kel<-par2+runif(1,min=-factor2,max=factor2)
                     while(kel<=0){
                        kel<-par2+runif(1,min=-factor2,max=factor2)}
                     k12<-par3+runif(1,min=-factor3,max=factor3)
                     while(k12<=0){
                        k12<-par3+runif(1,min=-factor3,max=factor3)}
                     k21<-par4+runif(1,min=-factor4,max=factor4)
                     while(k21<=0){
                        k21<-par4+runif(1,min=-factor4,max=factor4)}
                     Vd<-par5+runif(1,min=-factor5,max=factor5)
                     while(Vd<=0){
                        Vd<-par5+runif(1,min=-factor5,max=factor5)}},
                        
                     {ka<-par1*rnorm(1,mean=0,sd=factor1)+par1
                      while(ka<=0){
                         ka<-par1*rnorm(1,mean=0,sd=factor1)+par1}         
                      kel<-par2*rnorm(1,mean=0,sd=factor2)+par2
                      while(kel<=0){
                         kel<-par2*rnorm(1,mean=0,sd=factor2)+par2}
                      k12<-par3*rnorm(1,mean=0,sd=factor3)+par3
                      while(k12<=0){
                         k12<-par3*rnorm(1,mean=0,sd=factor3)+par3}
                      k21<-par4*rnorm(1,mean=0,sd=factor4)+par4
                      while(k21<=0){
                         k21<-par4*rnorm(1,mean=0,sd=factor4)+par4}
                      Vd<-par5*rnorm(1,mean=0,sd=factor5)+par5
                      while(Vd<=0){
                         Vd<-par5*rnorm(1,mean=0,sd=factor5)+par5}},
                         
                      {ka<-par1*runif(1,min=-factor1,max=factor1)+par1
                       while(ka<=0){
                          ka<-par1*runif(1,min=-factor1,max=factor1)+par1}
                       kel<-par2*runif(1,min=-factor2,max=factor2)+par2
                       while(kel<=0){
                          kel<-par2*runif(1,min=-factor2,max=factor2)+par2}
                       k12<-par3*runif(1,min=-factor3,max=factor3)+par3
                       while(k12<=0){
                          k12<-par3*runif(1,min=-factor3,max=factor3)+par3}
                       k21<-par4*runif(1,min=-factor4,max=factor4)+par4
                       while(k21<=0){
                          k21<-par4*runif(1,min=-factor4,max=factor4)+par4}
                       Vd<-par5*runif(1,min=-factor5,max=factor5)+par5
                       while(Vd<=0){
                          Vd<-par5*runif(1,min=-factor5,max=factor5)+par5}}
               )
               time1<-PKtime$time
               parms<-c(ka=ka,kel=kel,k12=k12,k21=k21,Vd=Vd)
               XX<-data.frame(lsoda(c(Dose,0,0),c(0,time1),defun,parms,rtol=1e-6,atol=1e-6))
               C1.lsoda[[j]]<-data.frame(XX[2:(length(time1)+1),1],XX[2:(length(time1)+1),3])
               colnames(C1.lsoda[[j]])<-list("time","concentration") 
          }      
    PKindex[[i]]<-montecarlo(C1.lsoda,time1,i,re) 
    }  
    PKindex<-as.data.frame(do.call("rbind",PKindex))
    rownames(PKindex)<-seq(nrow(PKindex)) 
    savefile(PKindex)
  }  
}   
