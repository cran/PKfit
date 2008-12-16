### PKindex is usually the Dataset.


### Normal fitting
### One exponential term

fmacro.one <- function(PKindex,
                       A=NULL,
                       a=NULL) 
{
   #options(warn=-1)
        
   ## Input and initial value for A and a
   
   if (is.null(A) || is.null(a) ) {
       par<-data.frame(Parameter=c("A","a"),Initial=c(0))
       par<-edit(par)
       repeat{
           if ( par[1,2] == 0 || par[2,2] ==0){
             cat("\n")
             cat("**********************************\n")
             cat(" Parameter value can not be zero. \n")
             cat(" Press enter to continue.         \n")
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
   }
   
   cat("\n")
   
   defun<-function(time,A,a){
     pred<- A*exp(-a*time)
   }
   
   ## Select weighting schemes
   file.menu <-c("equal weight",
                 "1/Cp",
                 "1/Cp^2")
   pick <- menu(file.menu, title="<< weighting schemes >>")

   with(entertitle(),{
   
   for( i in 1:length(unique(PKindex$Subject)))  {
      cat("\n\n               << Subject",i,">>\n\n" )  
      objfun <- function(par) {
         out<-defun(PKindex$time[PKindex$Subject==i],par[1],par[2])
         gift<- which( PKindex$conc[PKindex$Subject==i] != 0 )
         switch(pick,
                sum((PKindex$conc[PKindex$Subject==i][gift]-out[gift])^2),
                sum((PKindex$conc[PKindex$Subject==i][gift]-out[gift])^2/PKindex$conc[gift]),
                sum(((PKindex$conc[PKindex$Subject==i][gift]-out[gift])/PKindex$conc[gift])^2)
                )
      }
      
      gen<-genoud(objfun,nvars=2,max=FALSE,pop.size=30,max.generations=20,
           wait.generations=10,starting.value=c(par[1,2],par[2,2]),
           BFGS=FALSE,print.level=0,boundary.enforcement=2,
           Domains=matrix(c(0,0,100,10),2,2),MemoryMatrix=TRUE)
      cat("<< The value of parameter obtained from genetic algorithm >>\n\n")
      namegen<-c("A","a")
      outgen<-c(gen$par[1],gen$par[2])
      print(data.frame(Parameter=namegen,Value=outgen)) 
      F<-objfun(gen$par)
      
      opt<-optim(c(gen$par[1],gen$par[2]),objfun,method="Nelder-Mead")       
      nameopt<-c("A","a")
      outopt<-c(opt$par[1],opt$par[2])
      cat("\n<< The value of parameter fitted by Nelder-Mead Simplex algorithm >>\n\n")
      print(data.frame(Parameter=nameopt,Value=outopt))
      cat("\n<< Residual sum-of-squares and parameter values fitted by nls >>\n\n")
      fm<-nls(conc~defun(time,A,a),data=PKindex, algorithm="default",subset=Subject==i,
          start=list(A=opt$par[1],a=opt$par[2]),trace=TRUE,
          nls.control(maxiter=500,tol=1e-2,minFactor=1/1024))
      cat("\n")   
      coef<-data.frame(coef(fm)["a"])      
      plotting.lin(PKindex, fm, i, pick, coef, xaxis, yaxis)
      }
     })
     cat("\n") 
}      
      
      
### Two exponential term
fmacro.two<- function(PKindex,
                      A=NULL, 
                      a=NULL,
                      B=NULL,
                      b=NULL) 
{
   #options(warn=-1)
        
   ## Input initial value for A, a, B and b
   
   if (is.null(A) || is.null(a) || is.null(B) || is.null(b) ) {
       par<-data.frame(Parameter=c("A","a","B","b"),Initial=c(0))
       par<-edit(par)
       repeat{
           if ( par[1,2] == 0 || par[2,2] ==0 || par[3,2]==0 || par[4,2]==0){
             cat("\n")
             cat("**********************************\n")
             cat(" Parameter value can not be zero. \n")
             cat(" Press enter to continue.         \n")
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
      }
      
    cat("\n")
    
    defun<-function(time,A,a,B,b){
      pred<-A*exp(-a*time)+B*exp(-b*time)
    }
    
    ## Select weighting schemes
    file.menu <-c("equal weight",
                 "1/Cp",
                 "1/Cp^2")
    pick <- menu(file.menu, title="<< weighting schemes >>")

    with(entertitle(),{
    
    for( i in 1:length(unique(PKindex$Subject)))  {
      cat("\n\n               << Subject",i,">>\n\n" )  
      objfun<-function(par) {
         out<-defun(PKindex$time[PKindex$Subject==i],par[1],par[2],par[3],par[4])
         gift<- which( PKindex$conc[PKindex$Subject==i] != 0 )
         switch(pick,
                sum((PKindex$conc[PKindex$Subject==i][gift]-out[gift])^2),
                sum((PKindex$conc[PKindex$Subject==i][gift]-out[gift])^2/PKindex$conc[gift]),
                sum(((PKindex$conc[PKindex$Subject==i][gift] - out[gift])/PKindex$conc[gift])^2)
                )
      }
      
      gen<-genoud(objfun,nvars=4,max=FALSE,pop.size=30,max.generations=20,
           wait.generations=10,
           starting.value=c(par[1,2],par[2,2],par[3,2],par[4,2]),
           BFGS=FALSE,print.level=0,boundary.enforcement=2,
           Domains=matrix(c(1,0.01,0.1,0.01,100,10,50,1),4,2),
           MemoryMatrix=TRUE)     
           
      cat("<< The value of parameter obtained from genetic algorithm >>\n\n")   
      namegen<-c("A","a","B","b")
      outgen<-c(gen$par[1],gen$par[2],gen$par[3],gen$par[4])
      print(data.frame(Parameter=namegen,Value=outgen)) 
      F<-objfun(gen$par)
      
      opt<-optim(c(gen$par[1],gen$par[2],gen$par[3],gen$par[4]),objfun,method="Nelder-Mead")             
      nameopt<-c("A","a","B","b")
      outopt<-c(opt$par[1],opt$par[2],opt$par[3],opt$par[4])
      cat("\n<< The value of parameter fitted by Nelder-Mead Simplex algorithm >>\n\n")
      print(data.frame(Parameter=nameopt,Value=outopt))
      cat("\n<< Residuals sum-of-squares and parameter values fitted by nls >>\n\n")
      fm<-nls(conc~defun(time,A,a,B,b),data=PKindex, algorithm="default",subset=Subject==i,
          start=list(A=opt$par[1], a=opt$par[2], B=opt$par[3], b=opt$par[4]),trace=TRUE,
          nls.control(maxiter=500,tol=1e-2,minFactor=1/1024))
      cat("\n")        
      coef<-data.frame(coef(fm)["b"])     
      plotting.lin(PKindex, fm, i, pick, coef, xaxis, yaxis)
      }
     })
     cat("\n") 
}   


### Three exponential term
fmacro.three<- function(PKindex,
                        A=NULL, 
                        a=NULL,
                        B=NULL,
                        b=NULL,
                        C=NULL,
                        c=NULL) 
{
   #options(warn=-1)
        
   ## Input initial value for A, a, B, b, C and c
   
   if (is.null(A) || is.null(a) || is.null(B) || is.null(b) || is.null(C) || is.null(c) ) {
       par<-data.frame(Parameter=c("A","a","B","b","C","c"),Initial=c(0))
       par<-edit(par)
       repeat{
           if ( par[1,2] == 0 || par[2,2] ==0 || par[3,2]==0 || par[4,2]==0 || par[5,2]==0 || par[6,2]==0){
             cat("\n")
             cat("**********************************\n")
             cat(" Parameter value can not be zero. \n")
             cat(" Press enter to continue.         \n")
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
      }
      
    cat("\n")
    
    defun<-function(time,A,a,B,b,C,c){
      pred<- A*exp(-a*time)+B*exp(-b*time)+C*exp(-c*time)
    }  
    
    ## Select weighting schemes
    file.menu <-c("equal weight",
                 "1/Cp",
                 "1/Cp^2")
    pick <- menu(file.menu, title="<< weighting schemes >>")

    with(entertitle(),{
    
    for( i in 1:length(unique(PKindex$Subject)))  {
      cat("\n\n               << Subject",i,">>\n\n" )  
      objfun<-function(par) {
         out<-defun(PKindex$time[PKindex$Subject==i],par[1],par[2],par[3],par[4],par[5],par[6])
         gift<- which( PKindex$conc[PKindex$Subject==i] != 0 )
         switch(pick,
                sum((PKindex$conc[PKindex$Subject==i][gift]-out[gift])^2),
                sum((PKindex$conc[PKindex$Subject==i][gift]-out[gift])^2/PKindex$conc[gift]),
                sum(((PKindex$conc[PKindex$Subject==i][gift] - out[gift])/PKindex$conc[gift])^2)
                )
      }
      
      gen<-genoud(objfun,nvars=6,max=FALSE,pop.size=30,max.generations=20,
           wait.generations=10,
           starting.value=c(par[1,2],par[2,2],par[3,2],par[4,2],par[5,2],par[6,2]),
           BFGS=FALSE,print.level=0,boundary.enforcement=2,
           Domains=matrix(c(1,0.1,1,0.01,0.1,0.001,100,10,50,10,10,1),6,2),
           MemoryMatrix=TRUE)
           
      cat("<< The value of parameter obtained from genetic algorithm >>\n\n")   
      namegen<-c("A","a","B","b","C","c")
      outgen<-c(gen$par[1],gen$par[2],gen$par[3],gen$par[4],gen$par[5],gen$par[6])
      print(data.frame(Parameter=namegen,Value=outgen)) 
      F<-objfun(gen$par)
      
      opt<-optim(c(gen$par[1],gen$par[2],gen$par[3],gen$par[4],gen$par[5],gen$par[6]),objfun,method="Nelder-Mead")       
      nameopt<-c("A","a","B","b","C","c")
      outopt<-c(opt$par[1],opt$par[2],opt$par[3],opt$par[4],opt$par[5],opt$par[6])
      cat("\n<< The value of parameter fitted by Nelder-Mead Simplex algorithm >>\n\n")
      print(data.frame(Parameter=nameopt,Value=outopt))
      cat("\n<< Residuals sum-of-squares and parameter values fitted by nls >>\n\n")
      fm<-nls(conc~defun(time,A,a,B,b,C,c),data=PKindex, algorithm="default",subset=Subject==i,
          start=list(A=opt$par[1],a=opt$par[2],B=opt$par[3],b=opt$par[4],C=opt$par[5],c=opt$par[6]),
          trace=TRUE,nls.control(maxiter=500,tol=1,minFactor=1/2048))
      cat("\n")
      coef<-data.frame(coef(fm)["c"])     
      plotting.lin(PKindex, fm, i, pick, coef, xaxis, yaxis)
      }
     })
     cat("\n") 
}   

### Simulation
###One exponential term
smacro.one <-function(Subject=NULL,  # N Subj's 
                      PKtime=NULL,   # times for sampling
                      A=NULL,
                      a=NULL)
{
   #options(warn=-1)
   
   if (is.null(Subject) || !is.integer(Subject)) {
     cat("How many subject do you want?\n")
     Subject<-scan(nlines=1,quiet=TRUE)
   }

   if (is.null(PKtime)) { 
     ## need to verify is a numeric ordered vector.
     PKtime<-data.frame(time=c(0))
     PKtime<-edit(PKtime)
   }
  
   cat("\n")
   show(PKtime)
   
   if (is.null(A) || is.null(a)) {
       par<-data.frame(Parameter=c("A","a"),Initial=c(0))
       par<-edit(par)
       repeat{
           if ( par[1,2] == 0 || par[2,2] ==0){
             cat("\n")
             cat("**********************************\n")
             cat(" Parameter value can not be zero. \n")
             cat(" Press enter to continue.         \n")
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
    } 
    else {
       par1 <- A
       par2 <- a
    } 
    
    file.menu <- c("Simulation with Error",
                   "Monte Carlo Simulation")
    pick <- menu(file.menu, title = "<< Simulation Type >>") 
    
    if (pick ==1){
      cat("\n\n")
      file.menu <- c("no error",
                     "error=normal error", 
                     "error=uniform error",
                     "error=normal error*true value",
                     "error=uniform error*true value")
      pick <- menu(file.menu, title = "<< Error type >>")      
      
      type<-switch(pick, 
                   "No Error",
                   "Error=Normal Error", 
                   "Error=Uniform Error",
                   "Error=Normal Error*True Value",
                   "Error=Uniform Error*True Value")
      
      if (pick ==1){
          PKindex<-vector(Subject,mode="list")
          for( i in 1:Subject)  {
            cat("\n\n             << Subject",i,">>\n\n" ) 
            A<-par1
            a<-par2  
            PKindex[[i]]<-smacro.one.out(PKtime,A,a,defun,par1,par2,i,type)
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
                cat(" Press enter to continue.         \n")
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
                cat(" Press enter to continue.         \n")
                cat("**********************************\n\n")
                readline()
                cat("\n")
                factor2<-scan(nlines=1,quiet=TRUE)}   
           else{
                break
                return(factor2)}
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
                       a<-par2+rnorm(1,mean=0,sd=factor2)}},
                       
                    {A<-par1+runif(1,min=-factor1,max=factor1)
                     while(A<=0){
                        A<-par1+runif(1,min=-factor1,max=factor1)}
                     a<-par2+runif(1,min=-factor2,max=factor2)
                     while(a<=0){
                        a<-par2+runif(1,min=-factor2,max=factor2)}},
                        
                    {A<-par1*rnorm(1,mean=0,sd=factor1)+par1
                     while(A<=0){
                        A<-par1*rnorm(1,mean=0,sd=factor1)+par1}
                     a<-par2*rnorm(1,mean=0,sd=factor2)+par2
                     while(a<=0){
                        a<-par2*rnorm(1,mean=0,sd=factor2)+par2}},
                        
                    {A<-par1*runif(1,min=-factor1,max=factor1)+par1
                     while(A<=0){
                        A<-par1*runif(1,min=-factor1,max=factor1)+par1}
                     a<-par2*runif(1,min=-factor2,max=factor2)+par2
                     while(a<=0){
                        a<-par2*runif(1,min=-factor2,max=factor2)+par2}}   
                    )
            PKindex[[i]]<-smacro.one.out(PKtime,A,a,defun,par1,par2,i,type)    
            }  
         PKindex<- as.data.frame(do.call("rbind",PKindex))
         rownames(PKindex) <- seq(nrow(PKindex)) 
         savefile(PKindex)    
        }  
      }
      else if (pick ==2){ 
      cat("\n\n")
      file.menu <- c("error=normal error", 
                     "error=uniform error",
                     "error=normal error*true value",
                     "error=uniform error*true value")
      pick <- menu(file.menu, title = "<< Error type >>")
      
      type<-switch(pick, 
                  "Error=Normal Error", 
                  "Error=Uniform Error",
                  "Error=Normal Error*True Value",
                  "Error=Uniform Error*True Value")
      
      cat("\n\nHow many times do you want to iteration ?\n")
      re<-scan(nlines=1,quiet=TRUE)
      
      cat("\nEnter error factor for A\n")
      factor1<-scan(nlines=1,quiet=TRUE)
      repeat{
              if ( factor1 == 0 ){
                cat("\n")
                cat("**********************************\n")
                cat(" Parameter value can not be zero. \n")
                cat(" Press enter to continue.         \n")
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
                cat(" Press enter to continue.         \n")
                cat("**********************************\n\n")
                readline()
                cat("\n")
                factor2<-scan(nlines=1,quiet=TRUE)}   
           else{
                break
                return(factor2)}
           }
      
      cat("\n")
      cat("**********************************\n")
      cat("Summary Table                     \n")
      cat("Model: One-exponential Term Model \n") 
      cat("Subject #:", Subject,"            \n")
      cat("Error Type:", type,"              \n")
      cat("Simulation #:", re,"              \n\n")
      sim<-matrix(c(par1,par2,factor1,factor2),2,2)
      dimnames(sim)<-list(c("A","a"),c("Original","Error factor"))
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
                         a<-par2+rnorm(1,mean=0,sd=factor2)}},
                         
                     {A<-par1+runif(1,min=-factor1,max=factor1)
                      while(A<=0){
                         A<-par1+runif(1,min=-factor1,max=factor1)} 
                      a<-par2+runif(1,min=-factor2,max=factor2)
                      while(a<=0){
                      a<-par2+runif(1,min=-factor2,max=factor2)}},
                      
                      {A<-par1*rnorm(1,mean=0,sd=factor1)+par1
                       while(A<=0){
                          a<-par1*rnorm(1,mean=0,sd=factor1)+par1}
                       a<-par2*rnorm(1,mean=0,sd=factor2)+par2
                       while(a<=0){
                          a<-par2*rnorm(1,mean=0,sd=factor2)+par2}},
                          
                      {A<-par1*runif(1,min=-factor1,max=factor1)+par1
                       while(A<=0){
                          A<-par1*runif(1,min=-factor1,max=factor1)+par1}
                       a<-par2*runif(1,min=-factor2,max=factor2)+par2
                       while(a<=0){
                          a<-par2*runif(1,min=-factor2,max=factor2)+par2}}
              )
              time1<-PKtime$time
              defun<- A*exp(-a*time1) 
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
     cat("How many subject do you want?\n")
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
             cat(" Press enter to continue.         \n")
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
    pick <- menu(file.menu, title = "<< Simulation Type >>") 
    if (pick ==1){
      cat("\n\n")
      file.menu <- c("No Error",
                     "Error=Normal Error", 
                     "Error=Uniform Error",
                     "Error=Normal Error*True Value",
                     "Error=Uniform Error*True Value")
      pick <- menu(file.menu, title = "<< Error Type >>")      
      
      type<-switch(pick, 
                  "No Error",
                  "Error=Normal Error", 
                  "Error=Uniform Error",
                  "Error=Normal Error*True Value",
                  "Error=Uniform Error*True Value")
      
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
                cat(" Press enter to continue.         \n")
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
                cat(" Press enter to continue.         \n")
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
                cat(" Press enter to continue.         \n")
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
                cat(" Press enter to continue.         \n")
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
      file.menu <- c("Error=Normal Error", 
                     "Error=Uniform Error",
                     "Error=Normal Error*True Value",
                     "Error=Uniform Error*True Value")
      pick <- menu(file.menu, title = "<< Error Type >>")
      
      type<-switch(pick, 
                  "Error=Normal Error", 
                  "Error=Uniform Error",
                  "Error=Normal Error*True Value",
                  "Error=Uniform Error*True Value")
      
      cat("\n\nHow many times do you want to iteration ?\n")
      re<-scan(nlines=1,quiet=TRUE)
      
      cat("\nEnter error factor for A\n")
      factor1<-scan(nlines=1,quiet=TRUE)
      repeat{
              if ( factor1 == 0 ){
                cat("\n")
                cat("**********************************\n")
                cat(" Parameter value can not be zero. \n")
                cat(" Press enter to continue.         \n")
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
                cat(" Press enter to continue.         \n")
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
                cat(" Press enter to continue.         \n")
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
                cat(" Press enter to continue.         \n")
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
      dimnames(sim)<-list(c("A","a","B","b"),c("Original","Error factor"))
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



###Three exponential term
smacro.three <-function(Subject=NULL,  # N Subj's 
                        PKtime=NULL,   # times for sampling
                        A=NULL,
                        a=NULL,
                        B=NULL,
                        b=NULL,
                        C=NULL,
                        c=NULL)
{
   #options(warn=-1)
  
   if (is.null(Subject) || !is.integer(Subject)) {
     cat("How many subject do you want?\n")
     Subject<-scan(nlines=1,quiet=TRUE)
   }

   if (is.null(PKtime)) { 
     ## need to verify is a numeric ordered vector.
     PKtime<-data.frame(time=c(0))
     PKtime<-edit(PKtime)
   }
  
   cat("\n")
   show(PKtime)
   
   if (is.null(A) || is.null(a) || is.null(B) || is.null(b) || is.null(C) || is.null(c)) {
       par<-data.frame(Parameter=c("A","a","B","b","C","c"),Initial=c(0))
       par<-edit(par)
       repeat{
           if ( par[1,2] == 0 || par[2,2] ==0 || par[3,2]==0 || par[4,2]==0 || par[5,2]==0 || par[6,2]==0){
             cat("\n")
             cat("**********************************\n")
             cat(" Parameter value can not be zero. \n")
             cat(" Press enter to continue.         \n")
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
       par5<-par[5,2]
       par6<-par[6,2]
    } 
    else {
       par1 <- A
       par2 <- a
       par3 <- B
       par4 <- b
       par5 <- C
       par6 <- c
    } 
    
    file.menu <- c("Simulation with Error",
                   "Monte Carlo Simulation")
    pick <- menu(file.menu, title = "<< Simulation Type >>") 
    if (pick ==1){
      cat("\n\n")
      file.menu <- c("No Error",
                     "Error=Normal Error", 
                     "Error=Uniform Error",
                     "Error=Normal Error*True Value",
                     "Error=Uniform Error*True Value")
      pick <- menu(file.menu, title = "<< Error Type >>")      
      
      type<-switch(pick, 
                  "No Error",
                  "Error=Normal Error", 
                  "Error=Uniform Error",
                  "Error=Normal Error*True Value",
                  "Error=Uniform Error*True Value")
      
      if (pick ==1){
          PKindex<-vector(Subject,mode="list")
          for( i in 1:Subject)  {
            cat("\n\n             << Subject",i,">>\n\n" ) 
            A<-par1
            a<-par2  
            B<-par3
            b<-par4
            C<-par5
            c<-par6
            PKindex[[i]]<-smacro.three.out(PKtime,A,a,B,b,C,c,defun,par1,par2,par3,par4,par5,par6,i,type)  
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
                cat(" Press enter to continue.         \n")
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
                cat(" Press enter to continue.         \n")
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
                cat(" Press enter to continue.         \n")
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
                cat(" Press enter to continue.         \n")
                cat("**********************************\n\n")
                readline()
                cat("\n")
                factor4<-scan(nlines=1,quiet=TRUE)}   
           else{
                break
                return(factor4)}
           }
         
         cat("\nEnter error factor for C\n")
         factor5<-scan(nlines=1,quiet=TRUE)
         repeat{
              if ( factor5 == 0 ){
                cat("\n")
                cat("**********************************\n")
                cat(" Parameter value can not be zero. \n")
                cat(" Press enter to continue.         \n")
                cat("**********************************\n\n")
                readline()
                cat("\n")
                factor5<-scan(nlines=1,quiet=TRUE)}   
           else{
                break
                return(factor5)}
           }
         
         cat("\nEnter error factor for c\n")
         factor6<-scan(nlines=1,quiet=TRUE)
         repeat{
              if ( factor6 == 0 ){
                cat("\n")
                cat("**********************************\n")
                cat(" Parameter value can not be zero. \n")
                cat(" Press enter to continue.         \n")
                cat("**********************************\n\n")
                readline()
                cat("\n")
                factor6<-scan(nlines=1,quiet=TRUE)}   
           else{
                break
                return(factor6)}
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
                       b<-par4+rnorm(1,mean=0,sd=factor4)}
                    C<-par5+rnorm(1,mean=0,sd=factor5)
                    while(C<=0){
                       C<-par5+rnorm(1,mean=0,sd=factor5)}
                    c<-par6+rnorm(1,mean=0,sd=factor6)
                    while(c<=0){
                       c<-par6+rnorm(1,mean=0,sd=factor6)}},
                       
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
                        b<-par4+runif(1,min=-factor4,max=factor4)}
                     C<-par5+runif(1,min=-factor5,max=factor5)
                     while(C<=0){
                        C<-par5+runif(1,min=-factor5,max=factor5)}
                     c<-par6+runif(1,min=-factor6,max=factor6)
                     while(c<=0){
                        c<-par6+runif(1,min=-factor6,max=factor6)}},
                        
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
                        b<-par4*rnorm(1,mean=0,sd=factor4)+par4}
                     C<-par5*rnorm(1,mean=0,sd=factor5)+par5
                     while(C<=0){
                        C<-par5*rnorm(1,mean=0,sd=factor5)+par5}
                     c<-par6*rnorm(1,mean=0,sd=factor6)+par6
                     while(c<=0){
                        c<-par6*rnorm(1,mean=0,sd=factor6)+par6}},
                        
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
                        b<-par4*runif(1,min=-factor4,max=factor4)+par4}
                     C<-par5*runif(1,min=-factor5,max=factor5)+par5
                     while(C<=0){
                        C<-par5*runif(1,min=-factor5,max=factor5)+par5}
                     c<-par6*runif(1,min=-factor6,max=factor6)+par6
                     while(c<=0){
                        c<-par6*runif(1,min=-factor6,max=factor6)+par6}}   
                    )
            PKindex[[i]]<-smacro.three.out(PKtime,A,a,B,b,C,c,defun,par1,par2,par3,par4,par5,par6,i,type)          
            }  
         PKindex<- as.data.frame(do.call("rbind",PKindex))
         rownames(PKindex) <- seq(nrow(PKindex)) 
         savefile(PKindex)    
        }  
      }
      else if (pick ==2){ 
      cat("\n\n")
      file.menu <- c("Error=Normal Error", 
                     "Error=Uniform Error",
                     "Error=Normal Error*True Value",
                     "Error=Uniform Error*True Value")
      pick <- menu(file.menu, title = "<< Error Type >>")
      
      type<-switch(pick, 
                  "Error=Normal Error", 
                  "Error=Uniform Error",
                  "Error=Normal Error*True Value",
                  "Error=Uniform Error*True Value")
      
      cat("\n\nHow many times do you want to iteration ?\n")
      re<-scan(nlines=1,quiet=TRUE)
      
      cat("\nEnter error factor for A\n")
      factor1<-scan(nlines=1,quiet=TRUE)
      repeat{
              if ( factor1 == 0 ){
                cat("\n")
                cat("**********************************\n")
                cat(" Parameter value can not be zero. \n")
                cat(" Press enter to continue.         \n")
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
                cat(" Press enter to continue.         \n")
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
                cat(" Press enter to continue.         \n")
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
                cat(" Press enter to continue.         \n")
                cat("**********************************\n\n")
                readline()
                cat("\n")
                factor4<-scan(nlines=1,quiet=TRUE)}   
           else{
                break
                return(factor4)}
           }
      
      cat("\nEnter error factor for C\n")
      factor5<-scan(nlines=1,quiet=TRUE)
      repeat{
              if ( factor5 == 0 ){
                cat("\n")
                cat("**********************************\n")
                cat(" Parameter value can not be zero. \n")
                cat(" Press enter to continue.         \n")
                cat("**********************************\n\n")
                readline()
                cat("\n")
                factor5<-scan(nlines=1,quiet=TRUE)}   
           else{
                break
                return(factor5)}
           }
      
      cat("\nEnter error factor for c\n")
      factor6<-scan(nlines=1,quiet=TRUE)
      repeat{
              if ( factor6 == 0 ){
                cat("\n")
                cat("**********************************\n")
                cat(" Parameter value can not be zero. \n")
                cat(" Press enter to continue.         \n")
                cat("**********************************\n\n")
                readline()
                cat("\n")
                factor6<-scan(nlines=1,quiet=TRUE)}   
           else{
                break
                return(factor6)}
           }
      
      cat("\n")
      cat("************************************\n")
      cat("Summary Table                       \n")
      cat("Model: Three-exponential Term Model \n") 
      cat("Subject #:", Subject,"              \n")
      cat("Error Type:", type,"                \n")
      cat("Simulation #:", re,"                \n\n")
      sim<-matrix(c(par1,par2,par3,par4,par5,par6,factor1,factor2,factor3,factor4,factor5,factor6),6,2)
      dimnames(sim)<-list(c("A","a","B","b","C","c"),c("Original","Error factor"))
      show(sim)   
      cat("************************************\n\n")
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
                       b<-par4+rnorm(1,mean=0,sd=factor4)}
                    C<-par5+rnorm(1,mean=0,sd=factor5)
                    while(C<=0){
                       C<-par5+rnorm(1,mean=0,sd=factor5)}
                    c<-par6+rnorm(1,mean=0,sd=factor6)
                    while(c<=0){
                       c<-par6+rnorm(1,mean=0,sd=factor6)}},
                       
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
                        b<-par4+runif(1,min=-factor4,max=factor4)}
                     C<-par5+runif(1,min=-factor5,max=factor5)
                     while(C<=0){
                        C<-par5+runif(1,min=-factor5,max=factor5)}
                     c<-par6+runif(1,min=-factor6,max=factor6)
                     while(c<=0){
                        c<-par6+runif(1,min=-factor6,max=factor6)}},
                        
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
                        b<-par4*rnorm(1,mean=0,sd=factor4)+par4}
                     C<-par5*rnorm(1,mean=0,sd=factor5)+par5
                     while(C<=0){
                        C<-par5*rnorm(1,mean=0,sd=factor5)+par5}
                     c<-par6*rnorm(1,mean=0,sd=factor6)+par6
                     while(c<=0){
                        c<-par6*rnorm(1,mean=0,sd=factor6)+par6}},
                        
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
                        b<-par4*runif(1,min=-factor4,max=factor4)+par4}
                     C<-par5*runif(1,min=-factor5,max=factor5)+par5
                     while(C<=0){
                        C<-par5*runif(1,min=-factor5,max=factor5)+par5}
                     c<-par6*runif(1,min=-factor6,max=factor6)+par6
                     while(c<=0){
                        c<-par6*runif(1,min=-factor6,max=factor6)+par6}}   
              )
              time1<-PKtime$time
              defun<- A*exp(-a*time1)+B*exp(-b*time1)+C*exp(-c*time1)
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

