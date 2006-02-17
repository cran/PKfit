#PKfit main menu
PKmenu <- function() 
{
  
  ## Need to clean up the requirements.
  
  options(warn=-1)
  
  if (!require(odesolve)) {
    ## lsoda is belong to odesolve package
    stop("Package odesolve not found.")
  }

  #if (!require(stats4)) {
  #  ## BIC is belong to stats4 package
  #  stop("Package stats4 not found.")
  #}

  if (!require(rgenoud)) {
    ## genoud is belong to rgenoud package
    stop("Package rgenoud not found.")
  }  
  
  cat("\n")
  file.menu <- c("Normal Fitting", 
                 "Simulation", 
                 "Quit")
  pick <- menu(file.menu, title = "<< PK startup menu >>")  
  if (pick == 1){
     cat("\n\n")
     cat("**********************************\n")
     cat(" Please manipulating data first!! \n")
     cat("**********************************\n\n")
     nor.fit()
  } else {
    if (pick == 2){
      cat("\n\n") 
      PK.sim()
    } else {
      if (pick == 3){
        cat("\nQuit !!\n\n")     
      }
    }  
  }
}

check<-function(par)
{
   repeat{
     if ( par[1,2] == 0 ){
       cat("\n")
       cat("*********************************\n")
       cat(" Parameter value can not be zero \n")
       cat(" Press enter to continue !!      \n")
       cat("*********************************\n\n")
       readline()
       cat("\n")
       par<-edit(par)
       }   
     else{
       break
       return(edit(par))
      }
  } 
  cat("\n")       
  show(par)     
}


#Normal fitting
nor.fit <- function(PKindex)
{
  file.menu <- c("Data Manipulation",
                 "Selection of PK Model",
                 "Go Back One Upper Level")
  pick <- menu(file.menu, title = "<< Normal Fitting >>")
  if (pick == 1){
     cat("\n\n")      
     data.manipulate()
  } 
  else if (pick == 2){
     cat("\n\n") 
     if (missing(PKindex)){
        cat("***********************************\n")
        cat(" Please manipulating data first !! \n")
        cat("***********************************\n\n")
        return(nor.fit())
       } 
     else
        return(PK.fit(PKindex))
  }      
  else if (pick == 3){
         cat("\n\n") 
         PKmenu()
  }  
}

#Data manipulation
data.manipulate <- function()
{
  file.menu <- c("Input Path",
                 "Load Data Files (.CSV)", 
                 "Load Data Files (.RData)", 
                 "Key In Data", 
                 "Go Back One Upper Level")  
  pick <- menu(file.menu, title = "<< Data edit >>")
  
  if (pick == 1){
     cat("\n")
     cat("****************************\n")
     cat(" Enter path for '.csv' file \n")
     cat(" Without quote              \n")
     cat(" With extension             \n")
     cat("****************************\n\n")
     pk.path<-readline()
     pk.path<-paste(pk.path,".csv",sep="")
     PKindex<-read.csv(pk.path)
     cnames<-c("Subject", "time", "conc")
     PKindex<-edit(PKindex)
     cat("\n\n")
     show(PKindex)
     cat("\n\n")
     cat("***************************\n")
     cat(" Please select PK model !! \n")
     cat("***************************\n\n")  
     return(nor.fit(PKindex))   
  }
  
  else if (pick == 2){
     cat("\n")
     cat("*************************************************************\n")
     cat(" Enter data file name(.csv)                                  \n")
     cat(" Data should consist of subject no., time, and concnetration \n")
     cat("*************************************************************\n\n")
     PK.file <-readline()
     PK.file<-paste(PK.file,".csv",sep="")
     cnames<-c("Subject", "time", "conc")
     PKindex<-read.csv(PK.file,header=TRUE,sep=",",row.names=NULL,col.names=cnames)
     PKindex<-edit(PKindex)
     cat("\n\n")
     show(PKindex)
     cat("\n\n")
     cat("***************************\n")
     cat(" Please select PK model !! \n")
     cat("***************************\n\n")  
     return(nor.fit(PKindex))   
  } 
  else if (pick == 3){
     cat("\nEnter data file name\n") 
     PKname <-readline()
     PKname<-paste(PKname,".RData",sep="")
     load(PKname)
     PKindex<-edit(PKindex)
     colnames(PKindex)<-list("Subject", "time", "conc")
     cat("\n\n")
     show(PKindex)
     save(PKindex,file=PKname)
     cat("\n\n")
     nor.fit(PKindex)
  }   
  else if (pick == 4){
     cat("\n\n") 
     PKindex<-data.frame(Subject=c(1),time=c(0),conc=c(0))
     PKindex<-edit(PKindex)
     show(PKindex)     
     cat("\nSave data (y/n) ?\n")
     ans<-readline()
     cat("\n")
     if (ans == "n" | ans == "N"){
        return(nor.fit(PKindex))
        }
     else {
        cat("Enter name you want to call this data\n")
        PKname <-readline() 
        PKname<-paste(PKname,".RData",sep="")      
        if(file.exists(PKname)){
           cat("\n")
           cat("****************************************\n")
           cat(" The file name have been existed !!     \n")
           cat(" Would you want to overwrite it ? (y/n) \n")
           cat("****************************************\n")
           ans<-readline()
             if (ans == "y" | ans == "Y"){
                save(PKindex,file=PKname)
                cat("\n")
              }
              else{
                cat("\nEnter name you want to call this data\n")
                PKname <-readline() 
                PKname<-paste(PKname,".RData",sep="") 
                repeat{
                    if(file.exists(PKname)){
                      cat("\n")
                      cat("***********************************\n")
                      cat(" The file name have been existed **\n")
                      cat(" Enter name again, OK !!         **\n")
                      cat("***********************************\n")
                      PKname<-readline()
                      PKname<-paste(PKname,".RData",sep="") 
                      }
                     else{
                      break                       
                      }
                  }        
              }   
              save(PKindex,file=PKname)   
           }
        else{
           save(PKindex,file=PKname)
          }                            
        cat("\n")  
        return(nor.fit(PKindex))
      }      
  } 
  else if (pick == 5){
     cat("\n\n") 
     return(nor.fit(PKindex))
  } 
}

#PK model option
PK.fit <- function(PKindex)
{
  file.menu <- c("1-Compartment PK Model", 
                 "2-Compartment PK Model",
                 "Macroconstant Exponential Functions",
                 "Go Back One Upper Level")
  pick <- menu(file.menu, title = "<< Selection of PK Model >>")
  if (pick== 1){
     cat("\n\n")  
     one.list(PKindex)
  }     
  else if (pick == 2){
     cat("\n\n")
     two.list(PKindex)
  }
  else if (pick == 3){
     cat("\n\n")
     macro(PKindex)
  }
  else if (pick == 4){
     cat("\n\n")
     PKmenu()
  }        
}

#------------Normal fitting menu----------------
one.list <- function(PKindex)
{
  file.menu <- c("IV Route", 
                 "Non IV Route",
                 "Go Back One Upper Level")
  pick <- menu(file.menu, title = "<< 1-Compartment Model >>")
  if (pick == 1){
     cat("\n\n")  
     iv.route(PKindex)
  }
  else if (pick == 2){
     cat("\n\n")
     noniv.route(PKindex)
  }
  else if (pick == 3){
     cat("\n\n")
     PK.fit(PKindex)
  }              
}

two.list <- function(PKindex)
{
  cat("************************\n")
  cat(" SD: Single-Dose        \n")
  cat(" 1st-Ord: First-Ordered \n")
  cat(" Abs: Absorption        \n")
  cat(" w/o: without           \n")
  cat(" Tlag: Lag Time         \n")
  cat("************************\n\n")
  cat("\n")
  
  file.menu <- c("IV-Bolus, & SD", 
                 "IV-Infusion, & SD",
                 "Extravascular, SD, & 1st-Ord Abs w/o Tlag",
                 "Go Back One Upper Level")
  pick <- menu(file.menu, title = "<< 2-Compartment Model >>")
  if (pick == 1){
     cat("\n\n")  
     fbolus2(PKindex)
  }
  else if (pick == 2){
     cat("\n\n")
     finfu2(PKindex)
  }
  else if (pick == 3){
     cat("\n\n")
     ffirst2(PKindex)
  }
  else if (pick == 4){
     cat("\n\n")
     PK.fit(PKindex)
  }
}

PK.sim <- function()
{
  file.menu <- c("1-Compartment Model: IV Route", 
                 "1-Compartment Model: Non IV Route",
                 "2-Compartment Model",
                 "Macroconstant Exponential Functions",
                 "Go Back One Upper Level")
  pick <- menu(file.menu, title = "<< Selection of PK Model >>")
  if (pick == 1){
     cat("\n\n")  
     sone.iv.route()
  }
  else if (pick == 2){
     cat("\n\n")
     sone.noniv.route()
  }
  else if (pick == 3){
     cat("\n\n")
     stwo.all()
  }      
  else if (pick == 4){
     cat("\n\n")
     smacro()
  }
  else if (pick == 5){
     cat("\n\n")
     PKmenu()
  }
}

iv.route <- function(PKindex)
{
  cat("***************************************\n")
  cat(" SD: Single-Dose                       \n")
  cat(" MM Elim: Michaelis-Menten Elimination \n")
  cat("***************************************\n\n")
  cat("\n")
  
  file.menu <- c("IV-Bolus, & SD", 
                 "IV-Bolus, SD, & MM Elim", 
                 "IV-Infusion, & SD",
                 "IV-Infusion, SD, & MM Elim",
                 "Go Back One Upper Level")
  pick <- menu(file.menu, title = "<< IV Route >>")
  if (pick ==1){
     cat("\n\n")
     fbolus1(PKindex)
  }
  else if (pick == 2){
     cat("\n\n") 
     fbolus.mm(PKindex)
  }
  else if (pick == 3){
     cat("\n\n") 
     finfu1(PKindex)
  }
  else if (pick == 4){
     cat("\n\n") 
     finfu.mm(PKindex)
  }
  else if (pick == 5){
     cat("\n\n") 
     one.list(PKindex)
  }
}

noniv.route <- function(PKindex)
{
  cat("***************************************\n")
  cat(" SD: Single-Dose                       \n")
  cat(" 1st-Ord: First-Ordered                \n")
  cat(" Zero-Ord: Zero-Ordered                \n")
  cat(" Abs: Absorption                       \n")
  cat(" w: with                               \n")
  cat(" w/o: without                          \n")
  cat(" Tlag: Lag Time                        \n")
  cat(" MM Elim: Michaelis-Menten Elimination \n")
  cat("***************************************\n\n")
  cat("\n")

  file.menu <- c("Extravascular, SD, & 1-Ord Abs w Tlag",
                 "Extravascular, SD, & 1-Ord Abs w/o Tlag",
                 "Extravascular, SD, & Zero-Ord Abs w/o Tlag",
                 "Extravascular, SD, 1-Ord Abs, & MM Elim w Tlag ",
                 "Extravascular, SD, 1-Ord Abs, & MM Elim w/o Tlag",
                 "Extravascular, SD, Zero-Ord Abs, & MM Elim w/o Tlag",
                 "Go Back One Upper Level")
  pick <- menu(file.menu, title = "<< Non IV Route >>")
  if (pick == 1){
     cat("\n\n") 
     ffirst.lag(PKindex)
  }
  else if (pick == 2){
     cat("\n\n") 
     ffirst.nolag(PKindex)
  }
  else if (pick == 3){
     cat("\n\n") 
     fzero.nolag(PKindex)
  }
  else if (pick == 4){
     cat("\n\n") 
     ffirst.lagm(PKindex)
  }
  else if (pick == 5){
     cat("\n\n") 
     ffirst.nolagm(PKindex)
  }
  else if (pick == 6){
     cat("\n\n") 
     fzero.nolagm(PKindex)
  }
  else if (pick == 7){
     cat("\n\n") 
     one.list(PKindex)
  }
}

macro <- function(PKindex)
{
  file.menu <- c("1-Exponential Term", 
                 "2-Exponential Term",
                 "3-Exponential Term",
                 "Go Back One Upper Level")
  pick <- menu(file.menu, title = "<< Macroconstant Exponential Functions >>")
  if (pick ==1){
     cat("\n\n")
     fmacro.one(PKindex)
  }
  else if (pick == 2){
     cat("\n\n") 
     fmacro.two(PKindex)
  }
  else if (pick == 3){
     cat("\n\n") 
     fmacro.three(PKindex)
  }
  else if (pick == 4){
     cat("\n\n") 
     PK.fit(PKindex)
  }
}

#------------Simulation menu----------------
sone.iv.route <- function()
{
  cat("***************************************\n")
  cat(" SD: Single-Dose                       \n")
  cat(" MM Elim: Michaelis-Menten Elimination \n")
  cat("***************************************\n\n")
  cat("\n")
  file.menu <- c("IV-Bolus, & SD", 
                 "IV-Bolus, SD, & MM Elim", 
                 "IV-Infusion, & SD",
                 "IV-Infusion, SD, & MM Elim",
                 "Go Back One Upper Level")
  pick <- menu(file.menu, title = "<< IV (Bolus or Infusion) Route >>")
  if (pick ==1 ){
     cat("\n\n")
     sbolus1()
  }
  else if (pick == 2){
     cat("\n\n") 
     sbolus.mm()
  }
  else if (pick == 3){
     cat("\n\n") 
     sinfu1()
  }
  else if (pick == 4){
     cat("\n\n") 
     sinfu.mm()
  }
  else if (pick == 5){
     cat("\n\n") 
     PK.sim()
  }
}

sone.noniv.route <- function()
{
  cat("***************************************\n")
  cat(" SD: Single-Dose                       \n")
  cat(" 1st-Ord: First-Ordered                \n")
  cat(" Zero-Ord: Zero-Ordered                \n")
  cat(" Abs: Absorption                       \n")
  cat(" w: with                               \n")
  cat(" w/o: without                          \n")
  cat(" Tlag: Lag Time                        \n")
  cat(" MM Elim: Michaelis-Menten Elimination \n")
  cat("***************************************\n\n")
  cat("\n")
  
  file.menu <- c("Extravascular, SD, & 1-Ord Abs w Tlag",
                 "Extravascular, SD, & 1-Ord Abs w\o Tlag",
                 "Extravascular, SD, & Zero-Ord Abs w\o Tlag",
                 "Extravascular, SD, 1-Ord Abs, & MM Elim w Tlag",
                 "Extravascular, SD, 1-Ord Abs, & MM Elim w\o Tlag",
                 "Extravascular, SD, Zero-Ord Abs, & MM Elim w\o Tlag",
                 "Go Back One Upper Level")
  pick <- menu(file.menu, title = "<< Non IV Route >>")
  if (pick == 1){
     cat("\n\n") 
     sfirst.lag()
  }
  else if (pick == 2){
     cat("\n\n") 
     sfirst.nolag()
  }
  else if (pick == 3){
     cat("\n\n") 
     szero.nolag()
  }
  else if (pick == 4){
     cat("\n\n") 
     sfirst.lagm()
  }
  else if (pick == 5){
     cat("\n\n") 
     sfirst.nolagm()
  }
  else if (pick == 6){
     cat("\n\n") 
     szero.nolagm()
  }
  else if (pick == 7){
     cat("\n\n") 
     PK.sim()
  }
}

smacro <- function()
{
  file.menu <- c("1-Exponential Term", 
                 "2-Exponential Term", 
                 "3-Exponential Term",
                 "Go Back One Upper Level")
  pick <- menu(file.menu, title = "<< Macroconstant Exponential Functions >>")
  if (pick ==1){
     cat("\n\n")
     smacro.one()
  }
  else if (pick == 2){
     cat("\n\n") 
     smacro.two()
  }
  else if (pick == 3){
     cat("\n\n") 
     smacro.three()
  }
  else if (pick == 4){
     cat("\n\n") 
     PK.sim()
  }
}

stwo.all <- function()
{
  cat("************************\n")
  cat(" SD: Single-Dose        \n")
  cat(" 1st-Ord: First-Ordered \n")
  cat(" Abs: Absorption        \n")
  cat(" w/o: without           \n")
  cat(" Tlag: Lag Time         \n")
  cat("************************\n\n")
  cat("\n")
  file.menu <- c("IV-Bolus, & SD", 
                 "IV-Infusion, & SD",
                 "Extravascular, SD, & 1-Ord Abs w/o Tlag",
                 "Go Back One Upper Level")
  pick <- menu(file.menu, title = "<< 2-Compartment Model >>")
  if (pick ==1){
     cat("\n\n")
     sbolus2()
  }
  else if (pick == 2){
     cat("\n\n") 
     sinfu2()
  }
  else if (pick == 3){
     cat("\n\n") 
     sfirst2()
  }
  else if (pick == 4){
     cat("\n\n") 
     PK.sim()
  }  
}

#---------------plot for linear----------------
plotting.lin <- function (PKindex, fm, i, pick, coef, xaxis, yaxis,
                          separateWindows=TRUE)
{               
  #options(warn=-1)
  
  j<-1:length(PKindex$time[PKindex$Subject==i])
  x<-PKindex$time[PKindex$Subject==i]
  y<-PKindex$conc[PKindex$Subject==i]
        
 #Calculated concentration
  cal<-predict(fm,list(time=x))
    
 #Weighted residuals   
 if (!(pick %in% 1:3)) {
    stop("Pick is illegal")
  }
  
 wei <- switch(pick,
               ifelse(y[j]==0.0, 0, y[j]-cal[j]),
               ifelse(y[j]==0.0, 0, sqrt(1/(y[j]))*(y[j]-cal[j])),
               ifelse(y[j]==0.0, 0, sqrt(1/((y[j])^2))*(y[j]-cal[j])))
  
 #calculate AUC and AUMC   
  add<-function(time,conc){
     auc<-0 ; aumc<-0
     for(i in 2:length(time)) {
     auc[i]<-1/2*(time[i]-time[i-1])*(conc[i]+conc[i-1])
     auc[i]<-auc[i]+auc[i-1]
     aumc[i]<-1/2*(time[i]-time[i-1])*(conc[i]*time[i]+conc[i-1]*time[i-1])
     aumc[i]<-aumc[i]+aumc[i-1]
     }
     return(list(auc=auc,aumc=aumc))
  }
  add<-add(x,y)
  AUC<-add$auc
  AUMC<-add$aumc
        
 #Output   
  cat("<< Output >>\n")  
  output<-data.frame(x,y,cal,wei,AUC,AUMC)
  colnames(output)<-list("time","Observed","Calculated","Wtd Residuals","AUC","AUMC")
  show(output)  
  
 #AUC (0 to infinity)              
  cat("\n<< AUC (0 to infinity) computed by trapezoidal rule >>\n\n")
  auc.infinity<-y[length(y)]/coef[1,1]
  auc<-AUC[length(y)]+auc.infinity
  show(auc) 
  
 #AUMC (0 to infinity) 
  cat("\n<< AUMC (0 to infinity) computed by trapezoidal rule >>\n\n")
  aumc.infinity<-(x[length(x)]*y[length(y)])/coef[1,1]+x[length(x)]/((coef[1,1])^2)
  aumc<-AUMC[length(y)]+aumc.infinity
  show(aumc)   
  cat("\n")
         
  aicllsbc(fm)
  cat("\n\n")
    
 #Divide plot console into four parts
  if (separateWindows) {
    get(getOption("device"))()
  }
  par(mfrow=c(2,2))
  main<-paste(c("Subject", i ,"plot"),collapse=" ")
     
 #Linear plot
  plot(y~x,data=PKindex,type='p',main=main, 
       xlab=xaxis, ylab=yaxis,pch=15,col="black",bty="l",
       font.lab=2,cex.lab=1,cex.axis=1,cex.main=1) 
  lines(x,predict(fm,list(time=x)),type="l",lty=1,
        col="firebrick3",lwd="2")
  mtext("Linear",side=3,cex=0.88)
    
 #Semi-log plot
  plot(x,y,log="y",type='p',main=main,
       xlab=xaxis, ylab=yaxis,pch=15,col="black",bty="l",
       font.lab=2,cex.lab=1,cex.axis=1,cex.main=1) 
  lines(x,predict(fm,list(time=x)),type="l",lty=1,
        col="firebrick3",lwd="2")
  mtext("Semi-log",side=3,cex=0.88)
     
 #Residual plot, time vs weighted residual
  plot(x,wei,pch=15,col="blue",bty="l",xlab=xaxis,
       ylab="Weighted Residual",main="Residual Plots",cex.lab=1,
       cex.axis=1,cex.main=1,font.lab=2)
  abline(h=0,lwd=2,col="black",lty=2)
    
 #Residual plot, calculated concentration vs weigthed residual
  plot(cal,wei,pch=15,col="blue",bty="l",xlab="Calc Cp(i)",
       ylab="Weighted Residual",main="Residual Plots",cex.lab=1,
       cex.axis=1,cex.main=1,font.lab=2)
  abline(h=0,lwd=2,col="black",lty=2)          
}    

### ---------------plot for nonlinear----------------
plotting.non <- function (PKindex, fm, i, pick, xaxis, yaxis,
                          separateWindows=TRUE) 
{  

  #options(warn=-1)
              
  j <- 1:length(PKindex$time[PKindex$Subject==i])
  x <- PKindex$time[PKindex$Subject==i]
  y <- PKindex$conc[PKindex$Subject==i]
   
  ## Calculated concentration
  cal<-predict(fm,list(time=x))

  if (!(pick %in% 1:3)) {
    stop("Pick is illegal")
  }
  wei <- switch(pick,
                ifelse(y[j]==0.0, 0, y[j]-cal[j]),
                ifelse(y[j]==0.0, 0, sqrt(1/(y[j]))*(y[j]-cal[j])),
                ifelse(y[j]==0.0, 0, sqrt(1/((y[j])^2))*(y[j]-cal[j])))
  
 #calculate AUC and AUMC       
  add <- function(time,conc) {
    auc<-0 ; aumc<-0
    for(i in 2:length(time)) {
      auc[i]<-1/2*(time[i]-time[i-1])*(conc[i]+conc[i-1])
      auc[i]<-auc[i]+auc[i-1]
      aumc[i]<-1/2*(time[i]-time[i-1])*(conc[i]*time[i]+conc[i-1]*time[i-1])
      aumc[i]<-aumc[i]+aumc[i-1]
    }
    return(list(auc=auc,aumc=aumc))
  }
  
  add1 <- add(x,y)
  AUC <- add1$auc
  AUMC <- add1$aumc
              
  cat("<< Output >>\n\n")      
  output <- data.frame(x,y,cal,wei,AUC,AUMC)
  colnames(output) <- list("time","Observed","Calculated","Wtd Residuals","AUC","AUMC")
  show(output)
  cat("\n") 
        
  aicllsbc(fm)
  cat("\n\n")
    
  ## Divide plot console into four parts
  if (separateWindows) {
    get(getOption("device"))()
  }
  par(mfrow=c(2,2))
  main<-paste(c("Subject", i ,"plot"),collapse=" ")
     
 #Linear plot
  plot(y~x,data=PKindex,type='p',main=main, 
       xlab=xaxis, ylab=yaxis,pch=15,col="black",bty="l",
       font.lab=2,cex.lab=1,cex.axis=1,cex.main=1) 
  lines(x,predict(fm,list(time=x)),type="l",lty=1,
        col="firebrick3",lwd="2")
  mtext("Linear",side=3,cex=0.88)
  
 #Semi-log plot
  plot(x,y,log="y",type='p',main=main,
       xlab=xaxis, ylab=yaxis,pch=15,col="black",bty="l",
       font.lab=2,cex.lab=1,cex.axis=1,cex.main=1) 
  lines(x,predict(fm,list(time=x)),type="l",lty=1,
        col="firebrick3",lwd="2")
  mtext("Semi-log",side=3,cex=0.88)
  
  ## Residual plot, time vs weighted residual----- 
  plot(x,wei,pch=15,col="blue",bty="l",xlab=xaxis,
       ylab="Weighted Residual",main="Residual Plots",cex.lab=1,
       cex.axis=1,cex.main=1,font.lab=2)
  abline(h=0,lwd=2,col="black",lty=2)
  
  ## Residual plot, calcukated concentration vs weigthed residual-----
  plot(cal,wei,pch=15,col="blue",bty="l",xlab="Calc Cp(i)",
       ylab="Weighted Residual",main="Residual Plots",cex.lab=1,
       cex.axis=1,cex.main=1,font.lab=2)
  abline(h=0,lwd=2,col="black",lty=2)
}    

### ----------------plot for simulation----------------
### Draw 2 windows, and consider that 2 more windows are coming
plotting.sim <- function(i,x,y,separateWindows=TRUE) 
{
  #options(warn=-1)
  
  if (separateWindows) {
    get(getOption("device"))()
  }
  par(mfrow=c(2,2))
  main<-paste(c("Subject", i ,"plot"),collapse=" ")
   
  ## linear plot
  plot(y~x,type='p',main=main, 
       xlab="Time",ylab="Concentration",pch=15,col="black",bty="l",
       font.lab=2,cex.lab=1,cex.axis=1,cex.main=1) 
  text("Linear",side=3,cex=0.88)

  ## semi-log plot
  plot(x,y,log="y",type='p',main=main,
       xlab="Time",ylab="Concentration",pch=15,col="black",bty="l",
       font.lab=2,cex.lab=1,cex.axis=1,cex.main=1) 
  mtext("Semi-log",side=3,cex=0.88) 
  
  cat("\n\n")
}

#Estimate model fitting: AIC, Log likelihood, SBC
aicllsbc <- function(fm)
{   
  cat("\n") 
  cat("<< Akaike's Information Criterion (AIC) >>\n\n")
  show(AIC(fm))
    
  cat("\n<< Log likelihood >>\n\n")
  show(logLik(fm))
    
  if (!require(stats4)) {
    ## BIC is belong to stats4 package
    stop("Package stats4 not found.")
  }  
    
  cat("\n<< Schwarz's Bayesian Criterion (SBC) >>\n\n")
  show(BIC(fm))
  cat("\n")     
    
 #Summary the results of nls
  print(summary(fm))  
  cat("\n")   
  cat(date(),"\n\n")     
}  


savefile<-function(PKindex)  
{
  cat("\nSave data as a .RData file (y/n) ?\n")
  ans<-readline()
  cat("\n")
  if (ans == "n" | ans == "N"){
     cat("\nQuit !!\n")
     }
  else {
     cat("Enter name you want to call this data\n")
     PKname <-readline() 
     dataExt<- ".RData"
     PKname<-paste(PKname,dataExt,sep="")
     if(file.exists(PKname)){
           cat("\n")
           cat("*******************************************\n")
           cat("** The file name has been existed !!     **\n")
           cat("** Would you want to overwrite it ? (y/n)**\n")
           cat("*******************************************\n")
           ans<-readline()
             if (ans == "y" | ans == "Y"){
                save(PKindex,file=PKname)
              }
              else{
                cat("\nEnter name you want to call this data\n")
                PKname <-readline() 
                PKname<-paste(PKname,".RData",sep="") 
                repeat{
                  if (file.exists(PKname)){
                    cat("\n")
                    cat("***************************************\n")
                    cat("** The file name has been existed !! **\n")
                    cat("** Enter name again, OK !!           **\n")
                    cat("***************************************\n")
                    PKname<-readline()
                    PKname<-paste(PKname,".RData",sep="") 
                    }
                    else{
                     break
                     }
                   }
                 } 
                save(PKindex,file=PKname)
              }   
        else {
           save(PKindex,file=PKname)  
           }
    }
  cat("\n")  
  cat(date(),"\n")
}


#for iv bolus 
sbolus1.out<-function(PKtime,kel,Vd,defun,par1,par2,Dose,i,type)
{
  time<-PKtime$time
  parms<-c(kel=kel,Vd=Vd) 
  C1.lsoda<-data.frame(lsoda(Dose/Vd,c(0,time),defun,parms)) 
  cat("\n")
  cat("****************************************************\n")
  cat("Summary Table                                       \n")
  cat("Model: 1-Compartment, IV-Bolus, & Single-Dose Model \n")
  cat("Error Type:", type,"                                \n\n") 
  sim<-matrix(c(kel,Vd,par1,par2),2,2)
  dimnames(sim)<-list(c("kel","Vd"),c("Value","Original"))
  show(sim)
  cat("****************************************************\n\n")
  
  good <- ifelse(C1.lsoda[2:(length(time)+1),2]<=0,
                 0,
                 C1.lsoda[2:(length(time)+1),2])
  PKindex <- data.frame(i,
                        C1.lsoda[2:(length(time)+1),1],
                        good)
  colnames(PKindex)<-list("Subject","time","conc")
  show(PKindex)
  x <- C1.lsoda[2:(length(time)+1),1]
  y <- good
  plotting.sim(i,x,y)
  return(PKindex) 
}

## for iv bolus nonlinear elimination 
sbolus.mm.out<-function(PKtime,Vm,Km,Vd,defun,par1,par2,par3,Dose,i,type) 
{ 
  time<-PKtime$time
  parms<-c(Vm=Vm,Km=Km,Vd=Vd)  
  C1.lsoda<-data.frame(lsoda(Dose/Vd,c(0,time),defun,parms))     
  cat("\n")
  cat("*********************************************\n")
  cat("Summary Table                                \n")
  cat("Model: 1-Compartment, IV-Bolus, Single-Dose, \n") 
  cat("       & Michaelis-Menten Elimination Model  \n") 
  cat("Error Type:", type,"                         \n\n") 
  sim<-matrix(c(Vm,Km,Vd,par1,par2,par3),3,2)
  dimnames(sim)<-list(c("Vm","Km","Vd"),c("Value","Original"))
  show(sim)
  cat("*********************************************\n\n")
  
  good<-ifelse(C1.lsoda[2:(length(time)+1),2]<=0,
               0,
               C1.lsoda[2:(length(time)+1),2])
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

#for iv infusion 
sinfu1.out<-function(PKtime,kel,Vd,defun,par1,par2,Dose,i,type) 
{
  time<-PKtime$time
  parms<-c(kel=kel,Vd=Vd) 
  C1.lsoda<-data.frame(lsoda(0,c(0,time),defun,parms)) 
  cat("\n")
  cat("*******************************************************\n")
  cat("Summary Table                                          \n")
  cat("Model: 1-Compartment, IV-Infusion, & Single-Dose Model \n")
  cat("Error Type:", type,"                                   \n\n")
  sim<-matrix(c(kel,Vd,par1,par2),2,2)
  dimnames(sim)<-list(c("kel","Vd"),c("Value","Original"))
  show(sim)
  cat("*******************************************************\n\n")
  
  good<-ifelse(C1.lsoda[2:(length(time)+1),2]<=0,
               0,
               C1.lsoda[2:(length(time)+1),2])
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

#for iv infusion nonlinear elimination
sinfu.mm.out<-function(PKtime,Vm,Km,Vd,defun,par1,par2,par3,Dose,i,type) 
{
  time<-PKtime$time
  parms<-c(Vm=Vm,Km=Km,Vd=Vd)  
  C1.lsoda<-data.frame(lsoda(0,c(0,time),defun,parms))    
  cat("\n")
  cat("************************************************\n")
  cat("Summary Table                                   \n")
  cat("Model: 1-Compartment, IV-Infusion, Single-Dose, \n") 
  cat("       & Michaelis-Menten Elimination Model     \n")
  cat("Error Type:", type,"                            \n\n")
  sim<-matrix(c(Vm,Km,Vd,par1,par2,par3),3,2)
  dimnames(sim)<-list(c("Vm","Km","Vd"),c("Value","Original"))
  show(sim)    
  cat("************************************************\n\n")
  
  good<-ifelse(C1.lsoda[2:(length(time)+1),2]<=0,
               0,
               C1.lsoda[2:(length(time)+1),2])
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

#for extravascular first order absorption with/without lag time
sfirst1.out<-function(PKtime,ka,kel,Vd,defun,par1,par2,par3,Dose,i,type,Tlag)  
{
  time<-PKtime$time
  parms<-c(ka=ka,kel=kel,Vd=Vd)  
  C1.lsoda<-data.frame(lsoda(c(Dose,0),c(0,time), defun, parms))  
  
  if (missing(Tlag)){
  
  cat("\n\n")
  cat("*****************************************************\n")
  cat("Summary Table                                        \n")
  cat("Model: 1-Compartment, Extravascular, Single-Dose,    \n")
  cat("       & 1-Ordered Absorption without Lag Time Model \n")
  cat("Error Type:", type,"                                 \n\n")
  sim<-matrix(c(ka,kel,Vd,par1,par2,par3),3,2)
  dimnames(sim)<-list(c("ka","kel","Vd"),c("Value","Original"))
  show(sim)
  cat("*****************************************************\n\n")
  }
  
  else{
  
  cat("\n\n")
  cat("**************************************************\n")
  cat("Summary Table                                     \n")
  cat("Model: 1-Compartment, Extravascular, Single-Dose, \n")
  cat("       & 1-Ordered Absorption with Lag Time Model \n")
  cat("Error Type:", type,"                              \n\n")
  sim<-matrix(c(ka,kel,Vd,par1,par2,par3),3,2)
  dimnames(sim)<-list(c("ka","kel","Vd"),c("Value","Original"))
  show(sim)
  cat("**************************************************\n\n")
  }  
  
  good<-ifelse(C1.lsoda[2:(length(time)+1),3]<=1e-5,
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


#for extravascular first order absorption with/without lag time nonlinear elimination
sfirst.mm.out<-function(PKtime,ka,Vm,Km,Vd,defun,par1,par2,par3,par4,Dose,i,type,Tlag) 
{       
  time<-PKtime$time
  parms<-c(ka=ka,Vm=Vm,Km=Km,Vd=Vd) 
  C1.lsoda<-data.frame(lsoda(c(Dose,0),c(0,time), defun, parms)) 
  
  if (missing(Tlag)){
  
  cat("\n\n")
  cat("************************************************************************\n")
  cat("Summary Table                                                           \n")
  cat("Model: 1-Compartment, Extravascular, Single-Dose, 1-Ordered Absorption, \n")
  cat("       & Michaelis-Menten Elimination without Lag Time Model            \n")
  cat("Error Type:", type,"                                                    \n\n")
  sim<-matrix(c(ka,Vm,Km,Vd,par1,par2,par3,par4),4,2)
  dimnames(sim)<-list(c("ka","Vm","Km","Vd"),c("Value","Original"))
  show(sim)
  cat("************************************************************************\n\n")
  }
  
  else{
  cat("\n\n")
  cat("************************************************************************\n")
  cat("Summary Table                                                           \n")
  cat("Model: 1-Compartment, Extravascular, Single-Dose, 1-Ordered Absorption, \n")
  cat("       & Michaelis-Menten Elimination with Lag Time Model               \n")
  cat("Error Type:", type,"                                                    \n\n")
  sim<-matrix(c(ka,Vm,Km,Vd,par1,par2,par3,par4),4,2)
  dimnames(sim)<-list(c("ka","Vm","Km","Vd"),c("Value","Original"))
  show(sim)
  cat("************************************************************************\n\n")
  }
  
  
  good<-ifelse(C1.lsoda[2:(length(time)+1),3]<=1e-5,
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

#for extravascular zero order absorption without lag time
szero.out<-function(PKtime,Tabs,kel,Vd,defun,par1,par2,par3,Dose,i,type) 
{
  time<-PKtime$time
  parms<-c(Tabs=Tabs,kel=kel,Vd=Vd)  
  C1.lsoda<-data.frame(lsoda(0, c(0,time), defun, parms))   
  cat("\n")
  cat("********************************************************\n")
  cat("Summary Table                                           \n")
  cat("Model: 1-Compartment, Extravascular, Single-Dose,       \n")
  cat("       & Zero-Ordered Absorption without Lag Time Model \n")
  cat("Error Type:", type,"                                    \n\n")
  sim<-matrix(c(Tabs,kel,Vd,par1,par2,par3),3,2)
  dimnames(sim)<-list(c("Tabs","kel","Vd"),c("Value","Original"))
  show(sim)  
  cat("********************************************************\n\n")
  
  good<-ifelse(C1.lsoda[2:(length(time)+1),2]<=0,
               0,
               C1.lsoda[2:(length(time)+1),2])
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

#for extravascular zero order absorption without lag time nonlinear elimination
szero.mm.out<-function(PKtime,Tabs,Vm,Km,Vd,defun,par1,par2,par3,par4,Dose,i,type) 
{
  time<-PKtime$time
  parms<-c(Tabs=Tabs,Vm=Vm,Km=Km,Vd=Vd) 
  C1.lsoda<-data.frame(lsoda(0,c(0,time), defun, parms))   
  cat("\n\n")
  cat("***************************************************************************\n")
  cat("Summary Table                                                              \n")
  cat("Model: 1-Compartment, Extravascular, Single-Dose, Zero-Ordered Absorption, \n") 
  cat("       & Michaelis-Menten Elimination without Lag Time Model               \n")
  cat("Error Type:", type,"                                                       \n\n")
  sim<-matrix(c(Tabs,Vm,Km,Vd,par1,par2,par3,par4),4,2)
  dimnames(sim)<-list(c("Tabs","Vm","Km","Vd"),c("Value","Original"))
  show(sim)
  cat("***************************************************************************\n\n")
  
  good<-ifelse(C1.lsoda[2:(length(time)+1),2]<=0,
               0,
               C1.lsoda[2:(length(time)+1),2])
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

#for one exponential term
smacro.one.out<-function(PKtime,A,a,defun,par1,par2,Dose,i,type) 
{
  time<-PKtime$time
  defun<- A*exp(-a*time)
  output<-data.frame(time,defun) 
  cat("\n\n")
  cat("*******************************************\n")
  cat("Summary Table                              \n")
  cat("Model: One-exponential Term Model          \n") 
  cat("Error Type:", type,"                       \n\n")
  sim<-matrix(c(A,a,par1,par2),2,2)
  dimnames(sim)<-list(c("A","a"),c("Value","Original"))
  show(sim)
  cat("******************************************\n\n")
  
  PKindex<-data.frame(i,output)
  colnames(PKindex)<-list("Subject","time","conc")
  show(PKindex)
  x<-PKindex[,2]
  y<-PKindex[,3]
  plotting.sim(i,x,y)
  return(PKindex) 
}

#for two exponential term
smacro.two.out<-function(PKtime,A,a,B,b,defun,par1,par2,par3,par4,Dose,i,type)
{
  time<-PKtime$time
  defun<- A*exp(-a*time)+B*exp(-b*time)
  output<-data.frame(time,defun) 
  cat("\n\n")
  cat("*******************************************\n")
  cat("Summary Table                              \n")
  cat("Model: Two-exponential Term Model          \n") 
  cat("Error Type:", type,"                       \n\n")
  sim<-matrix(c(A,a,B,b,par1,par2,par3,par4),4,2)
  dimnames(sim)<-list(c("A","a","B","b"),c("Value","Original"))
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
  dimnames(sim)<-list(c("A","a","B","b","C","c"),c("Value","Original"))
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

#for two compartment iv bolus
sbolus2.out<-function(PKtime,kel,k12,k21,Vd,defun,par1,par2,par3,par4,Dose,i,type)
{ 
  time<-PKtime$time
  parms<-c(kel=kel,k12=k12,k21=k21,Vd=Vd) 
  C1.lsoda<-data.frame(lsoda(c(Dose/Vd,0),c(0,time),defun,parms)) 
  cat("\n\n")
  cat("****************************************************\n")
  cat("Summary Table                                       \n")
  cat("Model: 2-Compartment, IV-Bolus, & Single-Dose Model \n") 
  cat("Error Type:", type,"                                \n\n")
  sim<-matrix(c(kel,k12,k21,Vd,par1,par2,par3,par4),4,2)
  dimnames(sim)<-list(c("kel","k12","k21","Vd"),c("Value","Original"))
  show(sim)
  cat("****************************************************\n\n")
  good<-ifelse(C1.lsoda[2:(length(time)+1),2]<=0,
               0,
               C1.lsoda[2:(length(time)+1),2])
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

#for two compartment iv infusion
sinfu2.out<-function(PKtime,kel,k12,k21,Vd,defun,par1,par2,par3,par4,Dose,i,type)
{       
  time<-PKtime$time
  parms<-c(kel=kel,k12=k12,k21=k21,Vd=Vd) 
  C1.lsoda<-data.frame(lsoda(c(0,0),c(0,time),defun,parms)) 
  cat("\n\n")
  cat("*******************************************************\n")
  cat("Summary Table                                          \n")
  cat("Model: 2-Compartment, IV-Infusion, & Single-Dose Model \n") 
  cat("Error Type:", type,"                                   \n\n")
  sim<-matrix(c(kel,k12,k21,Vd,par1,par2,par3,par4),4,2)
  dimnames(sim)<-list(c("kel","k12","k21","Vd"),c("Value","Original"))
  show(sim)
  cat("*******************************************************\n\n")
  good<-ifelse(C1.lsoda[2:(length(time)+1),2]<=0,
               0,
               C1.lsoda[2:(length(time)+1),2])
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

#for two compartment extravascular first order absorption
sfirst2.out<-function(PKtime,ka,kel,k12,k21,Vd,defun,par1,par2,par3,par4,par5,Dose,i,type)
{                   
  time<-PKtime$time
  parms<-c(ka=ka,kel=kel,k12=k12,k21=k21,Vd=Vd) 
  C1.lsoda<-data.frame(lsoda(c(Dose,0,0),c(0,time),defun,parms)) 
  cat("\n\n")
  cat("*******************************************************\n")
  cat("Summary Table                                          \n")
  cat("Model: 2-Compartment, Extravascular,                   \n") 
  cat("       Single-Dose, & 1-Ordered without Lag Time Model \n") 
  cat("Error Type:", type,"                                   \n\n")
  sim<-matrix(c(ka,kel,k12,k21,Vd,par1,par2,par3,par4,par5),5,2)
  dimnames(sim)<-list(c("ka","kel","k12","k21","Vd"),c("Value","Original"))
  show(sim)
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

entertitle<-function()
{
  cat("\nEnter the title of x-axis(time)\n")
  cat("(or a blank line to use default)\n\n") 
  xaxis<-readline()
  if (substr(xaxis, 1, 1) == "")  xaxis<-"Time"  else xaxis<-xaxis
  cat("\nEnter the title of y-axis(Cp)\n")
  cat("(or a blank line to use default)\n\n") 
  yaxis<-readline()
  #cat("\n\n Please Wait !!  Data is Processing !! \n")
  if (substr(yaxis, 1, 1) == "")  yaxis<-"Concentration"  else yaxis<-yaxis
  #cat("\n\n Please Wait !!  Data is Processing !! \n")
  return(list(xaxis=xaxis,yaxis=yaxis))
}

montecarlo<-function(C1.lsoda,time,i,re)
{
  #options(warn=-1)
  conc<-rowMeans(as.data.frame(lapply(C1.lsoda,"[","concentration")))
  C1.lsoda<-do.call("rbind",C1.lsoda)
  rownames(C1.lsoda)<-seq(nrow(C1.lsoda))
  
  conc<-ifelse(conc<=0,
               0,
               conc)
  
  PKindex<-data.frame(i,time,conc)
  get(getOption("device"))()
  par(mfrow=c(2,2))
  x<-C1.lsoda$time
  #y<-C1.lsoda$concentration
  
  y<-ifelse(C1.lsoda$concentration<=0,
            0,
            C1.lsoda$concentration)
  
  main<-paste(c("Subject", i ,"plot"),collapse=" ")
  plot(C1.lsoda,type="p",main=main,xlab="Time",ylab="Concentration")
  lines(PKindex$time,PKindex$conc,type="l",lty=1,col="firebrick3",lwd="2")
  mtext("Linear",side=3,cex=0.8)
  
  plot(x,y,log="y",type='p',main=main,xlab="Time",ylab="Concentration")
  lines(PKindex$time,PKindex$conc,log="y",type="l",lty=1,col="firebrick3",lwd="2")
  mtext("Semi-log",side=3,cex=0.8) 
  
  len<-length(time)
  sub<-len*re
  a<-seq(0,sub,by=len)
  AA<-a[2:(length(a)-1)]
  colnames(PKindex)<-list("Subject","time","conc")
  show(PKindex)
  for(i in rev(AA))
      C1.lsoda<-C1.lsoda[append(1:(sub+(length(AA)-1)),NA,after=i),]   
  
  plot(C1.lsoda,type="l",main=main,xlab="Time",ylab="Concentration")  
  lines(PKindex$time,PKindex$conc,type="l",lty=1,col="firebrick3",lwd="2")
  mtext("Linear",side=3,cex=0.8)
  
  x<-C1.lsoda$time
  y<-C1.lsoda$concentration
  
  plot(x,y,log="y",type='l',main=main,xlab="Time",ylab="Concentration")
  lines(PKindex$time,PKindex$conc,log="y",type="l",lty=1,col="firebrick3",lwd="2")
  mtext("Semi-log",side=3,cex=0.8) 
  
  return(PKindex)
}