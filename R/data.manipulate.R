#Data manipulation
data.manipulate <- function()
{
  file.menu <- c("Load Data Files (.CSV)", 
                 "Load Data Files (.RData)", 
                 "Input data from keyboard", 
                 "Go Back One Upper Level")  
  pick <- menu(file.menu, title = "<< Data edit >>")
  
  if (pick == 1){
     cat("\n")
     cat("*************************************************************\n")
     cat(" Enter data file name(.csv)                                  \n")
     cat(" Data should consist of subject no., time, and drug conc. (Cp). \n")
     cat("*************************************************************\n\n")
     ### PK.file <-readline()
     ### PK.file<-paste(PK.file,".csv",sep="")
     cnames<-c("Subject", "time", "conc")
     PKindex<-read.csv(file.choose(),header=TRUE,sep=",",row.names=NULL,col.names=cnames)
     PKindex<-edit(PKindex)
     cat("\n\n")
     show(PKindex)
     cat("\n\n")
     cat("**********************************\n")
     cat(" Now you can select a PK model.   \n")
     cat("**********************************\n\n")   
     return(PK.fit(PKindex))
  }
     
  else if (pick == 2){
     ### cat("\nEnter data file name\n") 
     ### PKname <-readline()
     ### PKname<-paste(PKname,".RData",sep="")
     PKindex<-readRDS(file.choose())
     PKindex<-edit(PKindex)
     colnames(PKindex)<-list("Subject", "time", "conc")
     write.csv(PKindex,file=PKname,row.names=FALSE)    ###
     cat("\n\n")
     show(PKindex)
     cat("\n\n")
     PK.fit(PKindex)
  }   
  else if (pick == 3){
     cat("\n\n") 
     PKindex<-data.frame(Subject=c(1),time=c(0),conc=c(0))
     PKindex<-edit(PKindex)
     show(PKindex)
     ans<-readline("\nSave the data (y/n)?\n")
     cat("\n")
     if (ans == "n" | ans == "N"){
        return(nor.fit(PKindex))
        }
     else {
        PKname <-readline("Enter the file name for input data (no file extension!):\n") 
        PKnameRData<-paste(PKname,".RData",sep="")
        PKnameCSV<-paste(PKname,".csv",sep="")
        if(file.exists(PKnameRData)||file.exists(PKnameCSV)){
           cat("\n")
           cat("****************************************\n")
           cat(" The file has been existed.       \n")
           cat(" Would you like to overwrite it? (y/n) \n")
           cat("****************************************\n")
           ans<-readline()
             if (ans == "y" | ans == "Y"){
                saveRDS(PKindex,file=PKnameRData)
                write.csv(PKindex,file=PKnameCSV,row.names=FALSE)
                cat("\n")
              }
              else{
                PKname <-readline("Enter the file name for input data (no file extension!):\n") 
                PKnameRData<-paste(PKname,".RData",sep="")
                PKnameCSV<-paste(PKname,".csv",sep="") 
                repeat{
                    if(file.exists(PKnameRData)||file.exists(PKnameCSV)){
                      cat("\n")
                      cat("***********************************\n")
                      cat(" The file has been existed. \n")
                      cat(" Please enter the file name again.\n")
                      cat("***********************************\n")
                      PKname<-readline()
                      PKnameRData<-paste(PKname,".RData",sep="")
                      PKnameCSV<-paste(PKname,".csv",sep="") 
                      }
                     else{
                      break                       
                      }
                  }        
              }   
           }
           else{
              saveRDS(PKindex,file=PKnameRData)
              write.csv(PKindex,file=PKnameCSV,row.names=FALSE)  
           }
          }
        cat("\n")  
        return(PK.fit(PKindex))
      }
  else if (pick == 4){
     cat("\n\n") 
     return(nor.fit(PKindex))
  } 
}
