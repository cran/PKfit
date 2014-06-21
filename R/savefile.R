savefile<-function(PKindex)  
{
  cat("\nSave input data now (y/n)?\n")
  ans<-readline()
  cat("\n")
  if (ans == "n" | ans == "N"){
     ### run()  ### do nothing here
     }
  else {
     cat("Enter the data file name (no file extension!):\n")
     PKname0 <-readline() 
     PKname<-paste(PKname0,".RData",sep="")
     PKnameCSV<-paste(PKname0,".csv",sep="")
     if(file.exists(PKname)){
           cat("\n")
           cat("*******************************************\n")
           cat("** The file name has been existed.     **\n")
           cat("** Would you like to overwrite it? (y/n)**\n")
           cat("*******************************************\n")
           ans<-readline()
             if (ans == "y" | ans == "Y"){
                saveRDS(PKindex,file=PKname)
                write.csv(PKindex,file=PKnameCSV,row.names=FALSE)  ### now save as .csv with Header too!  -YJ
              }
              else{
                cat("Enter the data file name (no file extension!):\n")
                PKname0<-readline() 
                PKname<-paste(PKname0,".RData",sep="") 
                PKnameCSV<-paste(PKname0,".csv",sep="")
                repeat{
                  if (file.exists(PKname)){
                    cat("\n")
                    cat("***************************************\n")
                    cat("** The file name has been existed. **\n")
                    cat("** Please enter file name again.           **\n")
                    cat("***************************************\n")
                    PKname0<-readline()
                    PKname<-paste(PKname0,".RData",sep="")
                    PKnameCSV<-paste(PKname0,".csv",sep="") 
                    }
                    else{
                     break
                     }
                   }
                 } 
                saveRDS(PKindex,file=PKname)
                write.csv(PKindex,file=PKnameCSV,row.names=FALSE)  ### now save as .csv too!  -YJ
              }   
        else {
           saveRDS(PKindex,file=PKname)
           write.csv(PKindex,file=PKnameCSV,row.names=FALSE)  ### now save as .csv too!  -YJ  
           }
    }
  cat("\n")  
  cat(date(),"\n\n")
}
