### menu for entering PK model selection for normal fitting (nor.fit?)  --YJ
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
        cat("::::::::::::::::::::::::::::::::::::::::::::::::::::\n")
        cat(" Please enter the data or load your data file first.  \n")
        cat("::::::::::::::::::::::::::::::::::::::::::::::::::::\n\n")
        return(nor.fit())
       } 
     else
        return(PK.fit(PKindex))
  }      
  else if (pick == 3){
         cat("\n\n") 
         run()
  }  
}