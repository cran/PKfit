#------------Normal fitting menu----------------
one.list <- function(PKindex)
{
  file.menu <- c("IV Route", 
                 "Non IV Route",
                 "Go Back One Upper Level")
  pick <- menu(file.menu, title = "<< 1-Compartment PK Model >>")
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
