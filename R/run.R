#PKfit main menu
run <- function() 
{
  
  ## Need to clean up the requirements.
  options(warn=-1) 
  graphics.off()
  cat("\n")
  file.menu <- c("Normal Fitting", 
                 "Simulation", 
                 "Quit")
  pick <- menu(file.menu, title = "<< PKfit:- Top menu >>")  
  if (pick == 1){
     cat("\n\n")
     cat("****************************************************\n")
     cat(" Please enter the data or load your data file first.  \n")
     cat("****************************************************\n\n")
     nor.fit()
  } else {
    if (pick == 2){
      cat("\n\n") 
      PK.sim()
    } else {
      if (pick == 3){
        cat("\n Thank for using PKfit. Bye now.\n\n")     
      }
    }  
  }
}
