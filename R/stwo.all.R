stwo.all <- function()
{
  cat("************************\n")
  cat("      SD: single-dose   \n")
  cat(" 1st-Ord: first-ordered \n")
  cat("     Abs: absorption    \n")
  cat("     w/o: without       \n")
  cat("    Tlag: lag-time      \n")
  cat("************************\n\n")
  cat("\n")
  file.menu <- c("IV-Bolus, & SD", 
                 "IV-Infusion, & SD",
                 "Extravascular, SD, & 1-Ord Abs w/o Tlag",
                 "Go Back One Upper Level")
  pick <- menu(file.menu, title = "<< 2-Compartment PK Model >>")
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
