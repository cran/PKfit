stwo.MD.all <- function()
{
  OutputFilez() ### reset all output file names when running PK.sim()
  cat("************************\n")
  cat("      MD: multiplee-dose   \n")
  cat(" 1st-Ord: first-ordered \n")
  cat("     Abs: absorption    \n")
  cat("     w/o: without       \n")
  cat("    Tlag: lag-time      \n")
  cat("************************\n\n")
  cat("\n")
  file.menu <- c("IV-Bolus, MD", 
                 "IV-Infusion, MD",
                 "Extravascular, MD, 1st-Ord Abs w/o Tlag",
                 "Go Back One Upper Level",
                 "Go Back to Top Menu")
  pick <- menu(file.menu, title = "<< 2-Compartment PK Model >>")
  if (pick ==1){
     cat("\n\n")
     sbolus2(MD=TRUE)
  }
  else if (pick == 2){
     cat("\n\n") 
     sinfu2(MD=TRUE)
  }
  else if (pick == 3){
     cat("\n\n") 
     sfirst2(MD=TRUE)
  }
  else if (pick == 4){
     cat("\n\n") 
     PK.sim.MD()
  }  
  else if (pick == 5){
     cat("\n\n") 
     run()
  }  
}
