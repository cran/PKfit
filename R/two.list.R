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
  pick <- menu(file.menu, title = "<< 2-Compartment PK Model >>")
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
