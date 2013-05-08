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
