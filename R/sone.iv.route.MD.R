#------------Simulation menu----------------
sone.iv.route.MD <- function()
{
  OutputFilez() ### reset all output file names when running PK.sim()
  cat("***************************************\n")
  cat(" MD: Multiple-Dose                       \n")
  cat(" MM Elim: Michaelis-Menten Elimination \n")
  cat("***************************************\n\n")
  cat("\n")
  file.menu <- c("IV-Bolus & MD", 
                 "IV-Bolus, MD & MM Elim", 
                 "intermittent IV-Infusion & MD",           ### continuous iv infuion does not have so-called 'multiple-dose' at all.
                 "intermittent IV-Infusion, MD & MM Elim",  ### continuous iv infuion does not have so-called 'multiple-dose' at all.
                 "Go Back One Upper Level",
                 "Go Back to Top Menu")
  pick <- menu(file.menu, title = "<< IV (Bolus or Infusion) Route >>")
  if (pick ==1 ){
     cat("\n\n")
     sbolus1(MMe=FALSE,MD=TRUE)
  }
  else if (pick == 2){
     cat("\n\n") 
     sbolus1(MMe=TRUE,MD=TRUE)
  }
   else if (pick == 3){
     cat("\n\n") 
     sinfu1(MMe=FALSE,MD=TRUE)
  }
   else if (pick == 4){
     cat("\n\n") 
     sinfu1(MMe=TRUE,MD=TRUE)
  }
  else if (pick == 5){
     cat("\n\n") 
     PK.sim.MD()
  }
  else if (pick == 6){
     cat("\n\n") 
     run()
  }
}
