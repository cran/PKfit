#------------Simulation menu----------------
sone.iv.route.SD <- function()
{
  OutputFilez() ### reset all output file names when running PK.sim()
  cat("***************************************\n")
  cat(" SD: Single-Dose                       \n")
  cat(" MM Elim: Michaelis-Menten Elimination \n")
  cat("***************************************\n\n")
  cat("\n")
  file.menu <- c("IV-Bolus & SD", 
                 "IV-Bolus, SD & MM Elim", 
                 "IV-Infusion with 1st-ordered elim.",
                 "IV-Infusion with MM Elim",
                 "Go Back One Upper Level",
                 "Go Back to Top Menu")
  pick <- menu(file.menu, title = "<< IV (Bolus or Infusion) Route >>")
  if (pick ==1 ){
     cat("\n\n")
     sbolus1(MMe=FALSE,MD=FALSE)
  }
  else if (pick == 2){
     cat("\n\n") 
     sbolus1(MMe=TRUE,MD=FALSE)
  }
  else if (pick == 3){
     cat("\n\n") 
     sinfu1(MMe=FALSE,MD=FALSE)
  }
  else if (pick == 4){
     cat("\n\n") 
     sinfu1(MMe=TRUE,MD=FALSE)
  }
  else if (pick == 5){
     cat("\n\n") 
     PK.sim.SD()
  }
  else if (pick == 6){
     cat("\n\n") 
     run()
  }
}
