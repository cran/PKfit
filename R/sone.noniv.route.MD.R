sone.noniv.route.MD <- function()
{
  OutputFilez() ### reset all output file names when running PK.sim()
  
  cat("***************************************\n")
  cat(" MD: Multiple-Dose                     \n")
  cat(" 1st-Ord: First-Ordered                \n")
  cat(" Zero-Ord: Zero-Ordered                \n")
  cat(" Abs: Absorption                       \n")
  cat(" w: with                               \n")
  cat(" w/o: without                          \n")
  cat(" Tlag: Lag Time                        \n")
  cat(" MM Elim: Michaelis-Menten Elimination \n")
  cat("***************************************\n\n")
  cat("\n")
  
  file.menu <- c("Extravascular, MD & 1st-Ord Abs w Tlag",
                 "Extravascular, MD & 1st-Ord Abs w/o Tlag",
                 "Extravascular, MD, 1st-Ord Abs. & MM Elim w Tlag",
                 "Extravascular, MD, 1st-Ord Abs. & MM Elim w/o Tlag",
                 "Extravascular, MD, Zero-Ord Abs w/o Tlag",
                 "Extravascular, MD, Zero-Ord Abs. & MM Elim w/o Tlag",
                 "Go Back One Upper Level",
                 "Go Back to Top Menu")
  pick <- menu(file.menu, title = "<< Non IV Route - Multiple-Dose >>")
  if (pick == 1){
     cat("\n\n") 
     sfirst.nolag(Tlag=TRUE,MMe=FALSE,MD=TRUE)
  }
  else if (pick == 2){
     cat("\n\n") 
     sfirst.nolag(Tlag=FALSE,MMe=FALSE,MD=TRUE)
  }
  else if (pick == 3){
     cat("\n\n") 
     sfirst.nolag(Tlag=TRUE,MMe=TRUE,MD=TRUE)
  }
  else if (pick == 4){
     cat("\n\n") 
     sfirst.nolag(Tlag=FALSE,MMe=TRUE,MD=TRUE)
  }
  else if (pick == 5){
     cat("\n\n") 
     szero.nolag(MMe=FALSE,MD=TRUE)
  }
  else if (pick == 6){
     cat("\n\n") 
     szero.nolag(MMe=TRUE,MD=TRUE)
  }
  else if (pick == 7){
     cat("\n\n") 
     PK.sim.MD()
  }
  else if (pick == 8){
     cat("\n\n") 
     run()
  }
}
