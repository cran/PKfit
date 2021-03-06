noniv.route <- function(PKindex)
{
  cat("***************************************\n")
  cat(" SD: Single-Dose                       \n")
  cat(" 1st-Ord: First-Ordered                \n")
  cat(" Zero-Ord: Zero-Ordered                \n")
  cat(" Abs: Absorption                       \n")
  cat(" w: with                               \n")
  cat(" w/o: without                          \n")
  cat(" Tlag: Lag Time                        \n")
  cat(" MM Elim: Michaelis-Menten Elimination \n")
  cat("***************************************\n\n")
  cat("\n")

  file.menu <- c("Extravascular, SD & 1-Ord Abs w Tlag",
                 "Extravascular, SD & 1-Ord Abs w/o Tlag",
                 "Extravascular, SD & Zero-Ord Abs w/o Tlag",
                 "Extravascular, SD, 1-Ord Abs & MM Elim w Tlag ",
                 "Extravascular, SD, 1-Ord Abs & MM Elim w/o Tlag",
                 "Extravascular, SD, Zero-Ord Abs & MM Elim w/o Tlag",
                 "Go Back One Upper Level")
  pick <- menu(file.menu, title = "<< Non IV Route >>")
  if (pick == 1){
     cat("\n\n") 
     ffirst.lag(PKindex)
  }
  else if (pick == 2){
     cat("\n\n") 
     ffirst.nolag(PKindex)
  }
  else if (pick == 3){
     cat("\n\n") 
     fzero.nolag(PKindex)
  }
  else if (pick == 4){
     cat("\n\n") 
     ffirst.lagm(PKindex)
  }
  else if (pick == 5){
     cat("\n\n") 
     ffirst.nolagm(PKindex)
  }
  else if (pick == 6){
     cat("\n\n") 
     fzero.nolagm(PKindex)
  }
  else if (pick == 7){
     cat("\n\n") 
     one.list(PKindex)
  }
}
