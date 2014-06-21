PK.sim <- function()
{
  OutputFilez() ### reset all output file names when running PK.sim()
  file.menu <- c("Single-Dose", 
                 "Multiple-Dose",
                 "Macroconstant Exponential Functions",
                 "Go Back One Upper Level")
  pick <- menu(file.menu, title = "<< Selection of Dose Type >>")
  if (pick == 1){
     cat("\n\n")  
     PK.sim.SD()
  }
  else if (pick == 2){
     cat("\n\n")
     PK.sim.MD()
  }
  else if (pick == 3){
     cat("\n\n")
     smacro()
  }
  else if (pick == 4){
     cat("\n\n")
     run()
  }
}