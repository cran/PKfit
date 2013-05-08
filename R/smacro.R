smacro <- function()
{
  file.menu <- c("1-Exponential Term", 
                 "2-Exponential Term", 
                 "3-Exponential Term",
                 "Go Back One Upper Level")
  pick <- menu(file.menu, title = "<< Macroconstant Exponential Functions >>")
  if (pick ==1){
     cat("\n\n")
     smacro.one()
  }
  else if (pick == 2){
     cat("\n\n") 
     smacro.two()
  }
  else if (pick == 3){
     cat("\n\n") 
     smacro.three()
  }
  else if (pick == 4){
     cat("\n\n") 
     PK.sim()
  }
}
