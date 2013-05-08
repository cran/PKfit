PK.sim <- function()
{
  file.menu <- c("One-Compartment PK Model: IV Route", 
                 "One-Compartment PK Model: Non IV Route",
                 "Two-Compartment PK Model",
                 "Macroconstant Exponential Functions",
                 "Go Back One Upper Level")
  pick <- menu(file.menu, title = "<< Selection of PK Model >>")
  if (pick == 1){
     cat("\n\n")  
     sone.iv.route()
  }
  else if (pick == 2){
     cat("\n\n")
     sone.noniv.route()
  }
  else if (pick == 3){
     cat("\n\n")
     stwo.all()
  }      
  else if (pick == 4){
     cat("\n\n")
     smacro()
  }
  else if (pick == 5){
     cat("\n\n")
     run()
  }
}