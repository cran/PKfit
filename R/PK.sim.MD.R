PK.sim.MD <- function()
{
  OutputFilez() ### reset all output file names when running PK.sim()
  file.menu <- c("One-Compartment PK Model: IV Route", 
                 "One-Compartment PK Model: Non IV Route",
                 "Two-Compartment PK Model",
                 "*** Three-Compartment PK Model (not avaiable)",
                 "Go Back One Upper Level")
  pick <- menu(file.menu, title = "<< Selection of PK Model >>")
  if (pick == 1){
     cat("\n\n")  
     sone.iv.route.MD()
  }
  else if (pick == 2){
     cat("\n\n")
     sone.noniv.route.MD()
  }
  else if (pick == 3){
     cat("\n\n")
     stwo.MD.all()
  }      
  else if (pick == 4){
     cat("\n\n")
     ### sthree.MD.all()
     readline(" This function is not available yet. Press any key to retry.")
     PK.sim.MD()
  }
  else if (pick == 5){
     cat("\n\n")
     PK.sim()
  }
}