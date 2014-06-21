#PK model option
PK.fit <- function(PKindex)
{
  OutputFilez()   ### reset all output file names when running PK.fit()
  
  file.menu <- c("One-Compartment PK Model", 
                 "Two-Compartment PK Model",
                 "Three-Compartment PK Model (iv bolus with 1st-ordered elim.)",
                 "Three-Compartment PK Model (iv bolus with MM elim.)",
                 "Macroconstant Exponential Functions",
                 "Go Back One Upper Level")
  pick <- menu(file.menu, title = "<< Selection of a PK Model >>")
  if (pick== 1){
     cat("\n\n")  
     one.list(PKindex)
  }     
  else if (pick == 2){
     cat("\n\n")
     two.list(PKindex)
  }
  else if (pick == 3){
     cat("\n\n")
     fbolus3(PKindex)
  }
  else if (pick == 4){
     cat("\n\n")
     fbolus3.mm(PKindex)
  }
  else if (pick == 5){
     cat("\n\n")
     macro(PKindex)
  }
  else if (pick == 6){
     cat("\n\n")
     run()
  }        
}
