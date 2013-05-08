#PK model option
PK.fit <- function(PKindex)
{
  file.menu <- c("One-Compartment PK Models", 
                 "Two-Compartment PK Models",
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
     macro(PKindex)
  }
  else if (pick == 4){
     cat("\n\n")
     run()
  }        
}