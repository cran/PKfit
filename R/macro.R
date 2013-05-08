macro <- function(PKindex)
{
  file.menu <- c("1-Exponential Term", 
                 "2-Exponential Term",
                 "3-Exponential Term",
                 "Go Back One Upper Level")
  pick <- menu(file.menu, title = "<< Macroconstant Exponential Functions >>")
  if (pick ==1){
     cat("\n\n")
     fmacro.one(PKindex)
  }
  else if (pick == 2){
     cat("\n\n") 
     fmacro.two(PKindex)
  }
  else if (pick == 3){
     cat("\n\n") 
     fmacro.three(PKindex)
  }
  else if (pick == 4){
     cat("\n\n") 
     PK.fit(PKindex)
  }
}
