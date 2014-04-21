entertitle<-function()
{
  cat("\n Enter the title of x-axis (time): \n")
  cat("(or press Enter to use default)\n\n") 
  xaxis<-readline()
  if (substr(xaxis, 1, 1) == "")  xaxis<-"Time after dosing (hr)"  else xaxis<-xaxis
  cat("\n Enter the title of y-axis(Cp): \n")
  cat("(or press Enter to use default)\n\n") 
  yaxis<-readline()
  #cat("\n\n Please Wait.  Data is Processing. \n")
  if (substr(yaxis, 1, 1) == "")  yaxis<-"Drug plasma conc."  else yaxis<-yaxis
  #cat("\n\n Please Wait.  Data is Processing. \n")
  return(list(xaxis=xaxis,yaxis=yaxis))
}