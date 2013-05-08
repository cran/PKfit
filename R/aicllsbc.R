#Estimate model fitting: AIC, Log likelihood, SBC
aicllsbc <- function(fm)
{   
  cat("\n") 
  cat("<< Akaike's Information Criterion (AIC) >>\n")
  show(AIC(fm))
    
  cat("\n<< Log likelihood >>\n")
  show(logLik(fm))
    
  cat("\n<< Schwarz's Bayesian Criterion (SBC/BIC) >>\n")
  show(BIC(fm))
  cat("\n")     
    
 #Summary the results of nls
  print(summary(fm))
  cat("\n")   
  cat(date(),"\n")     
}  
