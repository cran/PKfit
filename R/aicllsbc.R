#Estimate model fitting: AIC, Log likelihood, SBC
aicllsbc <- function(fm)
{   

  #Summary the results of nls
  print(summary(fm))
  ModelSelect<-data.frame(Model_Select=c("AIC","Log Likelihood","SBC/BIS"),
                             Values=c(AIC(fm),logLik(fm),BIC(fm)))
  show(ModelSelect);cat("\n")
  
  cat("<< Variance-Covariance Matrix >>\n")
  print(vcov(fm))
    
  cat("\n<< weights >>\n")  ### for debugging purpose. -YJ
  print(weights(fm))
  
  ### cat("\n<< Akaike's Information Criterion (AIC) >>\n")
  ### show(AIC(fm))
  ###   
  ### cat("\n<< Log likelihood >>\n")
  ### show(logLik(fm))
  ###   
  ### cat("\n<< Schwarz's Bayesian Criterion (SBC/BIC) >>\n")
  ### show(BIC(fm));cat("\n")
}  
