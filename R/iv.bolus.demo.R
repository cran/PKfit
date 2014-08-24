iv.bolus.demo<-function(){
cat("\n\n")
options(warn=-1)
modfun<-NULL
conc<-NULL
###
### x and y axis labelling;
###
xaxis<- 'Time after dosing (hr)'
yaxis<- 'Drug X plasma conc. (ng.mL)'

PKindex<-data.frame(Subject=c(1),time=c(1,2,3,4,6,10,12),
                    conc=c(14.94,13.73,10.55,8.16,5.21,3.19,2.62))
Dose<-500
defun<- function(time, y, parms) { 
      dCpdt <- -parms["kel"] * y[1] 
      list(dCpdt) 
} 
    
modfun <<- function(time,kel, Vd) {  
      out <- lsoda(Dose/Vd,c(0,time),defun,parms=c(kel=kel,Vd=Vd),
                   rtol=1e-6,atol=1e-6) 
      out[-1,2] 
}

objfun <- function(par) {
        out <- modfun(PKindex$time, par[1], par[2])
        gift <- which( PKindex$conc != 0 )
        ### sum((PKindex$conc[gift]-out[gift])^2)
        sum(((PKindex$conc[gift]-out[gift])/PKindex$conc[gift])^2)
}        
     
opt<-optim(c(0.21,10),objfun,method="Nelder-Mead",control=list(maxit=5000))  
nameopt<-c("kel","Vd")
outopt<-c(opt$par[1],opt$par[2])


  if(opt$par[1]<0) {opt$par[1]<-0.01}
  if(opt$par[2]<0) {opt$par[2]<-0.01}

fm<-nlsLM(conc ~ modfun(time, kel, Vd),data=PKindex,start=list(kel=opt$par[1],Vd=opt$par[2]),
         control=nls.lm.control(maxiter=500),weights=(1/conc^2)) ### lower of Vd should not be zero due to Dose/Vd. --YJ
        
coef<-data.frame(coef(fm)["kel"])

### i = # of subj;
### pick = 1, 2, 3 --> equal, 1/conc, and 1/conc^2 weighting, respectively;
### change it if necessary.
i<-1; pick<- 3
cat(" ********************************\n\n")
cat("      --- Subject:- #",i,"---    \n\n")
cat(" ********************************\n\n")
cat("--- input data ---\n")
conc<-PKindex$conc[PKindex$Subject==i]
time<-PKindex$time[PKindex$Subject==i]
this_subj<-data.frame(time, conc)
show(this_subj);cat("\n")     # show input data 

### show weighting scheme   
cat("--- weighting scheme: ")
##  cat("equal weight");cat("\n\n")
##  cat("1/Cp");cat("\n\n")
cat("1/Cp^2");cat("\n\n")
cat("--- model selection: a one-compartment, iv bolus pk model\n    with 1st-ordered elim.\n\n") 
cat("<< PK parameter obtained from Nelder-Mead Simplex algorithm >>\n\n")
print(data.frame(Parameter=nameopt,Value=outopt));cat("\n")                        
plotting.lin(PKindex, fm, i, pick, coef, xaxis, yaxis)
}
       