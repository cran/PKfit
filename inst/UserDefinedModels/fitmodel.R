###
### please copy this file to your working path; use getwd() to find out where it is.
### v0.2.1
###
require(PKfit)
require(minpack.lm)
require(deSolve)
require(optimx)
graphics.off()
cat("\n\n")
options(warn=-1)
modfun<-NULL
conc<-NULL

###
### x and y axis labeling;
###
xaxis<- 'Time after dosing (hr)'
yaxis<- 'Drug X plasma conc. (ng.mL)'

###
### enter your data as time vs. conc.;
### default is only for single subject only;
### for multiple subjects, you have to change this.
### same as 'Dose'.
###
PKindex<-data.frame(Subject=c(1),time=c(1,2,3,4,6,10,12),
                    conc=c(14.94,13.73,10.55,8.16,5.21,3.19,2.62))
Dose<-500

### define your model with following ODEs;
### default is a iv, bolus, one-compartment model
### with 1st-order elimination rate constant 'kel';
###
defun<- function(time, y, parms) { 
      dCpdt <- -parms["kel"] * y[1] 
      list(dCpdt) 
} 

###
### define the numerical solution for ODEs here
###
modfun <<- function(time,kel, Vd) {  
      out <- lsoda(Dose/Vd,c(0,time),defun,parms=c(kel=kel,Vd=Vd),
                   rtol=1e-6,atol=1e-6) 
      out[-1,2] 
}

###
### define the objective function here to be optimized or minimized
### with general-purposed optim(); as well as weighting schemes
###
objfun <- function(par) {
        out <- modfun(PKindex$time, par[1], par[2])
        gift <- which( PKindex$conc != 0)
        sum((PKindex$conc[gift]-out[gift])^2)                       ### equal weight
        ### sum((PKindex$conc[gift]-out[gift])^2/PKindex$conc[gift])   ### 1/conc weight
        ### sum(((PKindex$conc[gift]-out[gift])/PKindex$conc[gift])^2) ### 1/conc^2 weight
}        

###
### applying a general-purposed optimization first;
### where c(0.2,10) contains the initial values for kel, Vd, respectively;
###     
opt<-optimx(c(2,20),objfun,method="Nelder-Mead",control=list(maxit=5000))
print(opt)
nameopt<-c("kel","Vd")
### outopt<-c(opt$par[1],opt$par[2])
outopt<-c(opt$p1,opt$p2)

###
### reset parameter values here if less than zero obtained from optim()...
###
###  if(opt$par[1]<0) {opt$par[1]<-0.001}
###  if(opt$par[2]<0) {opt$par[2]<-0.001}
 
  if(opt$p1<0) {opt$p1<-0.001}
  if(opt$p2<0) {opt$p2<-0.001}

###
### pass initial estimate obtained from optim() to final nonlinear regression algorithm, nlsLM()
###
### fm<-nlsLM(conc ~ modfun(time, kel, Vd),data=PKindex,start=list(kel=opt$par[1],Vd=opt$par[2]),
###         control=nls.lm.control(maxiter=500),weights=(1/conc^2))

fm<-nlsLM(conc ~ modfun(time, kel, Vd),data=PKindex,start=list(kel=opt$p1,Vd=opt$p2),
         control=nls.lm.control(maxiter=500),weights=(1/conc^2))

### obtain coefficient for nlsLM()
coef<-data.frame(coef(fm)["kel"])

### i = # of subj; do not change this unless you know what you are doing.
### pick = 1, 2, 3 -->(1) equal, (2) 1/conc, and (3) 1/conc^2 weighting, respectively; change it if necessary.
i<-1; pick<- 3
cat(" ********************************\n\n")
cat("      --- Subject:- #",i,"---    \n\n")
cat(" ********************************\n\n")
cat("--- input data ---\n")
conc<-PKindex$conc[PKindex$Subject==i]
time<-PKindex$time[PKindex$Subject==i]
this_subj<-data.frame(time, conc)
show(this_subj);cat("\n")     # show input data    
cat("--- weighting scheme: ")
### show weighting scheme
##  cat("equal weight");cat("\n\n")
##  cat("1/Cp");cat("\n\n")
cat("1/Cp^2");cat("\n\n")
cat("--- the model: a one-compartment, iv bolus with 1st-ordered elim.\n\n") ### <-  modify this if required.
cat("<< PK parameter obtained from Nelder-Mead Simplex algorithm >>\n\n")
print(data.frame(Parameter=nameopt,Value=outopt));cat("\n")                        
plotting.lin(PKindex, fm, i, pick, coef, xaxis, yaxis)
       