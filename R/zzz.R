# change the following line since R v.2.15.3 as .onAttach [2013/3/4 AM06:27:15]

.onAttach <- function(lib, pkg)  {

# echo output to screen
packageStartupMessage("

....................... PKfit v1.2.1 .....................

   Please type 'run()' to get started.                
                                                         
   If you want to see some demos, type                   
                                                         
   'demo(iv.bolus)' for Normal Fitting,                 
                                                         
   'demo(mcsim)' for Monte-Carlo Simulation, or
                                                         
   'demo(mmpk)' for Michaelis-Menten Model.              

..........................................................")
}
