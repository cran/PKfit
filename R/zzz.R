# change the following line since R v.2.15.3 as .onAttach [2013/3/4 AM06:27:15]

.onAttach <- function(lib, pkg)  {

# echo output to screen
packageStartupMessage("

.............. PKfit v1.2.5 ............

   Please type 'run()' to get started.                
                                                         
   To run demos, type 'demo(iv.bolus)', 
   'demo(mcsim)' or 'demo(mmpk)'

........................................\n")
}
