"# NakaEstim" 

This code is used for implementing the simulation part of 

   "Classical and Non-Informative Bayesian Inference of spmk for
    Nakagami Distribution Based On First-Failure Progressively" (under review)
    
   by  Sanku Dey, Devendra Kumar and Riaydh Al-Mosawi

This code is written by :Riyadh Al-Mosawi

To use this code first doanload the package using

library(devtools)
    
install_github("RiyadhAl-Mosawi/NakaEstim")

Then use the main function r "Estimate" included in "ParaEstim.r" file with identifying the following parameters:

    para,  # True parameter
    
    n,     # sample size
    
    m,     # No. of failure times
    
    k,     # failure time
    
    L,     # Lower specifiaction limit
    
    U,     # Upper specifiaction limit
    
    T,     # Target value
    
    MC.size  =11000,  # Size of MCMC sample   
    
    MC.burn  =1000,   # Burn-in sample
    
    q.val=-0.5, # genarlized entropy loss parameter
    
    c.val=0.5,  # linex loss parameter
    
    Sim.no,# No. pof replications
    
    randomseed=c(2021,6,30), #seed nomber
    
    direc   # directory saving the results

For example:

Estimate(c(1,1.5),20,10,2,0.4,10,1,1100,100,-0.5,0.5,2,c(2021,6,30),"C:/Users/hpPC/Desktop")


