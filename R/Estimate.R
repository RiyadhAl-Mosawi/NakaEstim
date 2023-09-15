Estimate=function(
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
    ){
  
  cat("\014")
  library(Rcpp)
  library(pracma)
  library(roptim)
  library(VGAM)
  library(coda)
  
  setwd(direc)

  MSE.fun=function(eta,eta.hat, type="SEL",q){
    if(type=="SEL") return(mean(apply(as.matrix(1:length(eta.hat)),1,function(i) (eta.hat[i]-eta)^2)))
    if(type=="GEL") return(mean(apply(as.matrix(1:length(eta.hat)),1,function(i) (eta.hat[i]/eta)^q-q*log(eta.hat[i]/eta)-1)))
    if(type=="LINEX") return(mean(apply(as.matrix(1:length(eta.hat)),1,function(i) exp(q*(eta.hat[i]-eta))-q*(eta.hat[i]-eta)-1)))
  }
  
  options(width=100,length=200,max.print = 10000,digits = 3, scipen=999)

  ## Censroing schemes
  Scheme=matrix(nr=3,nc=m,0)
  Scheme[1,]=c(rep(0,m-1),n-m)
  Scheme[2,]=c(n-m,rep(0,m-1))
  Scheme[3,]=c(rep(1,n-m),rep(0,2*m-n))
  
  eps=para[1]
  eta=para[2]
  
  spmk=spmk_fun(para,L,U,T,mu(para[1],para[2]),sig(para[1],para[2]))
  for(i1 in 1:3){
    R=Scheme[i1,]  
    it=1
    set.seed(randomseed)
    while (it<=Sim.no){
      if(it==1){
        MLE_Res=MPS_Res=NULL
        MCMC_MH_MLE_eps=MCMC_MH_MLE_eta=MCMC_MH_MLE_spmk=NULL
        MCMC_MH_MPS_eps=MCMC_MH_MPS_eta=MCMC_MH_MPS_spmk=NULL
        start.time=date()
      }  
      
      X=try(Gen_Data(para,R,k),silent=TRUE)
      if(is.character(X)) next
      
      
      
      #################################################################
      #------------- MLE using NR method
      ############s#####################################################
      mle.res=suppressWarnings(try(Estim(para,X,R,k,"LF","L-BFGS-B",c(0.5,0),c(3*para[1],3*para[2])),silent=TRUE))
      # Check convergence
      if(is.character(mle.res)) next
      
      
      mle.est=as.numeric(mle.res$par)
      mle.ese=sqrt(abs(diag(solve(mle.res$inform))))
      
      mle.aci=matrix(nr=2,nc=3,NA)
      mle.aci[1,1:2]=c(max(0,mle.est[1]-1.96*mle.ese[1]),mle.est[1]+1.96*mle.ese[1])
      mle.aci[2,1:2]=c(max(0,mle.est[2]-1.96*mle.ese[2]),mle.est[2]+1.96*mle.ese[2])
      mle.aci[1,3]=ifelse((mle.aci[1,2]>=para[1])&(mle.aci[1,1]<=para[1]),1,0)
      mle.aci[2,3]=ifelse((mle.aci[2,2]>=para[2])&(mle.aci[2,1]<=para[2]),1,0)
      
      mu.hat=try(mu(mle.est[1],mle.est[2]),silent=TRUE)
      if(is.character(mu.hat)) mu.hat=(gamma(mle.est[1]+0.5)/gamma(mle.est[1]))*sqrt(mle.est[2]/mle.est[1])
      sig.hat=try(sig(mle.est[1],mle.est[2]),silent=TRUE)
      if(is.character(sig.hat)) sig.hat=sqrt(mle.est[2]-((gamma(mle.est[1]+0.5)/gamma(mle.est[1]))*sqrt(mle.est[2]/mle.est[1]))^2)
      
      spmk.est.mle=spmk_fun(mle.est,L,U,T,mu.hat,sig.hat)
      spmk.ese.mle=sqrt(abs(as.vector(spmk_grad(mle.est,L,U,T,mu.hat,sig.hat))%*%solve(mle.res$inform)%*%as.vector(spmk_grad(mle.est,L,U,T,mu.hat,sig.hat))))
      spmk.aci.mle=numeric(3)
      spmk.aci.mle[1:2]=c(max(0,spmk.est.mle-1.96*spmk.ese.mle),spmk.est.mle+1.96*spmk.ese.mle)
      spmk.aci.mle[3]=ifelse((spmk.aci.mle[2]>=spmk)&(spmk.aci.mle[1]<=spmk),1,0)
      
      #################################################################
      #------------- MPS using NR method
      ############s#####################################################
      mps.res=suppressWarnings(try(Estim(para,X,R,k,"SPF","L-BFGS-B",c(0.5,0),c(3*para[1],3*para[2])),silent=TRUE))
      # Check convergence
      if(is.character(mps.res)) next
      
      mps.est=as.numeric(mps.res$par)
      mps.ese=sqrt(abs(diag(-solve(mps.res$inform))))
      mps.aci=matrix(nr=2,nc=3,NA)
      mps.aci[1,1:2]=c(max(0,mps.est[1]-1.96*mps.ese[1]),mps.est[1]+1.96*mps.ese[1])
      mps.aci[2,1:2]=c(max(0,mps.est[2]-1.96*mps.ese[2]),mps.est[2]+1.96*mps.ese[2])
      mps.aci[1,3]=ifelse((mps.aci[1,2]>=para[1])&(mps.aci[1,1]<=para[1]),1,0)
      mps.aci[2,3]=ifelse((mps.aci[2,2]>=para[2])&(mps.aci[2,1]<=para[2]),1,0)
      
      mu.hat=try(mu(mps.est[1],mps.est[2]),silent=TRUE)
      if(is.character(mu.hat)) mu.hat=(gamma(mps.est[1]+0.5)/gamma(mps.est[1]))*sqrt(mps.est[2]/mps.est[1])
      sig.hat=try(sig(mps.est[1],mps.est[2]),silent=TRUE)
      if(is.character(sig.hat)) sig.hat=sqrt(mps.est[2]-((gamma(mps.est[1]+0.5)/gamma(mps.est[1]))*sqrt(mps.est[2]/mps.est[1]))^2)
      
      spmk.est.mps=spmk_fun(mps.est,L,U,T,mu.hat,sig.hat)
      spmk.ese.mps=sqrt(abs(as.vector(spmk_grad(mps.est,L,U,T,mu.hat,sig.hat))%*%solve(mps.res$inform)%*%as.vector(spmk_grad(mps.est,L,U,T,mu.hat,sig.hat))))
      spmk.aci.mps=numeric(3)
      spmk.aci.mps[1:2]=c(max(0,spmk.est.mps-1.96*spmk.ese.mps),spmk.est.mps+1.96*spmk.ese.mps)
      spmk.aci.mps[3]=ifelse((spmk.aci.mps[2]>=spmk)&(spmk.aci.mps[1]<=spmk),1,0)
      

      cat("\n Generating MH samples based on likelihood function...\n");
      Bayes_MH_MLE=MH_sample(type="LF",mle.est,mle.ese,R,X,k,L,U,T,MC.size,MC.burn,q=q.val,c=c.val,verbose = F, display_progress=TRUE)
      
      if(is.infinite(Bayes_MH_MLE$est_SEL[1])) next
      if(is.infinite(Bayes_MH_MLE$est_SEL[2])) next
      if(is.infinite(Bayes_MH_MLE$est_SEL[3])) next
      if(is.infinite(Bayes_MH_MLE$est_GEL[1])) next
      if(is.infinite(Bayes_MH_MLE$est_GEL[2])) next
      if(is.infinite(Bayes_MH_MLE$est_GEL[3])) next
      if(is.infinite(Bayes_MH_MLE$est_LINEX[1])) next
      if(is.infinite(Bayes_MH_MLE$est_LINEX[2])) next
      if(is.infinite(Bayes_MH_MLE$est_LINEX[3])) next
      
      
      cat("\n Generating MH samples based on product of spacing function...\n");
      Bayes_MH_MPS=MH_sample(type="PSF",mle.est,mle.ese,R,X,k,L,U,T,MC.size,MC.burn,q=q.val,c=c.val,verbose = F, display_progress=TRUE)
      
      if(is.infinite(Bayes_MH_MPS$est_SEL[1])|is.nan(Bayes_MH_MPS$est_SEL[1])|is.na(Bayes_MH_MPS$est_SEL[1])|Bayes_MH_MPS$est_SEL[1]>3*eps) next
      if(is.infinite(Bayes_MH_MPS$est_SEL[2])|is.nan(Bayes_MH_MPS$est_SEL[2])|is.na(Bayes_MH_MPS$est_SEL[2])|Bayes_MH_MPS$est_SEL[2]>3*eta) next
      if(is.infinite(Bayes_MH_MPS$est_SEL[3])|is.nan(Bayes_MH_MPS$est_SEL[3])|is.na(Bayes_MH_MPS$est_SEL[3])|Bayes_MH_MPS$est_SEL[3]>3*spmk) next
      if(is.infinite(Bayes_MH_MPS$est_GEL[1])|is.nan(Bayes_MH_MPS$est_GEL[1])|is.na(Bayes_MH_MPS$est_GEL[1])|Bayes_MH_MPS$est_GEL[1]>3*eps) next 
      if(is.infinite(Bayes_MH_MPS$est_GEL[2])|is.nan(Bayes_MH_MPS$est_GEL[2])|is.na(Bayes_MH_MPS$est_GEL[2])|Bayes_MH_MPS$est_GEL[2]>3*eta) next 
      if(is.infinite(Bayes_MH_MPS$est_GEL[3])|is.nan(Bayes_MH_MPS$est_GEL[3])|is.na(Bayes_MH_MPS$est_GEL[3])|Bayes_MH_MPS$est_GEL[3]>3*spmk) next 
      if(is.infinite(Bayes_MH_MPS$est_LINEX[1])|is.nan(Bayes_MH_MPS$est_LINEX[1])|is.na(Bayes_MH_MPS$est_LINEX[1])|Bayes_MH_MPS$est_LINEX[1]>3*eps) next 
      if(is.infinite(Bayes_MH_MPS$est_LINEX[2])|is.nan(Bayes_MH_MPS$est_LINEX[2])|is.na(Bayes_MH_MPS$est_LINEX[2])|Bayes_MH_MPS$est_LINEX[2]>3*eta) next 
      if(is.infinite(Bayes_MH_MPS$est_LINEX[3])|is.nan(Bayes_MH_MPS$est_LINEX[3])|is.na(Bayes_MH_MPS$est_LINEX[3])|Bayes_MH_MPS$est_LINEX[3]>3*spmk) next 
      
      
      MLE_Res =rbind(MLE_Res,c(mle.est[1],mle.ese[1],mle.aci[1,],
                       mle.est[2],mle.ese[2],mle.aci[2,],
                       spmk.est.mle,spmk.ese.mle,spmk.aci.mle))
      MPS_Res =rbind(MPS_Res,c(mps.est[1],mps.ese[1],mps.aci[1,],
                               mps.est[2],mps.ese[2],mps.aci[2,],
                               spmk.est.mps,spmk.ese.mps,spmk.aci.mps))
      
      MCMC_MH_MLE_eps =rbind(MCMC_MH_MLE_eps ,c(Bayes_MH_MLE$est_SEL[1],
                             Bayes_MH_MLE$est_GEL[1],Bayes_MH_MLE$est_LINEX[1],
                             Bayes_MH_MLE$HPD_eps ,(Bayes_MH_MLE$HPD_eps[1]<eps)*(eps<Bayes_MH_MLE$HPD_eps[2])))
      MCMC_MH_MLE_eta =rbind(MCMC_MH_MLE_eta ,c(Bayes_MH_MLE$est_SEL[2],
                             Bayes_MH_MLE$est_GEL[2],Bayes_MH_MLE$est_LINEX[2],
                             Bayes_MH_MLE$HPD_eta ,(Bayes_MH_MLE$HPD_eta[1]<eta)*(eta<Bayes_MH_MLE$HPD_eta[2])))
      MCMC_MH_MLE_spmk=rbind(MCMC_MH_MLE_spmk,c(Bayes_MH_MLE$est_SEL[3],
                             Bayes_MH_MLE$est_GEL[3],Bayes_MH_MLE$est_LINEX[3],
                             Bayes_MH_MLE$HPD_spmk,(Bayes_MH_MLE$HPD_spmk[1]<spmk)*(spmk<Bayes_MH_MLE$HPD_spmk[2])))
      MCMC_MH_MPS_eps =rbind(MCMC_MH_MPS_eps ,c(Bayes_MH_MPS$est_SEL[1],Bayes_MH_MPS$est_GEL[1],
                             Bayes_MH_MPS$est_LINEX[1],Bayes_MH_MPS$HPD_eps ,
                             (Bayes_MH_MPS$HPD_eps[1]<eps)*(eps<Bayes_MH_MPS$HPD_eps[2])))
      MCMC_MH_MPS_eta =rbind(MCMC_MH_MPS_eta ,c(Bayes_MH_MPS$est_SEL[2],Bayes_MH_MPS$est_GEL[2],
                             Bayes_MH_MPS$est_LINEX[2],Bayes_MH_MPS$HPD_eta ,
                             (Bayes_MH_MPS$HPD_eta[1]<eta)*(eta<Bayes_MH_MPS$HPD_eta[2])))
      MCMC_MH_MPS_spmk=rbind(MCMC_MH_MPS_spmk,c(Bayes_MH_MPS$est_SEL[3],Bayes_MH_MPS$est_GEL[3],
                              Bayes_MH_MPS$est_LINEX[3],Bayes_MH_MPS$HPD_spmk,
                              (Bayes_MH_MPS$HPD_spmk[1]<spmk)*(spmk<Bayes_MH_MPS$HPD_spmk[2])))

      cat("============= ",it," ======================\n")
      MLE_res=Bay_res=NULL
      MLE_res=round(rbind(c(para[1],rep(NA,4),para[2],rep(NA,4),spmk,rep(NA,4)),
                          colMeans(MLE_Res,na.rm = TRUE),
                          colMeans(MPS_Res,na.rm = TRUE)),3)
      
      colnames(MLE_res)=c("est.eps","ese.eps","aci.eps"," ","cp.eps","est.eta","ese.eta",
                          "aci.eta"," ","cp.eta","est.spmk","ese.spmk","aci.spmk"," ","cp.spmk")
      rownames(MLE_res)=c("True", "MLE","MPS")
      cat("\n======= MLE and MPS Results from iteration 1 to iteration ",it,"============\n")
      print(MLE_res)
      
      eps_mh_res=round(c(colMeans(MCMC_MH_MLE_eps,na.rm=T),colMeans(MCMC_MH_MPS_eps,na.rm=T)),3)
      eta_mh_res=round(c(colMeans(MCMC_MH_MLE_eta,na.rm=T),colMeans(MCMC_MH_MPS_eta,na.rm=T)),3)
      spmk_mh_res=round(c(colMeans(MCMC_MH_MLE_spmk,na.rm=T),colMeans(MCMC_MH_MPS_spmk,na.rm=T)),3)
      Bay_res=data.frame(rbind(eps_mh_res,eta_mh_res,spmk_mh_res))
      cat("\n======= MCMC Results from iteration 1 to iteration ",it,"============\n")
      
      rownames(Bay_res)=c("eps(MH)","eta(MH)","spmk(MH)")
      colnames(Bay_res)=c("Est.MLE (S)","Est.MLE(G)","Est.MLE(L)","L.MLE","U.MLE","Len.MLE","CP.MLE",
                          "Est.MPS (S)","Est.MPS(G)","Est.MPS(L)","L.MPS","U.MPS","Len.MPS","CP.MPS")
      print(round(Bay_res,3))
      
      if(it<Sim.no) {it=it+1; next}
      else{
        mle_res=c(n,m,i1,
                  mean(MLE_Res[,1],na.rm=TRUE),
                  mean(MLE_Res[,1],na.rm=TRUE)-para[1],
                  mean((MLE_Res[,1]-para[1])^2,na.rm=TRUE),
                  sd(MLE_Res[,1],na.rm=TRUE),
                  mean(MLE_Res[,2],na.rm=TRUE),
                  mean(MLE_Res[,6],na.rm=TRUE),
                  mean(MLE_Res[,6],na.rm=TRUE)-para[2],
                  mean((MLE_Res[,6]-para[2])^2,na.rm=TRUE),
                  sd(MLE_Res[,6],na.rm=TRUE),
                  mean(MLE_Res[,7],na.rm=TRUE),
                  mean(MLE_Res[,11],na.rm=TRUE),
                  mean(MLE_Res[,11],na.rm=TRUE)-spmk,
                  mean((MLE_Res[,11]-spmk)^2,na.rm=TRUE),
                  sd(MLE_Res[,11],na.rm=TRUE),
                  mean(MLE_Res[,12],na.rm=TRUE))
        
        mps_res=c(n,m,i1,
                  mean(MPS_Res[,1],na.rm=TRUE),
                  mean(MPS_Res[,1],na.rm=TRUE)-para[1],
                  mean((MPS_Res[,1]-para[1])^2,na.rm=TRUE),
                  sd(MPS_Res[,1],na.rm=TRUE),
                  mean(MPS_Res[,2],na.rm=TRUE),
                  mean(MPS_Res[,6],na.rm=TRUE),
                  mean(MPS_Res[,6],na.rm=TRUE)-para[2],
                  mean((MPS_Res[,6]-para[2])^2,na.rm=TRUE),
                  sd(MPS_Res[,6],na.rm=TRUE),
                  mean(MPS_Res[,7],na.rm=TRUE),
                  mean(MPS_Res[,11],na.rm=TRUE),
                  mean(MPS_Res[,11],na.rm=TRUE)-spmk,
                  mean((MPS_Res[,11]-spmk)^2,na.rm=TRUE),
                  sd(MPS_Res[,11],na.rm=TRUE),
                  mean(MPS_Res[,12],na.rm=TRUE))
        
        Res_est=data.frame(rbind(mle_res,mps_res))
        colnames(Res_est)=c("n","m","I","Estim.eps","Bias.eps","MSE.eps","SSE.eps",
                            "ESE.eps","Estim.eta","Bias.eta","MSE.eta","SSE.eta",
                            "ESE.eta","Estim.spmk","Bias.spmk","MSE.spmk","SSE.spmk","ESE.spmk")
        rownames(Res_est)=c("MLE","MPS")
        
        mle_aci=c(n,m,i1,
                  max(0,mean(MLE_Res[,3],na.rm=TRUE)),
                  mean(MLE_Res[,4],na.rm=TRUE),
                  mean(MLE_Res[,4],na.rm=TRUE)-max(0,mean(MLE_Res[,3],na.rm=TRUE)),
                  mean(MLE_Res[,5],na.rm=TRUE),
                  max(0,mean(MLE_Res[,8],na.rm=TRUE)),
                  mean(MLE_Res[,9],na.rm=TRUE),
                  mean(MLE_Res[,9],na.rm=TRUE)-max(0,mean(MLE_Res[,8],na.rm=TRUE)),
                  mean(MLE_Res[,10],na.rm=TRUE),
                  max(0,mean(MLE_Res[,13],na.rm=TRUE)),
                  mean(MLE_Res[,14],na.rm=TRUE),
                  mean(MLE_Res[,14],na.rm=TRUE)-max(0,mean(MLE_Res[,13],na.rm=TRUE)),
                  mean(MLE_Res[,15],na.rm=TRUE))
        mps_aci=c(n,m,i1,
                  max(0,mean(MPS_Res[,3],na.rm=TRUE)),
                  mean(MPS_Res[,4],na.rm=TRUE),
                  mean(MPS_Res[,4],na.rm=TRUE)-max(0,mean(MPS_Res[,3],na.rm=TRUE)),
                  mean(MPS_Res[,5],na.rm=TRUE),
                  max(0,mean(MPS_Res[,8],na.rm=TRUE)),
                  mean(MPS_Res[,9],na.rm=TRUE),
                  mean(MPS_Res[,9],na.rm=TRUE)-max(0,mean(MPS_Res[,8],na.rm=TRUE)),
                  mean(MPS_Res[,10],na.rm=TRUE),
                  max(0,mean(MPS_Res[,13],na.rm=TRUE)),
                  mean(MPS_Res[,14],na.rm=TRUE),
                  mean(MPS_Res[,14],na.rm=TRUE)-max(0,mean(MPS_Res[,13],na.rm=TRUE)),
                  mean(MPS_Res[,15],na.rm=TRUE))
        
        Res_aci=data.frame(rbind(mle_aci,mps_aci))
        colnames(Res_aci)=c("n","m","R","Lower.eps","Upper.eps","Lenght.eps","ECP.eps",
                            "Lower.eta","Upper.eta","Lenght.eta","ECP.eta",
                            "Lower.spmk","Upper.spmk","Lenght.spmk","ECP.spmk")
        rownames(Res_aci)=c("MLE","MPS")
        
        sel_res=gel_res=linex_res=HPD_mh_res=NULL
        
        sel_res=round(c(mean(MCMC_MH_MLE_eps[,1]-eps),MSE.fun(eps,MCMC_MH_MLE_eps[,1],"SEL",0),
                        mean(MCMC_MH_MLE_eta[,1]-eta),MSE.fun(eta,MCMC_MH_MLE_eta[,1],"SEL",0),
                        mean(MCMC_MH_MLE_spmk[,1]),MSE.fun(spmk,MCMC_MH_MLE_spmk[,1],"SEL",0),
                        mean(MCMC_MH_MPS_eps[,1]-eps),MSE.fun(eps,MCMC_MH_MPS_eps[,1],"SEL",0),
                        mean(MCMC_MH_MPS_eta[,1]-eta),MSE.fun(eta,MCMC_MH_MPS_eta[,1],"SEL",0),
                        mean(MCMC_MH_MPS_spmk[,1]),MSE.fun(spmk,MCMC_MH_MPS_spmk[,1],"SEL",0)),3)
        
        gel_res=round(c(mean(MCMC_MH_MLE_eps[,2]-eps),MSE.fun(eps,MCMC_MH_MLE_eps[,2],"GEL",q.val),
                        mean(MCMC_MH_MLE_eta[,2]-eta),MSE.fun(eta,MCMC_MH_MLE_eta[,2],"GEL",q.val),
                        mean(MCMC_MH_MLE_spmk[,2]),MSE.fun(spmk,MCMC_MH_MLE_spmk[,2],"GEL",q.val),
                        mean(MCMC_MH_MPS_eps[,2]-eps),MSE.fun(eps,MCMC_MH_MPS_eps[,2],"GEL",q.val),
                        mean(MCMC_MH_MPS_eta[,2]-eta),MSE.fun(eta,MCMC_MH_MPS_eta[,2],"GEL",q.val),
                        mean(MCMC_MH_MPS_spmk[,2]),MSE.fun(spmk,MCMC_MH_MPS_spmk[,2],"GEL",q.val)),3)
        
        linex_res=round(c(mean(MCMC_MH_MLE_eps[,3]-eps),MSE.fun(eps,MCMC_MH_MLE_eps[,3],"LINEX",c.val),
                          mean(MCMC_MH_MLE_eta[,3]-eta),MSE.fun(eta,MCMC_MH_MLE_eta[,3],"LINEX",c.val),
                          mean(MCMC_MH_MLE_spmk[,3]),MSE.fun(spmk,MCMC_MH_MLE_spmk[,3],"LINEX",c.val),
                          mean(MCMC_MH_MPS_eps[,3]-eps),MSE.fun(eps,MCMC_MH_MPS_eps[,3],"LINEX",c.val),
                          mean(MCMC_MH_MPS_eta[,3]-eta),MSE.fun(eta,MCMC_MH_MPS_eta[,3],"LINEX",c.val),
                          mean(MCMC_MH_MPS_spmk[,3]),MSE.fun(spmk,MCMC_MH_MPS_spmk[,3],"LINEX",c.val)),3)
        
        Estim=data.frame(rbind(sel_res,gel_res,linex_res))
        colnames(Estim)=c("eps1.B","eps1.MSE","eta1.B","eta1.MSE","spmk1.B","spmk1.MSE","eps2.B","eps2.MSE","eta2.B","eta2.MSE","spmk2.B","spmk2.MSE")
        rownames(Estim)=c("SEL","GEL","LINEX")
        
        HPD_mh_res=unname(c(colMeans(MCMC_MH_MLE_eps[,6:7],na.rm=TRUE),
                            colMeans(MCMC_MH_MLE_eta[,6:7],na.rm=TRUE),
                            colMeans(MCMC_MH_MLE_spmk[,6:7],na.rm=TRUE),
                            colMeans(MCMC_MH_MPS_eps[,6:7],na.rm=TRUE),
                            colMeans(MCMC_MH_MPS_eta[,6:7],na.rm=TRUE),
                            colMeans(MCMC_MH_MPS_spmk[,6:7],na.rm=TRUE)))
        HPD_res=data.frame(rbind(HPD_mh_res))
        colnames(HPD_res)=c("eps.MLE.L","eps.MLE.CP","eta.MLE.L","eta.MLE.CP","spmk.MLE.L","spmk.MLE.CP",
                            "eps.MPS.L","eps.MPS.CP","eta.MPS.L","eta.MPS.CP","spmk.MPS.L","spmk.MPS.CP")
        rownames(HPD_res)=c("MH")
        
        
        sink(paste("ND_n",n,".txt",sep=""),append = TRUE)
        cat("====================================================","\n")
        cat("True Value of Epsilon      : ",para[1],"\n")
        cat("True Value of Eta          : ",para[2],"\n")
        cat("Number of groups (n)       : ",n,"\n")
        cat("Size of each group (k)     : ",k,"\n")
        cat("Sample size (m)            : ",m,"\n")
        cat("No. of iterations          : ",Sim.no,"\n")
        cat("MCMC sample size           : ",MC.size,"\n")
        cat("MCMC burn-in sample size   : ",MC.burn,"\n")
        cat("GEL loss function para     : ",q.val,"\n")
        cat("LINEX loss function para   : ",c.val,"\n")
        cat("Probability of censrong (I): ",i1,"\n")
        cat("The value of R             : ",R,"\n")
        cat("============================================","\n\n")
        print(Res_est[,1:13])
        print(Res_est[,c(1:3,14:18)])
        print(Res_aci[,1:11])
        print(Res_aci[,c(1:3,12:15)])
        cat("============================================","\n")
        print(round(Estim,3))
        print(round(HPD_res,3))
        cat("============================================","\n")
        sink()
        break
      }
    }
  }
}








  
  
  
  
