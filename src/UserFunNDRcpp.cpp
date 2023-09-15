#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <roptim.h>
// [[Rcpp::depends(roptim)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>
// [[Rcpp::depends(RcppProgress)]]

// Library functions of the paper:
// Classical and Bayesian Inference of Spmk for
//  Nakagami Distribution Based On First-Failure Progressively
//  Censored Samples

#include<iostream>
#include<algorithm>
#include"post.h"

using namespace Rcpp;
using namespace std;
using namespace Numer;
using namespace arma;
using namespace roptim;
// These codes for the progress bar inside the loop
#include <progress.hpp>
#include <progress_bar.hpp>

// Using qnaka(), pnaka(), dnaka() and rnaka() for the quantile, 
// cdf, pdf and generating fiunction of Nakgami distribution  using VGAM package


// gamma function

// [[Rcpp::export]]
double my_gam(double x){
  Function gam_rcpp("gamma");
  NumericVector res=gam_rcpp(x);
  return res[0];
}
// incomplete upper gamma function
double my_incgam(double x, double y){
  Function incgam_rcpp("incgam");
  NumericVector res=incgam_rcpp(x, y);
  return res[0];
}
// digamma function
double my_digam(double x){
  Function digam_rcpp("digamma");
  NumericVector res=digam_rcpp(x);
  return res[0];
}



// [[Rcpp::export]]
struct Dist{
  // The quantile of Nakagami Dist
  double qnt(double w,Rcpp::NumericVector para){
    Rcpp::Function qnaka_rcpp("qnaka");
    Rcpp::NumericVector res=qnaka_rcpp(w,para[1],para[0]);
    return res[0];
  }
  // The pdf of Nakagami Dist
  double pdf(double w,Rcpp::NumericVector para){
    Rcpp::Function dnaka_rcpp("dnaka");
    Rcpp::NumericVector res=dnaka_rcpp(w,para[1],para[0]);
    return res[0];
  }
  // The cdf of Nakagami Dist
  
  double cdf(double w,Rcpp::NumericVector para){
    Rcpp::Function pnaka_rcpp("pnaka");
    Rcpp::NumericVector res=pnaka_rcpp(w,para[1],para[0]);
    return res[0];
  }
  // The survival of Nakagami Dist
  
  double sur(double w,Rcpp::NumericVector para){
    Rcpp::Function pnaka_rcpp("pnaka");
    Rcpp::NumericVector res=pnaka_rcpp(w,para[1],para[0]);
    return 1.0-res[0];
  }
};

// The mean of Nakagami Dist
// [[Rcpp::export]]
double mu(double x, double y){
  return (my_gam(x+0.5)/my_gam(x))*std::sqrt(y/x);
}

// The standard deviation of Nakagami Dist
// [[Rcpp::export]]
double sig(double x, double y){
  return(std::sqrt(y-std::pow((my_gam(x+0.5)/my_gam(x))*std::sqrt(y/x),2)));
};


NumericVector my_as_mcmc(NumericVector x){
  Rcpp::Function as_mcmc_rcpp("as.mcmc");
  Rcpp::NumericVector res=as_mcmc_rcpp(x);
  return res;
}

// [[Rcpp::export]]
NumericVector my_HPD(NumericVector x){
  Rcpp::Function HPD_rcpp("HPDinterval");
  Rcpp::NumericVector res=HPD_rcpp(x);
  return res;
}


// likelihood function
double like(Rcpp::NumericVector para, 
            Rcpp::NumericVector X, 
            Rcpp::NumericVector R,
            int k){
  
  // Compute objective value
  double lk=1;
  int m=X.size();
  Dist dist;
  for(int i=0;i<m;i++){
    lk *= dist.pdf(X(i),para)*pow(dist.sur(X(i),para),k*(R(i)+1)-1);
  }
  return lk;
}  

// log-likelihood function
double loglike(Rcpp::NumericVector para, 
               Rcpp::NumericVector X, 
               Rcpp::NumericVector R,
               int k){
  
  return -std::log(like(para,X,R,k));
}

double mps(Rcpp::NumericVector para, 
           Rcpp::NumericVector X, 
           Rcpp::NumericVector R,
           int k){
  
  // Compute objective value
  int m=X.size();
  Dist dist;
  double mp=dist.sur(X(m-1),para)*dist.cdf(X(0),para)*pow(dist.sur(X(0),para),k*(R(0)+1)-1);
  for(int i=1;i<m;i++){
    mp *= (dist.cdf(X(i),para)-dist.cdf(X(i-1),para))*pow(dist.sur(X(i),para),k*(R(i)+1)-1);
  }
  return mp;
}


double logmps(Rcpp::NumericVector para, 
              Rcpp::NumericVector X, 
              Rcpp::NumericVector R,
              int k){
  
  return -std::log(mps(para,X,R,k));
}

class LogObjective : public Functor {
  private: arma::vec X;
    arma::vec R;
    int k;
    std::string Type;
public:
  LogObjective(arma::vec xx_, arma::vec rr_, int kk_, std::string type_) : X(xx_), R(rr_), k(kk_), Type(type_){}
  
  double operator()(const arma::vec &y) override {
    if(Type=="LF")
    return loglike(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),
                   Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(X)),
                   Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(R)),k);
    if(Type=="SPF")
      return logmps(Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(y)),
                     Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(X)),
                     Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(R)),k);
  }
};

struct Post{
  double postlike(NumericVector para, 
                  NumericVector X,
                  NumericVector R,
                  int k){
    return like(para,X,R,k)*std::sqrt(para[1]*R::trigamma(para[1])-1)/para[2];
  }  
  double postmps(NumericVector para, 
                 NumericVector X, 
                 NumericVector R,
                 int k){
    return mps(para,X,R,k)*std::sqrt(para[1]*R::trigamma(para[1])-1)/para[2];
  }
};


// These classes are used to compute derivative of incomplete gamma
class P1: public Func {
private:
  double eps;
  double eta;
public:
  P1(double a_, double b_) : eps(a_), eta(b_) {}
  
  double operator()(const double& x) const {
    return std::log(x)*std::pow(x,eps-1)*exp(-x);
  };
};

// Definition of Spmk function
// [[Rcpp::export]]
double spmk_fun(Rcpp::NumericVector para, double l, double u, double t, double m, double s){
  Dist dist;
  return(std::abs(R::qnorm(1.0-(1.0-(dist.cdf(u,para)-dist.cdf(l,para)))/2,0.0,1.0,1,0)/(3.0*std::sqrt(1.0+std::pow((m-t)/s,2)))));
}

// Gradient of Spmk function
// [[Rcpp::export]]
Rcpp::NumericVector spmk_grad(Rcpp::NumericVector para,
                              double l,
                              double u,
                              double t,
                              double m,
                              double s){
  
  double eps=para[0];
  double eta=para[1];
  
  double err_est1;
  int err_code1;
  double err_est2;
  int err_code2;
  
  P1 f1(eps, eta);
  double I_l = integrate(f1, eps*l*l/eta, 1000, err_est1, err_code1);
  double I_u = integrate(f1, eps*u*u/eta, 1000, err_est2, err_code2);
  
  double psi_l=std::log(my_incgam(eps*pow(l,2)/eta,eps));
  double psi_u=std::log(my_incgam(eps*pow(u,2)/eta,eps));
  double psi_eps_l=(-pow(eps,eps-1)*exp(-eps*l*l/eta)*pow(l,2*eps)*pow(eta,-eps)+I_l)/my_incgam(eps*l*l/eta,eps);
  double psi_eps_u=(-pow(eps,eps-1)*exp(-eps*u*u/eta)*pow(u,2*eps)*pow(eta,-eps)+I_u)/my_incgam(eps*u*u/eta,eps);
  double psi_eta_l=pow(eps,eps)*exp(-eps*l*l/eta)*pow(l,2*eps)*pow(eta,-eps-1)/my_incgam(eps*l*l/eta,eps);
  double psi_eta_u=pow(eps,eps)*exp(-eps*u*u/eta)*pow(u,2*eps)*pow(eta,-eps-1)/my_incgam(eps*u*u/eta,eps);
  
  double p=1-(my_incgam(eps*u*u/eta,eps)-my_incgam(eps*l*l/eta,eps))/my_gam(eps);
  double g1=1/(6*std::sqrt(1+pow(m-t,2)/pow(s,2))*R::dnorm(R::qnorm(1-p/2,0,1,0,0),0,1,0));
  double g2=my_incgam(eps*l*l/eta,eps)*(psi_eps_l-my_digam(eps))-my_incgam(eps*u*u/eta,eps)*(psi_eps_u-my_digam(eps));
  double g3=my_incgam(eps*l*l/eta,eps)*psi_eta_l-my_incgam(eps*u*u/eta,eps)*psi_eta_u;
  double G1=g1*g2/my_gam(eps);
  double G2=g1*g3/my_gam(eps);
  return {G1,G2};
} 

// [[Rcpp::export]]
Rcpp::NumericVector GenData(Rcpp::NumericVector para, Rcpp::NumericVector R){
  Rcpp::Function GenData_rcpp("Gen_Data");
  Rcpp::NumericVector res=GenData_rcpp(para, R);
  return(res);
}

// [[Rcpp::export]]
Rcpp::List Estim(arma::vec para, 
                 arma::vec X, 
                 arma::vec R, 
                 int k,
                 std::string type, 
                 std::string method,
                 arma::vec lw,
                 arma::vec up) {
  
    LogObjective rb(X,R,k,type);
    Roptim<LogObjective> opt(method);
    opt.control.trace = 0;
    opt.set_hessian(true);
    //arma::vec lw={0,0};
    //arma::vec up={10,10};
    opt.set_lower(lw);
    opt.set_upper(up); 
    opt.minimize(rb, para);
    return Rcpp::List::create(
      Rcpp::Named("par") = opt.par(),
      Rcpp::Named("value") = opt.value(),
      Rcpp::Named("inform") = opt.hessian(),
      Rcpp::Named("conv") = opt.convergence());
};  


// Function to compute Bayesian estimation and HPD using MCMC technique v=based on MH samples
// and LF (or PSF) functions
// [[Rcpp::export]]
Rcpp::List MH_sample(std::string type,
                     NumericVector para,NumericVector se,
                     NumericVector R,NumericVector X,
                     int k,double l, double u, double t, 
                     int MC_size,int MC_burn,
                     double q, double c,int verbose=0, bool display_progress=true){
  
  NumericVector MH_eps_sel(MC_size),MH_eta_sel(MC_size),MH_spmk_sel(MC_size);
  NumericVector MH_eps_gel(MC_size),MH_eta_gel(MC_size),MH_spmk_gel(MC_size);
  NumericVector MH_eps_linex(MC_size),MH_eta_linex(MC_size),MH_spmk_linex(MC_size);

  MH_eps_sel(0)=para[0];
  MH_eta_sel(0)=para[1];
  MH_spmk_sel(0)=spmk_fun(para,l,u,t,mu(para[0],para[1]),sig(para[0],para[1]));
  MH_eps_gel(0)=pow(MH_eps_sel(0),-q);
  MH_eta_gel(0)=pow(MH_eta_sel(0),-q);
  MH_spmk_gel(0)=pow(MH_spmk_sel(0),-q);
  MH_eps_linex(0)=std::exp(-c*MH_eps_sel(0));
  MH_eta_linex(0)=std::exp(-c*MH_eta_sel(0));
  MH_spmk_linex(0)=std::exp(-c*MH_spmk_sel(0));
  
  double eps1,eta1,del,dad;
  int rho=0;
  
  Post post;
  Functor_MH f(&post,&Post::postlike,&Post::postmps);
  
  Progress p(MC_size, display_progress);
  
  int i=1; 
  while(i<=MC_size){
    eps1=std::exp(rnorm(1,std::log(MH_eps_sel[i-1]),se[0])[0]);
    if(eps1==R_NaN) continue;
    if(eps1<=0.5) continue;
    if(eps1>3.0*para[0]) continue;
    eta1=std::exp(rnorm(1,std::log(MH_eta_sel[i-1]),se[1])[0]);
    if(eta1==R_NaN) continue;
    if(eta1>3.0*para[1]) continue;
    double per=f(type,{eps1,eta1},X,R,k)*eps1*eta1/(f(type,{MH_eps_sel[i-1],MH_eta_sel[i-1]},X,R,k)*MH_eps_sel[i-1]*MH_eta_sel[i-1]);
    if(per==R_NaN) continue;
    del=std::min(1.0,per);
    dad=runif(1)[0];
    rho=0;
    if(dad<del) rho=1;
    MH_eps_sel[i] =eps1*rho+(1-rho)*MH_eps_sel[i-1]; 
    MH_eta_sel[i] =eta1*rho+(1-rho)*MH_eta_sel[i-1]; 
    MH_spmk_sel[i]=spmk_fun({MH_eps_sel[i],MH_eta_sel[i]},l,u,t,mu(MH_eps_sel[i],MH_eta_sel[i]),sig(MH_eps_sel[i],MH_eta_sel[i]));
    
    MH_eps_gel[i] =pow(MH_eps_sel[i],-q);
    MH_eta_gel[i] =pow(MH_eta_sel[i],-q);
    MH_spmk_gel[i]=pow(MH_spmk_sel[i],-q);
    MH_eps_linex[i] =std::exp(-c*MH_eps_sel[i]);
    MH_eta_linex[i] =std::exp(-c*MH_eta_sel[i]);
    MH_spmk_linex[i]=std::exp(-c*MH_spmk_sel[i]);
    
    p.increment(); 
    if(verbose>0)
      std::cout<<"Estimate SEL= "<<MH_eps_sel[i]<<" "<<MH_eta_sel[i]<<" "<<MH_spmk_sel[i]<<std::endl; 
    i+=1;
  };
  
  MH_eps_sel.erase(0,MC_burn-1);
  MH_eta_sel.erase(0,MC_burn-1);
  MH_spmk_sel.erase(0,MC_burn-1);
  
  MH_eps_gel.erase(0,MC_burn-1);
  MH_eta_gel.erase(0,MC_burn-1);
  MH_spmk_gel.erase(0,MC_burn-1);
  
  MH_eps_linex.erase(0,MC_burn-1);
  MH_eta_linex.erase(0,MC_burn-1);
  MH_spmk_linex.erase(0,MC_burn-1);
  
  NumericVector est_sel={mean(MH_eps_sel),mean(MH_eta_sel),mean(MH_spmk_sel)};
  NumericVector est_gel={pow(mean(MH_eps_gel),-1/q),pow(mean(MH_eta_gel),-1/q),pow(mean(MH_spmk_gel),-1/q)};
  NumericVector est_linex={(-1/c)*log(mean(MH_eps_linex)),(-1/c)*log(mean(MH_eta_linex)),(-1/c)*log(mean(MH_spmk_linex))};
  NumericVector hpd_eps={my_HPD(my_as_mcmc(MH_eps_sel))[0],my_HPD(my_as_mcmc(MH_eps_sel))[1],my_HPD(my_as_mcmc(MH_eps_sel))[1]-my_HPD(my_as_mcmc(MH_eps_sel))[0]};
  NumericVector hpd_eta={my_HPD(my_as_mcmc(MH_eta_sel))[0],my_HPD(my_as_mcmc(MH_eta_sel))[1],my_HPD(my_as_mcmc(MH_eta_sel))[1]-my_HPD(my_as_mcmc(MH_eta_sel))[0]};
  NumericVector hpd_spmk={my_HPD(my_as_mcmc(MH_spmk_sel))[0],my_HPD(my_as_mcmc(MH_spmk_sel))[1],my_HPD(my_as_mcmc(MH_spmk_sel))[1]-my_HPD(my_as_mcmc(MH_spmk_sel))[0]};
  
  return Rcpp::List::create(
    Rcpp::Named("est_SEL")  =est_sel,
    Rcpp::Named("est_GEL")  =est_gel,
    Rcpp::Named("est_LINEX")=est_linex,
    Rcpp::Named("HPD_eps")=hpd_eps,
    Rcpp::Named("HPD_eta")=hpd_eta,
    Rcpp::Named("HPD_spmk")=hpd_spmk,
    Rcpp::Named("eps_sample")=MH_eps_sel,
    Rcpp::Named("eta_sample")=MH_eta_sel,
    Rcpp::Named("spmk_sample")=MH_spmk_sel
  );
};

