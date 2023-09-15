#ifndef POST_H
#define POST_H
using namespace std;
using namespace Rcpp;

struct Post;
struct Functor_MH{
  Post* pPost;
  double (Post::*pF1)(Rcpp::NumericVector para, 
          Rcpp::NumericVector X, 
          Rcpp::NumericVector R,
          int k);
  double (Post::*pF2)(Rcpp::NumericVector para, 
          Rcpp::NumericVector X, 
          Rcpp::NumericVector R,
          int k);
  
  Functor_MH(Post* ppost,double (Post::*pf1)(Rcpp::NumericVector para, 
                              Rcpp::NumericVector X, 
                              Rcpp::NumericVector R,
                              int k), 
                              double (Post::*pf2)(Rcpp::NumericVector para, 
                                      Rcpp::NumericVector X, 
                                      Rcpp::NumericVector R,
                                      int k)){
    pPost=ppost;
    pF1 = pf1;
    pF2 = pf2;
  }
  double operator()(std::string type,Rcpp::NumericVector para, 
                  Rcpp::NumericVector X, 
                  Rcpp::NumericVector R,
                  int k){
    if(type=="LF") return (pPost->*pF1)(para,X,R,k); else return (pPost->*pF2)(para,X,R,k);
  }
};
#endif
