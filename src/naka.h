#ifndef NAKA_H
#define NAKA_H

class Naka{
public:
  Naka(){};
  double qnt(double w,Rcpp::NumericVector para);
  double pdf(double w,Rcpp::NumericVector para);
  double cdf(double w,Rcpp::NumericVector para);
  double sur(double w,Rcpp::NumericVector para);
  // The mean of Nakagami Dist
  double mu(double x, double y);
  // The standard deviation of Nakagami Dist
  double sig(double x, double y);
};
#endif