#ifndef GAMMA_FUN_H
#define GAMMA_FUN_H

class Gamma_fun{
public:
  Gamma_fun(){};
  double my_gam(double x); // gamma function
  double my_incgam(double x, double y); // incomplete upper gamma function
  double my_digam(double x); // digamma function
};
#endif
