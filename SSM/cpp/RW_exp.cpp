//Movimiento no correlado para tiempo por una exponencial lambda

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

DataFrame cppRW_exp(double k, double nsteps, double mu, double w, double maxx){
  
  int ns = nsteps;
  double nk = k;
  double nmu = mu;
  double nw = w;
  double mm = maxx;
  
  NumericVector tm(nsteps);
  NumericVector tu(nsteps);
  NumericVector ls(nsteps);
  ls = NumericVector(ns)+1;
  
  NumericVector x(nsteps);
  NumericVector y(nsteps);
  NumericVector t(nsteps);
  NumericVector d(nsteps);
  
  RNGScope scp;
  Function rvonmises("rvonmises");
  Function rexp("rexp");
 
  
  tm = rexp(ns, nw);
  tu = rvonmises(ns, nmu, nk);
  
  d[0] = rand() * 2 *3.14159265358979323846;
  
  for(int i=1; i<ns & t[i]<mm; i++){
    t[i] = t[i-1] + tm[i];
    d[i] = d[i-1] + tu[i];
    x[i] = x[i-1] + cos(d[i]) * (ls[i]) * tm[i];
    y[i] = y[i-1] + sin(d[i]) * (ls[i]) * tm[i];
    }
  return DataFrame::create(Named("x")= x,Named("y")= y,Named("t")= t);
  //return DataFrame::create(Named("d")=d);
  
}

