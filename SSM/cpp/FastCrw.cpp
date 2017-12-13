#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

DataFrame cppFastCRW_exp(double k, double w, double ns, double maxx){

  RNGScope scp;
  Function rvonmises("rvonmises");
  Function rexp("rexp");
  Function runif("runif");
  Function cumsum("cumsum");
  
  double nk = k;
  double nw = w;
  double nsteps = ns;
  double mm = maxx;
  
  NumericVector vtm(nsteps);
  NumericVector vls(nsteps);
  
  
  
  NumericVector aux2(1);
  NumericVector aux(1);
  
  //NumericMatrix tu(nsteps,1);
  NumericVector tu(nsteps);
  NumericVector x(nsteps);
  NumericVector y(nsteps);
  NumericVector t(nsteps);
  NumericVector d(nsteps);
  NumericVector phi(nsteps);
  
  vtm=rexp(nsteps,nw);
  phi=rvonmises(nsteps, 0, nk);
  d=cumsum(phi);  
  
  for(int i=1; i<nsteps & t[i]<mm; i++){
  t[i] = t[i-1] + vtm[i];
    //d[i] = Rcpp::as<double>(rvonmises(1,d[i-1],nk));
  x[i] = x[i-1] + cos(d[i]) * vtm[i];
  y[i] = y[i-1] + sin(d[i]) * vtm[i]; 
    
    
  }
  
  return DataFrame::create(Named("x")= x,Named("y")= y,Named("t")= t);
  
}
  
  
  