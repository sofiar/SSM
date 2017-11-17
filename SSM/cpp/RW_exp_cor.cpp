#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

DataFrame cppRW_exp_cor(double k, double w, double ns, double maxx){
  RNGScope scp;
  Function rvonmises("rvonmises");
  Function rexp("rexp");
  Function runif("runif");
  
  double nk = k;
  double nw = w;
  double nsteps = ns;
  double mm = maxx;
  
  NumericVector vtm(nsteps);
  NumericVector vls(nsteps);
  
  vls = NumericVector(nsteps)+1;
  
  
  NumericVector aux2(1);
  NumericVector aux(1);
  
  //NumericMatrix tu(nsteps,1);
  NumericVector tu(nsteps);
  NumericVector x(nsteps);
  NumericVector y(nsteps);
  NumericVector t(nsteps);
  NumericVector d(nsteps);
  
  vtm=rexp(nsteps,nw);
  aux=rvonmises(1, 0, nk);
  
  //tu[0] = aux[0];// first angle
  
  d[0] = R::runif(0,1)* 2 *3.14159265358979323846; 
  for(int i=1; i<nsteps & t[i]<mm; i++){
    aux2=rvonmises(1,d[i-1],nk); 
    t[i] = t[i-1] + vtm[i];
    d[i] = aux2[0];
    x[i] = x[i-1] + cos(d[i]) * (vls[i]) * vtm[i];
    y[i] = y[i-1] + sin(d[i]) * (vls[i]) * vtm[i]; 
    
    
  }
  
  return DataFrame::create(Named("x")= x,Named("y")= y,Named("t")= t);
  //return DataFrame::create(Named("tu")= tu);
}
