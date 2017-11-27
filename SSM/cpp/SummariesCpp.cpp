// Calculo de los summaries en Cpp
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]


DataFrame SummariesCpp(NumericVector stps,NumericVector tur,NumericVector tt){

  double ms =1;
  double st =1;
  NumericVector tmp(1);
  
  NumericVector steps(stps);
  NumericVector turns(tur);
  NumericVector time(tt);
  int nt=time.size();
  
  ms=mean(steps);
  st=sd(turns);
  Function lm("lm");
  
  tmp = lm(steps~time+0);
 
  
  
  
  return DataFrame::create(Named("mean.steps")= ms,Named("sd.turns")=st,
                           Named("extra")=tmp);
  
}