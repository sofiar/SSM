#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

    
List PathelementsCpp(NumericVector xx, NumericVector yy){

  NumericVector x(xx);
  NumericVector y(yy);
  int n = x.size();
  int m = x.size()-1;
  int mm = x.size()-2;
  
  NumericVector adj(m);
  NumericVector co(m);
  NumericVector si(m);
  NumericVector op(m);
  NumericVector ang(m);
  NumericVector step(m);
  NumericVector sii(m);
  NumericVector adif(m-1);
  
  sii=sii-1;
 for(int ii=0; ii<m; ii++){
 
 adj[ii]=x[ii+1]-x[ii];
op[ii]=y[ii+1]-y[ii];
step[ii]=sqrt(pow(adj[ii], 2.0)+pow(op[ii], 2.0));
  if (step[ii]==0)
{
  step[ii]=0.0000001;      
}
  if (adj[ii]>0)
{
  sii[ii]=1;      
}
if (op[ii]<0)
  {
ang[ii]=sii[ii]*3.141593;      
}
  }
    
 ang = ang+atan(adj/op);
for(int j=0; j<(mm-1); j++){
adif[j] = ang[j+1]-ang[j];
//Corregimos para que quede entre -pi y pi
if (adif[j]<-3.141593)
  {adif[j]=adif[j]+2*(3.141593);}
  if (adif[j]>3.141593)
  {adif[j]=adif[j]-2*(3.141593);}
      }
  sii = adj/step;
  co = op/step;  
//  adif[mm]=  0.0 / 0.0;

  //List ret;
  //ret["steps"] = step;
  //ret["turns"] = adif;
  //ret["direction"] = ang;
  //ret["cosine"]=co;
  //ret["sine"]=sii;
  
  //return ret;  
  
  return List::create(Named("steps",step),
                      Named("turns", adif),
                      Named("direction",ang),
                      Named("cosine", co),
                      Named("sine", sii));
  
  
//return DataFrame::create(Named("steps")= step, Named("turns") = adif, Named("direction") = ang,Named("cosine")=co,Named("sine")=sii);
                                               

}