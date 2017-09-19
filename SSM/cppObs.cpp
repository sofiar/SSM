#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame cppObs(NumericVector x, NumericVector y, NumericVector t, double dt) {
NumericVector xs(x);
NumericVector ys(y);
NumericVector xt(t);
double dts = dt;
//int n = xt.size();
int nts = ceil(max(xt)/dts); //ceil(xt[n-1]/dts);
NumericVector st(nts);
NumericVector sx(nts);
NumericVector sy(nts);
NumericVector ite(nts); 

for(int i=1; i<nts; i++){
st[i] = st[i-1] + dts;
NumericVector::iterator it = std::upper_bound(xt.begin(), xt.end(), st[i]);
ite[i] = it -xt.begin()-1;
double dx = xs[ite[i]+1] - xs[ite[i]];
double dy = ys[ite[i]+1] - ys[ite[i]];
double dto = xt[ite[i]+1] - xt[ite[i]];
sx[i] = xs[ite[i]] + ((st[i]-xt[ite[i]])*dx)/dto;
sy[i] = ys[ite[i]] + ((st[i]-xt[ite[i]])*dy)/dto;
}
return DataFrame::create(Named("sx")= sx, Named("sy") = sy, Named("st") = st);
}
