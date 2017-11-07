// Si quiero llamar a otra funcon .cpp y usarla 
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

NumericMatrix cppCRW_exp(double k, double w, double ns, double maxx){
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
  int hm=1;
  for(int i=1; i<nsteps & t[i]<mm; i++){
    hm=hm+1;
    aux2=rvonmises(1,d[i-1],nk); 
    t[i] = t[i-1] + vtm[i];
    d[i] = aux2[0];
    x[i] = x[i-1] + cos(d[i]) * (vls[i]) * vtm[i];
    y[i] = y[i-1] + sin(d[i]) * (vls[i]) * vtm[i]; 
    
    
  }
  NumericMatrix sol(hm,3);
  
  sol(_,0)=x;
  sol(_,1)=y;
  sol(_,2)=t;
 
  return(sol);
  //return DataFrame::create(Named("tu")= tu);
}


NumericMatrix cppObsM(NumericVector x, NumericVector y, NumericVector t, double dt) {
  
  
  NumericVector xs(x);
  NumericVector ys(y);
  NumericVector xt(t);
  double dts = dt;
  int nts = ceil(max(xt)/dts); //ceil(xt[n-1]/dts);
  NumericVector st(nts);
  NumericVector sx(nts);
  NumericVector sy(nts);
  NumericVector ite(nts); 
  NumericMatrix sol(nts,3);
  
  
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
  
  sol(_,0)=sx;
  sol(_,1)=sy;
  sol(_,2)=st;
  
  
  //NumericMatrix::create(Named("sx")= sx, Named("sy") = sy, Named("st") = st)
  return sol;
}


// [[Rcpp::export]]
List ABC_CRW(int nsims,int nsams,int maxt,int nobs,double ddt){
  
  int ntot=nsims;
  int nobss=nobs;
  int nsam=nsams;
  int maxx=maxt;
  double dt=ddt;
  int sobs=1;
  NumericMatrix sX(nobss, ntot);
  NumericMatrix sY(nobss, ntot);
  NumericMatrix sT(nobss, ntot);
 
  // sample from priors
    NumericVector sw(ntot);
    NumericVector sk(ntot);
  
  sw = runif(nsims,0.1,10);
  sk = runif(nsims,5,100);

  //simulations
  for(int i=0; i<ntot; i++)
  {
  // observed data
  int dims=0;
  //NumericMatrix tray(dims,3);
  NumericMatrix tray(cppCRW_exp(sk[i],sw[i],nsam,maxx));
    NumericMatrix Obss(dims,3);
  Obss=cppObsM(tray(_,0),
               tray(_,1),
               tray(_,2),
               ddt);
  dims=Obss(_,0).size();
  int sv=0;
  sv=min(NumericVector::create(nobss,dims));
  
  if (sv<nobss){
    //NumericVector tail((nobss-sv),NA_REAL);
    sX(_,i)=Obss(_,0);
    sY(_,i)=Obss(_,1);
    sT(_,i)=Obss(_,2);
    
    for(int k=sv; k<nobss; k++)
    {
      sX(k,i)=NA_REAL;
      sY(k,i)=NA_REAL;
      sT(k,i)=NA_REAL;
    }
    
    //Rcout << sX(_,i) << "\n";
    //sX(_,i)= NumericVector::create(tail);
   //sY(_,i)= NumericVector::create(tail);
   //sT(_,i)= NumericVector::create(tail);
   
  
  }
  
  else{
    //Rcout << Obss << "\n";
    sX(_,i)= Obss(_,0);
    sY(_,i)= Obss(_,1);
    sT(_,i)= Obss(_,2);
  }
  
  //Rcout << dims << "\n";
  
   }
  return List::create(Named("sX")= sX,Named("sY")= sY,Named("sT")= sT,Named("w")=sw,
                            Named("k")=sk);
  }

