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


// [[Rcpp::export]]
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

double itcpp(NumericVector steps,NumericVector xx, NumericVector yy)
{
  Function sum("sum");
  NumericVector x(xx);
  NumericVector y(yy);
  NumericVector ste(steps);
  
  int N=x.size();
  double distancia2=0;
  double salida=0;
  
  distancia2= pow(x[0]-x[N-1],2)+pow(y[0]-y[N-1],2);
  salida=sqrt(distancia2)/Rcpp::as<double>(sum(ste));;
  
  return(salida);
  }


// [[Rcpp::export]]

double sicpp(double c,double s,double p,double b)
{
  double cc(c);
  double ss(s);
  double pp(p);
  double bb(b);
  double tmp=0;
  
  
  
  tmp = (2/sqrt((pp*((1-pow(cc,2)-pow(ss,2))/(pow(1-cc,2)+pow(ss,2))+pow(bb,2)))));
  return(tmp);
}




//double cdt2cpp (NumericVector sl,NumericVector t)
//{
  
//  NumericVector xx(sl);
//  NumericVector tt(t);
//  double tmp=0;
//  Function lm("lm");
  
//  lm(xx~tt + 0);
  
//  return(tmp);
//}





// [[Rcpp::export]]

NumericVector sumaries(NumericVector direction, NumericVector turns, 
                       NumericVector steps, NumericVector xx, 
                       NumericVector yy) {

  
  NumericVector dir(direction);
  NumericVector tur(turns);
  NumericVector ste(steps);
  NumericVector x(xx);
  NumericVector y(yy);
  
  
  NumericVector sol(8);
  
  Rcpp::Environment base("package:stats"); 
  Rcpp::Function acf_r = base["acf"];  
  Function acf("acf");
  Function mean("mean");
  Function sd("sd");
  Function sum("sum");
  
  double ct=0;
  double st=0;
  double r2=0;
  double t2=0;
  double bb=0;
  int dim=x.size();
  ct=Rcpp::as<double>(mean(cos(tur)));
  st=Rcpp::as<double>(mean(sin(tur)));
  r2=Rcpp::as<double>(sd(x))+Rcpp::as<double>(sd(y));
  
  
  IntegerVector donde = seq(0,dim-2);
  IntegerVector donde2 = seq(1,dim-1);
  t2=(Rcpp::as<double>(sum(x[donde]))-Rcpp::as<double>(sum(x[donde2])))/(dim-2)+
    (Rcpp::as<double>(sum(y[donde]))-Rcpp::as<double>(sum(y[donde2])))/(dim-2);


  sol[0]=Rcpp::as<double>(mean(ste));
  sol[1]=Rcpp::as<double>(sd(tur));
  sol[2]=Rcpp::as<double>(sd(ste));
  //sol[4]=bb;
  sol[3]= sqrt(pow(ct,2)+pow(st,2));
  sol[4]=r2;
  sol[5]=itcpp(ste,x,y);
  sol[6]=sicpp(ct,st,Rcpp::as<double>(mean(ste)),Rcpp::as<double>(sd(ste)));
  sol[7]=t2;
  return sol;
                          

  }

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
  co = adj/step;
  sii = op/step;  
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
  
}


// [[Rcpp::export]] 

double dNM(NumericVector mean,NumericMatrix mcov,NumericVector x)
{
  NumericMatrix Mcov(mcov);
  NumericVector xx(x);
  NumericVector Mean(mean);
  double n= xx.size();
  double sol(0);
  double sol2(0);
  double unidad(1);
  Function det("det");
  Function t("t");
  Function prod("%*%");
  Function solve("solve");
  NumericMatrix isigma(solve(Mcov));
  
  double detr(Rcpp::as<double> (det(Mcov)));
  
  NumericVector res(xx-Mean);
  
  sol2=(exp(prod(prod(-0.5*NumericVector(t(res)),isigma),res)))[0];
    
  sol=(unidad/pow(2*3.141593,n/2))*sol2;  
    
    
  
  
  return(sol);
  
}


// [[Rcpp::export]]

List slcpp(double ddt,double tk ,double tw, int nn, int nsim, int nstep, 
                    int snteps) {

  double dt(ddt);
  double t_k(tk);
  double t_w(tw);
  int n(nn);
  int nsims(nsim);
  int nsteps(nstep);
  int stpes(snteps);
  int cr(0);
  double inf=1.0/0.0;
  Function runif("runif");
  Function rtruncnorm("rtruncnorm");
  Function dtruncnorm("dtruncnorm");
  Function mean("mean");
  Function cov("cov");
  NumericVector mup(8);
  NumericMatrix sigmap(8,8); 
  
  // Real trajectorie 
  
  NumericMatrix Rtraj(cppCRW_exp(t_k,  t_w, nsteps, inf));
  NumericMatrix Robs(cppObsM(Rtraj(_,0),Rtraj(_,1),Rtraj(_,2),dt));  
  double maxt(max(Robs(_,2)));  
  int nobs (Robs(_,2).size());  
   
  // Summaries
  NumericVector SumTrue(8);
  Rcpp::List ps(PathelementsCpp(Robs(_,0),Robs(_,1)));
  SumTrue=sumaries(ps[2],ps[1],ps[0],Robs(_,0),Robs(_,1));
           
  //Simulations
  
  NumericVector sw(nsims);
  NumericVector sk(nsims);
  
  sw[0]=R::runif(0.1,10);
  sk[0]=R::runif(0.1,100);
   
  double ww=0; 
  double kk=0; 
  NumericVector muu(8);
  NumericMatrix sigmaa(8,8);
  
  
  //Firts
  ww=Rcpp::as<double>(rtruncnorm(1,0.1,10.0,sw[0],1.5));
  kk=Rcpp::as<double>(rtruncnorm(1,0.1,100.0,sk[0],20.0));
  
  NumericMatrix SumSim(n,8);
  for (int m=0;m<n;m++)
  {
  NumericMatrix Straj(cppCRW_exp(kk,ww, stpes, maxt));
  NumericMatrix Sobs(cppObsM(Straj(_,0),Straj(_,1),Straj(_,2),dt));  
  
  // Summaries
  Rcpp::List ps(PathelementsCpp(Sobs(_,0),Sobs(_,1)));
  SumSim(m,_)=sumaries(ps[2],ps[1],ps[0],Sobs(_,0),Sobs(_,1));  
  }
  for(int k=0; k<8; k++){
    muu(k)=Rcpp::as<double>(mean(SumSim(_,k)));  
  }
  sigmaa=cov(SumSim);
  
 //Rest
  
  for(int i=1; i<(nsims); i++){
   
  for (int m=0; m<n;m++ )
  {
  NumericMatrix Straj(cppCRW_exp(kk,ww, stpes, maxt));
  NumericMatrix Sobs(cppObsM(Straj(_,0),Straj(_,1),Straj(_,2),dt));  
    
  // Summaries
  Rcpp::List ps(PathelementsCpp(Sobs(_,0),Sobs(_,1)));
  SumSim(m,_)=sumaries(ps[2],ps[1],ps[0],Sobs(_,0),Sobs(_,1));  
  }
  
  for(int k=0; k<8; k++){
  mup[k]=Rcpp::as<double>(mean(SumSim(_,k)));  
    
  }
  sigmap=cov(SumSim);
  
  //Acepatamos?
  double coc(0);
  double num(0);
  double den(0);
  num=dNM(muu,sigmaa,SumTrue)*Rcpp::as<double>(dtruncnorm(sw[i-1],0.1,10.0,ww,1.5))*
    Rcpp::as<double>(dtruncnorm(sk[i-1],0.1,100.0,kk,20.0));
  den=dNM(muu,sigmaa,SumTrue)*Rcpp::as<double>(dtruncnorm(ww,0.1,10.0,sw[i-1],1.5))
    *Rcpp::as<double>(dtruncnorm(kk,0.1,100.0,sk[i-1],1.5));
    
    Rcout<<dNM(muu,sigmaa,SumTrue)<<'\n'; 
 
 
 coc=num/den;
 if (isnan(coc))
 {
   coc=0;
 }
  
 
 double u (R::runif(0.1,10));
 
 if (u<coc)
 {
   sk[i]=kk;
   sw[i]=ww;
   muu=mup;
   sigmaa=sigmap;
   cr=cr+1;
  // Rcout<<kk<<'\n'; 
  //Rcout<<ww<<'\n'; 
 }
 
 else
   {
   sk[i]=sk[i-1];
   sw[i]=sw[i-1];
   }

 }
  return List::create(Named("k",sk),
                      Named("w", sw));
  
  }

