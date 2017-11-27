
library(circular)
library(Rcpp)

library(truncnorm)
library(mvtnorm)
source('/home/sofia/proyecto_doctoral/pruebas/SSM/funaux.R')
setwd("/home/sofia/proyecto_doctoral/pruebas/SSM/cpp")


Rcpp::sourceCpp('PathelementsCpp.cpp')
Rcpp::sourceCpp('RW_exp_cor.cpp')
Rcpp::sourceCpp('/home/sofia/proyecto_doctoral/pruebas/SSM/cppObs.cpp')
Rcpp::sourceCpp('/home/sofia/proyecto_doctoral/pruebas/SSM/cpp/sl.cpp')

nsteps <- 800 # number of moves performed by the animal

# movement parameters:
t_w=6
t_k = 10
dt <- 3 # time interval for observations


true_sim <- cppRW_exp_cor(t_k,  t_w, nsteps, Inf) # simulate RW using the cpp function
oz <- cppObs(true_sim$x,true_sim$y,true_sim$t,dt) # these are the observed data
maxt=max(oz$st)
nobs=length(oz$st)

### sumaries

### for the real ones
ps=PathelementsCpp(oz$sx,oz$sy)
SumTrue=sumaries(ps$direction,ps$turns, ps$steps, oz$sx,oz$sy)

##########################   SIMULATIONS #############################

n=20
nsims=1000
snsteps=1000
sw=numeric(nsims)
sk=numeric(nsims)

sw[1]=runif(1,0.1,10)
sk[1]=runif(1,0.1,100)

mu=list()
sigma=list()

#### Firt
SumSim=matrix(NA,n,8)
for (j in 1:n)
{
  sim=cppRW_exp_cor(sk[1],sw[1],snsteps,maxt)
  obs=cppObs(sim$x,sim$y,sim$t,dt)
  pe=PathelementsCpp(obs$sx,obs$sy)
  SumSim[j,]=sumaries(pe$direction,pe$turns, pe$steps, obs$sx,obs$sy)

}
mu[[1]]=apply(SumSim,2,mean)
sigma[[1]]=cov(SumSim)
cc=0
nsims=10000
for (i in 2:nsims)
{
  ww=rtruncnorm(1,a=0.1,b=10,mean=sw[i-1],sd=1.5)
  kk=rtruncnorm(1,a=0.1,b=100,mean=sk[i-1],sd=20)
  
  for (j in 1:n)
  {
    sim=cppRW_exp_cor(kk,ww,snsteps,maxt)
    obs=cppObs(sim$x,sim$y,sim$t,dt)
    pe=PathelementsCpp(obs$sx,obs$sy)
    
    SumSim[j,]<-sumaries(pe$direction,pe$turns, pe$steps, obs$sx,obs$sy)
  }
  mup=apply(SumSim,2,mean)
  sigmap=cov(SumSim)
  
  ### Aceptamos?
  coc=(dmvnorm(SumTrue, mean=mup, sigma=sigmap, log=FALSE)*
         dtruncnorm(sw[i-1],a=0.1,b=10,mean=ww,sd=1.5)*
         dtruncnorm(sk[i-1],a=0.1,b=100,mean=kk,sd=20))/
    (dmvnorm(SumTrue, mean=mu[[i-1]], sigma=sigma[[i-1]], log=FALSE)*
       dtruncnorm(ww,a=0.1,b=10,mean=sw[i-1],sd=1.5)*
       dtruncnorm(kk,a=0.1,b=100,mean=sk[i-1],sd=20))
  
  
  if(is.na(coc))
  {
    coc=0
  }
  
  r=min(1,coc)
  u=runif(1)
  if(u<coc)
  {
    sk[i]=kk
    sw[i]=ww
    mu[[i]]=mup
    sigma[[i]]=sigmap
    print('acepto')
    cc=cc+1
    print(kk)
    print(ww)
    
    
  }
  else
  {
    sk[i]=sk[i-1]
    sw[i]=sw[i-1]
    mu[[i]]=mu[[i-1]]
    sigma[[i]]=sigma[[i-1]]
  }
  
}

#### Resultados plots

plot(sk[1:(i-1)],type='l',ylim=c(0,100))
abline(a=t_k,b=0,col="blue")

plot(sw[1:(i-1)],type='l',ylim=c(0,10))
abline(a=t_w,b=0,col="blue")


plot(density(sk[500:(i-1)],bw=.6),xlim=c(0,100),main='k')
abline(v=t_k,col='darkorange3',lwd=2)
abline(v=mean(sk[500:(i-1)]),col='dodgerblue3',lwd=2)


plot(density(sw[500:(i-1)],bw=.6),xlim=c(0,10),main=expression(lambda))
abline(v=t_w,col='darkorange3',lwd=2)
abline(v=mean(sw[500:(i-1)]),col='dodgerblue3',lwd=2)


