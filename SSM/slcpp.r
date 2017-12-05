
library(circular)
library(Rcpp)

library(truncnorm)
library(mvtnorm)
source('/home/sofia/proyecto_doctoral/pruebas/SSM/funaux.R')
setwd("/home/sofia/proyecto_doctoral/pruebas/SSM/cpp")


Rcpp::sourceCpp('PathelementsCpp.cpp')
Rcpp::sourceCpp('RW_exp_cor.cpp')
Rcpp::sourceCpp('/home/sofia/proyecto_doctoral/pruebas/SSM/cpp/cppObs.cpp')
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

n=25
n.par=2
nsims=1000
snsteps=1000
sw=numeric(nsims)
sk=numeric(nsims)

ac = matrix(0,n.par,2)  # for acceptance rates
sd.prop = c(1,3)       # sd for the proposal distributions
eig <- eigen(diag(n.par))
la = eig$values
vect = eig$vectors

params=matrix(0,nsims,n.par)
params[1,1]=runif(1,0.1,10)
params[1,2]=runif(1,0.1,100)


#### Firt
SumSim=matrix(NA,n,8)
for (j in 1:n)
{
  sim=cppRW_exp_cor(k=params[1,2],w=params[1,1],snsteps,maxt)
  obs=cppObs(sim$x,sim$y,sim$t,dt)
  pe=PathelementsCpp(obs$sx,obs$sy)
  SumSim[j,]=sumaries(pe$direction,pe$turns, pe$steps, obs$sx,obs$sy)

}
mu1=apply(SumSim,2,mean)
sigma1=cov(SumSim)

slik=dmvnorm(SumTrue, mean=mu1, sigma=sigma1, log=TRUE)


for (i in 2:nsims)
{
  params[i,]=params[i-1,]
  for (k in 1:n.par)
  {
  newp=params[i,]  
  newp[k]=rnorm(1,mean=params[i-1,k],sd=sd.prop[k])  
  #print(newp)
  print(i)
  if(newp[k]<0)
  {next}
  ac[k,1] = ac[k,1] + 1
 
  SumSim=sloop(dt,k=newp[2] ,w=newp[1], n, snsteps,maxt)
   
 # for (j in 1:n)
#  {
#    sim=cppRW_exp_cor(k=params[i,2],w=params[i,1],snsteps,maxt)
#    obs=cppObs(sim$x,sim$y,sim$t,dt)
#    pe=PathelementsCpp(obs$sx,obs$sy)
    
 #   SumSim[j,]<-sumaries(pe$direction,pe$turns, pe$steps, obs$sx,obs$sy)
#  }
 
  mup=apply(SumSim,2,mean)
  sigmap=cov(SumSim)
  
  n.slik=dmvnorm(SumTrue, mean=mup, sigma=sigmap, log=TRUE)
  
  ### Aceptamos?
  coc=exp(n.slik-slik)
  print(n.slik)
  
  r=min(1,coc)
  if(runif(1)<r) 
  {
    ac[k,2] = ac[k,2] + 1
    params[i,]=newp
    slik=n.slik
    #mu[[i]]=mup
    #sigma[[i]]=sigmap
    print('aceptamos')
    print(params[i,])
    
    
    
  }
  #else
  #{
  #  sk[i]=sk[i-1]
  #  sw[i]=sw[i-1]
  #  mu[[i]]=mu[[i-1]]
  #  sigma[[i]]=sigma[[i-1]]
  #}  
    
  }
  
}

#### Proporcion de Aceptacion
ac[,2]/ac[,1]


#### Resultados plots

plot(params[1:(i-1),1],type='l',ylim=c(0,10))
abline(a=t_w,b=0,col="blue")

plot(params[1:(i-1),2],type='l',ylim=c(0,100))
abline(a=t_k,b=0,col="blue")

plot(density(sk[500:(i-1)],bw=.6),xlim=c(0,100),main='k')
abline(v=t_k,col='darkorange3',lwd=2)
abline(v=mean(sk[500:(i-1)]),col='dodgerblue3',lwd=2)


plot(density(sw[500:(i-1)],bw=.6),xlim=c(0,10),main=expression(lambda))
abline(v=t_w,col='darkorange3',lwd=2)
abline(v=mean(sw[500:(i-1)]),col='dodgerblue3',lwd=2)


#################################################################################
library(coda)

n.b=nsims/3
w.chain=mcmc(params[n.b:nsims,1])  
summary(w.chain)
plot(w.chain)

k.chain=mcmc(params[n.b:nsims,2])  
summary(k.chain)
plot(k.chain)






