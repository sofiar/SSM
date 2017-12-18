##############################################################################
######################### Analisis de summaries ##############################
##############################################################################

#carga librerias y funciones 
library(circular)
library(Rcpp)
library(RcppArmadillo)
library(tidyverse)
library(ggplot2)
#library(adehabitatHR)
source('SSM/funaux.R')
setwd("/home/sofia/proyecto_doctoral/pruebas")

Rcpp::sourceCpp('SSM/cpp/RW_exp_cor.cpp')
Rcpp::sourceCpp('SSM/cpp/RW_exp.cpp')
Rcpp::sourceCpp('SSM/cpp/cppObs.cpp')
Rcpp::sourceCpp('SSM/cpp/PathelementsCpp.cpp')
Rcpp::sourceCpp('SSM/cpp/sl.cpp')


#############################################################################

################################### CRW #####################################


nsteps <- 2000 # number of moves performed by the animal

# movement parameters:
t_w=6
t_k = 20
dt <- 1.2# time interval for observations

# simulate true moves
NN=10# number of real trayectories to compare

TT=list()
maxt=NA
nobs=NA
for (i in 1:NN)
{
  true_sim <- cppRW_exp_cor(t_k,  t_w, nsteps, Inf) # simulate RW using the cpp function
  oz <- cppObs(true_sim$x,true_sim$y,true_sim$t,dt) # these are the observed data
  TT[[i]]=oz
  maxt=min(maxt,max(oz$st),na.rm=T)
  nobs=min(nobs,length(oz$st),na.rm=T)
}


nsims <- 1e3
nsam <- nsteps*3 # number of real steps in the simulations
maxt = max(oz$st)
nobs = length(oz$sx)
sX = matrix(NA,length(oz$sx),nsims)
sY = matrix(NA,length(oz$sx),nsims)
sT = matrix(NA,length(oz$sx),nsims)

# sample from priors

s_w <- runif(nsims,0.1,10)
#s_mu <- runif(nsims,-pi,pi)
s_k <- runif(nsims,5,100)

suppressWarnings(warning("as.circular"))
suppressWarnings(warning("conversion.circularmuradians0counter"))
suppressWarnings(warning("rvonmises"))

for( j in 1:nsims){
  
  sim=cppRW_exp_cor(s_k[j], s_w[j],nsam,maxt) # simulate RW using the cpp function. Ver aca el tema  de los tiempos
  # observe the movement process at time interval dt
  sz <- cppObs(sim$x,sim$y,sim$t,dt) # these are the observed data
  sobs = length(sz$sx)
  sv <- min(nobs,sobs) 
  
  sX[1:sv,j] <- sz$sx[1:sv]
  sY[1:sv,j] <- sz$sy[1:sv]
  sT[1:sv,j] <- sz$st[1:sv]
  
} 



####### Plot

df=bind_cols(Ssim)
names(df)<-c('Mean Steps','sd Turns', 'Sd Steps', 'SI', 'b. cosine', 
             'Sum(steps)/TF', 'MSD')
#---


#CRW
# simulated parameter against summary stat
data_frame(s_w, s_k) %>% bind_cols(df)%>%
  gather(stat, stat.val, -s_w, -s_k) %>%
  gather(param, par.sim, -stat, -stat.val) %>%
  ggplot() + geom_point(aes(par.sim, stat.val)) + 
  facet_wrap(stat~param, scales='free', ncol = 2, 
             labeller = label_wrap_gen(multi_line=FALSE))
#---




################################### RW #####################################


nsteps <- 2000 # number of moves performed by the animal

# movement parameters:
t_w=6
t_k = 20
t_mu=pi/20
dt <- 10# time interval for observations

# simulate true moves
NN=10# number of real trayectories to compare

TT=list()
maxt=NA
nobs=NA
for (i in 1:NN)
{
  true_sim <- cppRW_exp(t_k, nsteps, t_mu, t_w, Inf) # simulate RW using the cpp function
  oz <- cppObs(true_sim$x,true_sim$y,true_sim$t,dt) # these are the observed data
  TT[[i]]=oz
  maxt=min(maxt,max(oz$st),na.rm=T)
  nobs=min(nobs,length(oz$st),na.rm=T)
}


nsims <- 1e3
nsam <- nsteps*3 # number of real steps in the simulations
maxt = max(oz$st)
nobs = length(oz$sx)
sX = matrix(NA,length(oz$sx),nsims)
sY = matrix(NA,length(oz$sx),nsims)
sT = matrix(NA,length(oz$sx),nsims)

# sample from priors

s_w <- runif(nsims,0.1,10)
s_mu <- runif(nsims,-pi,pi)
s_k <- runif(nsims,5,100)

suppressWarnings(warning("as.circular"))
suppressWarnings(warning("conversion.circularmuradians0counter"))
suppressWarnings(warning("rvonmises"))

for( j in 1:nsims){
  
  sim=cppRW_exp(s_k[j],nsam, s_mu[j],s_w[j],maxt) # simulate RW using the cpp function. Ver aca el tema  de los tiempos
  # observe the movement process at time interval dt
  sz <- cppObs(sim$x,sim$y,sim$t,dt) # these are the observed data
  sobs = length(sz$sx)
  sv <- min(nobs,sobs) 
  
  sX[1:sv,j] <- sz$sx[1:sv]
  sY[1:sv,j] <- sz$sy[1:sv]
  sT[1:sv,j] <- sz$st[1:sv]
  
} 


##################### Calculo de Estadisticos Resumen ########################

Ssim<-data.frame(matrix(nrow=nsims,ncol=7,NA))

nonones=which(!is.na(sY[1,]))

for (m in 1:length(nonones))
{
  j=nonones[m]
  osx=na.omit(sX[,j])
  osy=na.omit(sY[,j])
  ost=na.omit(sT[,j])
  
  
  ps=PathelementsCpp(osx,osy)
  
  a=cdt(osx,osy,ost)
 b=fastLm(ps$cosine[2:length(ps$cosine)]~ps$cosine[1:(length(ps$cosine)-1)]+0)$coefficients[1]
#  bb=lm(ps$cosine[2:length(ps$cosine)]~ps$cosine[1:(length(ps$cosine)-1)]+0)$coefficients[1]

   ct=mean(ps$cosine)
  st=mean(ps$sine)
  
  Ssim[m,]<-c(mean(ps$steps),
              sd(ps$turns),
              sd(ps$steps),
              sicpp(ct,st,mean(ps$steps),sd(ps$steps)),
              b,
              sum(ps$steps)/ost[sv],
              a)

}

####### Plot

df=bind_cols(Ssim)
names(df)<-c('Mean Steps','sd Turns', 'Sd Steps', 'SI', 'b. cosine', 
             'Sum(steps)/TF', 'MSD')
#---

#RW
# simulated parameter against summary stat
data_frame(s_w, s_k,s_mu) %>% bind_cols(df)%>%
  gather(stat, stat.val, -s_w, -s_k,-s_mu) %>%
  gather(param, par.sim, -stat, -stat.val) %>%
  ggplot() + geom_point(aes(par.sim, stat.val)) + 
  facet_wrap(stat~param, scales='free', ncol = 3, 
             labeller = label_wrap_gen(multi_line=FALSE))
#---




