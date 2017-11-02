#### Modelo CRW con tiempo exponencial. 
#### Comparamos con NN trayectorias verdaderas para darle mayor robustez al ajuste

library(circular)
library(Rcpp)
source('/home/sofia/proyecto_doctoral/pruebas/SSM/funaux.R')
setwd("/home/sofia/proyecto_doctoral/pruebas/SSM/cpp")

Rcpp::sourceCpp('abc_crw.cpp')
Rcpp::sourceCpp('PathelementsCpp.cpp')
Rcpp::sourceCpp('RW_exp_cor.cpp')
Rcpp::sourceCpp('/home/sofia/proyecto_doctoral/pruebas/SSM/cppObs.cpp')

nsteps <- 800 # number of moves performed by the animal

# movement parameters:
t_w=3
t_k = 20
dt <- 3 # time interval for observations

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

#maxt = max(oz$st)
#nobs = length(oz$sx)

#########################################################################################
########################################## ABC ##########################################
#########################################################################################

# trayectories and observations

nsims=1e4
a=ABC_CRW(nsims,1000,maxt,nobs,dt)


########################################################################################
##################################### Summaries ########################################
########################################################################################

Summ=matrix(NA,NN,9)
### for the real ones

for (i in 1:NN)
{
  oz=TT[[i]][1:nobs,]
  ps=PathelementsCpp(oz$sx,oz$sy)
 
  
  #bb=acf(circular(ps$direction),plot=FALSE)
  #llmmm=lm(bb$lag[2:length(bb$lag)]~bb$acf[2:length(bb$lag)])
  
  aa=acf(ps$steps,plot=FALSE)
  #llmm=lm(aa$lag[2:length(aa$lag)]~aa$acf[2:length(aa$lag)])
  
 t2=(sum((oz$sx[1:(length(oz$sx)-1)]-oz$sx[2:(length(oz$sx))])^2))/(length(oz$sx)-1)/+
   (sum((oz$sy[1:(length(oz$sy)-1)]-oz$sy[2:(length(oz$sy))])^2))/(length(oz$sy)-1)
r2=sd(oz$sx)+sd(oz$sy)
  
  ct=mean(cos(ps$turns))
  st=mean(sin(ps$turns))
  bo=sd(ps$steps)#/abs(mean(ps$steps))
  
  Summ[i,]<-c(mean(ps$steps),
           sd(ps$turns),
           cdt2(ps$steps,oz$st[2:length(oz$st)]),
           sd(ps$steps),
           #mean(ps$turns),
           #mean(aa$acf),
           sqrt((mean(cos(ps$turns)))^2+(mean(sin(ps$turns)))^2),
         #  cdt2(nsd(ps$steps),oz$st[2:length(oz$st)]))
  it(ps$steps,oz$sx,oz$sy),
  #sd.circular(circular(ps$direction)),
  #angular.deviation(circular(ps$turns)),
  si(ct,st,mean(ps$steps),bo),
   t2,
     r2)
  
  
}

### for the simulation ones


SSum=matrix(NA,nsims,9)
for (j in 1:nsims)
{
  osx=a$sX[,j]
  osy=a$sY[,j]
  ost=a$sT[,j]
  
   ps=PathelementsCpp(osx,osy)
  
  
  #bb=acf(circular(ps$direction),plot=FALSE)
  #llmmm=lm(bb$lag[2:length(bb$lag)]~bb$acf[2:length(bb$lag)])
  
  aa=acf(ps$steps,plot=FALSE)
  #llmm=lm(aa$lag[2:length(aa$lag)]~aa$acf[2:length(aa$lag)])
  
  t2=(sum((osx[1:(length(osx)-1)]-osx[2:(length(osx))])^2))/(length(osx)-1)/+
    (sum((osy[1:(length(osy)-1)]-osy[2:(length(osy))])^2))/(length(osy)-1)
  r2=sd(na.omit(osx))+sd(na.omit(osy))
  
  ct=mean(cos(ps$turns))
  st=mean(sin(ps$turns))
  bo=sd(ps$steps)#/abs(mean(ps$steps))
  
  SSum[j,]<-c(mean(ps$steps),
               sd(ps$turns),
               cdt2(ps$steps,ost[2:length(ost)]),
               sd(ps$steps),
               #mean(ps$turns),
               #mean(aa$acf),
               sqrt((mean(cos(ps$turns)))^2+(mean(sin(ps$turns)))^2),
               #  cdt2(nsd(ps$steps),oz$st[2:length(oz$st)]))
               it(ps$steps,osx,osy),
               #sd.circular(circular(ps$direction)),
               #angular.deviation(circular(ps$turns)),
               si(ct,st,mean(ps$steps),bo),
               t2,
               r2)
 
}

View(SSum)

########################################################################################
##################################### Compare ##########################################
########################################################################################




