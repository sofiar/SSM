#### Probando cosas 

library(circular)
library(Rcpp)
source('/home/sofia/proyecto_doctoral/pruebas/SSM/funaux.R')
setwd("/home/sofia/proyecto_doctoral/pruebas/SSM/cpp")

Rcpp::sourceCpp('RW_exp.cpp')
Rcpp::sourceCpp('/home/sofia/proyecto_doctoral/pruebas/SSM/cppObs.cpp')
### Load ... WSdt
#load("~/proyecto_doctoral/pruebas/SSM/cpp/ws15.RData")


nsteps <- 800 # number of moves performed by the animal

# movement parameters:
t_w=2
t_k = 20
t_mu=pi/50

# time interval for observations
dt <- 1.5 

# simulate true moves
ntrue=10
SumObs=data.frame(matrix(nrow=ntrue,ncol=5,NA))

for (l in 1:ntrue)
{
  
true_sim <- cppRW_exp(t_k, nsteps, t_mu, t_w, Inf) # simulate RW using the cpp function
oz <- cppObs(true_sim$x,true_sim$y,true_sim$t,dt) # these are the observed data

## Sumaries
ps=pathelements(oz$sx,oz$sy)

SumObs[l,]<-c(mean(ps$steps),
         sd(ps$turns),
         cdt2(ps$steps,oz$st[2:length(oz$st)]),
         sd(ps$steps),
         sqrt((mean(cos(ps$turns)))^2+(mean(sin(ps$turns)))^2))



}

######## Veamos como cambia el ajuste dependiendo de la trayectoria Verdadera


nbest<-200
A=cov(na.omit(s))
# s son los summaries en "quienes" (que summaries) de las observaciones 

for (l in 1:ntrue)
{

Mahal=apply(s,1,mahalanobis,center=as.numeric(SumObs[l,]),cov=A)
which=numeric(nsims)
which[order(Mahal)[1:nbest]]=1 


par(mfrow=c(1,3),mar=c(8,5,4,1))


##### w
density_scale=density(s_w[which==1],bw=0.3,from=-1.0,to=10)
xscale <- seq(-0.5, 10.5, length=100)
y_scale <- dunif(xscale,min = 0,max=10)

plot(xscale,y_scale,type='l',main=expression(lambda),ylim=c(0,max(y_scale,max(density_scale$y)))
     ,col='red',xlab = '',ylab = 'Density',cex.axis=1.9,
     cex.main=2,cex.lab=2,lwd=2)

abline(v=mean(s_w[which==1]),col='dodgerblue3',lwd=2)
abline(v=t_w,col='darkorange3',lwd=2)
lines(density_scale,lwd=2)

##### mu
density_scale=density(s_mu[which==1],from=-pi,to=pi,bw=0.1)
xscale <- seq(-(pi)*5/4, pi*5/4, length=100)
y_scale <- dunif(xscale,min = -pi, max=pi)

plot(xscale,y_scale,type='l',main=expression(mu),
     ylim=c(0,max(y_scale,max(density_scale$y))),col='red',
     xlab = '',ylab = '',cex.main=2,cex.lab=2,lwd=2,cex.axis=1.9)
abline(v=mean(s_mu[which==1]),col='dodgerblue3',lwd=2)
abline(v=t_mu,col='darkorange3',lwd=2)
lines(density_scale,lwd=2)


##### k
density_scale=density(s_k[which==1],bw=4,from=-1,to=100)
xscale <- seq(0,105, length=100)
y_scale <- dunif(xscale,min = 5, max=100)

plot(xscale,y_scale,type='l',main=expression(k),
     ylim=c(0,max(y_scale,max(density_scale$y))),xlab='',col='red',ylab = '',
     cex.main=2,cex.lab=2,lwd=2,cex.axis=1.9)
abline(v=mean(s_k[which==1]),col='dodgerblue3',lwd=2)
abline(v=t_k,col='darkorange3',lwd=2)
lines(density_scale,lwd=2)

legend(-115,0.0008,xpd=NA,bty='n',
       horiz=TRUE,legend=c('Prior','True Value')
       ,fill=c('red','darkorange3'),cex=2,xjust=0.5,yjust=2)



legend(-115,-0.005,xpd=NA,bty='n',
       horiz=TRUE,legend=c('Estimated Mean', 'Estimated Posterior'),
       fill=c('dodgerblue3','black'),cex=2,xjust=0.5,yjust=2)



}


#################### Veamos como es la varianza de los summaries ############################
ntrue=1e4
SumObs=data.frame(matrix(nrow=ntrue,ncol=5,NA))
for (l in 1:ntrue)
{
  
  true_sim <- cppRW_exp(t_k, nsteps, t_mu, t_w, Inf) # simulate RW using the cpp function
  oz <- cppObs(true_sim$x,true_sim$y,true_sim$t,dt) # these are the observed data
  
  ## Sumaries
  ps=pathelements(oz$sx,oz$sy)
  
  SumObs[l,]<-c(mean(ps$steps),
                sd(ps$turns),
                cdt2(ps$steps,oz$st[2:length(oz$st)]),
                sd(ps$steps),
                sqrt((mean(cos(ps$turns)))^2+(mean(sin(ps$turns)))^2))
  
  
  
}


par(mfrow=c(2,5),mar=c(8,5,4,1))
hist(SumObs[,1])
hist(SumObs[,2])
hist(SumObs[,3])
hist(SumObs[,4])
hist(SumObs[,5])



hist(s[,1])
hist(s[,2])
hist(s[,3])
hist(s[,4])
hist(s[,5])


####### Si quisiera hacer una comparacion considerando las trayectorias originales de a grupo 
ngroup=100
SumObs_ng=apply(SumObs[1:ngroup,],2,mean)
SumObs_ng=apply(SumObs[501:600,],2,mean)

### Tengo que generar ngruop trayectorias con los mismos parametros para poder comparar 

Nnsim=4e4
nsam <- 1000 # number of real steps in the simulations
maxt = max(oz$st) # ver esto
SumSim=matrix(nrow=Nnsim,ncol=length(SumObs_ng),NA)


# sample from priors
s_w <- runif(Nnsim,0.1,10)
s_mu <- runif(Nnsim,-pi,pi)
s_k <- runif(Nnsim,0.1,100)

#-------------------------------------------------------------------------------
# simulate movement and observations
# note that there's only one simulation per parameter combination but we probably need several!
#-------------------------------------------------------------------------------

suppressWarnings(warning("as.circular"))
suppressWarnings(warning("conversion.circularmuradians0counter"))
suppressWarnings(warning("rvonmises"))

for( j in 1:Nnsim){
  ssum=matrix(ncol=length(s),nrow=ngroup)
  for(i in 1:ngroup)
  {
  sim=cppRW_exp(s_k[j], nsam, s_mu[j], s_w[j],maxt) # simulate RW using the cpp function. Ver aca el tema  de los tiempos
  # observe the movement process at time interval dt
  sz <- cppObs(sim$x,sim$y,sim$t,dt) # these are the observed data
  sobs = length(sz$sx)
  sv <- min(nobs,sobs) 
  
  if(sv>2) # if one simulation is too short to observe something sv<=1
  {
    
    pe=pathelements(sz$sx[1:sv],sz$sy[1:sv])
    cts=mean(cos(pe$turns))
    sts=mean(sin(pe$turns))
    
    
    ssum[i,]<-c(mean(pe$steps),
                sd(pe$turns),
                cdt2(pe$steps,sz$st[2:sv]),
                sd(pe$steps),
                sqrt((cts)^2+(sts)^2))
  }
  
      }

  SumSim[j,]=apply(ssum,2,mean)
  
  } 


########## Comparacion de a 100 #######

nbest<-100

#Mahalanobis
A=cov(na.omit(SumSim))
Mahal=apply(SumSim,1,mahalanobis,center=SumObs_ng,cov=A)

which=numeric(nsims)
which[order(Mahal)[1:nbest]]=1 



par(mfrow=c(1,3),mar=c(8,5,4,1))


##### w
density_scale=density(s_w[which==1],bw=0.3,from=-1.0,to=10)
xscale <- seq(-0.5, 10.5, length=100)
y_scale <- dunif(xscale,min = 0,max=10)

plot(xscale,y_scale,type='l',main=expression(lambda),ylim=c(0,max(y_scale,max(density_scale$y)))
     ,col='red',xlab = '',ylab = 'Density',cex.axis=1.9,
     cex.main=2,cex.lab=2,lwd=2)

abline(v=mean(s_w[which==1]),col='dodgerblue3',lwd=2)
abline(v=t_w,col='darkorange3',lwd=2)
lines(density_scale,lwd=2)

##### mu
density_scale=density(s_mu[which==1],from=-pi,to=pi,bw=0.1)
xscale <- seq(-(pi)*5/4, pi*5/4, length=100)
y_scale <- dunif(xscale,min = -pi, max=pi)

plot(xscale,y_scale,type='l',main=expression(mu),
     ylim=c(0,max(y_scale,max(density_scale$y))),col='red',
     xlab = '',ylab = '',cex.main=2,cex.lab=2,lwd=2,cex.axis=1.9)
abline(v=mean(s_mu[which==1]),col='dodgerblue3',lwd=2)
abline(v=t_mu,col='darkorange3',lwd=2)
lines(density_scale,lwd=2)


##### k
density_scale=density(s_k[which==1],bw=4,from=0,to=100)
xscale <- seq(0,105, length=100)
y_scale <- dunif(xscale,min = 5, max=100)

plot(xscale,y_scale,type='l',main=expression(k),
     ylim=c(0,max(y_scale,max(density_scale$y))),xlab='',col='red',ylab = '',
     cex.main=2,cex.lab=2,lwd=2,cex.axis=1.9)
abline(v=mean(s_k[which==1]),col='dodgerblue3',lwd=2)
abline(v=t_k,col='darkorange3',lwd=2)
lines(density_scale,lwd=2)

legend(-115,0.0008,xpd=NA,bty='n',
       horiz=TRUE,legend=c('Prior','True Value')
       ,fill=c('red','darkorange3'),cex=2,xjust=0.5,yjust=2)



legend(-115,-0.005,xpd=NA,bty='n',
       horiz=TRUE,legend=c('Estimated Mean', 'Estimated Posterior'),
       fill=c('dodgerblue3','black'),cex=2,xjust=0.5,yjust=2)






