#### Modelo exponencial con moviemiento No correlado (usando cpp)
#### Comparamos con NN trayectorias verdaderas para darle mayor robustez al ajuste

library(circular)
library(Rcpp)
source('/home/sofia/proyecto_doctoral/pruebas/SSM/funaux.R')
setwd("/home/sofia/proyecto_doctoral/pruebas/SSM/cpp")

Rcpp::sourceCpp('PathelementsCpp.cpp')
Rcpp::sourceCpp('RW_exp.cpp')
Rcpp::sourceCpp('cppObs.cpp')

nsteps <- 800 # number of moves performed by the animal

# movement parameters:
t_w=2
t_k = 20
t_mu=pi/50
dt <- 3 # time interval for observations

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


#######################################################################################################
##########################################       ABC       ############################################ 
#######################################################################################################


nsims <- 4e5#### see
nsam <- 1000 # number of real steps in the simulations
maxt = max(oz$st)
nobs = length(oz$sx)
sX = matrix(NA,length(oz$sx),nsims)
sY = matrix(NA,length(oz$sx),nsims)
sT = matrix(NA,length(oz$sx),nsims)

# sample from priors

s_w <- runif(nsims,0.1,10)
s_mu <- runif(nsims,-pi,pi)
s_k <- runif(nsims,5,100)

#-------------------------------------------------------------------------------
# simulate movement and observations
# note that there's only one simulation per parameter combination but we probably need several!
#-------------------------------------------------------------------------------

suppressWarnings(warning("as.circular"))
suppressWarnings(warning("conversion.circularmuradians0counter"))
suppressWarnings(warning("rvonmises"))

for( j in 1:nsims){
  
  sim=cppRW_exp(s_k[j], nsam, s_mu[j], s_w[j],maxt) # simulate RW using the cpp function. Ver aca el tema  de los tiempos
  # observe the movement process at time interval dt
  sz <- cppObs(sim$x,sim$y,sim$t,dt) # these are the observed data
  sobs = length(sz$sx)
  sv <- min(nobs,sobs) 
  
  sX[1:sv,j] <- sz$sx[1:sv]
  sY[1:sv,j] <- sz$sy[1:sv]
  sT[1:sv,j] <- sz$st[1:sv]
  
} 




############################### Summaries ######################################

## Calculation of the summaries for the TRUE observation


Stobs=matrix(NA,NN,8)
### for the real ones

for (i in 1:NN)
{
  oz=TT[[i]][1:nobs,]
  ps=PathelementsCpp(oz$sx,oz$sy)
  
  
  #bb=acf(circular(ps$direction),plot=FALSE)
  #llmmm=lm(bb$lag[2:length(bb$lag)]~bb$acf[2:length(bb$lag)])
  
  aa=acf(ps$steps,plot=FALSE)
  #llmm=lm(aa$lag[2:length(aa$lag)]~aa$acf[2:length(aa$lag)])
  
  #t2=(sum((oz$sx[1:(length(oz$sx)-1)]-oz$sx[2:(length(oz$sx))])^2))/(length(oz$sx)-1)/+
  #  (sum((oz$sy[1:(length(oz$sy)-1)]-oz$sy[2:(length(oz$sy))])^2))/(length(oz$sy)-1)
  #r2=sd(oz$sx)+sd(oz$sy)
  
  ct=mean(cos(ps$turns))
  st=mean(sin(ps$turns))
  bo=sd(ps$steps)#/abs(mean(ps$steps))
  
  Stobs[i,]<-c(mean(ps$steps),
           sd(ps$turns),
           cdt2(ps$steps,oz$st[2:length(oz$st)]),
           sd(ps$steps),
           #mean(ps$turns),
           #mean(aa$acf),
           sqrt((mean(cos(ps$turns)))^2+(mean(sin(ps$turns)))^2),
           cdt2(nsd(ps$steps),oz$st[2:length(oz$st)]),
  it(ps$steps,oz$sx,oz$sy),
  #sd.circular(circular(ps$direction)),
  #angular.deviation(circular(ps$turns)),
  si(ct,st,mean(ps$steps),bo))
  # t2,
  #   r2
  
}


Ssim<-data.frame(matrix(nrow=nsims,ncol=length(Stobs),NA))

for (i in 1:nsims)
{
  sv=length(na.omit(sX[,i]))
  if(sv>2) # if one simulation is too short to observe something sv<=1
  {
    
    # sv=length(na.omit(sX[,i]))
    pe=PathelementsCpp(sX[,i][1:sv],sY[,i][1:sv])
    
    # bbs=acf(circular(pe$direction),plot=FALSE)
    #llmmms=lm(bbs$lag[2:length(bbs$lag)]~bbs$acf[2:length(bbs$lag)])
    
    #aas=acf(pe$steps,plot=FALSE)
    #llmms=lm(aas$lag[2:length(aas$lag)]~aas$acf[2:length(aas$lag)])
    
    #t2s=(sum((sX[,i][1:(sv-1)]-sX[,i][2:sv])^2))/(sv-1)/+
    #  (sum((sY[,i][1:(sv-1)]-sY[,i][2:sv])^2))/(sv-1)
    #r2s=sd(sX[,i][1:sv])+sd(sY[,i][1:sv])
    
    cts=mean(cos(pe$turns))
    sts=mean(sin(pe$turns))
    #bos=sd(pe$steps)/abs(mean(pe$steps))
    
    
    Ssim[i,]<-c(mean(pe$steps),
                sd(pe$turns),
                cdt2(pe$steps,sT[,i][2:sv]),
                sd(pe$steps),
                #mean(pe$turns),
                #mean(aas$acf),
                sqrt((cts)^2+(sts)^2),
                cdt2(nsd(pe$steps),sT[,i][2:sv]),
    it(pe$steps,sX[,i][1:sv],sY[,i][1:sv]),
    #sd.circular(circular(pe$direction)),
    #angular.deviation(circular(pe$turns)),
    si(cts,sts,mean(pe$steps),bos))
    
    #t2s,
    #r2s)
    
  }
  
}


########################################################################################
##################################### Compare ##########################################
########################################################################################


#hist(Summ[,1])
#hist(SSum[,1])

nbest<-100
#quienes=c(sample(seq(10),5))

#quienes=c(1,2,3,4,5,6,7,8)
st=Stobs[,quienes]
ss=Ssim[,quienes]


# Via distancia de Mahalanobis 
Dists=matrix(NA,length(nonones),NN)
A=cov(na.omit(ss))

#Non=na.omit(SSum[,quienes])
#qq=which(is.finite((apply(Non,1,sum)))==T)

#indice=which(SSum[,1]==max(Non[qq,1]))

for (i in 1:NN)
{
  Dists[,i]=apply(ss,1,mahalanobis,center=st[i,],cov=A)
}


SumMahal=apply(Dists,1,sum)
which=numeric(length(nonones))
which[order(SumMahal)[1:nbest]]=1 

indices=nonones[which==1]


#####################################################################################
########################## Plot posterior and prior #################################
#####################################################################################


#png('/home/sofia/proyecto_doctoral/pruebas/SSM/abc_plots/p1.png', width=8, height=6.5, units="in", res=300)
par(mfrow=c(1,3),mar=c(8,5,4,1))


##### w
density_scale=density(a$w[indices],bw=0.4,from=-0.1,to=10)
xscale <- seq(-0.5, 10.5, length=100)
y_scale <- dunif(xscale,min = 0,max=10)

plot(xscale,y_scale,type='l',main=expression(lambda),ylim=c(0,max(y_scale,max(density_scale$y)))
     ,col='red',xlab = '',ylab = 'Density',cex.axis=1.5,
     cex.main=2,cex.lab=1.5,lwd=2)

abline(v=mean(a$w[indices]),col='dodgerblue3',lwd=2)
abline(v=t_w,col='darkorange3',lwd=2)
lines(density_scale,lwd=2)


##### k
density_scale=density(a$k[indices],bw=3,from=3,to=100)
xscale <- seq(0,105, length=100)
y_scale <- dunif(xscale,min = 3, max=100)

plot(xscale,y_scale,type='l',main=expression(k),
     ylim=c(0,max(y_scale,max(density_scale$y))),xlab='',col='red',ylab = '',
     cex.main=2,cex.lab=1.5,lwd=2,cex.axis=1.5)
abline(v=mean(a$k[indices]),col='dodgerblue3',lwd=2)
abline(v=t_k,col='darkorange3',lwd=2)
lines(density_scale,lwd=2)

##### joint 
plot(a$w[indices],a$k[indices],xlim=c(-2,12),ylim=c(-5,105),pch=3,
     xlab="w",ylab="k",main="Joint",cex.lab=1.5,cex.axis=1.5,cex.main=1.8)
rect(0,100,10,4,border='red',lwd=2)
points(mean(a$w[indices]),mean(a$k[indices]),col='dodgerblue3',pch=16,cex=1.5)
points(t_w,t_k,col='darkorange3',pch=16,cex=1.5)


#Legends
legend(-15,-2,xpd=NA,bty='n',
       horiz=TRUE,legend=c('Prior','True Value')
       ,fill=c('red','darkorange3'),cex=2,xjust=0.5,yjust=2)

legend(-45,-5,xpd=NA,bty='n',
       horiz=TRUE,legend=paste('dt = ',as.character(dt),sep='')
       ,cex=2,xjust=0.5,yjust=2)


legend(-15,-12,xpd=NA,bty='n',
       horiz=TRUE,legend=c('Estimated Mean', 'Estimated Posterior'),
       fill=c('dodgerblue3','black'),cex=2,xjust=0.5,yjust=2)



#dev.off()




