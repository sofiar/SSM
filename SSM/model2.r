
library(circular) 
library(mvtnorm)
source('/home/sofia/proyecto_doctoral/pruebas/SSM/funaux.R')

TT=450 # final time

# Parametrs of the model:
#t_scale = 2
t_mean=2
t_shape =2
t_mu = pi/50
t_k = 20

# simulate the movement
#tm <- rweibull(nsteps, scale=t_scale, shape=t_shape) # time per move

x <- c(0)
y <- c(0)
t <- c(0)

di <- c(runif(1) * 2*pi) # initial movement direction
tt=0
while(tt<TT){
  
  tm<-rgamma(1,shape=t_shape,scale=t_mean/t_shape)
  ls <- numeric(1) + 1 # assume constant speed for now
  tu <- rvonmises(1, mu=t_mu, kappa=t_k) # turning angle
  
  t<-c(t,t[length(t)]+tm)
  di<-c(di,di[length(di)]+tu)  
  x<-c(x,x[length(x)]+cos(di[length(di)])*ls * tm)
  y<-c(y,y[length(y)]+sin(di[length(di)])*ls * tm)
  
  tt=t[length(t)]
}

# we observe at time dt (pay attention to this parameter)

dt <- 5 # time interval for observations
oz = observe(x,y,t,dt) # these are the observed data


# plot true trajectory and the observed one

plot(x,y,type="p", asp=1)
points(oz$sx,oz$sy,col=2, pch=16)
lines(oz$sx,oz$sy,col=2)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# ABC
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


nsims <- 10e3#### see
nsam <- 1000 # number of real steps in the simulations
maxt = max(oz$st)
nobs = length(oz$sx)
h<-0.1
sX = matrix(NA,length(oz$sx),nsims)
sY = matrix(NA,length(oz$sx),nsims)
sT = matrix(NA,length(oz$sx),nsims)

# sample from priors
#s_scale <- runif(nsims,0,5)
s_mean<- runif(nsims,0,3)
s_shape <- runif(nsims,0.5,4)
s_mu <- runif(nsims,-pi,pi)
s_k <- runif(nsims,5,90)

#-------------------------------------------------------------------------------
# simulate movement and observations
# note that there's only one simulation per parameter combination but we probably need several!
#-------------------------------------------------------------------------------

for( j in 1:nsims){
  
  sx <- numeric(nsam)
  sy <- numeric(nsam)
  st <- numeric(nsam)
  sdi <- numeric(nsam)
  sdi[1] <- runif(1)*2*pi #angulo inicial
  #tm <- rweibull(nsam, scale=s_scale[j], shape=s_shape[j])
  tm<-rgamma(nsam,shape=s_shape[j],scale=s_mean[j]/s_shape[j])
  tur <- rvonmises(nsam, mu=circular(s_mu[j]), kappa=s_k[j]) 
  i=1
  ## simulation of the trajectory
  while(st[i]<maxt & i<nsam){
    i=i+1
    st[i] <- st[i-1] + tm[i]
    sdi[i] <- sdi[i-1]+tur[i]
    sx[i] <- sx[i-1] + cos(sdi[i]) * tm[i]
    sy[i] <- sy[i-1] + sin(sdi[i]) * tm[i]
  }
  # observe the movement process at time interval dt
  sz = observe(sx[1:i],sy[1:i],st[1:i],dt)
  sobs = length(sz$sx)
  sv <- min(nobs,sobs) ### ???? Para que ????
  
  sX[1:sv,j] <- sz$sx[1:sv]
  sY[1:sv,j] <- sz$sy[1:sv]
  sT[1:sv,j] <- sz$st[1:sv]
  
} 

#-------------------------------------------------------------------------------
## summaries
#-------------------------------------------------------------------------------

## Calculation of the summaries for the TRUE observation

Stobs<-c(sum(pathelements(oz$sx,oz$sy)$steps),mean(pathelements(oz$sx,oz$sy)$steps),sd(pathelements(oz$sx,oz$sy)$turns),
         cdt(oz$sx,oz$sy,oz$st),cha(oz$sx,oz$sy) ,sd(pathelements(oz$sx,oz$sy)$steps),mean(pathelements(oz$sx,oz$sy)$turns),
         sqrt(abs(max(oz$sx)-min(oz$sx))+abs(max(oz$sy)-min(oz$sy))),acf(pathelements(oz$sx,oz$sy)$turn,plot=FALSE)$acf[5],
         sd(ppp(oz$sx,oz$sy,oz$st)/dt),sd(pathelements(oz$sx,oz$sy)$direction),
         cdt2(pathelements(oz$sx,oz$sy)$steps,oz$st[2:length(oz$st)]),
         mean(ppp(oz$sx,oz$sy,oz$st)/dt),
         sum(abs(na.omit(pathelements(oz$sx,oz$sy)$turns)==0))/length(na.omit(oz$sx)))
Stobs<-Stobs[c(2,3,6,12,14)]


Ssim<-data.frame(matrix(nrow=nsims,ncol=length(Stobs),NA))

#names(Ssim)<-c('A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11','A12','A13','A14')
#rs=numeric(nsims)
names(Ssim)<-c('A2','A3','A6','A12','A14')

for (i in 1:nsims)
{
  sv=length(na.omit(sX[,i]))
  if(sv>2) # if one simulation is too short to observe something sv<=1
  {

    
    pe=pathelements(sX[,i][1:sv],sY[,i][1:sv])
    Ssim[i,]<-c(mean(pe$steps),
                sd(pe$turns),
                sd(pe$steps),
                cdt2(pe$steps,sT[,i][1:sv][2:length(sT[,i][1:sv])]),
                sum(abs(na.omit(pe$turns)==0))/length(na.omit(sX[,i]))
    )
    

  }
  
}

## Now we have to decide which indices we keep

s=Ssim
ss=Stobs

nbest<-30

Ssim=cbind(s$A2,s$A3,s$A6,s$A12)
#Stobs=c(ss[1:4])

#Mahalanobis
A=cov(na.omit(Ssim))
#XX=(as.matrix(Ssim[2,]-Stobs)%*%solve(A)%*%as.matrix(t(Ssim[2,]-Stobs)))
Mahal=apply(Ssim,1,mahalanobis,center=Stobs,cov=A)

which=numeric(nsims)
which[order(Mahal)[1:nbest]]=1

#which[order(abs(s$A10-ss[10]))[1:nbest]]=1
############################ PLOT MARGINS ##################################
b=layout(matrix(c(1,2,3,4,5,5), ncol=2, byrow=TRUE), heights=c(3,3,1))


### scale
hist_scale=hist(s_scale[which==1],breaks = 10,plot=FALSE)
density_scale=density(s_scale[which==1])
xscale <- seq(-1, 11, length=100)
y_scale <- dunif(xscale,min = 0,max=5)

plot(xscale,y_scale,type='l',main='Parameter: Scale',ylim=c(0,max(y_scale,max(density_scale$y))),col='darkseagreen4',xlab = 'scale',ylab = 'density')
abline(v=mean(s_scale[which==1]),col='dodgerblue3')
abline(v=t_scale,col='darkorange3')
lines(density(s_scale[which==1]))
#legend(x=1.2,y=2.2,legend=c('prior','true value','estimated value', 'estimated porterior')
#       ,fill=c('darkseagreen4','darkorange3','dodgerblue3','black'))


### mean
hist_mean=hist(s_mean[which==1],breaks = 10,plot=FALSE)
density_mean=density(s_mean[which==1])
xmean <- seq(-1, 11, length=100)
y_mean <- dunif(xmean,min = 0,max=5)

plot(xmean,y_mean,type='l',main='Parameter: Mean',ylim=c(0,max(y_mean,max(density_mean$y))),col='darkseagreen4',xlab = 'mean',ylab = 'density')
abline(v=mean(s_mean[which==1]),col='dodgerblue3')
abline(v=t_mean,col='darkorange3')
lines(density(s_mean[which==1]))



### shape
hist_shape=hist(s_shape[which==1],breaks = 10,plot=FALSE)
density_shape=density(s_shape[which==1])
xshape <- seq(0, 4.5, length=100)
y_shape <- dunif(xshape,min = 0.5,max=4)

plot(xshape,y_shape,type='l',main='Parameter: Shape',ylim=c(0,max(y_shape,max(density_shape$y))),col='darkseagreen4',xlab = 'shape',ylab = 'density')
abline(v=mean(s_shape[which==1]),col='dodgerblue3')
abline(v=t_shape,col='darkorange3')
lines(density(s_shape[which==1]))


### k
hist_k=hist(s_k[which==1],breaks = 10,plot=FALSE)
density_k=density(s_k[which==1])
xk <- seq(-1, 101, length=100)
y_k <- dunif(xk,min = 0,max=100)

plot(xk,y_k,type='l',main='Parameter: k',ylim=c(0,max(y_k,max(density_k$y))),col='darkseagreen4',xlab = 'k',ylab = 'density')
abline(v=mean(s_k[which==1]),col='dodgerblue3')
abline(v=t_k,col='darkorange3')
lines(density(s_k[which==1]))

### mu
hist_mu=hist(s_mu[which==1],breaks = 10,plot=FALSE)
density_mu=density(s_mu[which==1])
xmu <- seq(-3/2*pi,3/2*pi, length=100)
y_mu <- dunif(xmu,min =-pi,max=pi)

plot(xmu,y_mu,type='l',main='Parameter: mu',ylim=c(0,max(y_k,max(density_mu$y))),col='darkseagreen4',xlab = 'mu',ylab = 'density')
abline(v=mean(s_mu[which==1]),col='dodgerblue3')
abline(v=t_mu,col='darkorange3')
lines(density(s_mu[which==1]))

par(mai=c(0,0,0,0))
plot.new()
legend(x='center',ncol=4,inset=0,legend=c('prior','true value','estimated value', 'estimated porterior')
       ,fill=c('darkseagreen4','darkorange3','dodgerblue3','black'))
dev.off()


#### funciones de densidad de a dos 
op<-par(mfrow=c(3,2))

plot(s_mu[which==1],s_k[which==1],xlab='mu',ylab='k',xlim = c(-pi,pi),ylim=c(5,90))
points(t_mu,t_k,col='red',pch=17)

plot(s_mu[which==1],s_mean[which==1],xlab='mu',ylab='mean',xlim = c(-pi,pi),ylim=c(0,3))
points(t_mu,t_mean,col='red',pch=17)

plot(s_mu[which==1],s_shape[which==1],xlab='mu',ylab='shape',xlim = c(-pi,pi),ylim=c(0.5,4))
points(t_mu,t_shape,col='red',pch=17)

plot(s_shape[which==1],s_mean[which==1],xlab='shape',ylab='mean',xlim = c(0.5,4),ylim=c(0,3))
points(t_shape,t_mean,col='red',pch=17)

plot(s_shape[which==1],s_k[which==1],xlab='shape',ylab='k',xlim = c(0.5,4),ylim=c(5,90))
points(t_shape,t_k,col='red',pch=17)

plot(s_k[which==1],s_mean[which==1],xlab='k',ylab='mean',xlim=c(5,90),ylim = c(0,3))
points(t_k,t_mean,col='red',pch=17)

par(op)

### vamos a ver como son los summaries que tiene valores de scale cercanos al verdadero

aa=which(s_shape<(t_shape)+1/2 & s_shape>(t_shape)-1/2)
bb=which(s_shape>2)

#### prueba de nuevos summaries

val=rep(NA,nsims)

for (i in 1:nsims)
{
  sv=length(na.omit(sX[,i]))
  if(sv>2) # if one simulation is too short to observe something sv<=1
  {
    
    #st=pathelements(sX[,i],sY[,i])$steps
    #stt=numeric(length (st))
    #  for (m in 1:length(st))
    #  {
    #    stt[m]=sum(na.omit(st[1:m]))
    #  }
    #val[i]<-sd(na.omit(pathelements(sX[,i],sY[,i])$direction))
    #val[i]<-mean(na.omit(sX[,i]^2+sY[,i]^2))/sd(na.omit(sX[,i]^2+sY[,i]^2))        
    #rs[i]=rr(pathelements(oz$sx,oz$sy)$steps,pathelements(sX[,i][1:sv],sY[,i][1:sv])$steps)
    #val[i]=  sqrt(abs(max(sX[,i])-min(sX[,i]))+abs(max(sY[,i])-min(sY[,i])))
    #  val[i]=sd(stt/seq(1:length(stt)*dt))
    
    #val[i]=sum(abs(na.omit(pathelements(sX[,i],sY[,i])$turns)==0))/length(na.omit(sX[,i]))
    
    #val[i]=sd(pathelements(sX[,i][1:sv],sY[,i][1:sv])$turn)
    val[i]=sum(pathelements(sX[,i][1:sv],sY[,i][1:sv])$turn<mean(pathelements(sX[,i][1:sv],sY[,i][1:sv])$turn))/sv
    
  }
  
}


#val=s$val
par(mfrow=c(1,2))

vp=sum(abs(pathelements(oz$sx,oz$sy)$turns)==0)/length(na.omit(oz$sx))
hist(val[aa],freq=FALSE,xlim = c(min(na.omit(val)),max(na.omit(val))))
abline(v=mean(na.omit(val[aa])),col='red')
#abline(v=vp,col='blue')

hist(val[bb],freq=FALSE,xlim = c(min(na.omit(val)),max(na.omit(val))))
abline(v=mean(na.omit(val[-aa])),col='red')
abline(v=vp,col='blue')

hist(val[bb],freq=FALSE,xlim = c(min(na.omit(val)),max(na.omit(val))))
abline(v=mean(na.omit(val[-aa])),col='red')




