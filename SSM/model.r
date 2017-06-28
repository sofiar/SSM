#### PRUEBAS ssm

library(circular) 
library(mvtnorm)
source('/home/sofia/proyecto_doctoral/pruebas/SSM/funaux.R')

nsteps <- 500 # number of moves performed by the animal

# Parametrs of the model:
t_scale = 1
t_shape = 2
t_mu = pi/100
t_k = 20

# simulate the movement
tm <- rweibull(nsteps, scale=t_scale, shape=t_shape) # time per move
ls <- numeric(nsteps) + 1 # assume constant speed for now
tu <- rvonmises(nsteps, mu=t_mu, kappa=t_k) # turning angle

x <- numeric(nsteps)
y <- numeric(nsteps)
t <- numeric(nsteps)
di <- numeric(nsteps)
di[1] <- runif(1) * 2*pi # initial movement direction

for(i in 2:nsteps){
  t[i] <- t[i-1] + tm[i]
  di[i] <- di[i-1]+tu[i]
  x[i] <- x[i-1] + cos(di[i]) * (ls[i]) * tm[i]
  y[i] <- y[i-1] + sin(di[i]) * (ls[i]) * tm[i]
}

# we observe at time dt (pay attention to this parameter)

dt <- 10 # time interval for observations
oz = observe(x,y,t,dt) # these are the observed data


# plot true trajectory and the observed one

plot(x,y,type="l", asp=1)
points(oz$sx,oz$sy,col=2, pch=16)
lines(oz$sx,oz$sy,col=2)
sqrt(abs(max(oz$sx)-min(oz$sx))+abs(max(oz$sy)-min(oz$sy)))


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# ABC
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


nsims <- 10e4#### see
nsam <- 5000 # number of real steps in the simulations
maxt = max(oz$st)
nobs = length(oz$sx)
h<-0.1
sX = matrix(NA,length(oz$sx),nsims)
sY = matrix(NA,length(oz$sx),nsims)
sT = matrix(NA,length(oz$sx),nsims)

# sample from priors
s_scale <- runif(nsims,0,5)
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
  tm <- rweibull(nsam, scale=s_scale[j], shape=s_shape[j])
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
         sd(ppp(oz$sx,oz$sy,oz$st)/dt),sd(pathelements(oz$sx,oz$sy)$direction))



Ssim<-data.frame(matrix(nrow=nsims,ncol=length(Stobs),NA))
names(Ssim)<-c('A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11')
rs=numeric(nsims)
for (i in 1:nsims)
{
  sv=length(na.omit(sX[,i]))
  if(sv>2) # if one simulation is too short to observe something sv<=1
  {
    Ssim[i,]<-c(sum(pathelements(sX[,i][1:sv],sY[,i][1:sv])$steps), 
                mean(pathelements(sX[,i][1:sv],sY[,i][1:sv])$steps),
                sd(pathelements(sX[,i][1:sv],sY[,i][1:sv])$turns),
                cdt(sX[,i],sY[,i],sT[,i]),
                cha(sX[,i],sY[,i]),
                sd(pathelements(sX[,i][1:sv],sY[,i][1:sv])$steps),
                mean(pathelements(sX[,i][1:sv],sY[,i][1:sv])$turns),
                sqrt(abs(max(sX[,i])-min(sX[,i]))+abs(max(sY[,i])-min(sY[,i]))),
                acf(na.omit(pathelements(sX[,i],sY[,i])$turn),plot = FALSE)$acf[5],
                sd(ppp(sX[,i],sY[,i],sT[,i])/dt))
    
    rs[i]=rr(pathelements(oz$sx,oz$sy)$steps,pathelements(sX[,i][1:sv],sY[,i][1:sv])$steps)
  }
  
}

## Now we have to decide which indices we keep

s=Ssim
ss=Stobs

nbest<-30

Ssim=cbind(s$A1,s$A11)
Stobs=cbind(ss[1],ss[11])



#A (Maybe we migth use a different set of data to estimate this mean and sd)

NSsim=(Ssim-apply(na.omit(Ssim),2,mean))/apply(na.omit(Ssim),2,sd)
NStobs=(Stobs-apply(na.omit(Ssim),2,mean))/apply(na.omit(Ssim),2,sd)

#uniform
pp=(abs(NSsim-NStobs)<1)*1
which=(apply(pp,1,prod))
which[is.na(which)]=0 
sum(which==1)

## Dist
Dist=apply(abs(s/ss-1),1,max)
which=numeric(nsims)
which[order(Dist)[1:nbest]]=1

#Mahalanobis
A=cov(na.omit(Ssim))
#XX=(as.matrix(Ssim[2,]-Stobs)%*%solve(A)%*%as.matrix(t(Ssim[2,]-Stobs)))
Mahal=apply(Ssim,1,mahalanobis,center=Stobs,cov=A)

which=numeric(nsims)
which[order(Mahal)[1:nbest]]=1

#which[order(abs(s$A10-ss[10]))[1:nbest]]=1

## Dist 2
d2=sqrt(apply((s/sd(na.omit(s))-ss/sd(na.omit(s)))^2,1,sum))
which=numeric(nsims)
which[order(d2)[1:nbest]]=1

#Kmeans
sss=rbind(s,ss)
km=kmeans(na.omit(sss),20,nstart = 25)
quienes=apply(is.na(s),1,sum)
a=s_scale[-which(quienes==1)][which(km$cluster==km$cluster[1])]

hist_scale=hist(a,breaks = 10,plot=FALSE)
xscale <- seq(-1, 11, length=100)
y_scale <- dunif(xscale,min = 0,max=10)

plot(xscale,y_scale,type='l',main='Parameter: Scale',ylim=c(0,max(y_scale,max(density_scale$y))),col='darkseagreen4',xlab = 'scale',ylab = 'density')
abline(v=mean(a),col='dodgerblue3')
abline(v=t_scale,col='darkorange3')
lines(density(a))




### de a una variable

which=numeric(nsims)
which[order(abs(val-ss[11]))[1:nbest]]=1

############################ PLOT MARGINS ##################################
b=layout(matrix(c(1,2,3,4,5,5), ncol=2, byrow=TRUE), heights=c(3,3,1))


### scale
hist_scale=hist(s_scale[which==1],breaks = 10,plot=FALSE)
density_scale=density(s_scale[which==1])
xscale <- seq(-1, 11, length=100)
y_scale <- dunif(xscale,min = 0,max=10)

plot(xscale,y_scale,type='l',main='Parameter: Scale',ylim=c(0,max(y_scale,max(density_scale$y))),col='darkseagreen4',xlab = 'scale',ylab = 'density')
abline(v=mean(s_scale[which==1]),col='dodgerblue3')
abline(v=t_scale,col='darkorange3')
lines(density(s_scale[which==1]))
#legend(x=1.2,y=2.2,legend=c('prior','true value','estimated value', 'estimated porterior')
#       ,fill=c('darkseagreen4','darkorange3','dodgerblue3','black'))

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




### vamos a ver como son los summaries que tiene valores de scale cercanos al verdadero

aa=which(s_scale<(t_scale)+1/4 & s_scale>(t_scale)-1/4 & s_shape<(t_shape)+1/2 & s_shape>(t_shape)-1/2
         & s_mu<(t_mu)+1/2 & s_mu>(t_mu)-1/2)
bb=which(s_scale>4)


hist(s_scale[aa])
hist(s_mu[aa])


val=s$A9
val=rep(NA,nsims)
for (i in 1:nsims)
{
  sv=length(na.omit(sX[,i]))
  if(sv>2) # if one simulation is too short to observe something sv<=1
  {
    #val[i]<-sd(na.omit(pathelements(sX[,i],sY[,i])$direction))
    #val[i]<-mean(na.omit(sX[,i]^2+sY[,i]^2))/sd(na.omit(sX[,i]^2+sY[,i]^2))        
    #rs[i]=rr(pathelements(oz$sx,oz$sy)$steps,pathelements(sX[,i][1:sv],sY[,i][1:sv])$steps)
    val[i]=  sqrt(abs(max(sX[,i])-min(sX[,i]))+abs(max(sY[,i])-min(sY[,i])))
    }
  
}




for (k in 1:50)
{
plot(sX[,bb[k]],sY[,bb[k]],col=3, pch=16,xlab=as.character(s_scale[bb[k]]),
     main=as.character(sqrt(abs(max(sX[,bb[k]])-min(sX[,bb[k]]))+abs(max(sY[,bb[k]])-min(sY[,bb[k]])))))
lines(sX[,bb[k]],sY[,bb[k]],col=3)


plot(sX[,aa[k]],sY[,aa[k]],col=2, pch=16,xlab=as.character(s_scale[aa[k]]),
     main=as.character(sqrt(abs(max(sX[,aa[k]])-min(sX[,aa[k]]))+abs(max(sY[,aa[k]])-min(sY[,aa[k]])))))
lines(sX[,aa[k]],sY[,aa[k]],col=2)


}

par(mfrow=c(1,2))

hist(val[aa],freq=FALSE)
abline(v=mean(val[aa]),col='red')
hist(val,freq=FALSE)
abline(v=mean(na.omit(val)),col='red')



#sqrt(abs(max(oz$sx)-min(oz$sx))+abs(max(oz$sy)-min(oz$sy)))








