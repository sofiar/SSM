####### Pruebas
### Pruebo que los tiempos de espera tm sean modelados con una exponencial

library(circular) 
library(mvtnorm)
source('/home/sofia/proyecto_doctoral/pruebas/SSM/funaux.R')

nsteps <- 800 # number of moves performed by the animal

# Parametrs of the model:
#t_scale = 2
t_w=2
t_mu = pi/50
t_k = 20

# simulate the movement
#tm <- rweibull(nsteps, scale=t_scale, shape=t_shape) # time per move
tm<-rexp(nsteps,rate=t_w)
ls <- numeric(nsteps) + 1 # assume constant speed for now
tu <- rvonmises(nsteps, mu=circular(t_mu), kappa=t_k) # turning angle

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

dt <- 2 # time interval for observations
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

nsims <- 1e5#### see
nsam <- 1000 # number of real steps in the simulations
maxt = max(oz$st)
nobs = length(oz$sx)
sX = matrix(NA,length(oz$sx),nsims)
sY = matrix(NA,length(oz$sx),nsims)
sT = matrix(NA,length(oz$sx),nsims)

# sample from priors

#s_mean<- runif(nsims,0,3)
#s_shape <- runif(nsims,0.5,10)
s_w <- runif(nsims,0.1,5)
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
  #tm<-rgamma(nsam,shape=s_shape[j],scale=s_mean[j]/s_shape[j])
  tm<-rexp(nsam,rate=s_w[j])
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

ps=pathelements(oz$sx,oz$sy)

bb=acf(circular(ps$direction),plot=FALSE)
llmmm=lm(bb$lag[2:length(bb$lag)]~bb$acf[2:length(bb$lag)])

aa=acf(ps$steps,plot=FALSE)
llmm=lm(aa$lag[2:length(aa$lag)]~aa$acf[2:length(aa$lag)])

t2=(sum((oz$sx[1:(length(oz$sx)-1)]-oz$sx[2:(length(oz$sx))])^2))/(length(oz$sx)-1)/+
  (sum((oz$sy[1:(length(oz$sy)-1)]-oz$sy[2:(length(oz$sy))])^2))/(length(oz$sy)-1)
r2=sd(oz$sx)+sd(oz$sy)

ct=mean(cos(ps$turns))
st=mean(sin(ps$turns))
bo=sd(ps$steps)/abs(mean(ps$steps))

Stobs<-c(mean(ps$steps),
         sd(ps$turns),
         cdt2(ps$steps,oz$st[2:length(oz$st)]),
         sd(ps$steps),
         mean(ps$turns),
         mean(aa$acf),
         llmm$coefficients[1],
         llmm$coefficients[2],
         sqrt((mean(cos(ps$turns)))^2+(mean(sin(ps$turns)))^2),
         it(ps$steps,oz$sx,oz$sy),
         sd.circular(circular(ps$direction)),
         angular.deviation(circular(ps$turns)),
         llmmm$coefficients[1],
         llmmm$coefficients[2],
         si(ct,st,mean(ps$steps),bo),         
         sum(abs(na.omit(pathelements(oz$sx,oz$sy)$turns)==0))/length(na.omit(oz$sx)),
         t2,
         r2)


Ssim<-data.frame(matrix(nrow=nsims,ncol=length(Stobs),NA))

#names(Ssim)<-c('A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11','A12','A13','A14')
#rs=numeric(nsims)
#names(Ssim)<-c('A2','A3','A6','A12','A14')

for (i in 1:nsims)
{
  sv=length(na.omit(sX[,i]))
  if(sv>2) # if one simulation is too short to observe something sv<=1
  {
    pe=pathelements(sX[,i][1:sv],sY[,i][1:sv])
    
    bbs=acf(circular(pe$direction),plot=FALSE)
    llmmms=lm(bbs$lag[2:length(bbs$lag)]~bbs$acf[2:length(bbs$lag)])
    
    aas=acf(pe$steps,plot=FALSE)
    llmms=lm(aas$lag[2:length(aas$lag)]~aas$acf[2:length(aas$lag)])
    
    t2s=(sum((sX[,i][1:(sv-1)]-sX[,i][2:sv])^2))/(sv-1)/+
      (sum((sY[,i][1:(sv-1)]-sY[,i][2:sv])^2))/(sv-1)
    r2s=sd(sX[,i][1:sv])+sd(sY[,i][1:sv])
    
    cts=mean(cos(pe$turns))
    sts=mean(sin(pe$turns))
    bos=sd(pe$steps)/abs(mean(pe$steps))
    
    
    Ssim[i,]<-c(mean(pe$steps),
                sd(pe$turns),
                cdt2(pe$steps,sT[,i][2:sv]),
                sd(pe$steps),
                mean(pe$turns),
                mean(aas$acf),
                llmms$coefficients[1],
                llmms$coefficients[2],
                sqrt((cts)^2+(sts)^2),
                it(pe$steps,sX[,i][1:sv],sY[,i][1:sv]),
                sd.circular(circular(pe$direction)),
                angular.deviation(circular(pe$turns)),
                llmmms$coefficients[1],
                llmmms$coefficients[2],
                si(cts,sts,mean(pe$steps),bos),         
                sum(abs(na.omit(ps$turns)==0))/sv,
                t2s,
                r2s)
    
    }
  
}
par(mfrow=c(1,3))
## Now we have to decide which indices we keep

nbest<-20
quienes=c(1,2,3,4,9,12)
ss=Stobs[quienes]
s=Ssim[,quienes]

#Mahalanobis
A=cov(na.omit(s))
Mahal=apply(s,1,mahalanobis,center=ss,cov=A)

which=numeric(nsims)
which[order(Mahal)[1:nbest]]=1

par(mfrow=c(1,3))

plot(s_w[which==1],s_mu[which==1],xlim=c(0.1,5),ylim=c(-pi,pi))
points(t_w,t_mu,col='red',pch=21,bg='red')

plot(s_w[which==1],s_k[which==1],xlim=c(0.1,5),ylim=c(5,90))
points(t_w,t_k,col='red',pch=21,bg='red')

plot(s_k[which==1],s_mu[which==1],xlim=c(5,90),ylim=c(-pi,pi))
points(t_k,t_mu,col='red',pch=21,bg='red')

par(mfrow=c(1,2))
ind=12
pp=which(s_k<=30)
p=which(s_k>50)
hist(Ssim[pp,ind])
hist(Ssim[p,ind])







