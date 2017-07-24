####### Pruebas

#### verdadera

library(circular) 
library(mvtnorm)
source('/home/sofia/proyecto_doctoral/pruebas/SSM/funaux.R')

nsteps <- 500 # number of moves performed by the animal

# Parametrs of the model:
#t_scale = 2
t_mean=2
t_shape =2
t_mu = pi/30
t_k =20

# simulate the movement
tu<-numeric(nsteps)
tu[1]= rvonmises(1, mu=t_mu, kappa=t_k) # first angle
tm<-rgamma(nsteps,shape=t_shape,scale=t_mean/t_shape)
ls <- numeric(nsteps) + 1 # assume constant speed for now

#tm <- rweibull(nsteps, scale=t_scale, shape=t_shape) # time per move


x <- numeric(nsteps)
y <- numeric(nsteps)
t <- numeric(nsteps)
di <- numeric(nsteps)
di[1] <- runif(1) * 2*pi # initial movement direction

for(i in 2:nsteps){
  tu[i] <- rvonmises(1, mu=tu[i-1], kappa=t_k) # turning angle
  
  t[i] <- t[i-1] + tm[i]
  di[i] <- di[i-1]+tu[i]
  x[i] <- x[i-1] + cos(di[i]) * (ls[i]) * tm[i]
  y[i] <- y[i-1] + sin(di[i]) * (ls[i]) * tm[i]
}

# we observe at time dt (pay attention to this parameter)

dt <- 5 # time interval for observations
oz = observe(x,y,t,dt) # these are the observed data



plot(x,y,type="p", asp=1)
points(oz$sx,oz$sy,col=2, pch=16)
lines(oz$sx,oz$sy,col=2)

ps=pathelements(oz$sx,oz$sy)

bb=acf(as.circular(ps$turns))
llmm=lm(bb$lag[2:length(bb$lag)]~bb$acf[2:length(bb$lag)])
plot(bb$lag[2:length(bb$lag)],bb$acf[2:length(bb$lag)])


aa=acf(ps$steps)
llmm=lm(aa$lag[2:length(aa$lag)]~aa$acf[2:length(aa$lag)])



### funcion de autocorrelacion turns)
bb=acf(ps$turns,plot=FALSE)
llmmm=lm(bb$lag[2:length(bb$lag)]~bb$acf[2:length(bb$lag)])


posta=c(mean(ps$steps),
        sd(ps$turns),
        cdt2(ps$steps,oz$st[2:length(oz$st)]),
        sd(ps$steps),
        sum(abs(ps$turns-mean(ps$turns)))/length(na.omit(ps$turns)),
        mean(ps$turns),
        mean(aa$acf),
        llmm$coefficients[1],
        llmm$coefficients[2],
        cha(oz$sx,oz$sy),
        mean(bb$acf),
        llmmm$coefficients[1],
        llmmm$coefficients[2])


#################################   Simulaciones   #######################################
nsims=1e5
T_mean = runif(nsims,0.5,4)
T_shape = runif(nsims,1,10)
T_mu = runif(nsims,-pi,pi)
T_k = runif(nsims,1,100)
val1=numeric(nsims)
val2=numeric(nsims)
val3=numeric(nsims)
val4=numeric(nsims)
val5=numeric(nsims)
val6=numeric(nsims)
val7=numeric(nsims)
val8=numeric(nsims)
val9=numeric(nsims)
val10=numeric(nsims)
val11=numeric(nsims)
val12=numeric(nsims)
val13=numeric(nsims)
# simulate the movement
for (k in 1:nsims)
{
t_shape=T_shape[k]  
t_mean=T_mean[k]  
t_mu=T_mu[k]
t_k=T_k[k]

tm<-rgamma(nsteps,shape=t_shape,scale=t_mean/t_shape)
ls <- numeric(nsteps) + 1 # assume constant speed for now
tu=numeric(nsteps)
tu[1]= rvonmises(1, mu=t_mu, kappa=t_k) # first angle


x <- numeric(nsteps)
y <- numeric(nsteps)
t <- numeric(nsteps)
di <- numeric(nsteps)
di[1] <- runif(1) * 2*pi # initial movement direction

for(i in 2:nsteps){
  tu[i] <- rvonmises(1, mu=tu[i-1], kappa=t_k) # turning angle
  t[i] <- t[i-1] + tm[i]
  di[i] <- di[i-1]+tu[i]
  x[i] <- x[i-1] + cos(di[i]) * (ls[i]) * tm[i]
  y[i] <- y[i-1] + sin(di[i]) * (ls[i]) * tm[i]
}

# we observe at time dt (pay attention to this parameter)

oz = observe(x,y,t,dt) # these are the observed data
ps=pathelements(oz$sx,oz$sy)
val1[k]=mean(ps$steps)
val2[k]=sd(ps$turns)
val3[k]=cdt2(ps$steps,oz$st[2:length(oz$st)])
val4[k]=sd(ps$steps)
val5[k]=sum(abs(ps$turns-mean(ps$turns)))/length(na.omit(ps$turns))
val6[k]=mean(ps$turns)

### funcion de autocorrelacion steps (para el k)
aa=acf(ps$steps,plot=FALSE)
llmm=lm(aa$lag[2:length(aa$lag)]~aa$acf[2:length(aa$lag)])
val7[k]=mean(aa$acf)
val8[k]=llmm$coefficients[1]
val9[k]=llmm$coefficients[2]

val10[k]=cha(oz$sx,oz$sy)

### funcion de autocorrelacion turns)
bb=acf(ps$turns,plot=FALSE)
llmm=lm(bb$lag[2:length(bb$lag)]~bb$acf[2:length(bb$lag)])
val11[k]=mean(bb$acf)
val12[k]=llmm$coefficients[1]
val13[k]=llmm$coefficients[2]



}

############ Comparando summaries a ojo ###################

##### comparacion por distancia

Summaries=cbind(val1,val2,val3,val4,val6,val7,val8,val9)
A=cov(na.omit(Summaries))
Mahal=apply(Summaries,1,mahalanobis,center=posta[c(1,2,3,4,6,7,8,9)],cov=A)

which=numeric(nsims)
which[order(Mahal)[1:70]]=1

par(mfrow=c(1,3))

plot(T_shape[which==1],T_mean[which==1],xlim=c(1,10),ylim=c(0.5,4))
points(2,2,col='red',pch=21,bg='red')

plot(T_shape[which==1],T_k[which==1],xlim=c(1,10),ylim=c(1,100))
points(2,20,col='red',pch=21,bg='red')

plot(T_mean[which==1],T_k[which==1],xlim=c(0.5,4),ylim=c(1,100))
points(2,20,col='red',pch=21,bg='red')


######################### plot 3D ######################################
library(scatterplot3d)
par(mfrow=c(1,1),pty="s")

s3d=scatterplot3d(T_shape[which==1],T_mean[which==1],T_k[which==1],xlim=c(1,10)
              ,ylim=c(0.5,4),zlim=c(1,100),xlab='shape',ylab='                  mean',zlab='k',box=FALSE)
s3d$points3d(2,2,20, col="red", type="p", pch=16)


### Kmeans

q=which(kmeans(rbind(posta[c(1,2,3,4,6,7,8,9)],Summaries),1500)$cluster==1)
which=numeric(nsims)
which[q-1]=1

plot(T_shape[which==1],T_mean[which==1],xlim=c(1,10),ylim=c(0.5,4))
points(2,2,col='red')

plot(T_shape[which==1],T_k[which==1],xlim=c(1,10),ylim=c(1,100))
points(2,20,col='red')


