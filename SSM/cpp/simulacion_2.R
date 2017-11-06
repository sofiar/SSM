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
t_w=6
t_k = 10
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

Summ=matrix(NA,NN,10)
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
          mean(aa$acf),
          sqrt((mean(cos(ps$turns)))^2+(mean(sin(ps$turns)))^2),
          #cdt2(nsd(ps$steps),oz$st[2:length(oz$st)]))
          it(ps$steps,oz$sx,oz$sy),
          #sd.circular(circular(ps$direction)),
          #angular.deviation(circular(ps$turns)),
          si(ct,st,mean(ps$steps),bo),
          t2,
          r2)
  
  
}

### for the simulated ones


SSum=matrix(NA,nsims,10)
for (j in 1:nsims)
{
  osx=a$sX[,j]
  osy=a$sY[,j]
  ost=a$sT[,j]
  
   ps=PathelementsCpp(osx,osy)
  
  
  #bb=acf(circular(ps$direction),plot=FALSE)
  #llmmm=lm(bb$lag[2:length(bb$lag)]~bb$acf[2:length(bb$lag)])
  
  aa=acf(na.omit(ps$steps),plot=FALSE)
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
               mean(aa$acf),
               sqrt((mean(cos(ps$turns)))^2+(mean(sin(ps$turns)))^2),
               #  cdt2(nsd(ps$steps),oz$st[2:length(oz$st)]))
               it(ps$steps,osx,osy),
               #sd.circular(circular(ps$direction)),
               #angular.deviation(circular(ps$turns)),
               si(ct,st,mean(ps$steps),bo),
               t2,
               r2)
 
}

#View(SSum)

########################################################################################
##################################### Compare ##########################################
########################################################################################


#hist(Summ[,1])
#hist(SSum[,1])

nbest<-100
quienes=c(1,5)

st=Summ[,quienes]
ss=SSum[,quienes]


# Via distancia de Mahalanobis 
Dists=matrix(NA,nsims,NN)

Non=na.omit(SSum[,quienes])
qq=which(is.finite((apply(Non,1,sum)))==T)

indice=which(SSum[,1]==max(Non[qq,1]))


var(na.omit(SSum[,1]))

for (i in 1:NN)
{
  Dists[,i]=apply(ss,1,mahalanobis,center=st[i,],cov=A)
}


SumMahal=apply(Dists,1,sum)
which=numeric(nsims)
which[order(SumMahal)[1:nbest]]=1 

#####################################################################################
########################## Plot posterior and prior #################################
#####################################################################################

par(mfrow=c(1,3),mar=c(8,5,4,1))


##### w
density_scale=density(a$w[which==1],bw=0.4,from=-0.1,to=10)
xscale <- seq(-0.5, 10.5, length=100)
y_scale <- dunif(xscale,min = 0,max=10)

plot(xscale,y_scale,type='l',main=expression(lambda),ylim=c(0,max(y_scale,max(density_scale$y)))
     ,col='red',xlab = '',ylab = 'Density',cex.axis=1.5,
     cex.main=2,cex.lab=1.5,lwd=2)

abline(v=mean(a$w[which==1]),col='dodgerblue3',lwd=2)
abline(v=t_w,col='darkorange3',lwd=2)
lines(density_scale,lwd=2)


##### k
density_scale=density(a$k[which==1],bw=3,from=3,to=100)
xscale <- seq(0,105, length=100)
y_scale <- dunif(xscale,min = 3, max=100)

plot(xscale,y_scale,type='l',main=expression(k),
     ylim=c(0,max(y_scale,max(density_scale$y))),xlab='',col='red',ylab = '',
     cex.main=2,cex.lab=1.5,lwd=2,cex.axis=1.5)
abline(v=mean(a$k[which==1]),col='dodgerblue3',lwd=2)
abline(v=t_k,col='darkorange3',lwd=2)
lines(density_scale,lwd=2)

##### joint 
plot(a$w[which==1],a$k[which==1],xlim=c(-2,12),ylim=c(-5,105),pch=3,
     xlab="w",ylab="k",main="Joint",cex.lab=1.5,cex.axis=1.5,cex.main=1.8)
rect(0,100,10,4,border='red',lwd=2)
points(mean(a$w[which==1]),mean(a$k[which==1]),col='dodgerblue3',pch=16,cex=1.5)
points(t_w,t_k,col='darkorange3',pch=16,cex=1.5)






#Legends
legend(-115,0.0008,xpd=NA,bty='n',
       horiz=TRUE,legend=c('Prior','True Value')
       ,fill=c('red','darkorange3'),cex=2,xjust=0.5,yjust=2)



legend(-115,-0.005,xpd=NA,bty='n',
       horiz=TRUE,legend=c('Estimated Mean', 'Estimated Posterior'),
       fill=c('dodgerblue3','black'),cex=2,xjust=0.5,yjust=2)





