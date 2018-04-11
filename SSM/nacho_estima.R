
# estimamos con abc
library(circular)
library(Rcpp)
library(tidyverse)
library(mvtnorm)
library(abc)
library(randomForest)

source('SSM/funaux.R')
Rcpp::sourceCpp('SSM/cpp/abc_crw.cpp')
Rcpp::sourceCpp('SSM/cpp/PathelementsCpp.cpp')
Rcpp::sourceCpp('SSM/cpp/RW_exp_cor.cpp')
Rcpp::sourceCpp('SSM/cpp/cppObs.cpp')

# 1.2 funciones para calcular stat y simular datos
# funcion que simula datos a partir del valor de: w, k, dt, nsteps
simdata <- function(ws, ks, dt, nsteps) {
  cppRW_exp_cor(k=ks, w= ws, ns=nsteps, maxx= Inf) %>% # simulate RW using the cpp function
    with( cppObs(x=x, y=y, t=t, dt=dt) ) # these are the observed data
}
# funcion que calcula los estadisticos a un conjunto de datos
# no me queda claro para que se usa el 'nobs' 
stat_fun <- function( dd ) {
  ps = PathelementsCpp(dd$sx,dd$sy)
  bb = acf(circular(ps$direction),plot=FALSE)
  ct = mean(cos(ps$turns))
  st = mean(sin(ps$turns))
  bo = sd(ps$steps)#/abs(mean(ps$steps))
  # aa = acf(ps$steps,plot=FALSE)
  
  tr2 =  mean( diff(dd$sx)^2 ) /+ mean( diff(dd$sy)^2 )
  r2 = sd(dd$sx)+sd(dd$sy)
  mx <- which.max( abs(bb$acf[-1]) )
  
  sale <-c(mean(ps$steps),
           sd(ps$turns),
           cdt2(ps$steps,dd$st[2:length(dd$st)]),
           sd(ps$steps),
           mean(bb$acf[2:6]),
           bb$acf[mx+1],
           sqrt((mean(cos(ps$turns)))^2+(mean(sin(ps$turns)))^2),
           it(ps$steps,dd$sx,dd$sy),
           si(ct,st,mean(ps$steps),bo),
           tr2
  )
  # no uso r2 mean(bb$acf), 
  # incluyo corr en los valores iniciales y el lag de la max corr
  return(sale)
}

# un solo conjunto de datos 'reales'

nsteps <- 800 # number of moves performed by the animal
# movement parameters:
t_w <- 2
t_k <- 20
dt  <- 2 # time interval for observations
data.obs <- simdata(ws=t_w, ks=t_k, dt=dt, nsteps=nsteps)
stat.obs <- stat_fun(data.obs) %>% set_names(nm = paste('T', 1:10, sep=''))

# Usamos dos previas, para cada previa podemos usar, los summary stat o un rf de los summary stat
statsim.unif <- readRDS('SSM/statsim_unif.rds')
statsim.norm <- readRDS('SSM/statsim_norm.rds')

out.unif <- abc(target = stat.obs, 
              param = statsim.unif[, 1:2], 
              sumstat = statsim.unif[, 3:12] , tol = .5, 
              method='neuralnet')
summary(out.unif)
plot(out.unif, param = statsim.unif[, 1:2])


rf.unif.w <- randomForest( s_w ~ . , data = dplyr::select(statsim.unif, -s_k) )
rf.unif.k <- randomForest( s_k ~ . , data = dplyr::select(statsim.unif, -s_w) )
rf.statsim.unif <- data_frame(s_w = predict(rf.unif.w), s_k = predict(rf.unif.k) )

rf.stat.obs <- c( predict(rf.unif.w, newdata = t( as.data.frame(stat.obs) ) ),
                  predict(rf.unif.k, newdata = t( as.data.frame(stat.obs) ) ) )

out.unif.rf <- abc(target = rf.stat.obs, 
                param = statsim.unif[, 1:2], 
                sumstat = rf.statsim.unif , tol = .5, 
                method='neuralnet')

summary(out.unif.rf)
#-----


sd.norm <- statsim.norm %>% dplyr::select(T1:T10) %>% apply( 2, sd )
# scale(statsim.norm[, 3:12] ,center=FALSE)
out.norm <- abc(target = stat.obs, 
                param = statsim.norm[, 1:2], 
                sumstat = statsim.norm[,3:12] , tol = .2, 
                method='neuralnet')
summary(out.norm)
plot(out.norm, param = statsim.norm[, 1:2])

statsim.norm %>% 
  gather(stat, stat.val, -s_w, -s_k) %>% 
  gather(param, param.val, s_w, s_k) %>% 
  ggplot(aes(stat.val, param.val)) + geom_hex() +
  facet_wrap(param~stat, scales = 'free')



  

