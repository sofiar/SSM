# Intento de nacho de hacer una simulacion

# 1) Funciones

# 1.1 las funciones de cpp son las mismas 
#---------

library(circular)
library(Rcpp)
library(tidyverse)
library(mvtnorm)
library(abc)
# source('/home/sofia/proyecto_doctoral/pruebas/SSM/funaux.R')
# setwd("/home/sofia/proyecto_doctoral/pruebas/SSM/cpp")
# 
# Rcpp::sourceCpp('abc_crw.cpp')
# Rcpp::sourceCpp('PathelementsCpp.cpp')
# Rcpp::sourceCpp('RW_exp_cor.cpp')
# Rcpp::sourceCpp('cppObs.cpp')

source('SSM/funaux.R')

Rcpp::sourceCpp('SSM/cpp/abc_crw.cpp')
Rcpp::sourceCpp('SSM/cpp/PathelementsCpp.cpp')
Rcpp::sourceCpp('SSM/cpp/RW_exp_cor.cpp')
Rcpp::sourceCpp('SSM/cpp/cppObs.cpp')

# 1.2 funciones para calcular stat y simular datos
# funcion que simula datos a partir del valor de: w, k, dt, nsteps
simdata <- function(ws, ks, dt.list, nsteps) {
  xx <- cppRW_exp_cor(k=ks, w= ws, ns=nsteps, maxx= Inf) # simulate RW using the cpp function
  names(dt.list) <- paste('dt', dt.list, sep='_')
  lapply(dt.list, function(z) with(xx, cppObs(x=x, y=y, t=t, dt=z) ) ) # these are the observed data
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

#stat_fun(oz)

# funcion que simula un dato y calcula el stat
f1 <- function( par, dt.list, ns ) {
  simdata(ws=par[1], ks=par[2], dt.list = dt.list, nsteps=ns) %>%
   lapply( stat_fun )
}

# f1( par=c(2, 20), dt.list=list(1,2,4), ns=500 ) %>% str( max.level =1)

# 2) Simulacion de conjunta (theta, s)
#  Uso dos previas: 
#       pr1 = independientes con marginales uniformes
#       pr2 = normales con dependencia positiva

nsims <- 10e3
pr.unif <- data_frame(s_w = runif(nsims,0.1,10), s_k = runif(nsims,5,100) )

# una forma rapida de hacer las normales positivas
sds <- c(2, 15) 
rr  <- matrix( c(1, .4, .4, 1 ), ncol=2, byrow = TRUE)
pr.norm <- rmvnorm(nsims, c(5, 50), sigma = diag(sds) %*% rr %*% diag(sds) ) %>%
  as.data.frame() %>% mutate(s_w = abs(V1), s_k = abs(V2)) %>%
  dplyr::select(s_w, s_k)
# pr.norm %>% ggplot() + geom_point(aes(s_w, s_k) )

# obtener y salvar stats for simulated data (cada fila un stat, cada columna una sim)
apply( pr.unif, 1,  
             function(pp) {
               f1(par=pp, dt.list = list( .05, .5, 5 ) , ns = 800) %>%
                 bind_rows() %>%
                 mutate(s_w = pp[1], s_k=pp[2], stat.nm = paste('T', 1:10, sep=''))
             }) %>% 
  bind_rows() %>%
  gather(dt, stat.val, starts_with('dt') ) %>%
  spread(stat.nm, stat.val) %>%
  saveRDS(file='SSM/statsim_unif.rds')

apply( pr.norm, 1,  
       function(pp) {
         f1(par=pp, dt.list = list( .05, .5, 5 ) , ns = 800) %>%
           bind_rows() %>%
           mutate(s_w = pp[1], s_k=pp[2], stat.nm = paste('T', 1:10, sep=''))
       }) %>% 
  bind_rows() %>%
  gather(dt, stat.val, starts_with('dt') ) %>%
  spread(stat.nm, stat.val) %>%
  saveRDS(file='SSM/statsim_norm.rds')
