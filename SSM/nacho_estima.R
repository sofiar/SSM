
# estimamos con abc
library(circular)
library(Rcpp)
library(tidyverse)
library(mvtnorm)
library(abc)
library(randomForest)

# ==============================================================
# FUNCIONES
# ==============================================================
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
# ==============================================================

# un solo conjunto de datos 'reales', con 3 dt = .05, .5, 5 
# tener cuidado que los nsteps, y dt.list coincidan con los que se 
# usaron para obtener las simulaciones

t_w <- 2
t_k <- 20

data.obs <- simdata(ws=t_w, ks=t_k, dt.list = list(.05, .5, 5), nsteps=800)

stat.obs <- lapply(data.obs, stat_fun) %>%
  bind_rows() %>%
  mutate( stat.nm = paste('T', 1:10, sep='') ) %>%
  gather(dt, stat.val, starts_with('dt') ) %>%
  spread(stat.nm, stat.val)

# funcion que estima dos ajuste por abc: 
#   usando los summary stat
#   usando un modelo RF para obterner una transformacion de los summary stat
abc_fn <- function(d) {
  # primero hago los modelos RF
  rf.w <- randomForest( s_w ~ . , data = dplyr::select(d, -s_k, -dt) )
  rf.k <- randomForest( s_k ~ . , data = dplyr::select(d, -s_w,-dt) )
  rf.stsim <- data_frame(s_w = predict(rf.w), s_k = predict(rf.k) )
  
  stobs = stat.obs %>% filter(dt == d$dt[1])
  rf.stat.obs <- c( predict(rf.w, newdata = stobs  ), predict(rf.k, newdata = stobs ) ) 
  
  # abc usando los summary stat
  out <- abc(target = as.numeric(stobs[-1]), 
             param = d[, 1:2], 
             sumstat = d[, 4:13] , tol = .5, 
             method='neuralnet')
  
  # abc usando los RF de summary stat
  out.rf <- abc(target = rf.stat.obs, 
                param = d[, 1:2], 
                sumstat = rf.stsim , tol = .5, 
                method='neuralnet')
  list(out, out.rf)
}


# Usamos dos previas, para cada previa podemos usar, los summary stat o un rf de los summary stat
statsim.unif <- readRDS('SSM/statsim_unif.rds')
statsim.norm <- readRDS('SSM/statsim_norm.rds')

out.unif <- statsim.unif %>% 
  split.data.frame( factor(statsim.unif$dt)  ) %>%
  lapply( abc_fn )

out.norm <- statsim.norm %>% 
  split.data.frame( factor(statsim.norm$dt)  ) %>%
  lapply( abc_fn )

saveRDS(out.unif, file='SSM/out_unif.rds')
saveRDS(out.norm, file='SSM/out_norm.rds')


  

