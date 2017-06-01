########################### Ajuste de HMM ##################################

## Implementacion de LOOCV. Ajusta HMM con 49 individuos y testea con el otro

#Estima las probas iniciales como proporcion de la cantidad de veces que el estado k es inicial
#para k=1,2,3,4
delta_est <- list()
for (i in 1:50) {
  #### ACA NO ES AL REVEZ???
  
  #asi estaba
  #delta_est[[i]] <- table(states[-i, 1])/49  
  #asi lo deje
  delta_est[[i]] <- table(states[1, -i])/49
  }

# Estimacion de las medias para las normales multivariadas
mu_est <- list()
for (i in 1:50) {
  
  #se queda con los valores de data que NO correspondan al indv i y los agrupa por estado. Pide
  #el promedio de cada una de las variables X1, X2 y X3 
  mu_est[[i]] <- data %>% filter(timeseries != tracks[i]) %>% group_by(state) %>%
    summarize(x1m = mean(X1), x2m = mean(X2), x3m = mean(X3))
}
# Estimacion de las matrices de covarianzas
cov_est <- list()
for (j in 1:50) {
  temp1 <- filter(data, timeseries != tracks[i] & state == 1)
  c1 <- cov(temp1[, c("X1", "X2", "X3")])
  temp2 <- filter(data, timeseries != tracks[i] & state == 2)
  c2 <- cov(temp2[, c("X1", "X2", "X3")])
  temp3 <- filter(data, timeseries != tracks[i] & state == 3)
  c3 <- cov(temp3[, c("X1", "X2", "X3")])
  temp4 <- filter(data, timeseries != tracks[i] & state == 4)
  c4 <- cov(temp4[, c("X1", "X2", "X3")])
  cov_est[[j]] <- data.frame(state = rep(1:4, each = 3), rbind(c1, c2, c3,
                                                               c4))
}

# estimacion de la matriz de probabilidad de transicion

tpm_est <- list()
for (j in 1:50) {
  tpm_est[[j]] <- diag(4)
  n <- length(states[, j])
  smat <- matrix(0, nrow = 4, ncol = 4)
  for (k in c(1:50)[-j]) {
    w1 <- which(states[1:(n - 1), k] == 1) #quienes tienen estado 1
    s1 <- states[w1, k] - states[w1 + 1, k] # resta estos estados menos el consecutivo
    w2 <- which(states[1:(n - 1), k] == 2)
    s2 <- states[w2, k] - states[w2 + 1, k]
    w3 <- which(states[1:(n - 1), k] == 3)
    s3 <- states[w3, k] - states[w3 + 1, k]
    w4 <- which(states[1:(n - 1), k] == 4)
    s4 <- states[w4, k] - states[w4 + 1, k]
    
    for (m in 1:4)
    {
      #cuenta la cantidad de veces que se pasa del estado m a uno l (l=1,2,3,4)
      smat[1,m]<-smat[1,m]+length(which(s1==(1-m))) 
      smat[2,m]<-smat[2,m]+length(which(s2==(2-m)))
      smat[3,m]<-smat[3,m]+length(which(s3==(3-m)))
      smat[4,m]<-smat[4,m]+length(which(s4==(4-m)))
    }
    
  }
  ##promedia
  wsm<-rowSums(smat)
  
  tpm_est[[j]][1,]<-smat[1,]/wsm[1] 
  tpm_est[[j]][2,]<-smat[2,]/wsm[2]
  tpm_est[[j]][3,]<-smat[3,]/wsm[3]
  tpm_est[[j]][4,]<-smat[4,]/wsm[4]
}
    

########## Estimacion del error de clasificacion

# para cada indviduo i se predice la secuencia de estados mas probables a partir del valor de los 
# parametros ajustados por todos los demas individuos. Para encontrar la secuencia de estados mas probable
# utiliza el algoritmo recursivo de Viterbi

source("funAux.R")
data$predstate <- NA
for(u in 1:50){
  ## Datos del individuo u
  tempdat <- filter(data, timeseries==tracks[u]) 
  ## Matriz para guardar las probabilidades de cada uno de los estados en cada tiempo
  allprobs <- matrix(NA, nrow=500, 4)
  for(j in 1:500) { 
    for(k in 1:4) {
      #proba de estado k a tiempo j
      allprobs[j,k] <- dmvnorm(x=tempdat[j,c("X1", "X2", "X3")],mean = unlist(mu_est[[u]][k,2:4]),
                               sigma = as.matrix(filter(cov_est[[u]],state==k)[,c("X1", "X2", "X3")]))
  }
      
  }
  
          
  
  ## calcula la secuencia de estados mas probables por viterbi
  data$predstate[which(data$timeseries==tracks[u])] <- HMM.viterbi(x=tempdat, m = 4,
                                                                   gamma = tpm_est[[u]],
                                                                   allprobs = allprobs,
                                                                   delta = delta_est[[u]])
}


plot(data[which(data$timeseries==tracks[3]),]$predstate,ylab='',lab=c(4,4,5),ylim=c(1,4))
points(data[which(data$timeseries==tracks[3]),]$state,col='red')


#prediction error
prederror <- length(which(data$state-data$predstate!=0))/dim(data)[1]
prederror

      
#confusion matrix
mat <- acast(data, state~predstate)
mat



      
      
      
  
                               



