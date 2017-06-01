##Simulacion 2

########## Simulating data from 2-state HMM with covariate included for wind speed (9 days)
#2 states
N <- 2
#simulating values of wind
wind.mat <- matrix(NA, nrow=9, ncol=200)
#wind matrix
#valores iniciales
for(j in 1:9){
  wind.mat[j,1] <- runif(1, min=0, max=7)
}
#completa valores dependen del anterior
for(j in 1:9){
  for(i in 2:200){
    wind.mat[j,i] <- runif(n=1, min = max(0, wind.mat[j,i-1]-2),
                           max = min(7, wind.mat[j,i-1]+2))
  }
}
#coefficients
beta.mat <- matrix(c(-2.5, -.1, -3, .1), nrow=2, byrow=T)
#initial distribution
delta <- c(.5, .5)
#states
# gamma es la matriz de transición (depende de los coeficientes beta y del viento)
# gamma define entonces los estados se samplea con las probas de gamma
state.mat <- matrix(NA, nrow=9, ncol=200)
for(j in 1:9){
  state.mat[j,1] <- sample(x = c(1,2), size = 1, prob = delta)
  for(i in 2:200){
    gamma <- diag(N)
    gamma[!gamma] <- exp(beta.mat[,1] + beta.mat[,2]*wind.mat[j,i]) ## nu del manuscrito
    gamma <- gamma/apply(gamma,1,sum)#psi del manuscrito
    state.mat[j,i] <- sample(x = c(1,2), size=1, prob=gamma[state.mat[j,i-1],])
  }
}

#means and variances of the state-depenent gamma distributions
#X en el estado uno se distribuye como una mezcla de dos gamas (en propocion mp1). En el estado dos se 
# distribuye como una unica gamma

mu.2st <- list()
mu.2st[[1]] <- c(.1, .3)
mu.2st[[2]] <- .7
var.2st <- list()
var.2st[[1]] <- c(.2, .3)
var.2st[[2]] <- 1
#mixing proportion for mixture of gammas in state 1
mp1 <- .5
#simulating observations
msa.mat <- matrix(NA, nrow=9, ncol=200)
for(j in 1:9){
  for(i in 1:200){
    if(state.mat[j,i]==1){
      msa.mat[j,i] <- sample(x=c(rgamma(n = 1,
                                        shape=mu.2st[[1]][1]^2/var.2st[[1]][1]^2,
                                        scale=var.2st[[1]][1]^2/mu.2st[[1]][1]),
                                 rgamma(n = 1, shape=mu.2st[[1]][2]^2/var.2st[[1]][2]^2,
                                        scale=var.2st[[1]][2]^2/mu.2st[[1]][2])),
                             size=1, prob = c(1-mp1, mp1))
    } else{
      msa.mat[j,i] <- rgamma(n = 1, shape=mu.2st[[2]][1]^2/var.2st[[2]][1]^2,
                             scale=var.2st[[2]][1]^2/mu.2st[[2]][1])
    }
  }
}


# formatting the data
trcks <- paste("T", 1:9, sep="")
msafull <- list()
for(k in 1:9) msafull[[k]] <- data.frame(msa=msa.mat[k,], wind=wind.mat[k,], track=trcks[k])