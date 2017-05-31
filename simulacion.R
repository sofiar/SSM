library(dplyr)
library(mvtnorm)
library(ggplot2)
library(MCMCpack)
library(reshape2)


### Simulacion para cuatro estados 50 individuos 500 pasos de tiempo

# probabilidad inicial de los estados 
del <- c(0.2, 0.1, 0.6, 0.1)
# matriz de probabilidad de transición
tpm <- matrix(c(0.9, 0.05, 0.03, 0.02, 0.1, 0.8, 0.07, 0.03, 0.01, 0.15, 0.75,
                  0.09, 0.04, 0.01, 0.11, 0.86), byrow = T, nrow = 4)

#secuencia de estados para cada individuos en cada tiempo
states <- matrix(NA, nrow = 500, ncol = 50)
for (i in 1:50) {
  for (j in 1:500) {
    if (j == 1) {
      states[1, i] <- sample(x = 1:4, size = 1, prob = del)
    } else {
      states[j, i] <- sample(x = 1:4, size = 1, prob = tpm[states[j -1, i], ])
    }
  }
}


# medias y matrices de varianzas para cada estado. Cada estado se considera distribuido a partir de una
# normal multivariada (n=3)
mus <- list()
vars <- list()
# state 1
mus[[1]] <- c(0, -7, 10)
vars[[1]] <- riwish(3, matrix(c(3, 0.3, 2, 0.3, 1, 0.5, 2, 0.1, 2), 3, 3))
# state 2
mus[[2]] <- c(-5, 1, 1)
vars[[2]] <- riwish(3, matrix(c(3, 0.3, 2, 0.3, 1, 0.5, 0.4, 0.1, 3), 3, 3))
# state 3
mus[[3]] <- c(-2, -2, 3)
vars[[3]] <- riwish(3, matrix(c(3, 0.1, 2, 0.3, 1, 0.5, 2, 0.1, 2), 3, 3))
# state 4
mus[[4]] <- c(5, 5, 1)
vars[[4]] <- riwish(3, matrix(c(1, 0.3, 0.2, 0.3, 2, 0.1, 0.5, 0.1, 2), 3, 3))

## Simulacion de las observaciones

obs <- list()
for (k in 1:50) obs[[k]] <- matrix(NA, nrow = 500, ncol = 3)
for (i in 1:50) {
  for (j in 1:500) {
    obs[[i]][j, ] <- rmvnorm(n = 1, mean = mus[[states[j, i]]], sigma = vars[[states[j,
                                                                                     i]]])
  }
}

## formato y estructura de los datos
tracks <- paste("T", 1:50, sep = "")
data <- data.frame(obs[[1]], state = states[, 1], timeseries = rep(tracks[1],
                                                                 500))
for (j in 2:50) data <- rbind(data, data.frame(obs[[j]], state = states[, j],
                                               timeseries = rep(tracks[j], 500)))

### plot de la serie de estados de indiv 1
data.1=data[data$timeseries=='T1',]

plot(which(data.1$state==1),rep ('1',length(which(data.1$state==1))),xlab = 'Indice',ylab='Estados',ylim = c(1,4),xlim=c(1,500),col='lightblue',lab=c(16,4,5),pch=16)
points(which(data.1$state==2),rep ('2',length(which(data.1$state==2))),pch=16,col='darkcyan')
points(which(data.1$state==3),rep ('3',length(which(data.1$state==3))),pch=16,col='coral')
points(which(data.1$state==4),rep ('4',length(which(data.1$state==4))),pch=16,col='gold2')








