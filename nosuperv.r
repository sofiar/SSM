source('funAux.R')

N <- 2

mu0 <- c(0.1, 0.3, 0.7)
sigma0 <- c(0.2, 0.3, 1.0)
beta0 <- matrix(c(-2.5, -.1, -3, .1), nrow=2, byrow=T)
delta0 <- c(.5, .5)
pi0 <- .5
x1 <- "wind"

##inicializa de los valores reales?
mod <- mle(obs=msafull, mu0, sigma0, beta0, delta0, pi0, N, x1)

names(mod)
mod$delta


#### PLOTS
## VALORES MSA SIMULADOS PARA LOS 9 DIAS
ggplot(data = data.frame(rbind_all(msafull), x = rep(1:200, 9)), aes(x = x,msa)) + geom_line() + facet_grid(track ~ .) + ylab("Simulated MSA") + xlab("")

### VALORES DE MSA SIMULADOS PARA EL DIA 2 Y LOS DIFERENTES ESTADOS 
temp.dat <- data.frame(filter(data.frame(rbind_all(msafull), x = rep(1:200,
                                                                     9)), track == "T2"), state = state.mat[2, ])
tdm <- melt(temp.dat[, c("msa", "x", "state")], id.vars = "x")
ggplot(data = tdm, aes(x = x, value)) + geom_line() + facet_grid(variable ~
                                                                   ., scales = "free_y") + xlab("")

                                        
## Code for the marginal and state-dependent densities
xr <- seq(0.0001, 7, length.out=100)

# cheating a bit here, true proportion of observations
# in each state
dels <- table(state.mat)/1800
d1 <- dels[1]*(mod$pi*dgamma(xr,shape=mod$mu[1]^2/mod$sigma[1]^2,scale=mod$sigma[1]^2/mod$mu[1])+
                 (1-mod$pi)*dgamma(xr,shape=mod$mu[2]^2/mod$sigma[2]^2,scale=mod$sigma[2]^2/mod$mu[2]))
d2 <- dels[2]*dgamma(xr,shape=mod$mu[3]^2/mod$sigma[3]^2,scale=mod$sigma[3]^2/mod$mu[3])
dfres <- data.frame(xr, d1, d2, d1+d2)
colnames(dfres) <- c("MSA", "State 1", "State 2", "Marginal")
dfres.m <- melt(dfres, id.vars="MSA")
ggplot(data=rbind_all(msafull), aes(x=msa)) +
  geom_histogram(aes(y=..density..), binwidth = .1, alpha=.5, origin=0) +
  geom_line(data=dfres.m, aes(MSA, value, color=variable), size=1) +
  ylim(0, 1.5) + xlim(0, 4) + xlab("MSA") + ylab("Density") +
  scale_color_manual(name="Densities",
                     values=c("State 1"="blue", "State 2"="purple", "Marginal"="black")) +
  ggtitle("State-dependent and Marginal Densities for MSA")




