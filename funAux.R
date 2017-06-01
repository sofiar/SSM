################ Funciones auxiliares ##############################


## Algoritmo recursivo de viterbi para calcular la secuencia de estados mas probable

# x = vector valores de las tres variables para cada tiempo del individuo en cuestion
# m = cantidad total de posibles estados (en este caso son cuatro)
# gamma = matriz con valores de las medias (estimados)
# allprobs = matriz con probabilidad de cada estado en cada tiempo
# delta = probas iniciales para cada estado

HMM.viterbi <- function(x, m, gamma, allprobs, delta=NULL, ...)
{
  if(is.null(delta)) delta <- solve(t(diag(m) - gamma + 1), rep(1, m))
  n <- dim(x)[1]
  xi <- matrix(0, n, m)
  foo <- delta*allprobs[1,]
  xi[1,] <- foo/sum(foo)
  for(i in 2:n)
  {
    foo <- apply(xi[i-1,]*gamma, 2, max)*allprobs[i,]
    xi[i,] <- foo/sum(foo)
  }
  iv <- numeric(n)
  iv[n] <- which.max(xi[n,])
  for(i in (n-1):1)
    iv[i] <- which.max(gamma[, iv[i+1]]*xi[i,])
  return(iv)
}

########### Calculating the negative log-likelihood
## function that converts 'natural' parameters (possibly constrained)
## to 'working' parameters (all of which are real-valued) - necessary
## to use the unconstrained optimizer nlm()
pn2pw <- function(mu,sigma,beta,delta,pi,N){
  tmu <- log(mu) # N means for the
  tsigma <- log(sigma)
  tbeta<-as.vector(beta)
  teta<-log(delta[-1]/delta[1])
  tpi<-qlogis(pi)
  parvect <- c(tmu,tsigma,tbeta,teta,tpi)
  return(parvect)
}
## function that performs the inverse transformation
pw2pn <- function(parvect,N){
  mu <- exp(parvect[1:(N+1)])
  sigma <- exp(parvect[(N+2):(2*N+2)])
  beta <- matrix(parvect[(2*N+3):(2*N^2+2)],nrow=N*(N-1))
  delta <- exp(c(0,parvect[(2*N^2+3):(2*N^2+N+1)]))
  delta <- delta/sum(delta)
  pi <- plogis(parvect[2*N^2+N+2])
  return(list(mu=mu,sigma=sigma,beta=beta,delta=delta,pi=pi))
}

## function that computes the negative log-likelihood of the HMM
mllk <- function(parvect,obs,N,x1){ #parvec: vector de los parametros
  lpn <- pw2pn(parvect,N)
  K <- length(obs)
  mllk.all <- 0
  for(k in 1:K){
    n <- dim(obs[[k]])[1]
    allprobs <- matrix(rep(1, N*n), nrow=n)
    ind.step <- which(!is.na(obs[[k]][,"msa"]))## lugares donde la variable exite
    step.prob <- rep(1,n)
    allprobs[ind.step,1] <- lpn$pi*dgamma(obs[[k]][ind.step,"msa"],
                                          shape=lpn$mu[1]^2/lpn$sigma[1]^2,
                                          scale=lpn$sigma[1]^2/lpn$mu[1])+
      (1-lpn$pi)*dgamma(obs[[k]][ind.step, "msa"],
                        shape=lpn$mu[2]^2/lpn$sigma[2]^2,scale=lpn$sigma[2]^2/lpn$mu[2]) #densidad para estado1
    allprobs[ind.step,2] <- dgamma(obs[[k]][ind.step,"msa"],
                                   shape=lpn$mu[3]^2/lpn$sigma[3]^2,scale=lpn$sigma[3]^2/lpn$mu[3])
    lscale <- 0
    # implementing the forward algorithm
    foo <- lpn$delta
    foo <- foo*allprobs[1,]
    sumfoo <- sum(foo); lscale <- lscale + log(sumfoo); foo <- foo/sumfoo #scaling
    for (i in 2:n){
      gamma <- diag(N)
      gamma[!gamma] <- exp(lpn$beta[,1]+lpn$beta[,2]*obs[[k]][i,x1])
      gamma <- gamma/apply(gamma,1,sum)
      foo <- foo%*%gamma*allprobs[i,]
      sumfoo <- sum(foo); lscale <- lscale+log(sumfoo); foo <- foo/sumfoo # scaling
    }
    mllk.all <- mllk.all + lscale
  }

return(-1*mllk.all)
}
                                   

## function that runs the numerical minimization of 'mllk' (i.e. tries to find the MLE)
mle <- function(obs,mu0,sigma0,beta0,delta0,pi0,N,x1){
  parvect <- pn2pw(mu0,sigma0,beta0,delta0,pi0,N)
  mod <- nlm(mllk,parvect,obs,N,x1,print.level=2,iterlim=1000,stepmax=5)#minimiza mllk por algoritmo d Newton
  pn <- pw2pn(mod$estimate,N)
  return(list(mu=pn$mu,sigma=pn$sigma,beta=pn$beta,
              delta=pn$delta,pi=pn$pi,mllk=mod$minimum))
}

