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

