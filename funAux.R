################ Funciones auxiliares ##############################

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

