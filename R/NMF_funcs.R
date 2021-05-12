#defining functions to peform non-negative matrix factorisation
#based on ESLII textbook, KL updates formulae

nmf_likelihood <- function(counts, W, H){
  WH <- W %*% H
  return(sum(counts * log(WH) - WH))
}

nmf <- function(counts, K, max_iter=50, thresh=NULL, seed=NULL){
  N <- nrow(counts)
  V <- ncol(counts)

  # initialise matrix factors W and H
  # if a seed is given use it, but otherwise it will be random each run
  if(!is.null(seed)) set.seed(seed)
  W <- matrix(runif(N*K), N, K)
  H <- matrix(runif(K*V), K, V)

  loglik <- rep(NA, max_iter)
  conv <- F

  for(iter in 1:max_iter){
    #update all elements of W
    WH <- W %*% H
    for(i in 1:N){
      for(k in 1:K){
        W[i,k] <- W[i,k] * ((sum(H[k,] * counts[i,] / WH[i,])) / sum(H[k,]))
      }
    }

    #update all elements of H
    WH <- W %*% H
    for(k in 1:K){
      for(j in 1:V){
        H[k,j] <- H[k,j] * ((sum(W[,k] * counts[,j] / WH[,j])) / sum(W[,k]))
      }
    }

    loglik[iter] <- nmf_likelihood(counts, W, H)
    #by default set the convergence threshold to a % of the last value
    if(iter > 5){
      if(is.null(thresh)) thresh <- abs(1e-5 * loglik[iter-1])
      if(abs(loglik[iter] - loglik[iter-1]) < thresh){
        conv <- T
        break
      }
    }
  }

  return(list("W"=W, "H"=H, "ll"=loglik[iter], lls=loglik[1:iter],
              "iterations" = iter, "converged"=conv))
}



