### Implementing non-negative matrix factorisation from ESLII textbook, KL updates

nmf_likelihood <- function(counts, W, H){
  WH <- W %*% H
  return(sum(counts * log(WH) - WH))
}

#' Run NMF using a count matrix
#'
#' @param counts matrix of word counts
#' @param K internal dimension of matrix factors
#' @param max_iter the maximum number of iterations to run
#' @param thresh threshold for L convergence, (L_i - L_{i-1})/L_i < thresh
#' @param seed for the random initalisation of factors W and H
#' @return A list of all parameters
#' @export
nmf <- function(counts, K, max_iter=50, thresh=1e-4, seed=NULL){
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

    #check for convergence
    loglik[iter] <- nmf_likelihood(counts, W, H)
    if(L_converged(loglik, iter, thresh)){
      conv <- T
      break
    }
  }

  return(list("W"=W, "H"=H, "L"=loglik[iter], "Ls"=loglik[1:iter],
              "iterations" = iter, "converged"=conv))
}



