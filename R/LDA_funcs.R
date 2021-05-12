#Making algorithm 1 from Hoffman 2010 into functions so I can source the script
#https://www.di.ens.fr/~fbach/mdhnips2010.pdf

#we need input parameters:
#a matrix of counts, rows are samples, columns are words
#number of topics to try and find

##### Useful functions
library(MCMCpack)

E_logtheta <- function(gammas, d, k){
  return(digamma(gammas[d, k]) - digamma(sum(gammas[d,])))
}

E_logbeta <- function(lambda, k, w){
  return(digamma(lambda[k, w]) - digamma(sum(lambda[k,])))
}

opt_phi_d <- function(d, phis, gammas, lambda, W, K){
  for(w in 1:W){
    for(k in 1:K){
      phis[[d]][w, k] <- exp(E_logtheta(gammas, d, k) + E_logbeta(lambda, k, w))
    }
  }
  phis[[d]] <- phis[[d]] / rowSums(phis[[d]])
  return(phis)
}

opt_gamma_d <- function(d, gammas, alpha, counts, phis, W, K){
  for(k in 1:K){
    gammas[d, k] <- alpha
    for(w in 1:W){
      gammas[d, k] <-  gammas[d, k] + (counts[d, w] * phis[[d]][w, k])
    }
  }
  return(gammas)
}

opt_lambda <- function(lambda, eta, counts, phis, W, K, D){
  for(k in 1:K){
    for(w in 1:W){
      lambda[k, w] <- eta
      for(d in 1:D){
        lambda[k, w] <- lambda[k, w] + (counts[d, w] * phis[[d]][w, k])
      }
    }
  }
  return(lambda)
}

gamma_converged <- function(d, gammas, gammas_old, tol=0.0001){
  diff <- mean(abs(gammas[d,] - gammas_old[d,]))
  return(diff < tol)
}

likelihood_doc <- function(d, phis, gammas, lambda, counts, alpha, eta, W, K, D){
  L <- 0
  for(w in 1:W){
    for(k in 1:K){
      L <- L + counts[d, w] * phis[[d]][w, k] * (E_logtheta(gammas, d, k) + E_logbeta(lambda, k, w) - log(phis[[d]][w, k]))
    }
  }
  L <- L - lgamma(sum(gammas[d,]))

  for(k in 1:K){
    L <- L + (alpha - gammas[d, k]) * E_logtheta(gammas, d, k) +
      lgamma(gammas[d, k])

    L <- L - lgamma(sum(lambda[k,]))/D

    for(w in 1:W){
      L <- L + ((eta - lambda[k, w]) * E_logbeta(lambda, k, w) + lgamma(lambda[k, w]))/D
    }
  }

  L <- L + lgamma(K*alpha) - K * lgamma(alpha) + (lgamma(W*eta) - W * lgamma(eta))/D

  return(L)
}

likelihood_corp <- function(phis, gammas, lambda, counts, alpha, eta, W, K, D){
  L <- 0
  for(d in 1:D){
    L <- L + likelihood_doc(d, phis, gammas, lambda, counts, alpha, eta, W, K, D)
  }
  return(L)
}



###### Function to run the LDA
lda <- function(counts, K, max_iter=50, thresh=NULL, messages=TRUE){
  #get number of samples and words from counts
  W <- ncol(counts)
  D <- nrow(counts)

  #initalise parameters
  alpha <- 1/K #still not sure about why this isn't a vector
  eta <- 1/K
  lambda <- matrix(sample(1:2, K*W, T), K, W)

  loglik <- rep(NA, max_iter)
  conv <- F

  for(iter in 1:max_iter){
    # E-step
    if(messages) print(paste("E-Step iteration", iter))

    #reinitialise phis and gammas
    phis <- vector("list", D)
    gammas <- matrix(NA, D, K)
    for(d in 1:D){
      N <- sum(docs[[d]])
      gammas[d,] <- alpha + N/K
      phis[[d]] <- matrix(1/K, W, K)
    }

    #update phi and gamma for each document
    for(d in 1:D){
      for(iter2 in 1:20){
        gammas_old <- gammas

        #update phi and gamma for this document
        phis <- opt_phi_d(d, phis, gammas, lambda, W, K)
        gammas <- opt_gamma_d(d, gammas, alpha, counts, phis, W, K)

        #check for convergence
        if(gamma_converged(d, gammas, gammas_old)) break
      }
    }

    # M-step
    if(messages) print(paste("M-Step iteration", iter))
    lambda <- opt_lambda(lambda, eta, counts, phis, W, K, D)

    #check for convergence
    loglik[iter] <- likelihood_corp(phis, gammas, lambda, counts, alpha, eta, W, K, D)
    #by default set the convergence threshold to 0.1% of the first likelihood
    if(iter == 1){
      if(is.null(thresh)) thresh=abs(0.001*loglik[iter])
    } else {
      if(abs(loglik[iter] - loglik[iter-1]) < thresh){
        conv <- T
        break
      }
    }
  }

  #retrieve estimates for theta and beta
  beta <- matrix(NA, K, W)
  for(k in 1:K){
    beta[k, ] <- rdirichlet(n=1, alpha=lambda[k,])
  }

  thetas <- matrix(NA, D, K)
  for(d in 1:D){
    thetas[d, ] <- rdirichlet(n=1, alpha=gammas[d,])
  }

  return(list("beta"=beta, "thetas"=thetas,
              "phis"=phis, "gammas"=gammas,
              "lambda"=lambda, "loglik"=loglik[iter],
              "lls"=loglik[1:iter],
              "iterations"=iter,
              "converged"=conv))
}
