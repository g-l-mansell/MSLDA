#Making algorithm 1 from Hoffman 2010 into functions so I can source the script
#https://www.di.ens.fr/~fbach/mdhnips2010.pdf

#we need input parameters:
#a matrix of counts, rows are samples, columns are words
#number of topics to try and find

##### Useful functions
library(MCMCpack)

E_logtheta <- function(gamm, k){
  return(digamma(gamm[k]) - digamma(sum(gamm)))
}

E_logbeta <- function(lambda, k, w){
  return(digamma(lambda[k, w]) - digamma(sum(lambda[k,])))
}

opt_phi_d <- function(phi, gamm, lambda, W, K){
  a <- matrix(digamma(gamm), W, K, byrow=T)
  b <- digamma(sum(gamm))
  c <- t(digamma(lambda))
  d <- matrix(digamma(rowSums(lambda)), W, K, byrow=T)

  phi <- exp(a - b + c - d)
  phi <- phi / rowSums(phi)
  return(phi)
}

opt_gamma_d <- function(gamm, alpha, count, phi){
  temp <- matrix(count, nrow(phi), ncol(phi)) * phi
  gamm <- colSums(temp) + alpha
  return(gamm)
}

gamma_converged <- function(gamm, gamm_old, tol=0.0001){
  diff <- mean(abs(gamm - gamm_old))
  return(diff < tol)
}

e_step_d <- function(gamm, phi, alpha, lambda, count, W, K){
  for(iter in 1:20){
    gamm_old <- gamm

    #update phi and gamma for this document
    phi <- opt_phi_d(phi, gamm, lambda, W, K)
    gamm <- opt_gamma_d(gamm, alpha, count, phi)

    #check for convergence
    if(gamma_converged(gamm, gamm_old)) break
  }
  return(list("gamm"=gamm, "phi"=phi))
}


opt_lambda <- function(lambda, eta, counts, phis, W, K, D){
  #make a temp matrix with D rows and W columns (same shape as counts) per each k
  for(k in 1:K){
    phis_k <- t(sapply(phis, function(phi) phi[,k]))
    counts_phis <- phis_k * counts
    lambda[k,] <- eta + colSums(counts_phis)
  }
  return(lambda)
}

likelihood_corp <- function(phis, gammas, lambda, counts, alpha, eta, W, K, D){
  #line 1 of likelihood (eq 4)
  #for each k make a D * W matrix
  #then sum over ks, element multiply by counts, then sum
  mat <- matrix(0, D, W)
  for(k in 1:K){
    phis_k <- t(sapply(phis, function(phi) phi[,k]))
    mat_a <- matrix(digamma(gammas[,k]), D, W)
    mat_b <- matrix(digamma(rowSums(gammas)), D, W)
    mat_c <- matrix(digamma(lambda[k,]), D, W, byrow=T)
    mat_d <- digamma(sum(lambda[k,]))
    mat_e <- -log(phis_k)

    mat_k <- phis_k * (mat_a - mat_b + mat_c - mat_d - mat_e)
    mat <- mat + mat_k
  }
  L <- sum(counts * mat)

  #line 2, first term
  L <- L - sum(lgamma(rowSums(gammas)))

  #line 2, sum over k term
  #make a D*K matrix for each term and sum
  mat_a <- alpha - gammas
  mat_b <- digamma(gammas)
  mat_c <- matrix(digamma(rowSums(gammas)), D, K)
  mat_d <- lgamma(gammas)
  L <- L + sum(mat_a * (mat_b - mat_c) - mat_d)

  #line 3, first term
  # (sum over k applies to every term) (skip multiplying then dividing by D)
  L <- L - sum(lgamma(rowSums(lambda)))

  #line 3, sum over w term
  #making a K*W matrix and summing
  mat_a <- eta - lambda
  mat_b <- digamma(lambda)
  mat_c <- matrix(digamma(rowSums(lambda)), K, W)
  mat_d <- lgamma(lambda)
  L <- L + sum(mat_a * (mat_b - mat_c) + mat_d)

  #line 4
  L <- L + D * (lgamma(K*alpha) - K*lgamma(alpha))
  L <- L + lgamma(W*eta) - W*lgamma(eta)

  return(L)
}

initalise_lambda <- function(counts, K, nmf){
  if(nmf){
    print("Intialising lambda with NMF")
    lambda <- nmf(counts, K)$H + 1
  } else{
    idx <- sample(1:nrow(counts), K)
    lambda <- counts[idx,] + 1
  }
  return(lambda)
}


library(foreach)
library(doParallel)
numCores <- detectCores()
registerDoParallel(numCores)

###### Function to run the LDA
lda <- function(counts, K, max_iter=50, thresh=NULL, alpha=NULL, eta=NULL, nmf=T, messages=TRUE){
  #get number of samples and words from counts
  W <- ncol(counts)
  D <- nrow(counts)

  #initalise parameters
  if(is.null(alpha)) alpha <- 1/K
  if(is.null(eta)) eta <- 1/K

  #initalise lambda using K random documents or with NMF
  lambda <- initalise_lambda(counts, K, nmf)

  loglik <- rep(NA, max_iter)
  conv <- F

  for(iter in 1:max_iter){
    # E-step
    if(messages) print(paste("E-Step iteration", iter))

    #initialise phis and gammas
    phis <- vector("list", D)
    gammas <- matrix(0, D, K)
    for(d in 1:D){
      gammas[d,] <- 1
      phis[[d]] <- matrix(1/K, W, K)
    }

    #update phi and gamma for each document
    res_lists <- foreach (d=1:D) %dopar% {
      e_step_d(gamm=gammas[d,], phi=phis[[d]], alpha, lambda, count=counts[d,], W, K)
    }
    #combine results
    for(d in 1:D){
      phis[[d]] <- res_lists[[d]]$phi
      gammas[d,] <- res_lists[[d]]$gamm
    }

    # M-step
    if(messages) print(paste("M-Step iteration", iter))
    lambda <- opt_lambda(lambda, eta, counts, phis, W, K, D)

    # lambda <- foreach(k=1:K, .combine=rbind) %dopar% {
    #   opt_lambda_k(k, eta, counts, phis, K)
    # }

    #check for convergence
    if(messages) print("Calculating likelihood")
    loglik[iter] <- likelihood_corp(phis, gammas, lambda, counts, alpha, eta, W, K, D)
    #by default set the convergence threshold to a % of the last value
    if(iter > 5){
      if(is.null(thresh)) thresh <- abs(1e-5 * loglik[iter-1])
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


