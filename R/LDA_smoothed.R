### Implementing the batch LDA algorithm from Hoffman 2010 paper

opt_phi_d <- function(phi, gamm, lambda, V, K){
  t <- digamma(gamm) - digamma(sum(gamm))
  B <- digamma(lambda) - matrix(digamma(rowSums(lambda)), K, V)

  phi <- t(exp(matrix(t, K, V) + B))
  phi <- phi / rowSums(phi)
  return(phi)
}

opt_gamma_d <- function(gamm, alpha, n, phi, V, K){
  temp <- matrix(n, V, K) * phi
  gamm <- colSums(temp) + alpha
  return(gamm)
}

gamma_converged <- function(gamm, gamm_old, tol=0.0001){
  diff <- mean(abs(gamm - gamm_old))
  return(diff < tol)
}

e_step_d <- function(gamm, phi, alpha, lambda, n, V, K){
  for(iter in 1:30){
    gamm_old <- gamm

    #update phi and gamma for this document
    phi <- opt_phi_d(phi, gamm, lambda, V, K)
    gamm <- opt_gamma_d(gamm, alpha, n, phi, V, K)

    #check for convergence
    if(gamma_converged(gamm, gamm_old)) break
  }
  return(list("gamm"=gamm, "phi"=phi))
}

opt_lambda <- function(lambda, eta, N, phis, V, K, D){
  #make a temp matrix with D rows and V columns (same shape as N) per each k
  for(k in 1:K){
    phis_k <- t(phis[,k,])
    counts_phis <- phis_k * N
    lambda[k,] <- eta + colSums(counts_phis)
  }
  return(lambda)
}

likelihood_corp <- function(phis, gammas, lambda, N, alpha, eta, V, K, D){
  t <- digamma(gammas) - digamma(rowSums(gammas)) #would be a capital T
  B <- digamma(lambda) - matrix(digamma(rowSums(lambda)), K, V)

  Np <- array(apply(N, 1, function(n) rep(n, K)), c(V, K, D))
  Tp <- array(apply(t, c(2, 1), function(t) rep(t, V)), c(V, K, D))
  Bp <- array(rep(t(B), D), c(V, K, D))

  L <- D * (lgamma(K*alpha) - K*lgamma(alpha)) +
    K*(lgamma(V*eta) - V*lgamma(eta)) +
    sum(Np * phis * (Tp + Bp - log(phis))) +
    sum(lgamma(gammas) + (alpha - gammas) * t) +
    sum(lgamma(lambda) + (eta - lambda) * B) -
    sum(lgamma(rowSums(gammas))) -
    sum(lgamma(rowSums(lambda)))

  return(L)
}

initalise_lambda <- function(N, K, nmf){
  if(nmf){
    print("Intialising lambda with NMF")
    lambda <- nmf(N, K)$H + 1
  } else{
    idx <- sample(1:nrow(N), K)
    lambda <- N[idx,] + 1
  }
  return(lambda)
}


library(foreach)
library(doParallel)
numCores <- detectCores()
registerDoParallel(numCores)

###### Function to run the LDA
lda <- function(N, K, max_iter=50, thresh=NULL, alpha=NULL, eta=NULL, NMF=T, messages=TRUE){
  #get number of samples and words from N
  V <- ncol(N)
  D <- nrow(N)

  #initialise parameters
  if(is.null(alpha)) alpha <- 1/K
  if(is.null(eta)) eta <- 1/K

  #initialise lambda using K random documents or with NMF
  lambda <- initalise_lambda(N, K, NMF)

  loglik <- rep(NA, max_iter)
  conv <- F

  for(iter in 1:max_iter){
    # E-step
    if(messages) print(paste("E-Step iteration", iter))

    #initialise phis and gammas
    phis <- array(1/K, c(V, K, D))
    gammas <- matrix(1, D, K)

    #update phi and gamma for each document
    res_lists <- foreach (d=1:D) %dopar% {
      e_step_d(gamm=gammas[d,], phi=phis[,,d], alpha, lambda, n=N[d,], V, K)
    }
    #combine results
    for(d in 1:D){
      phis[,,d] <- res_lists[[d]]$phi
      gammas[d,] <- res_lists[[d]]$gamm
    }

    # M-step
    if(messages) print(paste("M-Step iteration", iter))
    lambda <- opt_lambda(lambda, eta, N, phis, V, K, D)

    # lambda <- foreach(k=1:K, .combine=rbind) %dopar% {
    #   opt_lambda_k(k, eta, N, phis, K)
    # }

    #check for convergence
    if(messages) print("Calculating likelihood")
    loglik[iter] <- likelihood_corp(phis, gammas, lambda, N, alpha, eta, V, K, D)
    #by default set the convergence threshold to a % of the last value
    if(iter > 5){
      if(is.null(thresh)) thresh <- abs(1e-5 * loglik[iter-1])
      if(abs(loglik[iter] - loglik[iter-1]) < thresh){
        conv <- T
        break
      }
    }
  }

  #retrieve estimates for thetas and beta
  thetas <- gammas / rowSums(gammas)
  beta <- lambda / rowSums(lambda)

  return(list("beta"=beta, "thetas"=thetas,
              "phis"=phis, "gammas"=gammas,
              "lambda"=lambda,
              "L"=loglik[iter],
              "Ls"=loglik[1:iter],
              "alpha"=alpha, "K"=K,
              "iterations"=iter,
              "converged"=conv))
}


