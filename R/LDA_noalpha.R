### Adapting LDA_original to treat alpha as a univariate hyperparameter

opt_phi_d <- function(phi, gamm, beta){
  t <- digamma(gamm) - digamma(sum(gamm))
  phi <- t(beta * exp(t - 1))

  #normalise rows of phi to sum to 1
  phi[rowSums(phi)==0, ] <- 1/K #to prevent dividing by 0
  phi <- phi / rowSums(phi)
  phi[phi==0] <- 1e-16 #and write over any 0s

  return(phi)
}

opt_gamma_d <- function(gamm, phi, n, alpha, V, K){
  p <- colSums(matrix(n, V, K) * phi)
  gamm <- p + alpha
  return(gamm)
}

gamma_converged <- function(gamm, gamm_old, tol=0.0001){
  diff <- mean(abs(gamm - gamm_old))
  return(diff < tol)
}

e_step_d <- function(gamm, phi, alpha, beta, n, V, K){
  for(iter in 1:20){
    gamm_old <- gamm

    #update phi and gamma for this document
    phi <- opt_phi_d(phi, gamm, beta)
    gamm <- opt_gamma_d(gamm, phi, n, alpha, V, K)

    #check for convergence
    if(gamma_converged(gamm, gamm_old)) break
  }
  return(list("gamm"=gamm, "phi"=phi))
}

update_beta <- function(beta, phis, N, V, K){
  #this step could be optimised, but keeping it long for now to try avoid errors
  for(k in 1:K){
    for(v in 1:V){
      beta[k,v] <- sum(N[,v] * phis[v,k,])
    }
  }

  #normalise rows of beta to sum to 1
  beta <- beta / rowSums(beta)

  #write over any 0s
  beta[beta==0] <- 1e-16

  return(beta)
}


likelihood_corp <- function(phis, gammas, alpha, beta, N, V, K, D){
  t <- digamma(gammas) - digamma(rowSums(gammas))  #would be a capital T
  N_prime <- array(apply(N, 1, function(n) rep(n, K)), c(V, K, D))
  T_prime <- array(apply(t, c(2, 1), function(t) rep(t, V)), c(V, K, D))
  beta_prime <- array(rep(t(beta), D), c(V, K, D))

  L <- D * (lgamma(K*alpha) - K*lgamma(alpha)) -
    sum(lgamma(rowSums(gammas))) +
    sum(lgamma(gammas) + (alpha - gammas)*t) +
    sum(N_prime * phis * (T_prime + log(beta_prime) - log(phis)))

  return(L)
}

initalise_beta <- function(N, V, K, D){
  #take K random samples, and use their word counts as the initial topic distributions
  idx <- sample(1:D, K)
  beta <- N[idx, ]
  #normalise so rows sum to 1
  beta <- beta / rowSums(beta)
  #write over any 0s
  beta[beta==0] <- 1e-16
  return(beta)
}


library(foreach)
library(doParallel)
numCores <- detectCores()
registerDoParallel(numCores)

###### Function to run the LDA
lda_counts <- function(N, K, max_iter=50, thresh=NULL, alpha=NULL, seed=NULL, messages=TRUE){
  #get number of samples and words from N
  V <- ncol(N)
  D <- nrow(N)

  #initialise parameters
  if(is.null(alpha)) alpha <- 1/K

  # if a seed is given use it, but otherwise it will be random each run
  if(!is.null(seed)) set.seed(seed)
  beta <- initalise_beta(N, V, K, D)

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
      e_step_d(gamm=gammas[d,], phi=phis[,,d], alpha, beta, n=N[d,], V, K)
    }
    #combine results
    for(d in 1:D){
      phis[,,d] <- res_lists[[d]]$phi
      gammas[d,] <- res_lists[[d]]$gamm
    }

    # M-step
    if(messages) print(paste("M-Step iteration", iter))
    beta <- update_beta(beta, phis, N, V, K)

    #check for convergence
    if(messages) print("Calculating likelihood")
    loglik[iter] <- likelihood_corp(phis, gammas, alpha, beta, N, V, K, D)
    #by default set the convergence threshold to a % of the last value
    if(iter > 5){
      if(is.null(thresh)) thresh <- abs(1e-5 * loglik[iter-1])
      if(abs(loglik[iter] - loglik[iter-1]) < thresh){
        conv <- T
        break
      }
    }
  }

  #retrieve estimates for thetas (proportional to gammas, or officially theta <- Dir(gamma))
  thetas <- gammas / rowSums(gammas)
  #library(MCMCpack)
  #thetas <- matrix(NA, D, K)
  #for(d in 1:D){
  #  thetas[d, ] <- rdirichlet(n=1, alpha=gammas[d,])
  #}

  return(list("beta"=beta, "thetas"=thetas,
              "phis"=phis, "gammas"=gammas,
              "loglik"=loglik[iter],
              "lls"=loglik[1:iter],
              "alpha"=alpha, "K"=K,
              "iterations"=iter,
              "converged"=conv))
}


