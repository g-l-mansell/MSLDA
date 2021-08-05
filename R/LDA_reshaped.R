### Adapting LDA_original to use a count matrix

update_phi_d2 <- function(phi, gamm, beta, K){
  t <- digamma(gamm) - digamma(sum(gamm))
  phi <- t(beta * exp(t - 1))

  #normalise rows of phi to sum to 1
  phi[rowSums(phi)==0, ] <- 1/K #to prevent dividing by 0
  phi <- phi / rowSums(phi)
  phi[phi==0] <- 1e-16 #and write over any 0s

  return(phi)
}

update_gamma_d2 <- function(gamm, phi, n, alpha, V, K){
  p <- colSums(matrix(n, V, K) * phi)
  gamm <- p + alpha
  return(gamm)
}

e_step_d2 <- function(gamm, phi, alpha, beta, n, V, K){
  gamm <- rep(1/K, K)
  for(iter in 1:20){
    gamm_old <- gamm

    #update phi and gamma for this document
    phi <- update_phi_d2(phi, gamm, beta, K)
    gamm <- update_gamma_d2(gamm, phi, n, alpha, V, K)

    #check for convergence
    if(gamma_converged(gamm, gamm_old)) break
  }
  return(list("gamm"=gamm, "phi"=phi))
}

update_beta2 <- function(beta, phis, N, V, K){
  #this step could be optimised, but keeping it long for now to try avoid errors
  for(k in 1:K){
    for(v in 1:V){
      beta[k,v] <- sum(N[,v] * phis[v,k,])
    }
  }

  #normalise rows of beta to sum to 1
  beta <- beta / rowSums(beta)
  beta[beta==0] <- 1e-16
  return(beta)
}


loglik_corp2 <- function(phis, gammas, alpha, beta, N, V, K, D){
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

initalise_beta2 <- function(N, V, K, D){
  #take K random samples, and use their word counts as the initial topic distributions
  idx <- sample(1:D, K)
  beta <- N[idx, ]
  #normalise so rows sum to 1
  beta <- beta / rowSums(beta)
  #write over any 0s
  beta[beta==0] <- 1e-16
  return(beta)
}


#' Run LDA adapted to use a count matrix
#'
#' @inheritParams lda_noalpha
#' @param N matrix of word counts
#' @return A list of all parameters
#' @import foreach
#' @import doParallel
#' @importFrom parallel detectCores
#' @export
lda_reshaped <- function(N, K, max_iter=50, thresh=1e-4, seed=NULL, cores=NULL, alpha=NULL){

  #define parameters
  V <- ncol(N)
  D <- nrow(N)

  loglik <- rep(NA, max_iter) #actually the lower bound on the log likelihood
  conv <- F

  if(is.null(alpha)) alpha <- 1/K

  if(is.null(cores)) cores <- detectCores()
  registerDoParallel(cores)

  #initialise variables (phi and gamma are reinitialised each E step)
  phis <- array(NA, c(V, K, D))
  gammas <- matrix(NA, D, K)

  #initialise beta using a random K documents (using a seed if given)
  if(!is.null(seed)) set.seed(seed)
  beta <- initalise_beta2(N, V, K, D)

  for(iter in 1:max_iter){
    message("Iteration", iter)

    #E-step
    res_lists <- foreach (d=1:D) %dopar% {
      e_step_d2(gamm=gammas[d,], phi=phis[,,d], alpha, beta, n=N[d,], V, K)
    }
    #combine results
    for(d in 1:D){
      phis[,,d] <- res_lists[[d]]$phi
      gammas[d,] <- res_lists[[d]]$gamm
    }

    # M-step
    beta <- update_beta2(beta, phis, N, V, K)

    #Check for convergence
    loglik[iter] <- loglik_corp2(phis, gammas, alpha, beta, N, V, K, D)
    if(L_converged(loglik, iter, thresh)){
      conv <- T
      break
    }
  }

  #retrieve estimates for thetas
  thetas <- gammas / rowSums(gammas)

  return(list("beta"=beta, "thetas"=thetas,
              "phis"=phis, "gammas"=gammas,
              "L"=loglik[iter],
              "Ls"=loglik[1:iter],
              "alpha"=alpha, "K"=K,
              "iterations"=iter,
              "converged"=conv))
}


