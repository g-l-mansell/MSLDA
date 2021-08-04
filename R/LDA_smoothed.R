### Implementing the batch LDA algorithm from Hoffman 2010 paper

update_phi_d3 <- function(phi, gamm, lambda, V, K){
  t <- digamma(gamm) - digamma(sum(gamm))
  B <- digamma(lambda) - matrix(digamma(rowSums(lambda)), K, V)

  phi <- t(exp(matrix(t, K, V) + B))
  phi <- phi / rowSums(phi)
  return(phi)
}

update_gamma_d3 <- function(gamm, alpha, n, phi, V, K){
  temp <- matrix(n, V, K) * phi
  gamm <- colSums(temp) + alpha
  return(gamm)
}

e_step_d3 <- function(gamm, phi, alpha, lambda, n, V, K){
  gamm <- rep(1/K, K)
  for(iter in 1:30){
    gamm_old <- gamm

    #update phi and gamma for this document
    phi <- update_phi_d3(phi, gamm, lambda, V, K)
    gamm <- update_gamma_d3(gamm, alpha, n, phi, V, K)

    #check for convergence
    if(gamma_converged(gamm, gamm_old)) break
  }
  return(list("gamm"=gamm, "phi"=phi))
}

update_lambda <- function(lambda, eta, N, phis, V, K, D){
  #make a temp matrix with D rows and V columns (same shape as N) per each k
  for(k in 1:K){
    phis_k <- t(phis[,k,])
    counts_phis <- phis_k * N
    lambda[k,] <- eta + colSums(counts_phis)
  }
  return(lambda)
}

loglik_corp3 <- function(phis, gammas, lambda, N, alpha, eta, V, K, D){
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

initalise_lambda <- function(N, K, NMF){
  if(NMF){
    print("Intialising lambda with NMF")
    lambda <- nmf(N, K)$H + 1
  } else{
    idx <- sample(1:nrow(N), K)
    lambda <- N[idx,] + 1
  }
  return(lambda)
}


#' Run LDA adapted to use a count matrix
#'
#' @inheritParams lda_reshaped
#' @param NMF logical indicating if lambda should be initialised using non-negative matrix factorisation,
#' if FALSE it is generated using K random documents
#' @param eta if you want to set the exchangeable Dirichlet parameter for beta,
#' if NULL a default value of 1/K is used
#' @import foreach
#' @import doParallel
#' @importFrom parallel detectCores
#' @export
lda_smoothed <- function(N, K, max_iter=50, thresh=NULL, seed=NULL, cores=NULL, alpha=NULL, eta=NULL, NMF=FALSE){

  #get number of samples and words from N
  #define parameters
  V <- ncol(N)
  D <- nrow(N)
  loglik <- rep(NA, max_iter) #actually the lower bound on the log likelihood
  conv <- F
  if(is.null(alpha)) alpha <- 1/K
  if(is.null(eta)) eta <- 1/K

  #initialise variables (phi and gamma are reinitialised each E step)
  phis <- array(NA, c(V, K, D))
  gammas <- matrix(NA, D, K)

  #initialise lambda using K random documents or with NMF
  if(!is.null(seed)) set.seed(seed)
  lambda <- initalise_lambda(N, K, NMF)

  for(iter in 1:max_iter){
    print(paste("Iteration", iter))

    # E-step
    res_lists <- foreach (d=1:D) %dopar% {
      e_step_d3(gamm=gammas[d,], phi=phis[,,d], alpha, lambda, n=N[d,], V, K)
    }
    #combine results
    for(d in 1:D){
      phis[,,d] <- res_lists[[d]]$phi
      gammas[d,] <- res_lists[[d]]$gamm
    }

    # M-step
    lambda <- update_lambda(lambda, eta, N, phis, V, K, D)

    #Check for convergence
    loglik[iter] <- loglik_corp3(phis, gammas, lambda, N, alpha, eta, V, K, D)

    #if no convergence threshold is given, use a default % of the last value
    if(iter==1 & is.null(thresh)) default_thresh <- T

    if(iter > 5){
      if(default_thresh) thresh <- abs(1e-4 * loglik[iter-1])
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


