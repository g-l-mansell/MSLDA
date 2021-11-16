loglik_corp3 <- function(phis, gammas, alpha, beta, N, V, K, D){
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


#' @describeIn lda_reshaped Alpha is fixed
#' @inheritParams lda_noalpha
#' @param N matrix of word counts
#' @return A list of all parameters
#' @import foreach
#' @import doParallel
#' @importFrom parallel detectCores
#' @export
#' @order 2
lda_reshaped_noalpha <- function(N, K, max_iter=50, thresh=1e-4, seed=NULL, cores=NULL, alpha=NULL){

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
    loglik[iter] <- loglik_corp3(phis, gammas, alpha, beta, N, V, K, D)
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


