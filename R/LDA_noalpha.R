### Adapting LDA_original_par to treat alpha as a univariate hyperparameter

#find the overall log-likelihood - only minor adjustments to the loglik_corp now alpha is univaraite not K-length
loglik_corp_a <- function(gammas, phis, alpha, beta, docs, D, K){
  t <- digamma(gammas) - digamma(rowSums(gammas)) #would be capital T

  L <- D*(lgamma(K*alpha) - K*lgamma(alpha)) -
    sum(lgamma(rowSums(gammas))) +
    sum((alpha - gammas) * t + lgamma(gammas))

  for(d in 1:D){
    w <- docs[[d]]
    N <- length(w)
    t_prime <- matrix(t[d,], N, K, byrow=T)
    beta_prime <- beta[, w]
    L <- L + sum(phis[[d]] * (t_prime + log(t(beta_prime)) - log(phis[[d]])))
  }
  return(L)
}

#' @describeIn lda_original Alpha is fixed
#' @param alpha if you want to set the exchangeable Dirichlet parameter for theta,
#' if NULL a default value of 1/K is used
#' @import foreach
#' @import doParallel
#' @importFrom parallel detectCores
#' @export
#' @order 2
lda_noalpha <- function(docs, K, max_iter=50, thresh=NULL, seed=NULL, cores=NULL, alpha=NULL){

  #define parameters
  D <- length(docs)
  V <- length(unique(unlist(docs)))
  loglik <- rep(NA, max_iter) #actually the lower bound on the log likelihood
  conv <- F
  if(is.null(cores)) cores <- detectCores()
  registerDoParallel(cores)
  if(is.null(alpha)) alpha <- 1/K

  #initialise variables (phi and gamma are reinitialised each E step)
  phis <- vector("list", D)
  gammas <- matrix(NA, D, K)

  #initialise beta using a random K documents (using a seed if given)
  if(!is.null(seed)) set.seed(seed)
  beta <- initalise_beta(docs, V, K, D)

  for(iter in 1:max_iter){
    print(paste("Iteration", iter))

    #E-step
    res_lists <- foreach (d=1:D) %dopar% {
      e_step_d(gamm=gammas[d,], phi=phis[[d]], alpha, beta, w=docs[[d]], V, K)
    }
    #combine results
    for(d in 1:D){
      phis[[d]] <- res_lists[[d]]$phi
      gammas[d,] <- res_lists[[d]]$gamm
    }

    #M-step
    beta <- update_beta(beta, phis, docs, V, K, D)

    #Check for convergence
    loglik[iter] <- loglik_corp_a(gammas, phis, alpha, beta, docs, D, K)

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


