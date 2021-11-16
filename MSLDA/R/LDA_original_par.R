### Implementing the LDA algorithm from Blei 2003 paper, but this time running the E-step in parallel

#' @describeIn lda_original Runs E-step in parallel
#' @param cores number of cores to run the E-step in parallel,
#' if NULL all detected cores are used
#' @import foreach
#' @import doParallel
#' @importFrom parallel detectCores
#' @export
#' @order 1
lda_original_par <- function(docs, K, max_iter=50, thresh=1e-4, seed=NULL, cores=NULL){

  #define parameters
  D <- length(docs)
  V <- length(unique(unlist(docs)))
  loglik <- rep(NA, max_iter) #actually the lower bound on the log likelihood
  conv <- F
  if(is.null(cores)) cores <- detectCores()
  registerDoParallel(cores)

  #initialise variables (phi and gamma are reinitialised each E step)
  alpha <- rep(1/K, K)
  phis <- vector("list", D)
  gammas <- matrix(NA, D, K)

  #initialise beta using a random K documents (using a seed if given)
  if(!is.null(seed)) set.seed(seed)
  beta <- initalise_beta(docs, V, K, D)

  for(iter in 1:max_iter){
    message("Iteration", iter)

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
    alpha <- update_alpha(alpha, gammas, max_iter=20, thresh=0.1, K)

    #Check for convergence
    loglik[iter] <- loglik_corp(gammas, phis, alpha, beta, docs, D, K)
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


