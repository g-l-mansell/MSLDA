### Implementing the LDA algorithm from Blei 2003 paper

update_phi_d <- function(phi, gamm, beta, w, K){
  beta_prime <- beta[, w]
  phi <- t(beta_prime * exp(digamma(gamm)))

  #normalise rows of phi to sum to 1
  phi[rowSums(phi)==0, ] <- 1/K #to prevent dividing by 0
  phi <- phi / rowSums(phi)
  phi[phi==0] <- 1e-16 #and write over any 0s
  return(phi)
}

update_gamma_d <- function(phi, alpha){
  gamm <- alpha + colSums(phi)
  return(gamm)
}

gamma_converged <- function(gamm, gamm_old, tol=0.0001){
  diff <- mean(abs(gamm - gamm_old))
  return(diff < tol)
}

e_step_d <- function(gamm, phi, alpha, beta, w, V, K){
  gamm <- rep(1/K, K) #alpha + length(w)/K

  for(iter in 1:20){
    gamm_old <- gamm

    #update phi and gamma for this document
    phi <- update_phi_d(phi, gamm, beta, w, K)
    gamm <- update_gamma_d(phi, alpha)

    #check for convergence
    if(gamma_converged(gamm, gamm_old)) break
  }
  return(list("gamm"=gamm, "phi"=phi))
}

#testing new update beta method based on
#https://github.com/akashii99/Topic-Modelling-with-Latent-Dirichlet-Allocation/blob/master/LDA_blei_implement.ipynb
update_beta <- function(beta, phis, docs, V, K, D){
  for(v in 1:V){
    beta[,v] <- 0
    for(d in 1:D){
      w <- docs[[d]]
      phi <- phis[[d]]
      wnj <- matrix(w == v, nrow(phi), K)
      beta[,v] <- beta[,v] + colSums(phi * as.numeric(wnj))
    }
  }
  #normalise rows of beta to sum to 1
  beta <- beta / rowSums(beta)

  #write over any 0s
  beta[beta==0] <- 1e-16

  return(beta)
}

#only terms of log-likelihood that contain alpha - to check alpha convergence
loglik_alpha <- function(alpha, gammas, D){
  L <- D * (lgamma(sum(alpha)) - sum(lgamma(alpha)))
  for(d in 1:D){
    gamm <- gammas[d,]
    dig <- digamma(gamm) - digamma(sum(gamm))
    L <- L + sum((alpha-1) * dig)
  }
  return(L)
}

update_alpha <- function(alpha, gammas, max_iter=50, thresh=1e-4){
  D <- nrow(gammas)
  loglik <- rep(NA, max_iter)
  for(iter in 1:max_iter){
    #calculate g using dL/da in section A.4.2
    term1 <- digamma(sum(alpha)) - digamma(alpha)
    term2 <- colSums(digamma(gammas) - digamma(rowSums(gammas)))
    g <- D*term1 + term2

    #calculate H as diag(h) + z (z just a new constant)
    z <- D *trigamma(sum(alpha))
    h <- - D * trigamma(alpha)

    #use linear time Newton Raphson method in section A.2
    c <- sum(g/h) / (z^-1 + sum(h^-1))
    alpha <- alpha - (g-c)/h

    #check for convergence
    loglik[iter] <- loglik_alpha(alpha, gammas, D)
    if(iter > 5){
      if(abs(loglik[iter] - loglik[iter-1]) < thresh) break
    }
  }
  return(alpha)
}

#find the log likelihood of a document using eq 15
#not currently used
loglik_doc <- function(gamm, phi, alpha, beta, w, K){
  N <- length(w)
  t <- digamma(gamm) - digamma(sum(gamm))
  t_prime <- matrix(t, N, K, byrow=T)
  beta_prime <- beta[, w]

  L <- lgamma(sum(alpha)) - lgamma(sum(gamm)) +
    sum(lgamma(gamm) - lgamma(alpha) + (alpha-gamm)*t) +
    sum(phi * (t_prime + log(t(beta_prime)) - log(phi)))
  return(L)
}


#find the overall log-likelihood (sum of documents)
loglik_corp <- function(gammas, phis, alpha, beta, docs, D, K){
  t <- digamma(gammas) - digamma(rowSums(gammas)) #would be capital T
  alpha_prime <- matrix(alpha, D, K, byrow=T)

  L <- D*(lgamma(sum(alpha)) - sum(lgamma(alpha))) -
    sum(lgamma(rowSums(gammas))) +
    sum((alpha_prime - gammas) * t + lgamma(gammas))

  for(d in 1:D){
    w <- docs[[d]]
    N <- length(w)
    t_prime <- matrix(t[d,], N, K, byrow=T)
    beta_prime <- beta[, w]
    L <- L + sum(phis[[d]] * (t_prime + log(t(beta_prime)) - log(phis[[d]])))
  }
  return(L)
}

initalise_beta <- function(docs, V, K, D){
  beta <- matrix(0, K, V)
  idx <- sample(1:D, K)
  for(k in 1:K){
    d <- idx[k]
    for(v in 1:V){
      beta[k, v] <- sum(docs[[d]] == v)
    }
  }
  beta <- beta / rowSums(beta)
  return(beta)
}


#' Run LDA as in the Blei 2003 paper
#'
#' @param docs a list containing all the documents, with the vocabulary encoded
#' e.g. docs[[1]] = c(1, 5, 2)  would represent the word indices from a pre-defined vocabulary
#' @param K the number of topics to look for
#' @param max_iter the maximum number of EM iterations to run
#' @param thresh if you want to set a specific threshold for L convergence,
#' if NULL it's defined as 0.01 percent of the previous value
#' @param seed if you want a reproducible result you can set a seed,
#' if NULL the topic distributions are initialised from K random documents,
#' @return A list of all parameters
#' @export
#' @order 1
lda_original <- function(docs, K, max_iter=50, thresh=NULL, seed=NULL){

  #define parameters
  D <- length(docs)
  V <- length(unique(unlist(docs)))
  loglik <- rep(NA, max_iter) #actually the lower bound on the log likelihood
  conv <- F

  #initialise variables (phi and gamma are reinitialised each E step)
  alpha <- rep(1/K, K)
  phis <- vector("list", D)
  gammas <- matrix(NA, D, K)

  #initialise beta using a random K documents (using a seed if given)
  if(!is.null(seed)) set.seed(seed)
  beta <- initalise_beta(docs, V, K, D)

  for(iter in 1:max_iter){
    print(paste("Iteration", iter))

    #E-step
    for(d in 1:D){
      temp <- e_step_d(gammas[d,], phis[[d]], alpha, beta, docs[[d]], V, K)
      gammas[d,] <- temp$gamm
      phis[[d]] <- temp$phi
    }

    #M-step
    beta <- update_beta(beta, phis, docs, V, K, D)
    alpha <- update_alpha(alpha, gammas, max_iter=20, thresh=0.1)

    #Check for convergence
    loglik[iter] <- loglik_corp(gammas, phis, alpha, beta, docs, D, K)

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

  #retrieve estimates for thetas (proportional to gammas, or officially theta <- Dir(gamma))
  thetas <- gammas / rowSums(gammas)

  return(list("beta"=beta, "thetas"=thetas,
              "phis"=phis, "gammas"=gammas,
              "L"=loglik[iter],
              "Ls"=loglik[1:iter],
              "alpha"=alpha, "K"=K,
              "iterations"=iter,
              "converged"=conv))
}


