### Implementing the LDA algorithm from Blei 2003 paper
# Does not work!!

#here we have M documents, V words, K topics
#each document has N_d words which we observe
#the topic each word w_dn comes from z_dn is latent (unobserved)
#what we want to find are theta_d (the topic proportions),
#alpha (the Dirichlet parameter), and beta (the topic distributions over words)

#we introduce variational parameters phi_d and gamma_d per document
#then use an EM like algorithm
#E-step: estimates of alpha and beta are fixed, update phis and gammas
#M-step: estimates of phis and gammas are fixed, update alpha and beta

#most of these updates have a simple formula
#but alpha needs to be optimised using Newton-Raphson
#we continue the updates until the likelihood has converged

### Defining useful update function
e_step_d <- function(gamm, phi, alpha, beta, w, thresh=1e-4){
  N <- length(w)
  K <- length(alpha)
  gamm <- alpha + N/K
  phi <- matrix(1/K, N, K)

  #loop until gamma converges (but using a limited for loop to prevent getting stuck in a while loop)
  for(iter in 1:20){
    gamm_old <- gamm

    #update each element of phi - could vectorise to speed up?
    for(n in 1:N){
      for(k in 1:K){
        phi[n, k] <- beta[k, w[n]] * exp(digamma(gamm[k]))
      }
    }

    #normalise rows of phi to sum to 1
    phi[rowSums(phi)==0, ] <- 1/K #to prevent dividing by 0!
    phi <- phi / rowSums(phi)

    #update gamma
    gamm <- alpha + colSums(phi)
    #gamm <- gamm / sum(gamm) #test test test

    #check gamma (more efficient than checking loglik)
    if(max(abs(gamm - gamm_old)) < thresh) break
  }
  return(list(phi, gamm))
}

e_step <- function(gammas, phis, alpha, beta, docs, thresh){
  #update the variational parameters for each document in turn
  #could update each doc in parallel
  for(d in 1:length(docs)){
    res <- e_step_d(gamm=gammas[d,], phi=phis[[d]], alpha, beta, w=docs[[d]], thresh=thresh)
    phis[[d]] <- res[[1]]
    gammas[d,] <- res[[2]]
  }
  return(list(gammas, phis))
}

#testing new update beta method based on
#https://github.com/akashii99/Topic-Modelling-with-Latent-Dirichlet-Allocation/blob/master/LDA_blei_implement.ipynb
update_beta <- function(beta, phis, docs){
  for(v in 1:ncol(beta)){
    beta[,v] <- 0
    for(d in 1:length(docs)){
      w <- docs[[d]]
      phi <- phis[[d]]
      wnj <- matrix(w == v, nrow(phi), ncol(phi))
      beta[,v] <- beta[,v] + colSums(phi * as.numeric(wnj))
    }
  }
  #normalise rows of beta to sum to 1
  beta <- beta / rowSums(beta)
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

nr_alpha <- function(alpha, gammas, max_iter=20, thresh=1e-4){
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
loglik_doc <- function(gamm, phi, alpha, beta, w, K){

  #calculate the third line
  term3 <- 0
  for(n in 1:length(w)){
    for(i in 1:K){
      term3 <- term3 + phi[n, i] * log(beta[i,w[n]])
    }
  }

  #store common variable as a K length vector
  dig <- digamma(gamm) - digamma(sum(gamm))

  L <- lgamma(sum(alpha))
  - sum(lgamma(alpha))
  + sum((alpha - 1) * dig)
  + sum(phi * dig)
  + term3
  - lgamma(sum(gamm))
  + sum(lgamma(gamm))
  - sum((gamm-1) * dig)
  - sum(phi * log(phi))

  return(L)
}

#find the overall log-likelihood (sum of documents)
loglik_corp <- function(gammas, phis, alpha, beta, docs, D, K){
  L <- 0
  for(d in 1:D){
    L <- L + loglik_doc(gammas[d,], phis[[d]], alpha, beta, docs[[d]], K)
  }
  return(L)
}



initalise_beta <- function(docs, K, V){
  beta <- matrix(0, K, V)
  idx <- sample(1:length(docs), K)
  for(k in 1:K){
    d <- idx[k]
    for(v in 1:V){
      beta[k, v] <- sum(docs[[d]] == v)
    }
  }
  beta <- beta / rowSums(beta)
  return(beta)
}


#function to perform lda
lda_blei03 <- function(docs, K, max_iter=20, thresh=NULL, seed=NULL){

  #find parameters
  D <- length(docs)
  V <- length(unique(unlist(docs)))
  Ns <- sapply(docs, length)
  loglik <- rep(NA, max_iter) #actually the lower bound on the log likelihood
  conv <- F

  #initialise variables
  alpha <- rep(1/K, K)
  #beta <- matrix(1/V, K, V)
  phis <- vector("list", D)
  gammas <- matrix(NA, D, K)

  #phi and gamma are reinitialised and tuned on every iteration
  #initialise beta using a random K documents
  # if a seed is given use it, but otherwise it will be random each run
  if(!is.null(seed)) set.seed(seed)
  beta <- initalise_beta(docs, K, V)

  for(iter in 1:max_iter){
    print(paste("Iteration", iter))

    e_res <- e_step(gammas, phis, alpha, beta, docs, thresh=0.1)
    gammas <- e_res[[1]]
    phis <- e_res[[2]]

    #M-step
    beta <- update_beta(beta, phis, docs)
    alpha <- update_alpha(alpha, gammas, max_iter=20, thresh=0.1)

    loglik[iter] <- loglik_corp(gammas, phis, alpha, beta, docs, D, K)
    #by default set the convergence threshold to 0.1% of the last value
    if(iter > 5){
      if(is.null(thresh)) thresh <- abs(0.001* loglik[iter-1])
      if(abs(loglik[iter] - loglik[iter-1]) < thresh){
        conv <- T
        break
      }
    }

  }

  return(list("beta"=beta, "alpha"=alpha,
              "phis"=phis, "gammas"=gammas,
              "loglik" = loglik[iter],
              "lls"=loglik[1:iter],
              "iterations"=iter,
              "converged"=conv))
}

#still not sure alpha update is correct, but at least gammas look reasonable

