## Generate dataset
library(MCMCpack) #for the dirchlet function - implement myself later
set.seed(21)

# Use the generative model from Blei 2003 to generate documents
# note the Hoffman 2010 paper has a slightly different model

D <- 100 #number of documents

themes <- read.csv("themes.csv", header = F)
K <- ncol(themes)
vocab <- unique(as.vector(as.matrix(themes)))
vocab <- vocab[!vocab==""] #not sure why blanks were not read as NAs
V <- length(vocab)

# define alpha
alpha <- rep(0.1, K)

# define beta using the themes csv
beta <- matrix(NA, K, V)
colnames(beta) <- vocab
for(i in 1:V){
  beta[,i] <- colSums(themes == vocab[i])
}
#normalise so probabilities sum to 1
beta <- beta / rowSums(beta)


# For each document: draw theta_i ~ Dirichlet(alpha)
thetas <- rdirichlet(D, alpha)
Ns <- rpois(D, 500) #N[d] is number of words in document d
docs <- vector("list", D) #w[[d]] is the vector of N words

for(d in 1:D){
  #draw topics z_n ~ Multinomial(theta_i)
  N <- Ns[d]
  z <- sample(1:K, N, replace=T, prob=thetas[d,])
  w <- rep(NA, N)
  for(n in 1:N){
    #draw words from topics
    w[n] <- sample(1:V, 1, prob=beta[z[n],])
  }
  docs[[d]] <- w
}

#look at an example 'document'
#paste(vocab[docs[[5]]], collapse = " ")

#matrix of word counts
counts <- matrix(NA, D, V)
for(d in 1:D){
  for(v in 1:V){
    counts[d,v] <- sum(docs[[d]] == v)
  }
}

#store the true parameter values
alpha_true <- alpha
beta_true <- beta
thetas_true <- thetas

vocab <- data.frame(idx = 1:V, word = vocab)

#save the dataset so we can use it in other scripts
save(docs, counts, alpha_true, beta_true, thetas_true, vocab, file="MyCorpus.Rdata")

#also the words as a csv for easy viewing
# docs <- sapply(docs, function(w) paste(vocab[w,2], collapse = " "))
# docs <- cbind(round(thetas_true,2), docs)
# colnames(docs) <- c("breakfast", "pets", "books", "document")
# write.csv(docs, "MyCorpus.csv", row.names = F)
