### Load the text dataset
setwd("/home/an20830/Documents/COMPASS/TB2/Mini Project/MSLDA/analysis/text")
load("MyCorpus.Rdata")

source("../../R/plot_funcs.R")
source("../../R/NMF_funcs.R")
source("../../R/LDA_original.R")
source("../../R/LDA_funcs_par.R")

#Run 3 analyses
K <- 3
max_iter <- 200

#take the orginal best run from analysis.R
load("RandomRunsBlei03.Rdata")
res2 <- res[[9]]

#try the other two models three times
lls <- rep(-Inf, 3)
for(i in 1:3){
  res <- nmf(counts, K, max_iter)
  if(res$ll > lls[1]) res1 <- res

  res <- lda(counts, K, max_iter, alpha=0.1)
  if(res$loglik > lls[3]) res3 <- res
}

save(res1, res2, res3, file="CompareAlgorithms.Rdata")

#Plots
cols <- brewer.pal(6, "Paired")

(p1 <- plot_mixture(thetas_true, width=0.6) +
  scale_fill_manual(values=cols[c(2,4,6)]) +
  labs(title="True Generating Mixture"))

(p2 <- plot_mixture(res1$W, topic_label = LETTERS[1:3], width=0.6) +
  scale_fill_manual(values=cols[c(1,3,5)]) +
  labs(title="NMF Estimated Mixture"))

(p3 <- plot_mixture(res2$gammas, topic_label = LETTERS[1:3], width=0.6) +
  scale_fill_manual(values=cols[c(1,3,5)]) +
  labs(title="LDA Estimated Mixture"))

(p4 <- plot_mixture(res3$thetas, topic_label = LETTERS[1:3], width=0.6) +
  scale_fill_manual(values=cols[c(1,3,5)]) +
  labs(title="Smoothed LDA Estimated Mixture"))

png("Comparison.png", width=800, height=500)
grid.arrange(p1, p2, p3, p4, ncol=2)
dev.off()


#Try to match up the topics (should be easy?) and evaluate some sort of error
res_list <- list(res1$W, res2$gammas, res3$thetas)
true_max <- apply(thetas_true, 1, which.max)
mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

error <- rep(NA, 3)

for(i in 1:3){
  res <- res_list[[i]]
  res <- res / rowSums(res)

  model_max <- apply(res, 1, which.max)
  order <- rep(NA, 3)
  for(j in 1:K){
    samples <- which(true_max == j)
    order[j] <- mode(model_max[samples])
  }
  res <- res[,order]

  error[i] <- mean(abs(thetas_true-res))
}

print(error)
plot(error)

