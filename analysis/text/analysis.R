#applying smoothed lda to the text dataset
load("MyCorpus.Rdata")
source("/home/an20830/Documents/COMPASS/TB2/Mini Project/MSLDA/R/LDA_funcs_par.R")
D <- length(docs) #number of documents
Ns <- sapply(docs, function(w) length(w)) #number of words in each document
W <- nrow(vocab) #size of vocabulary
K <- 3 #number of topics

results <- lda(counts, K, thresh=20)
#View(round(results$thetas, 2))

#applying nmf to the text dataset
source("/home/an20830/Documents/COMPASS/TB2/Mini Project/MSLDA/R/NMF_funcs.R")
test <- nmf(counts, K=3)
#View(round(test$W, 2))


#plots
source("/home/an20830/Documents/COMPASS/TB2/Mini Project/MSLDA/R/plot_funcs.R")
cols <- brewer.pal(6, "Paired")

p1 <- plot_mixture(test$W, topic_label = letters[1:3]) +
  scale_fill_manual(values=cols[c(1,3,5)]) +
  labs(title="NMF Estimated Mixture")

p2 <- plot_mixture(results$thetas, topic_label = LETTERS[1:3]) +
  scale_fill_manual(values=cols[c(1,3,5)]) +
  labs(title="LDA Estimated Mixture")

p3 <- plot_mixture(thetas_true) +
  scale_fill_manual(values=cols[c(2,4,6)]) +
  labs(title="True Generating Mixture")

#jpeg("/home/an20830/Documents/COMPASS/TB2/Mini Project/Report/text_proportions.jpg", width=1000, height=350)
grid.arrange(p3, p1, p2, ncol=3)
#dev.off()



#testing/comparing my different lda implementations

#applying smoothed lda to the text dataset
load("MyCorpus.Rdata")

#originial lda
source("/home/an20830/Documents/COMPASS/TB2/Mini Project/MSLDA/R/LDA_original2.R")

#smoothed lda (basic implementation)
source("/home/an20830/Documents/COMPASS/TB2/Mini Project/MSLDA/R/LDA_funcs.R")
res2 <- lda(counts, K=3, max_iter=20)
plot(res2$lls)

#smoothed lda (optimised implementation)
source("/home/an20830/Documents/COMPASS/TB2/Mini Project/MSLDA/R/LDA_funcs_par.R")
res3 <- lda(counts, K=3, max_iter=10)
plot(res3$lls)
