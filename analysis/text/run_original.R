#applying original lda method to the text dataset
load("MyCorpus.Rdata")
library(RColorBrewer)

#originial lda
source("/home/an20830/Documents/COMPASS/TB2/Mini Project/MSLDA/R/LDA_original.R")

#running it n times and plotting the lines together
n <- 10

#issues with my implementation - some initalisations work some dont
#so putting in a loop with try catch
res <- vector("list", n)
for(i in 1:n){

  #basically keep trying until it doesnt fail
  worked <- F
  while(worked == F){
    temp <- try(lda_blei03(docs, K=3, max_iter=30), silent=TRUE)
    if(class(temp) == "try-error") {
      print("trying again")
    } else {
      res[[i]] <- temp
      worked <- T
    }
  }
}

save(res, file="RandomRunsBlei03.Rdata")
load("RandomRunsBlei03.Rdata")
lls <- lapply(res, function(r) r$lls)


### Plot L over iterations
cols <- brewer.pal(n+1, "Set3")[-2] #get rid of the yellow
pd <- as.data.frame(matrix(NA, nrow=10, ncol=31))
colnames(pd) <- c("run", 1:30)
pd[,1] <- names(cols) <- paste("run", 1:n)
for(i in 1:n){
  iters <- length(lls[[i]])
  pd[i, 2:(iters+1)] <- lls[[i]]
}
pd <- pivot_longer(pd, cols=-"run", values_to = "ll", names_to="iter") %>%
  filter(!is.na(ll))
class(pd$iter) <- "numeric"
ggplot(pd, aes(x=iter, y=ll, color=run)) +
  geom_line() +
  labs(x="Iteration", y="L") +
  scale_color_manual(values = cols) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave("LDA_runs_lls.jpg")


#index 3 differnt results
ll <- sapply(res, function(r) r$loglik)
#ps <- c(which.min(ll),which.min(abs(ll-40)),which.max(ll))
ps <- c(3, 6, 9)


### Plot alphas
alphas <- lapply(res, function(r) r$alpha)
pd <- as.data.frame(matrix(NA, nrow=4, ncol=4))
colnames(pd) <- c("run", 1:3) #"a1", "a2", "a3")
pd[1,] <- c("True", alpha_true)
pd[-1,1] <- c("Low L", "Mid L", "High L")
for(i in 1:3){
  pd[(i+1),2:4] <- alphas[[ps[i]]]
}
pd <- pivot_longer(pd, cols=-"run", names_to="a", values_to = "values")
class(pd$values) <- class(pd$a) <- "numeric"
pd$run <- factor(pd$run, levels=unique(pd$run))

cols_ps <- c(1, cols[ps])
names(cols_ps) <- levels(pd$run)

ggplot(pd, aes(x=a, y=values, colour=run)) +
  geom_point() +
  labs(x="", y="alpha_k") +
  xlim(c(-1, 5)) +
  ylim(c(0, 0.4)) +
  scale_color_manual(values = cols_ps) +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x=element_blank(), panel.grid.minor.x = element_blank()) +
  facet_grid(~run, switch = "x")

ggsave("LDA_runs_alphas.jpg")

### Plot mixtures
source("/home/an20830/Documents/COMPASS/TB2/Mini Project/MSLDA/R/plot_funcs.R")
cols <- brewer.pal(6, "Paired")
(p0 <- plot_mixture(thetas_true, width=0.6) +
    scale_fill_manual(values=cols[c(6, 4, 2)]) +
    labs(title="True mixture"))
(p1 <- plot_mixture(res[[ps[3]]]$gammas, topic_label = LETTERS[1:3], width=0.6) +
    scale_fill_manual(values=cols[c(5, 3, 1)]) +
    labs(title="Estimated mixture - best run"))

jpeg("LDA_runs_mixtures.jpg", width=1000, height=350)
grid.arrange(p0, p1, ncol=2)
dev.off()


#running with different k
k <- 5
n <- 5
pd <- matrix(data=NA, k, n)
rownames(pd) <- 1:k+1
colnames(pd) <- paste("rep", 1:n)

for(k in 1:k){
  for(i in 1:n){
    worked <- F
    while(worked == F){
      temp <- try(lda_blei03(docs, K=(k+1), max_iter=30), silent=TRUE)
      if(class(temp) == "try-error") {
        print("trying again")
      } else {
        worked <- T
        pd[k, i] <- temp$loglik
      }
    }
  }
}
save(pd, file="Repeats.Rdata")
pd2 <- as.data.frame(pd)
pd2$K <- rownames(pd)
pd2 <- pivot_longer(pd2, cols=-"K", values_to="value", names_to="repeat")
ggplot(pd2, aes(x=K, y=value, group=1))+
  geom_smooth(method="loess", se=F, colour="grey", lwd=0.8) +
  geom_point() +
  labs(x="number of topics", y="L") +
  theme_minimal()

ggsave("LDA_vs_K.jpg")
