#applying original lda method to the text dataset
setwd("/home/an20830/Documents/COMPASS/Research Project/MSLDA/analysis/text")
load("data/MyCorpus.Rdata")
library(RColorBrewer)
library(tidyverse)

#originial lda
source("/home/an20830/Documents/COMPASS/Research Project/MSLDA/R/LDA_original.R")

#running it n times and plotting the lines together
n <- 5
K <- 3

#issues with my implementation - some initalisations work some dont
#so putting in a loop with try catch
res <- vector("list", n)
for(i in 1:n){

  #basically keep trying until it doesnt fail
  worked <- F
  while(worked == F){
    temp <- try(lda_blei03(docs, K, max_iter=30, thresh=10), silent=TRUE)
    if(class(temp) == "try-error") {
      print("trying again")
    } else {
      res[[i]] <- temp
      worked <- T
    }
  }
}

save(res, file="res/5RunsLDAOriginal.Rdata")
load("res/5RunsLDAOriginal.Rdata")
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

ggsave("res/5RunsLDAOriginal_L.jpg", height=4, width=5)


#index 3 different results
ll <- sapply(res, function(r) r$loglik)
#ps <- c(which.min(ll),which.min(abs(ll-40)),which.max(ll))
#ps <- c(3, 6, 9)
ps <- 1:5


### Plot alphas
alphas <- lapply(res, function(r) r$alpha)
pd <- as.data.frame(matrix(NA, nrow=n+1, ncol=K+1))
colnames(pd) <- c("run", 1:K) #"a1", "a2", "a3")
pd[1,] <- c("True", alpha_true)
pd[-1,1] <- paste("run", 1:n) #c("Low L", "Mid L", "High L")
for(i in 1:n){
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
  ylim(c(0, 0.3)) +
  scale_color_manual(values = cols_ps) +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x=element_blank(), panel.grid.minor.x = element_blank()) +
  facet_grid(~run, switch = "x")

ggsave("res/5RunsLDAOriginal_a.jpg", height=4, width=5)

### Plot mixtures
best <- which.max(sapply(lls, max))
source("/home/an20830/Documents/COMPASS/TB2/Mini Project/MSLDA/R/plot_funcs.R")
cols <- brewer.pal(6, "Paired")
(p0 <- plot_mixture(thetas_true, width=0.6) +
    scale_fill_manual(values=cols[c(6, 4, 2)]) +
    labs(title="True mixture"))
(p1 <- plot_mixture(res[[ps[best]]]$gammas, topic_label = LETTERS[1:3], width=0.6) +
    scale_fill_manual(values=cols[c(5, 3, 1)]) +
    labs(title="Estimated mixture - best run"))

jpeg("res/LDAOriginal_mix.jpg", width=1000, height=350)
grid.arrange(p0, p1, ncol=2)
dev.off()


#running with different k
ks <- 2:6
n <- 2
res <- matrix(data=NA, length(ks), n)
rownames(res) <- ks
colnames(res) <- paste("rep", 1:n)

for(j in 1:length(ks)){
  for(i in 1:n){
    worked <- F
    while(worked == F){
      temp <- try(lda_blei03(docs, K=ks[j], max_iter=30, thresh=10), silent=TRUE)
      if(class(temp) == "try-error") {
        print("trying again")
      } else {
        worked <- T
        res[j, i] <- temp$loglik
      }
    }
  }
}
save(res, file="res/LDAOriginal_Ks.Rdata")

load("res/LDAOriginal_Ks.Rdata")
res <- as.data.frame(cbind(ks, res))
colnames(res)[1] <- "K"
res_lls <- pivot_longer(res, cols=-"K", names_to="Rep", values_to="L")
max_line <- data.frame(K=ks, Max=apply(res[,2:(n+1)], 1, max))

ggplot(res_lls, aes(x=K, y=L, group=1))+
  geom_line(data=max_line, aes(x=K, y=Max), colour="grey") +
  geom_point() +
  labs(x="number of topics", y="L") +
  theme_minimal()

ggsave("res/LDAOriginal_LvK.jpg")
