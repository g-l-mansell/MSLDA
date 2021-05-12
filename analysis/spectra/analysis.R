### Load the main dataset
setwd("/home/an20830/Documents/COMPASS/TB2/Mini Project/MSLDA/analysis/spectra")
load("Spectra.Rdata")

source("../../R/plot_funcs.R")
source("../../R/NMF_funcs.R")
source("../../R/LDA_funcs_par.R")

#Experiment 1 - run K=9 NMF, LDA, LDA with NMF initialisation
#we want at least 3 runs of each, and a line plot of ll over iterations
#plot - best run topic mixtures / max topics heatmap

#Experiment 2 - run LDA for different values of K
#again we want at least 3 runs of each, and a plot of final ll vs number of clusters
#plot best run topic mixtures / max topics heatmap

#Experiment 3 - run a small grid search of alpha and eta?


#Not sure how to measure accuracy of the mixtures without just a plot?
#Note to also run an experiment comparing the accuracy of the differnt methods on the text data


### Experiment 1
reps <- 3
K <- 9
max_iter <- 200
res_list1 <- vector("list", reps*3)

for(i in 1:reps){
  #nmf
  res_list1[[i*3-2]] <- nmf(counts, K, max_iter)
  #lda
  res_list1[[i*3-1]] <- lda(counts, K, max_iter, nmf=F)
  #lda with nmf
  res_list1[[i*3]] <- lda(counts, K, max_iter, nmf=T)
}
save(res_list1, file="LDA_runs.Rdata")

#make 3 plots over iterations with replicates
cols <- brewer.pal(reps+1, "Set3")[-2] #get rid of yellow
titles <- c("NMF", "LDA", "LDA with NMF initalisation")

load("LDA_runs.Rdata")
res_lls <- as.data.frame(matrix(NA, nrow=reps*3, ncol=max_iter+2))
colnames(res_lls) <- c("Model", "Run", 1:max_iter)
for(i in 1:reps){
  for(j in 1:3){
    id <- i*3-(3-j)
    res_lls[id,1] <- titles[j]
    res_lls[id,2] <- i
    res <- res_list1[[id]]
    lls <- res$lls
    res_lls[id,3:(length(lls)+2)] <- lls
  }
}
res_lls <- pivot_longer(res_lls, cols=-c("Model", "Run"), values_to="L", names_to="Iteration") %>%
  filter(!is.na(L)) %>%
  mutate(Iteration=as.numeric(Iteration),
         Model=factor(Model, levels=titles))

ggplot(res_lls, aes(x=Iteration, y=L, color=Model, group=Run)) +
  geom_line() +
  facet_wrap(~Model, scales = "free") +
  theme_minimal() +
  scale_color_manual(values=cols) +
  #xlim(c(0, max_iter)) +
  theme(legend.position = "none")
ggsave("LDA_runs.png", width=9, height=2.5)


#plotting the results of the best run
cols <- brewer.pal(11, "Set3")[-c(2, 9)]
res <- res_list1[[6]]
(p1 <- plot_mixture(res$thetas[seq(1, 72, 8),], nsamples=9, sample_label = idx[seq(1, 72, 8)], width=0.6) +
    scale_fill_manual(values=cols))
(p2 <- plot_clusters(res$thetas, idx))

png("LDA_results.png", width=800, height=300)
grid.arrange(p1, p2, ncol=2)
dev.off()


### Experiment 2
ks <- 2:15
max_iter <- 100
res_list2 <- vector("list", length(ks))

for(j in 1:length(ks)){
  res_list2[[j]] <- lda(counts, K=ks[j], max_iter, nmf=T)
}
#save(res_list2, file="L_vs_K2.RData")

#ran the above twice and saved the results as (too big to do in 1 go)
files <- c("L_vs_K.RData", "L_vs_K2.RData")
reps <- length(files)
res_lls <- matrix(NA, nrow=length(ks), ncol=reps)

for(i in 1:2){
  load(files[i])
  res_lls[, i] <- sapply(res_list2, function(r) r$loglik)
}

res_lls <- as.data.frame(cbind(ks, res_lls))
colnames(res_lls) <- c("K", paste("rep", 1:2))
res_lls <- pivot_longer(res_lls, cols=-"K", names_to="Rep", values_to="L")

ggplot(res_lls, aes(x=K, y=L, group=1))+
  geom_smooth(method="loess", se=F, colour="grey", lwd=0.8) +
  geom_point() +
  labs(x="number of topics", y="L") +
  theme_minimal()
ggsave("LDA_vs_K.jpg")

# res <- res_list2[[13]]
# (p1 <- plot_mixture(res$thetas[seq(1, 72, 8),], nsamples=9, sample_label = idx[seq(1, 72, 8)], width=0.6))
# (p2 <- plot_clusters(res$thetas, idx))
# grid.arrange(p1, p2, ncol=2)


### Plot the topic-distributions of the best model from experiment 1
beta <- res_list1[[6]]$beta
for(i in 1:nrow(beta)){
  plot(mz_locations, beta[i,], type="l")
}

#the best cluster looks like topic 3 - Psm
png("Topic3.png", width=600, height=300)
plot(mz_locations, beta[3,], type="l", xlab="m/z", ylab="P(w | z = 3)", col="darkgrey", lwd=1.3)
dev.off()

#plot of all samples and the spectra both zoomed in at m/z=3000
pd <- as.data.frame(counts)
colnames(pd) <- mz_locations
pd$Strain <- idx
pd$Sample <- factor(1:nrow(pd))
pd <- pivot_longer(pd, cols=-c("Strain", "Sample"), values_to="Intensity", names_to="mz")
class(pd$mz) <- "numeric"

p1 <- ggplot(pd, aes(x=mz, y=Intensity, color=Strain, group=Sample)) +
  geom_line(lwd=0.2) +
  labs(x="m/z") +
  theme_minimal() +
  guides(color = guide_legend(override.aes = list(size = 0.6))) +
  xlim(c(2900, 3100))

pd2 <- as.data.frame(beta)
colnames(pd2) <- mz_locations
pd2$Topic <- factor(paste("Topic", 1:9))
pd2 <- pivot_longer(pd2, cols=-"Topic", values_to="Intensity", names_to="mz")
class(pd2$mz) <- "numeric"

cols <- brewer.pal(11, "Set3")[-c(2, 9)]
p2 <- ggplot(pd2, aes(x=mz, y=Intensity, color=Topic, group=Topic)) +
  geom_line(lwd=0.5) +
  labs(x="m/z") +
  theme_minimal() +
  guides(color = guide_legend(override.aes = list(size = 0.6))) +
  xlim(c(2900, 3100)) +
  scale_color_manual(values=cols)

png("TopicsVsSpectra_3000.png", width=1000, height=400)
grid.arrange(p1, p2, ncol=2)
dev.off()
