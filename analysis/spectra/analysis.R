### Load the spectra dataset
setwd("/home/an20830/Documents/COMPASS/Research Project/MSLDA/analysis/spectra")
load("data/Spectra.Rdata")

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
save(res_list1, file="res/LDA_runs.Rdata")

#make 3 plots over iterations with replicates
cols <- brewer.pal(reps+1, "Set3")[-2] #get rid of yellow
titles <- c("NMF", "LDA", "LDA with NMF initalisation")

load("res/LDA_runs.Rdata")
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
ggsave("res/LDA_runs.png", width=9, height=2.5)


#plotting the results of the best run
cols <- brewer.pal(11, "Set3")[-c(2, 9)]
best <- which.max(sapply(res_list1, function(x) ifelse(max(x$lls)>0, NA, max(x$lls))))
res <- res_list1[[best]]
(p1 <- plot_mixture(res$thetas[seq(1, 72, 8),], nsamples=9, sample_label = idx[seq(1, 72, 8)], width=0.6) +
    scale_fill_manual(values=cols))
(p2 <- plot_clusters(res$thetas, idx))

png("res/LDA_results.png", width=800, height=300)
grid.arrange(p1, p2, ncol=2, widths=c(1.1, 1))
dev.off()


### Experiment 2
ks <- 5:10
max_iter <- 100
reps <- 1
res <- matrix(NA, length(ks), reps)
res <- cbind(ks, res)
colnames(res) <- c("K", paste("rep", 1:reps))

for(i in 1:reps){
  for(j in 1:length(ks)){
    temp <- lda(counts, K=ks[j], max_iter, NMF=F)
    res[j, i+1] <- temp$loglik
    plot(temp$lls)
  }
}

save(res, file="res/LDA_runs_Ks.Rdata")

#we want the scatter but with a line connecting the maximum points
max_line <- data.frame(K=ks, Max=apply(res[,2:(reps+1)], 1, max))
res_lls <- pivot_longer(as.data.frame(res), cols=-"K", names_to="Rep", values_to="L")

ggplot(res_lls, aes(x=K, y=L, group=1))+
  geom_line(data=max_line, aes(x=K, y=Max), colour="grey") +
  geom_point() +
  labs(x="number of topics", y="L") +
  theme_minimal()
ggsave("LDA_vs_K.jpg")


#replot the mixtures with the k=7 model
res7 <- lda(counts, K=7, max_iter=200, nmf=T) #rerunning until res7$loglik > mean(res[3, 2:4])
save(res7, file="res/LDA_results_K7.Rdata")
cols <- brewer.pal(9, "Set3")[-c(2, 9)]

(p1 <- plot_mixture(res7$thetas[seq(1, 72, 8),], nsamples=9, sample_label = idx[seq(1, 72, 8)], width=0.6) +
  scale_fill_manual(values=cols) +
  labs(title="Estimated Mixtures"))

(p2 <- prcomp(res7$theta)$x[, 1:2] %>%
  as.data.frame %>%
  mutate(Strain=idx) %>%
  ggplot(aes(x=PC1, y=PC2, colour=Strain)) +
    geom_point() +
    theme_minimal() +
    labs(title="PCA of Mixtures")+
    theme(text = element_text(size = 14)))

(p3 <- prcomp(counts)$x[, 1:2] %>%
  as.data.frame %>%
  mutate(Strain=idx) %>%
  ggplot(aes(x=PC1, y=PC2, colour=Strain)) +
    geom_point() +
    theme_minimal() +
    labs(title="PCA of Spectra")+
    theme(text = element_text(size = 14)))

png("res/LDA_results_K7.png", width=1300, height=400)
grid.arrange(p1, p2, p3, ncol=3)
dev.off()



### Plot the topic-distributions of the run above
beta <- res7$beta
png("res/BestTopics.png", width=600, height=600)
par(mfrow=c(2, 1))
plot(mz_locations, beta[6,], type="l", xlab="m/z", ylab="P(m/z | topic 6)",
     col="darkgrey", lwd=1.3, main="Topic 6 asssociated with Psm")
plot(mz_locations, beta[2,], type="l", xlab="m/z",  ylab="P(m/z | topic 2)",
     col="darkgrey", lwd=1.3, main="Topic 2 associated with RO340")
dev.off()

#plot of all samples and the spectra both zoomed in at m/z=3000
# pd <- as.data.frame(counts)
# colnames(pd) <- mz_locations
# pd$Strain <- idx
# pd$Sample <- factor(1:nrow(pd))
# pd <- pivot_longer(pd, cols=-c("Strain", "Sample"), values_to="Intensity", names_to="mz")
# class(pd$mz) <- "numeric"
#
# p1 <- ggplot(pd, aes(x=mz, y=Intensity, color=Strain, group=Sample)) +
#   geom_line(lwd=0.2) +
#   labs(x="m/z") +
#   theme_minimal() +
#   guides(color = guide_legend(override.aes = list(size = 0.6))) +
#   xlim(c(2900, 3100))
#
# pd2 <- as.data.frame(beta)
# colnames(pd2) <- mz_locations
# pd2$Topic <- factor(paste("Topic", 1:9))
# pd2 <- pivot_longer(pd2, cols=-"Topic", values_to="Intensity", names_to="mz")
# class(pd2$mz) <- "numeric"
#
# cols <- brewer.pal(11, "Set3")[-c(2, 9)]
# p2 <- ggplot(pd2, aes(x=mz, y=Intensity, color=Topic, group=Topic)) +
#   geom_line(lwd=0.5) +
#   labs(x="m/z") +
#   theme_minimal() +
#   guides(color = guide_legend(override.aes = list(size = 0.6))) +
#   xlim(c(2900, 3100)) +
#   scale_color_manual(values=cols)
#
# png("TopicsVsSpectra_3000.png", width=800, height=250)
# grid.arrange(p1, p2, ncol=2)
# dev.off()


### Experiment 2
alphas <- etas <- c(0.1, 0.5, 1, 10)
n <- length(alphas)
max_iter <- 100
res <- matrix(NA, n, n)

for(i in 1:n){
  for(j in 1:n){
    temp <- lda(counts, K=7, max_iter, nmf=F, alpha=alphas[i], eta=etas[j])
    res[i, j] <- temp$loglik
    plot(temp$lls)
  }
}

