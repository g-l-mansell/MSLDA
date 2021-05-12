#applying smoothed lda method to the text dataset
load("MyCorpus.Rdata")

#smoothed parellelised lda
#still a question about L
#and comparison to NMF - e.g. can this be a starting point - or can we start with random samples like the original?
source("../../R/LDA_funcs_par.R")
source("../../R/NMF_funcs.R")
source("../../R/plot_funcs.R")

#orginial had 4 plots
#LDA vs K
#10 runs L vs iter
#scatter of alpha of different maxima
#example estimated proportions
test <- lda(counts, K=3, nmf=T)
