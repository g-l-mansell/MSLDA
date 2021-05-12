library(ggplot2)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(gridExtra)

#writing a function to plot mixing proportions as a horizontal ggplot barplot
plot_mixture <- function(dat, nsamples=10, sample_labels=NULL, topic_label=1:ncol(dat), width=0.9){

  if(is.null(sample_labels)) sample_labels <- paste("doc", 1:nsamples)
  Sample <- factor(sample_labels, levels=rev(sample_labels))

  plot_dat <- as.data.frame(dat)
  plot_dat <- plot_dat[1:nsamples,]
  colnames(plot_dat) <- paste("Topic", topic_label)
  plot_dat <- cbind(plot_dat, Sample)

  plot_dat <- plot_dat %>%
    pivot_longer(cols=-"Sample", names_to = "Topic", values_to="Proportion")

  p <- ggplot(plot_dat, aes(fill=Topic, x=Sample, y=Proportion)) +
      geom_bar(position="fill", stat="identity", width=width) +
      coord_flip() +
      labs(x="", fill="") +
      theme_minimal() +
      theme(text = element_text(size = 14))

  return(p)
}

plot_clusters <- function(dat, lables){
  #create a heatmap of groupings
  #for each sample find its max topic
  max_topics <- apply(dat, 1, which.max)
  hmd <- data.frame("Strain" = factor(lables, levels=rev(unique(lables))), Max_Topic=as.factor(max_topics))

  hmd <- hmd %>%
    group_by(Strain, Max_Topic) %>%
    count() %>%
    ungroup() %>%
    complete(Strain, Max_Topic, fill = list(n = 0))

  p <- ggplot(hmd, aes(x=Max_Topic, y=Strain, fill=n))+
    geom_tile() +
    scale_fill_distiller(palette = "Blues", direction=+1) +
    theme_minimal() +
    labs(y="") +
    theme(text = element_text(size = 14))

  return(p)
}
