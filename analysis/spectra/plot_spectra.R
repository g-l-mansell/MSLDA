### Load the main dataset
setwd("/home/an20830/Documents/COMPASS/TB2/Mini Project/MSLDA/analysis/spectra")

dat <- read.csv("/home/an20830/Documents/COMPASS/TB2/Mini Project/MALDI/Y_New_Raw_all_S_Aureus_Strains_Rescale_Median.csv", header=F)
idx <- dat[-1,1]
mz_locations <- as.numeric(gsub("X", "", dat[1,])[-1])
counts <- dat[-1, -1]

#remove first two columns which are all 0 and long tails of the spectra
n <- 3:which(mz_locations > 13000)[1]
counts <- counts[,n]
mz_locations <- mz_locations[n]

### Plot using ggplot facetwrap
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(gridExtra)

lookups <- data.frame("Name" = c("je2", "Ne1532", "Agr", "Hla", "LAC", "Psm", "RN4220", "RO340", "Hld"),
                      "ID" = c(2072, 1791532, 314, 309, 208, 336, 1604220, 161340, 312))

pd <- counts
colnames(pd) <- mz_locations
pd$Strain <- factor(lookups$Name[match(idx, lookups$ID)])
pd$Sample <- factor(1:nrow(pd))

#to just plot 1 sample of each strain
pd1 <- pd[match(lookups$Name, pd$Strain),]

pd1 <- pivot_longer(pd1, cols=-c("Strain", "Sample"), values_to="Intensity", names_to="mz")
class(pd1$mz) <- "numeric"

(p1 <- ggplot(pd1, aes(x=mz, y=Intensity, color=Strain)) +
    geom_line(lwd=0.4) +
    labs(x="m/z") +
    facet_wrap(~Strain) +
    theme_minimal() +
    theme(legend.position = "none"))

(p2 <- ggplot(pd1, aes(x=mz, y=Intensity, color=Strain)) +
    geom_line(lwd=0.2) +
    labs(x="m/z") +
    theme_minimal() +
    guides(color = guide_legend(override.aes = list(size = 0.6))))

#plot all samples of 1 strain
pd2 <- filter(pd, Strain=="Agr") %>%
  pivot_longer(cols=-c("Strain", "Sample"), values_to="Intensity", names_to="mz")
class(pd2$mz) <- "numeric"
pd2$Sample <- paste("Agr", pd2$Sample, "   ")

(p3 <- ggplot(pd2, aes(x=mz, y=Intensity, color=Sample)) +
    geom_line(lwd=0.2) +
    labs(x="m/z") +
    scale_color_manual(values = magma(12)[2:10]) +
    theme_minimal() +
    guides(color = guide_legend(override.aes = list(size = 0.6))))

png("Spectra.png", width=600, height=800)
grid.arrange(p1, p2, p3, ncol=1, heights = c(4, 2.4, 2.4))
dev.off()

counts <- as.matrix(counts)
idx <- as.character(pd$Strain)
save(counts, mz_locations, idx, file="Spectra.Rdata")

