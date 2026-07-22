# Purpose: Estimate ancestral recombination graph statistics while accounting for non-monophyletic relationships
# Author: Keaka Farleigh, Ph.D.
# Email: keakafarleigh@gmail.com
# Date: July 21st, 2026


### Load your packages
library(phytools)
library(ape)
library(tidyverse)
library(phangorn)
library(parallel)


### Read in data
arg.bed <- read.delim("./pyrrhus_b1_con_argweaver.bed", header = FALSE)

arg.bed <- arg.bed[,c(1,2,3,5)]

colnames(arg.bed) <- c("chromosome", "start", "end", "tree")


### Set populations 
pops <- read.delim("./b1_con.popmap",header = FALSE)

pops$species <- "continental"


# Set the individuals from baja1
pops$species[c(12:15,17)] <- "baja1"

colnames(pops)[1] <- "sample"

pops_hap1 <- pops
pops_hap2 <- pops

# Paste 1 and 2 on each individual to account for different haplotypes
pops_hap1$sample <- paste(pops_hap1$sample, "_1", sep = "")
pops_hap2$sample <- paste(pops_hap2$sample, "_2", sep = "")

pops_final <- rbind(pops_hap1,pops_hap2)

## Separate into lists of individuals from the two popuations
con <- pops_final[which(pops_final$species == "continental"),1]
b1 <- pops_final[which(pops_final$species == "baja1"),1]


source("./sliding_argstats.R")

### Set the number of cores we want to use 
# Split the data into a list so that we can parallelize
arg.list <- split(arg.bed, seq_len(nrow(arg.bed)))

pyr_b1_con_slidingarg_stats <- mclapply(arg.list, sliding_argstats, pop1 = con, pop2 = b1, pop1.name = "continental", pop2.name = "baja1", mc.cores = 6, mc.silent = TRUE)


save.image("pyrrhus_b1_con_slidingwindow_wholegenome_argstats.Rdata")
