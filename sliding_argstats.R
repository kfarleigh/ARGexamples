require(phytools)
require(ape)
require(tidyverse)
require(phangorn)
require(parallel)

### Write a custom function in case species are not monophyletic
# dat is a list element
# pop1 is a vector of individuals in one population
# pop2 is a vector of individuals in the other populations
# tip_idx is a vector of tip labels, this is used internally in the sliding_argstats function, do not modify
identify_clade <- function(dat, pop1, pop2, tip_idx){
  
  prop.df <- data.frame(node = NA, prop.pop1 = NA, prop.pop2 = NA, n.ind = NA, n.pop1 = NA, n.pop2 = NA)
  
  pop1.inds <- pop1
  pop2.inds <- pop2
  
  dat <- dat
  
  n.pop1 <- length(which(tip_idx[dat] %in% pop1.inds))
  n.pop2 <- length(which(tip_idx[dat] %in% pop2.inds))
  
  prop.pop1 <- n.pop1/length(pop1.inds)
  prop.pop2 <- n.pop2/length(pop2.inds)
  
  n.ind <- length(dat)
  
  prop.df[1,2] <- prop.pop1
  prop.df[1,3] <- prop.pop2
  prop.df[1,4] <- n.ind
  prop.df[1,5] <- n.pop1
  prop.df[1,6] <- n.pop2
  
  remove(pop1.inds, pop2.inds, n.pop1, n.pop2, prop.pop1, prop.pop2, n.ind, dat)
  
  return(prop.df)
  
  
}

### Sliding argstats arguments
# arg.dat is a list element that contains a dataframe. The data frame columns should be chromosome, start, end, and the phylogeney estimated in ARG analysis
# pop1 is a vector of individuals in one population
# pop2 is a vector of individuals in the other populations
# pop1.name is a character string, that tells us the name of population 1
# pop2.name is a character string, that tells us the name of population 2

sliding_argstats <- function(arg.dat, pop1, pop2, pop1.name, pop2.name){
  
  ### Create a data frame to store results 
  arg.stats.df <- data.frame(chromosome = arg.dat$chromosome, 
                             start = arg.dat$start, end = arg.dat$end,
                             tmrca = NA,
                             pop1.meantmrcaw = NA,
                             pop1.mediantmrcaw = NA,
                             pop1.mintmrcaw = NA,
                             pop1.maxtmrcaw = NA,
                             pop1.mono = NA,
                             pop2.meantmrcaw = NA,
                             pop2.mediantmrcaw = NA,
                             pop2.mintmrcaw = NA,
                             pop2.maxtmrcaw = NA,
                             pop2.mono = NA,
                             pop1 = pop1.name,
                             pop2 = pop2.name)
    
    i <- 1
    
    nwk <- arg.dat$tree[i]
    
    tree <- read.tree(text = nwk)
    
    # Get the TMRCA for the tree 
    tmrca <- max(nodeHeights(tree))
    
    # Get cross-coalescent events for each population/species 
    pop1_dist <- findMRCA(tree, tips = pop1, type = "height")
    
    pop2_dist <- findMRCA(tree, tips = pop2, type = "height")
    
    arg.stats.df[i,4] <- tmrca
    
    if(is.monophyletic(tree, tips= pop1)){
      
      
      # This is height from the root, which gives us the opposite of what we want; testing showed this to be equivalent to extracting a clade and working up from there
      pop1.dist.fix <- tmrca - pop1_dist
      
      arg.stats.df[i,5] <- pop1.dist.fix
      arg.stats.df[i,6] <- pop1.dist.fix
      arg.stats.df[i,7] <- pop1.dist.fix
      arg.stats.df[i,8] <- pop1.dist.fix
      arg.stats.df[i,9] <- TRUE
      
    } else {
      
      # Get a list of nodes and their descendants
      
      all_nodes <- Descendants(tree, type = "tips")
      tip_idx <- tree$tip.label
      
      test <- lapply(all_nodes, identify_clade, pop1 = pop1, pop2 = pop2, tip_idx = tip_idx)
      
      test.df <- do.call("rbind", test)
      
      test.df$node <- paste(1:length(all_nodes))
      
      # Remove nodes with only 1 individual and filter for only nodes where there are no contiental individuals
      pop1.df.filt <- test.df %>% filter(n.pop1 > 1, n.pop2 == 0)
      
      if(nrow(pop1.df.filt) > 0){
        
        # Determine if nodes have a parent/child relationship
        prop.child <- data.frame(node = NA, children = NA, is.child = NA)
        for(j in 1:nrow(pop1.df.filt)){
          
          nodes <- pop1.df.filt[,1]
          
          which(nodes %in% getDescendants(tree, node = pop1.df.filt$node[j],))
          
          prop.child[j,1] <- pop1.df.filt$node[j]
          prop.child[j,2] <- length(which(nodes %in% getDescendants(tree, node = pop1.df.filt$node[j],)))
          
          if(any(Ancestors(tree, pop1.df.filt$node[j], type = "all") %in% nodes)){
            
            prop.child[j,3] <- TRUE
            
          } else{
            
            prop.child[j,3] <- FALSE
            
          }
          
        }
        
        # Only select nodes where is.child == FALSE for calculations
        prop.child.filt <- prop.child %>% filter(is.child == FALSE)
        
        # Calculate tmrca-within
        tmrca_w <- c()
        for(k in 1:nrow(prop.child.filt)){
          
          tips <- getDescendants(tree, node = prop.child.filt$node[k])
          
          pop1_dist <- findMRCA(tree, tips = tips, type = "height")
          pop1.fix <- tmrca - pop1_dist
          
          tmrca_w <- c(tmrca_w, pop1.fix)
          
          remove(tips, pop1_dist, pop1.fix)
          
        }
        
        arg.stats.df[i,5] <- mean(tmrca_w)
        arg.stats.df[i,6] <- median(tmrca_w)
        arg.stats.df[i,7] <- min(tmrca_w)
        arg.stats.df[i,8] <- max(tmrca_w)
        arg.stats.df[i,9] <- FALSE
        
        remove(tmrca_w)
        
      } else{
        
        arg.stats.df[i,5] <- NA
        arg.stats.df[i,6] <- NA
        arg.stats.df[i,7] <- NA
        arg.stats.df[i,8] <- NA
        arg.stats.df[i,9] <- FALSE
        
      }
      
    }
    
    if(is.monophyletic(tree, tips = pop2)){
      
      pop2.fix <- tmrca - pop2_dist
      arg.stats.df[i,10] <- pop2.fix
      arg.stats.df[i,11] <- pop2.fix
      arg.stats.df[i,12] <- pop2.fix
      arg.stats.df[i,13] <- pop2.fix
      arg.stats.df[i,14] <- TRUE
      
    } else {
      
      # Get a list of nodes and their descendants
      
      all_nodes <- Descendants(tree, type = "tips")
      tip_idx <- tree$tip.label
      
      test <- lapply(all_nodes, identify_clade, pop1 = pop1, pop2 = pop2, tip_idx = tip_idx)
      
      test.df <- do.call("rbind", test)
      
      test.df$node <- paste(1:length(all_nodes))
      
      # Remove nodes with only 1 individual and filter for only nodes where there are no contiental individuals
      pop2.df.filt <- test.df %>% filter(n.pop2 > 1, n.pop1 == 0)
      
      if(nrow(pop2.df.filt) > 0){
        
        # Determine if nodes have a parent/child relationship
        prop.child <- data.frame(node = NA, children = NA, is.child = NA)
        for(l in 1:nrow(pop2.df.filt)){
          
          nodes <- pop2.df.filt[,1]
          
          which(nodes %in% getDescendants(tree, node = pop2.df.filt$node[l],))
          
          prop.child[l,1] <- pop2.df.filt$node[l]
          prop.child[l,2] <- length(which(nodes %in% getDescendants(tree, node = pop2.df.filt$node[l],)))
          
          if(any(Ancestors(tree, pop2.df.filt$node[l], type = "all") %in% nodes)){
            
            prop.child[l,3] <- TRUE
            
          } else{
            
            prop.child[l,3] <- FALSE
            
          }
          
        }
        
        # Only select nodes where is.child == FALSE for calculations
        prop.child.filt <- prop.child %>% filter(is.child == FALSE)
        
        # Calculate tmrca-within
        tmrca_w <- c()
        for(m in 1:nrow(prop.child.filt)){
          
          tips <- getDescendants(tree, node = prop.child.filt$node[m])
          
          con_dist <- findMRCA(tree, tips = tips, type = "height")
          pop2.fix <- tmrca - con_dist
          
          tmrca_w <- c(tmrca_w, pop2.fix)
          
          remove(tips, con_dist, pop2.fix)
          
        }
        
        arg.stats.df[i,10] <- mean(tmrca_w)
        arg.stats.df[i,11] <- median(tmrca_w)
        arg.stats.df[i,12] <- min(tmrca_w)
        arg.stats.df[i,13] <- max(tmrca_w)
        arg.stats.df[i,14] <- FALSE
        
        remove(tmrca_w)
        
        
      } else{
        
        arg.stats.df[i,10] <- NA
        arg.stats.df[i,11] <- NA
        arg.stats.df[i,12] <- NA
        arg.stats.df[i,13] <- NA
        arg.stats.df[i,14] <- FALSE
        
        
      }
      
    }
    
    remove(tmrca, tree, nwk)
    
    return(arg.stats.df)
    
}