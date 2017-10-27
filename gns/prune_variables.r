library(MASS)
library(igraph)
if (!require(gnsutils)){
  source("fastCor.R")
  #source("~/Documents/gns-001723/Shared/fastCor.R")
}

compute_vif <- function(mat){ #function for computing vif
  require(MASS)
  vif = diag(ginv(cor(mat)))
  return(vif)
}

compute_dvif <- function(mat, iremove){ #function for computing delta VIF after removing variable(s) iremove
  vif1 = compute_vif(mat)
  vif2 = compute_vif(mat[,-iremove])
  dvif = log(vif1[-iremove]/vif2)
  return(dvif)
}

shuffle_mat <- function(mat, iremove){ #shuffle columns of matrix except for variable iremove
  allcolumns = 1:dim(mat)[2]
  columns_for_shuffle = allcolumns[-iremove]
  shuffled_selection = sapply(columns_for_shuffle, function(col) sample(as.numeric(mat[,col])))
  shuffled_mat = mat
  shuffled_mat[,columns_for_shuffle] = shuffled_selection
  return(shuffled_mat)
}

compute_dvif_thresh <- function(mat, iremove, nshuffle, prctle){ #estimate threshold delta vif for removing iremove
  #by shuffling data (mat) nshuffle times and taking the prctle percentile value of the resulting distribution of delta vifs
  #actually we want nshuffle samples but each shuffle of mat gives us ncol(mat) - 1 values
  #so really we only need to shuffle ceiling(nshuffle/(ncol(mat) - 1)) times
  nreps = ceiling(nshuffle/(ncol(mat) - 1))
  dvif_shuff = c(replicate(nreps, compute_dvif(shuffle_mat(mat,iremove),iremove)))
  thresh = as.numeric(quantile(dvif_shuff,prctle))
  return(thresh)
}

prune_variables_iter <- function(mat, keep, nshuffle = 1000, prctle = 0.95, cor.thresh = 0.9) { #do one variable pruning iteration inside cluster with data mat. 
  #input keep is vector with indicator for each variable in mat. 
  #2: variable is kept. 1: variable has not been dealt with yet. 0: variable has been eliminated.
  #function returns keep2, which sets eliminated variables to 0 and kept variables to 2
  keep2 = keep
  this.mat = mat[,which(keep==1)]
  this.vars = which(keep==1)
  nvar = length(this.vars)
  if (length(this.vars)==2){ #if only two variables then keep one at random
    if (abs(cor(this.mat[,1],this.mat[,2]))>cor.thresh) {
      keep2[this.vars[1]] = 2
      keep2[this.vars[2]] = 0
    } else {
      keep2[this.vars[1]] = 2
      keep2[this.vars[2]] = 2
    }
  } else if(length(this.vars)>2){
    vif = compute_vif(this.mat)
    top_vif = which.max(vif)
    top_vif_abs = this.vars[top_vif]
    variables_considered = this.vars[-top_vif]
    dvif = compute_dvif(this.mat,top_vif)
    dvif_thresh = compute_dvif_thresh(this.mat, top_vif, nshuffle, prctle)
    corr_vars = variables_considered[which(dvif>dvif_thresh)]
    #plot(variables_considered,dvif,xlim=c(1,nvar))
    #abline(h=dvif_thresh)
    keep2[top_vif_abs] = 2
    keep2[corr_vars] = 0
  }
  return(keep2)
}

keep_top <- function(mat) {
  vif = compute_vif(mat)
  top_vif = which.max(vif)
  keep = rep(0,dim(mat)[2])
  keep[top_vif] = 2
  return(keep)
}


prune_cluster <- function(mat, cor.thresh = 0.9, toponly = FALSE, nshuffle = 1000, prctle = 0.95, verbose = TRUE) {
  if (verbose) cat(dim(mat)[2], " variables in cluster; starting pruning.\n")
  if (toponly) {
    if (verbose) cat("Keeping only variable with top VIF.\n")
    keep = keep_top(mat)
  } else {
    keep = rep(0,dim(mat)[2])
    keep2 = rep(1,dim(mat)[2])
    while(!identical(keep,keep2) && sum(keep2==1)>0){ 
      keep = keep2
      keep2 = prune_variables_iter(mat, keep, nshuffle, prctle, cor.thresh)
      print(keep2) # for debugging
      if(verbose) cat("Pruning iteratively, ", sum(keep2>0), " variables remaining in cluster.\n")
    }
    keep = keep2
  }
  if(verbose) cat("Finished pruning, ", sum(keep>0), " variables remaining in cluster.\n")
  if(verbose){
    vifmat = compute_vif(mat)
    plot(vifmat)
    points(which(keep>0), vifmat[which(keep>0)], pch=0, col='red')
  }
  names_to_keep = colnames(mat)[which(keep>0)]
  return(names_to_keep)
}

prune_dataframe <- function(df, cor.thresh = 0.9, toponly = FALSE, nshuffle = 1000, prctle = 0.95, verbose = TRUE){
  cor.table = fastCorTable(df, cor.thresh, 1, nblocks=ceiling(dim(df)[2]/1500), verbose=F)
  if (dim(cor.table)[1] > 0){
    #create igraph graph
    g = graph_from_edgelist(as.matrix(cor.table[,-3]), directed=F)
    # extract the clusters from the graph
    cluster.size = components(g)$csize
    cluster.membership = components(g)$membership
    if(verbose) cat("Building graph with ", length(cluster.membership), " nodes and ", length(cluster.size), "components.\n")
    # loop through each cluster and select representative genes
    tokeep = c()
    for(cluster in 1:length(cluster.size)){
      cat("Pruning cluster ", cluster, " out of ", length(cluster.size), ".\n")
      cluster.genes = names(which(cluster.membership==cluster))
      cluster.mat = df[,colnames(df)%in%cluster.genes]
      cluster.tokeep = prune_cluster(cluster.mat, cor.thresh, toponly, nshuffle, prctle, verbose)
      tokeep = append(tokeep,cluster.tokeep)
    }
    toremove = names(cluster.membership)[!names(cluster.membership)%in%tokeep]
    if (verbose) cat("Removing a total of ", length(toremove), " variables from dataframe.\n")
    pruned.df = df[,!colnames(df)%in%toremove]
  } else {
    cat("No two variables correlated above threshold!\n")
    pruned.df = df
  }

  return(pruned.df)
}

