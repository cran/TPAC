#
# TPAC.R
#
# Implementation of the Tissue-adjusted Pathway Analysis of Cancer (TPAC) method which
# performs single sample pathway analysis of cancer gene expression data while adjusting for
# gene activity in the associated normal tissue.
#
# @author rob.frost@dartmouth.edu
#

library(MASS)

#
# Computes a tissue-adjusted Mahalanobis distance and one-sided CDF values for
# all tumor samples in the specified gene expression matrix. This matrix should reflect the subset of the full 
# expression profile that corresponds to a single gene set. A gamma distribution is 
# fit on the squared tissue-adjusted Mahalanobis distances computed on data that is permuted to
# match to the null that the cancer expression data is uncorrelated with mean values matching those
# found in the associated normal tissue. 
#
# Inputs:
#
# -gene.expr: An n x p matrix of gene expression values for n samples and p genes. 
# -mean.expr: Vector of mean expression values for the p genes from which distances will be measured
#          These mean expression values must be computed using the same normalization method used for 
#          tumor expression data in gene.expr.
# -tissue.specificity: Optional vector of tissue-specificity values for the p genes, 
#          i.e, fold change between mean expression in the associated tissue and the average of 
#          mean expression values in other cancer-associated normal tissues. If specified, these
#          are used to modify the diagonal elements of the covariance matrix used in the Mahalanobis statistic.
#
# Output:
#  
#   A data.frame with the following elements (row names will match row names from gene.expr):
#     -s.pos: vector of TPAC scores (i.e., CDF values) computed from the positive squared adjusted Mahalanobis distances.
#     -s.neg: vector of TPAC scores computed from the negative squared adjusted Mahalanobis distances.
#     -s: vector of TPAC scores computed from the sum of the positive and negative squared adjusted Mahalanobis distances.
#
tpac = function(gene.expr, mean.expr, tissue.specificity) {
    
  if (missing(gene.expr)) {
    stop("Missing gene expression matrix!")
  }    
  
  if (missing(mean.expr)) {
    stop("Missing mean gene expression vector!")
  }    
  
  n = nrow(gene.expr)
  p = ncol(gene.expr)      
     
  #----------------------------------------------------------------------------------------------------  
  # Compute sample variance
  #----------------------------------------------------------------------------------------------------    
  gene.var = apply(gene.expr, 2, var)

  #----------------------------------------------------------------------------------------------------  
  # Check for NA variances
  #----------------------------------------------------------------------------------------------------    
  if (any(is.na(gene.var))) {
    stop("Computed ", length(which(is.na(gene.var))), " NA gene expression variances!")
  }  
  
  #----------------------------------------------------------------------------------------------------  
  # Check for 0 variances
  #----------------------------------------------------------------------------------------------------    
  genes.with.zero.var = which(gene.var == 0)
  if (length(genes.with.zero.var) > 0) {
    warning("Removing ", length(genes.with.zero.var), " genes with 0 variance")
    genes.to.keep = which(gene.var > 0)
    mean.expr = mean.expr[genes.to.keep]
    gene.expr = gene.expr[,genes.to.keep]
    gene.var = gene.var[genes.to.keep]    
    if (!missing(tissue.specificity)) {
      tissue.specificity = tissue.specificity[genes.to.keep]
    }
  }    
  num.genes = length(gene.var)
  
  #----------------------------------------------------------------------------------------------------
  # If tissue-specific weights are not used, the identity matrix is used for the sample covariance, i.e.,
  # assume no covariance between genes and unit variance. If tissue-specific weights are specified,
  # the diagonal elements for positive distances are set to the reciprocal of the tissue-specific weight, 
  # i.e., variation for genes with normal tissue-specific expression is down weighted, while the
  # diagonal elements for negative distances are set to the tissue-specific weights, i.e.,
  # variation for genes with normal tissue-specific expression is up weighted.
  #----------------------------------------------------------------------------------------------------
  
  pos.weights = rep(1, num.genes)          
  neg.weights = rep(1, num.genes)          
  if (!missing(tissue.specificity)) {
    neg.weights = tissue.specificity  
    pos.weights = 1/tissue.specificity      
  }   
    
  #----------------------------------------------------------------------------------------------------
  # Compute the squared Mahalanobis distance. Use custom logic
  # rather than standard R mahalanobis() to support distances from origin and custom
  # covariance matrix.
  #----------------------------------------------------------------------------------------------------
  mahalanobis.sq = computeMahalanobis(X=gene.expr, mean.vals=mean.expr, gene.var=gene.var,
      pos.weights=pos.weights, neg.weights=neg.weights)
  combined.mahalanobis.sq = unlist(mahalanobis.sq)
  
  # Check if any of the distances are infinite
  if (any(is.infinite(combined.mahalanobis.sq))) {
    stop("Computed ", length(which(is.infinite(combined.mahalanobis.sq))), " infinite distances")
  }  

  # Check if any of the distances are NaN
  if (any(is.nan(combined.mahalanobis.sq))) {
    stop("Computed ", length(which(is.nan(combined.mahalanobis.sq))), " NaN distances")
  }      
  
  # Check if any of the distances are NA
  if (any(is.na(combined.mahalanobis.sq))) {
    stop("Computed ", length(which(is.na(combined.mahalanobis.sq))), " NA distances")
  }  
  
  #----------------------------------------------------------------------------------------------------
  # Compute one-sided p-values using the squared Mahalanobis distances using a
  # gamma null distribution that is fit on permuted data.
  #----------------------------------------------------------------------------------------------------
  
  # To determine the null distribution, compute adjusted Mahalanobis distances
  # with sample labels permuted for each gene. This will break any correlation between
  # genes and should remove outlier samples.
  gene.expr.null = apply(gene.expr, 2, function(c) {
        return (sample(c, length(c), replace=F))
      })
  
  # Compute the squared distances on the permuted data 
  mahalanobis.sq.null = computeMahalanobis(X=gene.expr.null, mean.vals=mean.expr, gene.var=gene.var, 
      pos.weights=pos.weights, neg.weights=neg.weights)

  # Compute CDF values according to gamma distribution estimated on the non-zero null distances. Do this
  # separately for the negative distances, positive distances and sum of the absolute value of
  # the negative and positive distances
  pos.cdf.values = computeCDFValues(mahalanobis.sq=mahalanobis.sq$pos.dist.sq,
    mahalanobis.sq.null=mahalanobis.sq.null$pos.dist.sq)$cdf.values         
  neg.cdf.values = computeCDFValues(mahalanobis.sq=mahalanobis.sq$neg.dist.sq,
    mahalanobis.sq.null=mahalanobis.sq.null$neg.dist.sq)$cdf.values         
  total.cdf.values = computeCDFValues(mahalanobis.sq=mahalanobis.sq$total.dist.sq,
    mahalanobis.sq.null=mahalanobis.sq.null$total.dist.sq)$cdf.values         
  
  # Create a data.frame to hold the results
  results = data.frame(s.pos = pos.cdf.values,
                       s.neg = neg.cdf.values,
                       s = total.cdf.values)
                   
  rownames(results) = rownames(gene.expr)                   
  
  return (results)
}

computeMahalanobis = function(X, mean.vals, gene.var, neg.weights, pos.weights) {
  n = nrow(X)
  results = data.frame(pos.dist.sq = rep(0, n),
      neg.dist.sq = rep(0, n),
      total.dist.sq = rep(0, n))
  if (missing(neg.weights)) {
    neg.weights = 1
  }
  if (missing(pos.weights)) {
    pos.weights = 1
  }  
  for (i in 1:n) {
    r = X[i,]
    diffs = r-mean.vals
    neg.diffs = diffs
    neg.diffs[which(diffs > 0)] = 0
    results$neg.dist.sq[i] = sum(neg.weights * neg.diffs^2/gene.var)
    pos.diffs = diffs
    pos.diffs[which(diffs < 0)] = 0
    results$pos.dist.sq[i] = sum(pos.weights * pos.diffs^2/gene.var)

  }
  results$total.dist.sq = results$neg.dist.sq + results$pos.dist.sq
  
  return (results)
}

computeCDFValues = function(mahalanobis.sq, mahalanobis.sq.null) {
  
  cdf.values = rep(0, length(mahalanobis.sq))
  gamma.fit = NA  
  
  # Remove any 0 null squared distances since fitdistr will fail in this case.
  # If there are fewer than 2 non-zero entries, silently default p-values to 1.
  non.zero.entries = which(mahalanobis.sq.null !=0)
  if (length(non.zero.entries) >= 2) {
    # Fit a gamma distribution to the non-zero squared distances from the row permuated data
    # using maximum likelihood estimation as implemented by fitdistr() in the MASS package. 
    # Let fitdistr pick initial values (silence NaN warnings).  
    mahalanobis.sq.null.nonzero = mahalanobis.sq.null[non.zero.entries]
    gamma.fit = try(fitdistr(mahalanobis.sq.null.nonzero, "gamma", lower=0.01))
    if (inherits(gamma.fit, "try-error")) {
      warning("Estimation of gamma distribution failed, defaulting p-values to 1")
      gamma.fit=NA
    } else {
      # compute the one-sided p-values
      cdf.values = 1-pgamma(mahalanobis.sq, 
                            shape=gamma.fit$estimate[1], rate=gamma.fit$estimate[2], lower.tail=F)      
    }
  } 
  
  results = list()
  results$cdf.values = cdf.values
  results$gamma.fit = gamma.fit
  
  return (results)
}

#
# Calls the tpac() method for multiple gene sets.
# 
# Inputs:
#
# -gene.expr: A n x p matrix of gene expression values for n tumor samples and p genes.
# -mean.expr: See description in tpac(). 
# -tissue.specificity: See description in tpac(). 
# -gene.set.collection: List of m gene sets for which scores are computed.
#              Each element in the list corresponds to a gene set and the list element is a vector
#              of indices for the genes in the set. The index value is defined relative to the
#              order of genes in the gene.exprs matrix. Gene set names should be specified as list names.
#              See createGeneSetCollection() for utility function that can be used to 
#              help generate this list of indices.
#
# Output:
#  
#   A list containing three elements:
#     -S.pos: n x m matrix of TPAC scores computed using the positive distances
#     -S.neg: n x m matrix of TPAC scores computed using the negative distances
#     -S: n x m matrix of TPAC scores computed using both positive and negative distances
#
tpacForCollection = function(gene.expr, mean.expr, tissue.specificity, gene.set.collection) {
    
  if (missing(gene.expr)) {
    stop("Missing gene expression matrix!")
  }
  if (missing(mean.expr)) {
    stop("Missing mean expression!")
  }
  if (missing(gene.set.collection)) {
    stop("Missing gene set collection list!")
  }  
  
  p = ncol(gene.expr)
  n = nrow(gene.expr)
  
  sample.ids = rownames(gene.expr)
  num.sets = length(gene.set.collection)
  set.names = names(gene.set.collection)
  gene.ids = colnames(gene.expr)
  set.sizes = unlist(lapply(gene.set.collection, length))
  min.set.size = min(set.sizes)
  median.set.size = median(set.sizes)  
    
  # Prepare the result matrices
  results = list()
  results$S.pos = matrix(0, nrow=n, ncol=num.sets,
      dimnames=list(sample.ids, set.names))
  results$S.neg = matrix(0, nrow=n, ncol=num.sets,
      dimnames=list(sample.ids, set.names))  
  results$S = matrix(0, nrow=n, ncol=num.sets,
      dimnames=list(sample.ids, set.names))  
  
  message("Computing TPAC distances for ", num.sets, " gene sets, ", n, " samples and ", p, " genes.")
  message("Min set size: ", min.set.size, ", median size: ", median.set.size)

  # Process all gene sets in the collection
  for (i in 1:num.sets) {
    set.members = gene.set.collection[[i]]
    set.size = set.sizes[i]
    set.name = set.names[i]
    if (i %% 50 == 0) {
      message("Computing for gene set ", set.name, " of size ", set.size)
      #message("Set members: ", paste0(set.members, collapse=","))
    }
    set.gene.expr = gene.expr[,set.members]   
    set.mean.expr = mean.expr[set.members]
    if (set.size == 1) {
      # Force vector to matrix
      warning("Gene set ", i, " has just a single member!")
      set.exprs = as.matrix(set.exprs)
    }
    if (!missing(tissue.specificity)) {
      set.tissue.specificity = tissue.specificity[set.members]
      tpac.results = tpac(gene.expr=set.gene.expr, 
            mean.expr = set.mean.expr, 
            tissue.specificity = set.tissue.specificity)     
    } else {
      tpac.results = tpac(gene.expr=set.gene.expr, 
          mean.expr = set.mean.expr)
    }
    results$S.pos[,i] = tpac.results$s.pos
    results$S.neg[,i] = tpac.results$s.neg
    results$S[,i] = tpac.results$s
  }    

  return (results)
}

#
# Utility function that creates a gene set collection list in the format required
# by tpacForCollection() given the gene IDs measured in the expression matrix and a 
# list of gene sets as defined by the IDs of the member genes.
#
# Inputs:
#
# -gene.ids: Vector of gene IDs. This should correspond to the genes measured in the 
#            gene expression data.
# -gene.set.collection: List of m gene sets where each element in the list corresponds to
#            a gene set and the list element is a vector gene IDs. List names are gene set names.
# -min.size: Minimum gene set size after filtering out genes not in the gene.ids vector. 
#            Gene sets whose post-filtering size is below this are removed from the final
#            collection list. Default is 1 and cannot be set to less than 1.
# -max.size: Maximum gene set size after filtering out genes not in the gene.ids vector. 
#            Gene sets whose post-filtering size is above this are removed from the final
#            collection list. If not specified, no filtering is performed.
#
# Output:
#  
#   Version of the input gene.set.collection list where gene IDs have been replaced by position indices,
#   genes not present in the gene.ids vector have been removed and gene sets failing the min/max size
#   constraints have been removed.
#
createGeneSetCollection = function(gene.ids, gene.set.collection, min.size=1, max.size) {
  
  # min.size must be at least 1
  if (min.size < 1) {
     stop("Invalid min.size value! Must be 1 or greater.")
  }
  # If max size is set, make sure it is not less than min size
  if (!missing(max.size)) {
    if (max.size < min.size) {
      stop("max.size cannot be less than min.size!")
    }          
  }    
  
  num.genes = length(gene.ids)
  if (num.genes < 1) {
    stop("gene.ids must contain at least one genes!")
  }
  
  num.gene.sets = length(gene.set.collection)   
  if (num.gene.sets < 1) {
    stop("gene.set.collection must contain at least one gene set!")
  }  
  
  num.sets = length(gene.set.collection)
  set.names = names(gene.set.collection)
  gene.set.indices = list()
  for (i in 1:num.sets) {
    set.ids = gene.set.collection[[i]]
    # map IDs to indices
    set.indices = unlist(sapply(set.ids, function(x){which(gene.ids == x)}))
    set.size = length(set.indices)
    if (set.size < min.size) {
      next
    }
    if (!missing(max.size)) {
        if (set.size > max.size) {
          next
        }
    }
    current.index = length(gene.set.indices)+1
    gene.set.indices[[current.index]] = set.indices
    names(gene.set.indices)[current.index] = set.names[i]
  }
    
  return (gene.set.indices)
}



