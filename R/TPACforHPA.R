#
# TPACforHPA.R
#
# Wrapper around the Tissue-adjusted Pathway Analysis of Cancer (TPAC) method 
# that uses Human Protein Atlas (HPA) gene expression data to compute the
# gene-level normal tissue mean and tissue-specificity values.
#
# @author rob.frost@dartmouth.edu
#

#
# Calls the tpac() method for multiple gene sets using normal tissue information from the HPA.
# 
# Inputs:
#
# -cancer.gene.expr: A n x p matrix of gene expression values (FPKM+1 format) for n tumor samples and p genes. Column names must be gene Ensembl IDs.
# -cancer.type: The human cancer type associated with the gene.expr data. Must be one of the
#.      cancer types returned by getSupportedCancerTypes().
# -gene.set.collection: List of m gene sets for which scores are computed.
#              Each element in the list corresponds to a gene set and the list element is a vector
#              of Ensembl IDs. Note that this is the structure used as input to createGeneSetCollection() rather
#              than the output of that method.
# -min.set.size
# -max.set.size
#
# Output:
#  
#   A list containing three elements:
#     -S.pos: n x m matrix of TPAC scores computed using the positive distances
#     -S.neg: n x m matrix of TPAC scores computed using the negative distances
#     -S: n x m matrix of TPAC scores computed using both positive and negative distances
#
tpacForCancer = function(cancer.gene.expr, cancer.type, gene.set.collection, min.set.size=1, max.set.size) {
    
  if (missing(cancer.gene.expr)) {
    stop("Missing gene expression matrix!")
  }
  if (missing(cancer.type)) {
    stop("Missing cancer type!")
  }
  if (missing(gene.set.collection)) {
    stop("Missing gene set collection list!")
  }  
  
  # Retrieve the associated HPA data, find the intersection of genes and return the
  # cancer expression subset, normal tissue mean expression values and normal tissue specificity values
  reconcile.out = reconcileCancerAndHPAData(cancer.expr=cancer.gene.expr, cancer.type=cancer.type)
  
  # Create a gene set collection structure that replaces IDs with indices
  if (missing(max.set.size)) {
    collection.with.indices = createGeneSetCollection(gene.ids = reconcile.out$gene.ids, 
                                                      gene.set.collection=gene.set.collection,
                                                      min.size=min.set.size)
  } else {
    collection.with.indices = createGeneSetCollection(gene.ids = reconcile.out$gene.ids, 
                                                      gene.set.collection=gene.set.collection,
                                                      min.size=min.set.size,
                                                      max.size=max.set.size)
  }

  # Call the underlying TPAC method
  tpac.out = tpacForCollection(gene.expr = reconcile.out$cancer.expr,
                               mean.expr = reconcile.out$normal.mean.expr,
                               tissue.specificity = reconcile.out$tissue.specificity,
                               gene.set.collection = collection.with.indices)
  return (tpac.out)
}

#
# Returns the human cancer types that have associated normal tissue data in the HPA.
#
getSupportedCancerTypes = function() {
  return (getHPACancerToNormalMappings()[,1]) 
}

#
# Retrieves the HPA tissue type for the specified cancer type
#
getHPATissueForCancerType = function(cancer.type) {
  cancer.tissue.map = getHPACancerToNormalMappings()
  cancer.index = which(cancer.tissue.map[,1] == cancer.type)
  if (length(cancer.index) == 0) {
    stop("Cancer type ", cancer.type, " not supported!")
  }
  normal.tissue = cancer.tissue.map[cancer.index,2]
  return (normal.tissue)
}

# 
# Returns a data.frame mapping HPA cancer types to normal tissue types.
#
getHPACancerToNormalMappings = function() {
  mappings = matrix(nrow=18, ncol=2)
  mappings[1,] = c("urothelial cancer", "urinarybladder")
  mappings[2,] = c("breast cancer", "breast")
  mappings[3,] = c("cervical cancer", "cervix")
  mappings[4,] = c("colorectal cancer", "colon")
  mappings[5,] = c("glioma", "brain")
  mappings[6,] = c("head and neck cancer", "tonsil")
  mappings[7,] = c("renal cancer", "kidney")
  mappings[8,] = c("liver cancer", "liver")
  mappings[9,] = c("lung cancer", "lung")
  mappings[10,] = c("ovarian cancer", "ovary")
  mappings[11,] = c("pancreatic cancer", "pancreas")
  mappings[12,] = c("prostate cancer", "prostate")
  mappings[13,] = c("colorectal cancer", "rectum")
  mappings[14,] = c("melanoma", "skin")
  mappings[15,] = c("stomach cancer", "stomach")
  mappings[16,] = c("testis cancer", "testis")
  mappings[17,] = c("thyroid cancer", "thyroid")
  mappings[18,] = c("endometrial cancer", "endometrium")
  mappings = as.data.frame(mappings, stringsAsFactors = F)
  colnames(mappings) = c("cancer type", "normal tissue")
  return (mappings)
}

# 
# Retrieve the HPA data for the associated normal tissue, identify the common genes and 
# subset and reorder the associated data structures.
# 
reconcileCancerAndHPAData = function(cancer.expr, cancer.type) {
  tissue.type = getHPATissueForCancerType(cancer.type)
  # cosmetic use of is.data.table to avoid an empty importFrom("data.table") in NAMESPACE
  is.data.table(TPACData::hpa.data)
  hpa.for.tissue = TPACData::hpa.data[list(tissue.type)]
  #hpa.for.tissue = TPACData::hpa.data[.(tissue.type)]
  normal.mean.expr = hpa.for.tissue$FPKM
  tissue.specificity = hpa.for.tissue$TissueSpecificity
  hpa.gene.ids = hpa.for.tissue$Gene
  cancer.gene.ids = colnames(cancer.expr)
  intersection.indices = getIntersectionIndices(cancer.gene.ids, hpa.gene.ids)
  results = list()
  results$cancer.expr = cancer.expr[,intersection.indices$x.indices]
  results$normal.mean.expr = normal.mean.expr[intersection.indices$y.indices]
  results$tissue.specificity = tissue.specificity[intersection.indices$y.indices]
  results$gene.ids = colnames(results$cancer.expr)  
  return (results)
}

#
# Utility function that finds the common elements in the two lists and returns
# indices for a reordering of each original list according to these common
# elements.
#
getIntersectionIndices = function(x, y) {
  common.elements = intersect(x,y)
  intersect.x.indices = unlist(sapply(common.elements, function(z) {which(x == z)}))
  intersect.y.indices = unlist(sapply(common.elements, function(z) {which(y == z)}))
  return (list(x.indices=intersect.x.indices, y.indices=intersect.y.indices))
}

