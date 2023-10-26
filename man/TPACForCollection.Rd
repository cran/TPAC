\name{tpacForCollection}
\alias{tpacForCollection}
\title{TPAC method for multiple gene sets}
\description{
Executes the TPAC (tissue-adjusted pathway analysis for cancer) method (\code{\link{tpac}}) on multiple gene sets, i.e., a gene set collection.
}
\usage{
    tpacForCollection(gene.expr, mean.expr, tissue.specificity, gene.set.collection)
}
\arguments{
\item{gene.expr}{See description in \code{\link{tpac}}} 
\item{mean.expr}{See description in \code{\link{tpac}}} 
\item{tissue.specificity}{See description in \code{\link{tpac}}} 
\item{gene.set.collection}{List of m gene sets for which scores are computed.
    Each element in the list corresponds to a gene set and the list element is a vector
    of indices for the genes in the set. The index value is defined relative to the
    order of genes in the \code{gene.expr} matrix. Gene set names should be specified as list names.}
}
\value{
A list containing two elements:
    \itemize{
    \item \code{s.pos}: n x m matrix of TPAC scores computed using the positive squared adjusted Mahalanobis distances.
    \item \code{s.neg}: n x m matrix of TPAC scores computed using the negative squared adjusted Mahalanobis distances.
    \item \code{s}: n x m matrix of TPAC scores computed using the sum of the positive and negative squared adjusted Mahalanobis distances.
     }
}
\seealso{\code{\link{tpac}},\code{\link{tpacForCancer}}} 
\keyword{file}
\examples{
    # Simulate Gaussian expression data for 10 genes and 10 samples
    gene.expr=matrix(rnorm(100), nrow=10)
    # Use 0 as mean.expr 
    mean.expr=rep(0,10)
    # Simulate tissue-specific weights
    tissue.specificity = runif(10, min=0.5, max=1.5)
    # Define a collection with two disjoint sets that span the 10 genes
    collection=list(set1=1:5, set2=6:10)    
    # Execute TPAC on both sets
    tpacForCollection(gene.expr=gene.expr, mean.expr=mean.expr,
        tissue.specificity=tissue.specificity, gene.set.collection=collection)
}