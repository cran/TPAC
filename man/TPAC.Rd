\name{tpac}
\alias{tpac}
\title{Tissue-adjusted pathway analysis of cancer (TPAC) algorithm}
\description{
Implementation of the tissue-adjusted pathway analysis of cancer (TPAC) method, which performs single sample gene set scoring of
cancer transcriptomic data. The TPAC method computes a tissue-adjusted squared Mahalanobis distance and one-sided CDF values for
all samples in the specified cancer gene expression matrix. Importantly, these distances are measured from the
normal tissue mean and normal tissue specificity (i.e., fold-change in mean expression in the corresponding normal tissue relative to
the mean expression in other tissue types) is used to weight genes when computing the distance metric. The distances are decomposed into
a positive component, negative component and overall component, i.e., the positive component only captures expression values above the mean.
A gamma distribution is fit on the squared tissue-adjusted Mahalanobis distances computed on data that is permuted to match to the null 
that the cancer expression data is uncorrelated with mean values matching those found in the associated normal tissue. This null
gamma distribution is used to calculate the CDF values for the unpermuted distances.}
\usage{
    tpac(gene.expr, mean.expr, tissue.specificity)
}
\arguments{
\item{gene.expr}{An n x p matrix of gene expression values for n samples and p genes. This matrix should reflect the subset of the 
  full expression profile that corresponds to a single gene set.} 
\item{mean.expr}{Vector of mean expression values for the p genes from which distances will be measured, i.e., mean expression in 
  corresponding normal tissue.These mean expression values must be computed using the same normalization method used for 
  tumor expression data in \code{gene.expr}}
\item{tissue.specificity}{Optional vector of tissue-specificity values for the p genes, 
  i.e, fold change between mean expression in the associated normal tissue and the average of 
  mean expression values in other cancer-associated normal tissues. If specified, these specificity values are used to modify the 
  diagonal elements of the covariance matrix used in the Mahalanobis statistic.}
}
\value{
A \code{data.frame} with the following elements (row names will match row names from \code{gene.expr}):
    \itemize{
    \item \code{s.pos}: vector of TPAC scores (i.e., CDF values) computed from the positive squared adjusted Mahalanobis distances.
    \item \code{s.neg}: vector of TPAC scores computed from the negative squared adjusted Mahalanobis distances.
    \item \code{s}: vector of TPAC scores computed from the sum of the positive and negative squared adjusted Mahalanobis distances.
    }
}
\seealso{\code{\link{tpacForCollection}}}
%\code{\link{tpacForTCGA}}} 
\keyword{file}
\examples{
    # Simulate Gaussian expression data for 10 genes and 10 samples
    gene.expr=matrix(rnorm(100), nrow=10)
    # Use 0 as mean.expr 
    mean.expr=rep(0,10)
    # Simulate tissue-specific weights
    tissue.specificity = runif(10, min=0.5, max=1.5)
    # Execute TPAC to compute positive, negative and overall scores 
    tpac(gene.expr=gene.expr, mean.expr=mean.expr, tissue.specificity=tissue.specificity)
}
