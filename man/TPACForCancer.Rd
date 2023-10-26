\name{tpacForCancer}
\alias{tpacForCancer}
\title{Executes the TPAC method on cancer gene expression data}
\description{
Executes the TPAC (tissue-adjusted pathway analysis for cancer) method (\code{\link{tpacForCollection}}) on cancer gene expression data
using normal tissue expression data from the Human Protein Atlas (HPA) that is included in the package as \code{\link{hpa.data}}.
This HPA normal tissue data was specially processed by the HPA group as FPKM values using a pipeline similar to that employed by GDC for the TCGA data. For consistency with this HPA normal tissue data, the provided \code{cancer.gene.expr} data must be specified as FPKM+1 values. Please see the vignette for an example of calling this function using appropriately normalized TCGA gene expression data.
}
\usage{
    tpacForCancer(cancer.gene.expr, cancer.type, gene.set.collection, 
                  min.set.size=1, max.set.size)
}
\arguments{
\item{cancer.gene.expr}{An n x p matrix of gene expression values for n tumors of the specified tumor type and p genes. 
                        The data should be normalized as FPKM+1 values, row names should be sample ID, and column names should be Ensembl gene IDs.}
\item{cancer.type}{Cancer type of the expression data. Must be one of the supported cancer types as per \code{\link{getSupportedCancerTypes}}.}
\item{gene.set.collection}{List of m gene sets for which scores are computed.
    Each element in the list corresponds to a gene set and the list element is a vector
    of Ensembl IDs for genes in the set. Gene set names should be specified as list names.}
\item{min.set.size}{See description of \code{min.size} in \code{\link{createGeneSetCollection}}}
\item{max.set.size}{See description of \code{max.size} in \code{\link{createGeneSetCollection}}}
}
\value{
A list containing two elements:
    \itemize{
    \item \code{S.pos}: n x m matrix of TPAC scores computed using the positive squared adjusted Mahalanobis distances.
    \item \code{S.neg}: n x m matrix of TPAC scores computed using the negative squared adjusted Mahalanobis distances.
    \item \code{S}: n x m matrix of TPAC scorescomputed using the sum of the positive and negative squared adjusted Mahalanobis distances.
     }
}
\seealso{\code{\link{tpac}}, \code{\link{hpa.data}}, \code{\link{tpacForCollection}}, \code{\link{getSupportedCancerTypes}}, \code{\link{createGeneSetCollection}}}
\keyword{file}
\examples{
    # Simulate Gaussian expression data for 10 genes and 10 samples
    # (Note: cancer expression should be FPKM+1 for real applications)
    cancer.gene.expr=matrix(rnorm(200), nrow=20)
    # Create arbitrary Ensembl IDs
    gene.ids = c("ENSG00000000003","ENSG00000000005","ENSG00000000419",
                 "ENSG00000000457","ENSG00000000460","ENSG00000000938",
                 "ENSG00000000971","ENSG00000001036","ENSG00000001084",
                 "ENSG00000001167")
    colnames(cancer.gene.expr) = gene.ids
    # Define a collection with two disjoint sets that span the 10 genes
    collection=list(set1=gene.ids[1:5], set2=gene.ids[6:10])    
    # Execute TPAC on both sets
    tpacForCancer(cancer.gene.expr=cancer.gene.expr, cancer.type="glioma", 
                  gene.set.collection=collection)
}