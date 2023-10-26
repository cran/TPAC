\name{createGeneSetCollection}
\alias{createGeneSetCollection}
\title{Utility function to help create gene set collection list object}
\description{
Utility function that creates a gene set collection list in the format required by \code{\link{tpacForCollection}} given the gene IDs measured in the expression matrix and a list of gene sets as defined by the IDs of the member genes.
}
\usage{
    createGeneSetCollection(gene.ids, gene.set.collection, min.size=1, max.size)
}
\arguments{
\item{gene.ids}{Vector of gene IDs. This should correspond to the genes measured in the gene expression data.}
\item{gene.set.collection}{List of gene sets where each element in the list corresponds to a gene set and the list element is a vector of gene IDs. List names are gene set names. Must contain at least one gene set.}
\item{min.size}{Minimum gene set size after filtering out genes not in the gene.ids vector. Gene sets whose post-filtering size is below this are removed from the final collection list. Default is 1 and cannot be set to less than 1.}
\item{max.size}{Maximum gene set size after filtering out genes not in the gene.ids vector. Gene sets whose post-filtering size is above this are removed from the final collection list. If not specified, no filtering is performed.}
}
\value{
Version of the input gene.set.collection list where gene IDs have been replaced by position indices, genes not present in the gene.ids vector have been removed and gene sets failing the min/max size constraints have been removed.
}
\seealso{\code{\link{tpacForCollection}}}
\keyword{file}
\examples{
    # Create a collection with two sets defined over 3 genes
    createGeneSetCollection(gene.ids=c("A", "B", "C"),
        gene.set.collection = list(set1=c("A", "B"), set2=c("B", "C")),
        min.size=2, max.size=3)
}
