\name{XBSeq}
\alias{XBSeq}
\title{
   Express function to carry out XBSeq analysis
}
\description{
   A wrapper function to carry out XBSeq analysis procedure
}
\usage{
XBSeq(counts, bgcounts, conditions, method = "pooled", 
   sharingMode = "maximum", fitType = "local", pvals_only = FALSE)
}
\arguments{
  \item{counts}{
   A data.frame or matrix contains the observed signal 
}
  \item{bgcounts}{
   A data.frame or matrix contains the background noise
}
  \item{conditions}{
   A factor to specify the experimental design 
}
  \item{method}{
   Method used to estimate SCV
}
  \item{sharingMode}{
   Mode of sharing of information
}
  \item{fitType}{
   Option to fit mean-SCV relation
}
  \item{pvals_only}{
   Logical; Specify whether to extract pvalues only
}
}
\value{
A data.frame with following columns:
\item{id}{
   rownames of XBSeqDataSet
}
\item{baseMean}{
   The basemean for all genes
}
\item{baseMeanA}{
   The basemean for condition 'A'
}
\item{baseMeanB}{
   The basemean for condition 'B'
}
\item{foldChange}{
   The fold change compare condition 'B' to 'A'
}
\item{log2FoldChange}{
   The log2 fold change
}
\item{pval}{
   The p value for all genes
}
\item{padj}{
   The adjusted p value for all genes
}
}
\author{
   Yuanhang Liu
}
\references{
   Hung-I Harry Chen, Yuanhang Liu, Yi Zou, Zhao Lai, Devanand Sarkar, Yufei Huang, Yidong Chen
   doi: \url{http://dx.doi.org/10.1101/016196}
}
\seealso{
\code{\link{estimateRealCount}},  \code{\link{XBSeqDataSet}}, \code{\link{estimateSCV}}, \code{\link{XBSeqTest}}
}
\examples{
   conditions <- c(rep('C1', 3), rep('C2', 3))
   data(ExampleData)
   Stats <- XBSeq(Observed, Background, conditions)
}