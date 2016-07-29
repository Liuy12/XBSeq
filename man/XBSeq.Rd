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
   sharingMode = "maximum", fitType = "local", pvals_only = FALSE, paraMethod='NP', big_count = 900)
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
  \item{paraMethod}{
   Method to use for estimation of distribution parameters, 'NP' or 'MLE'. See details section for details 
}
  \item{big_count}{
   An integer specify a count number above which where be considerred as 'big' and beta approximation will be used instead for testing differential expression
}
}
\details{
   This is the express function for carry out differential expression analysis. Two methods can be choosen from for \code{paraMethod}. 'NP' stands for non-parametric method. 'MLE' stands for maximum liklihood estimation method. 'NP' is generally recommended for experiments with replicates smaller than 5.    
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
   H. I. Chen, Y. Liu, Y. Zou, Z. Lai, D. Sarkar, Y. Huang, et al.,
   "Differential expression analysis of RNA sequencing data by
   incorporating non-exonic mapped reads," BMC Genomics, vol. 16
   Suppl 7, p. S14, Jun 11 2015.
}
\seealso{
\code{\link{estimateRealCount}},  \code{\link{XBSeqDataSet}}, \code{\link{estimateSCV}}, \code{\link{XBSeqTest}}
}
\examples{
   conditions <- c(rep('C1', 3), rep('C2', 3))
   data(ExampleData)
   Stats <- XBSeq(Observed, Background, conditions)
}