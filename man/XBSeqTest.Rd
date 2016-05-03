\name{XBSeqTest}
\alias{XBSeqTest}
\title{
   XBSeq test for differential expression
}
\description{
   The same method is adopted from DESeq for testing differential expression 
}
\usage{
XBSeqTest(XB, condA, condB, pvals_only = FALSE, method = c("NP", "MLE"))
}
\arguments{
  \item{XB}{
   A XBSeqDataSet object
}
  \item{condA}{
   Factor level specified for condition A
}
  \item{condB}{
   Factor level specified for condition B
}
  \item{pvals_only}{
   Logical;whether or not only extract p values
}
  \item{method}{
   method to use for estimation of distribution parameters, 'NP' or 'MLE'. See details section for details
}
}
\details{
   Differential expression analysis based on statistical methods proposed for DESeq. Details about the method can be found in DESeq manual page. Two methods can be choosen from for \code{method}. 'NP' stands for non-parametric method. 'MLE' stands for maximum liklihood estimation method. 'NP' is generally recommended for experiments with replicates smaller than 5.    
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
\references{
   H. I. Chen, Y. Liu, Y. Zou, Z. Lai, D. Sarkar, Y. Huang, et al.,
   "Differential expression analysis of RNA sequencing data by
   incorporating non-exonic mapped reads," BMC Genomics, vol. 16
   Suppl 7, p. S14, Jun 11 2015.
}
\author{
   Yuanhang Liu
}
\seealso{
\code{\link{XBSeq}}, \code{\link{estimateSCV}}
}
\examples{
   data(ExampleData)
   conditions <- factor(c(rep('C1', 3), rep('C2', 3)))
   XB <- XBSeqDataSet(Observed, Background, conditions)
   XB <- estimateRealCount(XB)
   XB <- estimateSizeFactors(XB)
   XB <- estimateSCV(XB)
   Teststas <- XBSeqTest(XB, levels(conditions)[1L], levels(conditions)[2L])
   str(Teststas)
   }