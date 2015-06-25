\name{XBSeqTest}
\alias{XBSeqTest}
\title{
   XBSeq test for differential expression
}
\description{
   The same method is adopted from DESeq for testing differential expression 
}
\usage{
XBSeqTest(XB, condA, condB, pvals_only = FALSE)
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
\details{
   Differential expression analysis based on statistical methods proposed for DESeq. Details about the method can be found in DESeq manual page. 
}
   
\references{
   Hung-I Harry Chen, Yuanhang Liu, Yi Zou, Zhao Lai, Devanand Sarkar, Yufei Huang, Yidong Chen
   doi: \url{http://dx.doi.org/10.1101/016196}
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