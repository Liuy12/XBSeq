\name{estimateRealCount}
\docType{methods}
\alias{estimateRealCount}
\alias{estimateRealCount,XBSeqDataSet-method}
\title{
   Preliminary step to estimate the true signal based on observed signal
   and background noise
}
\description{
   Based on the observed signal as well as the background noise, estimate the true signal 
   for each gene. 
}
\usage{
   estimateRealCount(object)
}
\arguments{
  \item{object}{ A XBSeqDataSet object}
}
\details{
   The observed signal can be achieved by using HTSeq to count the reads map to exonic regions.
   The background noise can be extracted by using HTSeq the second time to count the reads map
   to non-exonic regions, the regions we defined by excluding potential functional elements. 
   The the underneath true signal is estimated by the simple subtraction of observed signal and
   background noise. The true signal of genes with background noise larger than observed signal
   will be assigned as 0.  
}
\value{
   A matrix contains the estimated true signal for each gene with the same length as observed
   signal.
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
\seealso{\code{\link{counts}}}
\examples{
   conditions <- factor(c(rep('C1', 3), rep('C2', 3)))
   data(ExampleData)
   XB <- XBSeqDataSet(Observed, Background, conditions)
   XB <- estimateRealCount(XB)
   str(counts(XB, 3))
}
