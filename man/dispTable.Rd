\name{dispTable}
\alias{dispTable}
\alias{dispTable,XBSeqDataSet-method}

\title{
    Access the dispersion table information for a XBSeqDataSet object
}
\description{
   A method adopted from DESeq to examine the dispersion table information for a XBSeqDataSet object
}
\usage{
dispTable(object, ...)
}
\arguments{
  \item{object}{
      a XBSeqDataSet
}
  \item{...}{
      further argumnts are ignored
}  
}
\value{
   Dispersion table information for a XBSeqDataSet object
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
\code{\link{estimateSCV}}, \code{\link{dispEst}}, \code{\link{XBSeqDataSet}}
}
\examples{
   conditions <- factor(c(rep('C1', 3), rep('C2', 3)))
   data(ExampleData)
   XB <- XBSeqDataSet(Observed, Background, conditions)
   XB <- estimateRealCount(XB)
   XB <- estimateSizeFactors(XB)
   XB <- estimateSCV(XB, fitType='local')
   dispTable(XB)
}   
