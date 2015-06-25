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
   Hung-I Harry Chen, Yuanhang Liu, Yi Zou, Zhao Lai, Devanand Sarkar, Yufei Huang, Yidong Chen
   doi: \url{http://dx.doi.org/10.1101/016196}
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
