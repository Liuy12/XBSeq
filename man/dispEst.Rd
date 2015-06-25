\name{dispEst}
\docType{methods}
\alias{dispEst}
\alias{dispEst,XBSeqDataSet-method}
\alias{dispEst<-}
\alias{dispEst<-,XBSeqDataSet-method}
\title{
   Function to access the dispersion estimation for each gene 
}
\description{
   The dispersion estimated for each gene are stored as a data.frame after user called \code{\link{estimateSCV}}
}
\usage{
dispEst(object, varname = NA)
dispEst(object, varname = NA) <- value
}
\arguments{
  \item{object}{
   XBSeqDataSet object
}
  \item{varname}{
   variable name of dispersion estimates
}
  \item{value}{
   The dispersion estimates for each gene 
}
}
\value{
   A data.frame which contains the dispersion estimates for each gene
}
\author{
   Yuanhang Liu
}
\references{
   Hung-I Harry Chen, Yuanhang Liu, Yi Zou, Zhao Lai, Devanand Sarkar, Yufei Huang, Yidong Chen
   doi: \url{http://dx.doi.org/10.1101/016196}
}
\seealso{
\code{\link{estimateSCV}}, \code{\link{dispTable}}, \code{\link{XBSeqDataSet}}
}
\examples{
   conditions <- factor(c(rep('C1', 3), rep('C2', 3)))
   data(ExampleData)
   XB <- XBSeqDataSet(Observed, Background, conditions)
   XB <- estimateRealCount(XB)
   XB <- estimateSizeFactors(XB)
   XB <- estimateSCV(XB, fitType='local')
   str(dispEst(XB))
}  
