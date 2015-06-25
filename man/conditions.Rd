\name{conditions}
\docType{methods}
\alias{conditions}
\alias{conditions<-,XBSeqDataSet-method}
\alias{conditions,XBSeqDataSet-method}
\title{
   Accessor functions for the 'conditions' information in a
   XBSeqDataSet object.}
\description{
   Conditions extract the experimental design information similar as used in DESeq. 
}
\usage{
\S4method{conditions}{XBSeqDataSet}(object,...)
\S4method{conditions}{XBSeqDataSet}(object,...) <- value
}
\arguments{
  \item{object}{a XBSeqDataSet}
  \item{value}{experimental design information}
  \item{...}{Further arguments will be ignored}
  }
\value{
   The experimental design information for a XBSeqDataSet object
}
\references{
   Hung-I Harry Chen, Yuanhang Liu, Yi Zou, Zhao Lai, Devanand Sarkar, Yufei Huang, Yidong Chen
   doi: \url{http://dx.doi.org/10.1101/016196}
}
\seealso{\code{\link{XBSeqDataSet}}}
\examples{
   data(ExampleData)
   conditions <- factor(c(rep('C1', 3), rep('C2', 3)))
   XB <- XBSeqDataSet(Observed, Background, conditions)
   conditions(XB)
}  