\name{counts}
\docType{methods}
\alias{counts}
\alias{counts,XBSeqDataSet-method}
\title{
   Accessor functions for the 'counts' information in a
   XBSeqDataSet object.}
\description{
   The 'counts' function extract a certian assay element from XBSeqDataSet object. The normalized assay element can be extracted by specifying 'normalized = TRUE'. 
}
\usage{
\S4method{counts}{XBSeqDataSet}(object,slot = 3, normalized = FALSE)
}
\arguments{
  \item{object}{a XBSeqDataSet}
  \item{slot}{a integer value to specify which assay element to extract (default to 3)}
  \item{normalized}{whether the normalized assay element should be returned}
  }
\details{
   counts is a function to access an array elemen which is specified by the end user. The difference between this function and the counts function for \code{DESeqDataSet} is that this function can be used to access a specific array elememt rather than a pre-defined array element "counts" in the case of \code{DESeqDataSet}.
   By default, the first array element contains information of observed signal. The second array element contains information of background noise. The third array element contains information of estimated true signal after calling the function \code{estimateRealCount}. 
}
\value{
   Either normalized or un-normalized assay element
}
\references{
   Hung-I Harry Chen, Yuanhang Liu, Yi Zou, Zhao Lai, Devanand Sarkar, Yufei Huang, Yidong Chen
   doi: \url{http://dx.doi.org/10.1101/016196}
}
\seealso{
\code{\link{XBSeqDataSet}}, \code{\linkS4class{DESeqDataSet}}
}
\examples{
   data(ExampleData)
   conditions <- factor(c(rep('C1', 3), rep('C2', 3)))
   XB <- XBSeqDataSet(Observed, Background, conditions)
   str(counts(XB, 1))
}