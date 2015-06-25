\name{XBSeqDataSet-class}
\Rdversion{1.1}
\docType{class}
\alias{XBSeqDataSet-class}
\alias{XBSeqDataSet}

\title{Class \code{"XBSeqDataSet"}}
\description{
\code{XBSeqDataSet} is a subclass of \code{"\linkS4class{DESeqDataSet}"}, used to store the input values, intermediate calculations and results of an
analysis of differential expression.  Different from the original \code{DESeqDataSet} class, \code{XBSeqDataSet} has some extra slots including:
\itemize{
\item{fitInfo:} {An object of \code{environment} class which contains the scv fitting information for a \code{XBSeqDataSet} object}
\item{dispTable:} {An object of \code{character} class which indicates  method used for scv fitting. Details can be found in \code{\link{estimateSCV}}.}
\item{conditions:} {An object of \code{factor} class which contains the experimental design information for a \code{XBSeqDataSet} object}
\item{dispEst:} {An object of \code{list} class which contains the final dispersion estimates for each gene. Details can be found in \code{\link{dispEst}}}
}}
\usage{
XBSeqDataSet(counts, bgcounts, conditions, sizeFactors=NULL, ...)
}
\arguments{
   \item{counts}{A data frame or matrix which contains the observed signal for each gene across all the samples. Rows represent genes and columns represent samples.}
   \item{bgcounts}{A data frame or matrix which contains the background noise for each gene across all the samples. Rows represent genes and columns represent samples.}
   \item{conditions}{Object of class \code{"character"}. The conditions for the experimental design. }
   \item{sizeFactors}{Numeric vector which contains normalizing factors for the data matrix. In most cases, it is recommended that you calculate
	sizeFactors by \code{estimateSizeFactors}. You are also able to provide sizeFactors yourself.}
	\item{...}{Further arguments provided will be ignored }
}
\section{Methods}{
  \describe{
    \item{conditions}{\code{signature(object = "XBSeqDataSet")}: ... }
    \item{conditions<-}{\code{signature(object = "XBSeqDataSet")}: ... }
    \item{counts}{\code{signature(object = "XBSeqDataSet")}: ... }
    \item{dispTable}{\code{signature(object = "XBSeqDataSet")}: ... }
    \item{dispEst}{\code{signature(object = "XBSeqDataSet")}: ... }
    \item{dispEst<-}{\code{signature(object = "XBSeqDataSet")}: ... }
	  \item{fitInfo}{\code{signature(object = "XBSeqDataSet")}: ... }
    \item{fitInfo<-}{\code{signature(object = "XBSeqDataSet")}: ... }
    \item{estimateSCV}{\code{signature(object = "XBSeqDataSet")}: ... }
    \item{estimateSizeFactors}{\code{signature(object = "XBSeqDataSet")}: ... }
    \item{estimateRealCount}{\code{signature(object = "XBSeqDataSet")}: ... }
	 }
}

\value{
A XBSeqDataSet object.
}

\author{
   Yuanhang Liu
}
\references{
   Hung-I Harry Chen, Yuanhang Liu, Yi Zou, Zhao Lai, Devanand Sarkar, Yufei Huang, Yidong Chen
   doi: \url{http://dx.doi.org/10.1101/016196}
}
\seealso{
\code{\link{estimateSCV}}, \code{\link{conditions}}, \code{\link{dispEst}}, \code{\link{dispTable}}, \code{\link{fitInfo}}, \code{\linkS4class{DESeqDataSet}}, \code{\link{counts}}, \code{\link{estimateRealCount}}
}
\examples{
   data(ExampleData)
   conditions <- factor(c(rep('C1', 3), rep('C2', 3)))
   XB <- XBSeqDataSet(Observed, Background, conditions)
   str(XB)
}