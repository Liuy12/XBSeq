\name{fitInfo}
\docType{methods}
\alias{fitInfo}
\alias{fitInfo,XBSeqDataSet-method}
\title{
    Accessor function for the fitInfo objects in a XBSeqDataSet
}
\description{
   Same method is adopted from DESeq to access the fit information from a XBSeqDataSet
}
\usage{
fitInfo( object, name)
}
\arguments{
  \item{object}{
      a XBSeqDataSet
}
  \item{name}{
      if \code{estimateSCV} was called with \code{method="per-condition"}
      a name hasd to specified. Try \code{ls(XB@fitInfo}.
}
}
\value{
   A list containing fitting information for a XBSeqDataSet object:
   \item{perGeneSCVEsts}{ SCV estimates for each gene, which has the same length as the number of rows as the assay elements in an object}
   \item{SCVFunc}{ The function used to predict the fitted SCV}
   \item{fittedSCVEsts}{ The fitted SCV estimates for each gene, which is of the same length as perGeneSCVEsts}
   \item{df}{ Integer value indicating the degree of freedom}
   \item{sharingMode}{ The sharing mode argument specified by the user}
}
\seealso{\code{\link{estimateSCV}}, \code{\link{XBSeqDataSet}}}
\author{
   Yuanhang Liu
}
\references{
   H. I. Chen, Y. Liu, Y. Zou, Z. Lai, D. Sarkar, Y. Huang, et al.,
   "Differential expression analysis of RNA sequencing data by
   incorporating non-exonic mapped reads," BMC Genomics, vol. 16
   Suppl 7, p. S14, Jun 11 2015.
}
\examples{
   conditions <- factor(c(rep('C1', 3), rep('C2', 3)))
   data(ExampleData)
   XB <- XBSeqDataSet(Observed, Background, conditions)
   XB <- estimateRealCount(XB)
   XB <- estimateSizeFactors(XB)
   XB <- estimateSCV(XB, fitType='local')
   str(fitInfo(XB))
}   

