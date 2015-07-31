\name{plotSCVEsts}
\alias{plotSCVEsts}
\title{
   Plot estimated squared coefficient of variation  
}
\description{
   Plot estimated SCV based on ggplot2
}
\usage{
plotSCVEsts(XB, name = NULL, ymin, linecol = "red3", 
   xlab = "mean of normalized counts", ylab = "SCV")
}
\arguments{
  \item{XB}{
   A XBSeqDataSet object
}
  \item{name}{
   The name of the fit information. Only specify this if you choose \code{method="per-condition"}
}
  \item{ymin}{
   The limit of y axis
}
  \item{linecol}{
   The linecolour of the SCV-mean trend
}
  \item{xlab}{
   The lable of x axis
}
  \item{ylab}{
   The lable of y axis
}
}
\value{
   Summary plot for the fitting and estimation of scv
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
\code{\link{estimateSCV}}
}
\examples{
   conditions <- factor(c(rep('C1', 3), rep('C2', 3)))
   data(ExampleData)
   XB <- XBSeqDataSet(Observed, Background, conditions)
   XB <- estimateRealCount(XB)
   XB <- estimateSizeFactors(XB)
   XB <- estimateSCV(XB, fitType='local')
   plotSCVEsts(XB)
}