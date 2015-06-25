\name{MAplot}
\alias{MAplot}
\title{
   Generate maplot after differential expression test
}
\description{
   Generate maplot after differential expression test based on ggplot2
}
\usage{
MAplot(stats, ylim, padj = TRUE, pcuff = 0.1, lfccuff = 1, 
   linecol = "red3", xlab = "mean of normalized counts",
   ylab = expression(log[2] ~ fold ~ change), shape) 
}
\arguments{
  \item{stats}{
   The output of \code{XBSeqTest}
}
  \item{ylim}{
   Range of limit for y axis
}
  \item{padj}{
   Whether to use adjusted p value or not
}
  \item{pcuff}{
   Threshold for pvalue
}
  \item{lfccuff}{
   Log fold change cutoff
}
  \item{linecol}{
   Colour of horizontal line 
}
  \item{xlab}{
   Lable for x axis
}
  \item{ylab}{
   Lable for y axis
}
  \item{shape}{
   The shape of the points used
}  
}
\details{
   Generate classic MAplot for DE analysis using ggplot2, where A and M are from slot \code{baseMean} and slot \code{log2FoldChange} of 
   the test statistics aftering calling \code{XBSeqTest}. The ggplot2 package generally generate figures of better quality as well as give 
   user better control of the plotting system compared with the base plotting system. 
}
\value{
   MAplot of test statistics
}
\author{
   Yuanhang Liu
}
\references{
   Hung-I Harry Chen, Yuanhang Liu, Yi Zou, Zhao Lai, Devanand Sarkar, Yufei Huang, Yidong Chen
   doi: \url{http://dx.doi.org/10.1101/016196}
}
\seealso{\code{\link{XBSeqTest}}}
\examples{
   conditions <- c(rep('C1', 3), rep('C2', 3))
   data(ExampleData)
   Stats <- XBSeq(Observed, Background, conditions)
   MAplot(Stats)
}