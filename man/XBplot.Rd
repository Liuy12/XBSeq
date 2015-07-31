\name{XBplot}
\alias{XBplot}
\title{
   Examine the distribution of observed signal and background noise
}
\description{
   Function to viewlize the distribution of observed signal X and background noise B across all genes for one specified sample
}
\usage{
XBplot(XB, Samplenum = NULL, unit = c('counts', 'LogRPKM'), Libsize = NULL, Genelength = NULL, xlab = 'log2 RPKM', ylab = 'Frequencies', col = c('blue', 'red'), alpha =c(1, 0.6))
}
\arguments{
  \item{XB}{
   An XBSeqDataSet object
}
  \item{Samplenum}{
   An integer number to specify which sample to examine
}
  \item{unit}{
   Whether to examine the distrbution in 'counts' unit or 'LogRPKM' unit. 'LogRPKM' is generally recommended
}
  \item{Libsize}{
   A single integer indicating the library size of the sample. By default, the sum of all reads mapped to exonic regions are used.
}
  \item{Genelength}{
   A numeric vector containing genelength information. Please make sure the length and order of the gene length information is the same as arrays in the XB object. 
}
  \item{xlab}{
   lab for x axis
}
  \item{ylab}{
   Lable for y axis
}
  \item{col}{
   A vector of two colours for observed signal and background noise
}
  \item{alpha}{
   A vector of two numeric numbers indicating transparency
}
}
\details{
   We strongly recommended users to apply XBplot to their datasets before differential expression analysis. According to our experience, for XBplot in 'logRPKM' unit, the peak of distribution of background noise generally coinsides with the left hump of distribution of observed signal. 
}
\value{
   Plot of distribution of observed signal and background noise.
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
\seealso{\code{\link{ExampleData}}}
\examples{
   conditions <- c(rep('C1', 3), rep('C2', 3))
   data(ExampleData)
   XB <- XBSeqDataSet(Observed, Background, conditions)
   XBplot(XB, Samplenum = 1, unit = "LogRPKM", Genelength = genelength[,2])
}