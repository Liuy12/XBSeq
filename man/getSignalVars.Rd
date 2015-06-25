\name{getSignalVars}
\alias{getSignalVars}
\title{
   Estimate variance of the signal based on variance summation law
}
\description{
   Based on variance of observed signal as well as background noise, estimate the variance of the true signal
}
\usage{
getSignalVars(counts, bgcounts)
}
\arguments{
  \item{counts}{
   A data frame or matrix which contains the observed signal (expression level) information for an experiment. Rows represent genes and 
   Columns represent samples. Please refer details for more information. 
}
  \item{bgcounts}{
   A data frame or matrix which contains the background noise information for an experiment. Rows represent genes and 
   Columns represent samples. Please refer details for more information. 
}
}
\details{
   Observed signal are the reads mapped to the exonic regions which can be obtained by applying HTSeq procedure with GTF files of 
   exonic regions. Background noise are the reads mapped to the non-exonic regions which can be obtained by applying HTSeq procedure 
   with GTF files of non-exonic regions we defined by certain criteria. Details regarding how to carry out the HTSeq procedure for observed
   signal as well as background noise can be found in the vignette of XBSeq. One example dataset is provided in \code{\link{ExampleData}}. 
   
   By assuming that the true signal and background noise are independent, the variance of the underneath signal (\eqn{\sigma_s^2})
   can be estimated by applying variance summation law:
   \deqn{\sigma_s^2 = \sigma_x^2 + \sigma_b^2 - 2{\rho}{\sigma_x}{\sigma_b}}
   where \eqn{\sigma_x^2} and \eqn{\sigma_b^2} are variance for observed signal and background noise respectively. 
}
\value{
   A matrix with the same number of rows as \code{counts}. Rows represent the estimated variance of true signal for each gene. 
}
\author{
   Yuanhang Liu
}
\references{
   Hung-I Harry Chen, Yuanhang Liu, Yi Zou, Zhao Lai, Devanand Sarkar, Yufei Huang, Yidong Chen
   doi: \url{http://dx.doi.org/10.1101/016196}
}
\seealso{
\code{\link{estimateSCV}}
}
\examples{
   conditions <- factor(c(rep('C1', 3), rep('C2', 3)))
   data(ExampleData)
   data_var <- getSignalVars(Observed, Background)  
}