\name{XBSeq-package}
\alias{XBSeq-package}
\docType{package}
\title{
   Differential expression analysis of RNA sequencing data by incorporating non-exonic mapped reads
}
\description{
   We developed a novel algorithm, XBSeq, where a statistical model was established based on the 
   assumption that observed signals are the convolution of true expression signals and sequencing 
   noises. The mapped reads in non-exonic regions are considered as sequencing noises, which follows 
   a Poisson distribution. Given measureable observed signal and background noise from RNA-seq data, true 
   expression signals, assuming governed by the negative binomial distribution, can be delineated 
   and thus the accurate detection of differential expressed genes.
}
\details{
\tabular{ll}{
Package: \tab XBSeq\cr
Type: \tab Package\cr
Version: \tab 1.1.0\cr
Date: \tab 2015-12-02\cr
License: \tab >=GPL3\cr
Imports: \tab DESeq2, Biobase, pracma, matrixStats, ggplot2,locfit, methods, BiocGenerics, dplyr, plotly, MetricsGraphics\cr
}
}
\author{
   Yuanhang Liu

Maintainer: Yuanhang Liu <liuy12@uthscsa.edu>
}
\references{
   H. I. Chen, Y. Liu, Y. Zou, Z. Lai, D. Sarkar, Y. Huang, et al.,
   "Differential expression analysis of RNA sequencing data by
   incorporating non-exonic mapped reads," BMC Genomics, vol. 16
   Suppl 7, p. S14, Jun 11 2015.
}