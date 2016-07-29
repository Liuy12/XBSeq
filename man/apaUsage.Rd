\name{apaUsage}
\docType{methods}
\alias{apaUsage}
\title{
   Testing differential alternative polyadenylation (apa) usage by using roar
}
\description{
   XBSeq now offers testing differential apa usage via funcitionality provided by roar package
}
\usage{
   apaUsage(bamTreatment, bamControl, apaAnno, paired = NULL)
}
\arguments{
  \item{bamTreatment}{ A list of full path of filenames of bam alignments with data for the treatment condition (by convention it is considered the “treated” condition: }
   \item{bamControl}{  A list of full path of filenames of bam alignments with data for the control condition to be considered. } 
   \item{apaAnno}{full path of apa annotation used by roar package. APA annotation for several organisms of various genome build can be downloaded from [here](https://github.com/Liuy12/XBSeq_files). For details regarding how to construct APA annotation, please refer to vignette}
   \item{paired}{ a vector indicating how samples are paired, leave to NULL if the experiment is not paired.
   }
}
\details{
Bioconductor package roar is used to detect preferential usage of shorter isoforms via alternative poly-adenylation from RNA-seq data. The approach is based on Fisher test to detect dis-equilibriums in the number of reads falling over the 3’UTRs when comparing two biological conditions. 
}
\value{
   The resulting data frame will have the "gene_id" of the initial annotation as row names and as columns the m/M ratio for the treatment and control conditions, the roar value and the Fisher test pvalue, expression value of treatment group, expression value of control groups (respectively: mM_treatment, mM_control, roar, pval, treatmentValue, controlValue). If more than one sample has been given for a condition the "pval" column will contain the multiplication of all the comparisons pvalue and there will be other columns containing the pvalues resulting from all the pairwise treatment vs control contrasts, with names "pvalue_X_Y" where X represent the position of the sample in the treatment list of bam files and Y the position for the control list
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