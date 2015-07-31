---
title: "Differential expression analysis of count data using XBSeq package"
author: "Yuanhang Liu"
date: "`r Sys.Date()`"
vignette: >
  %\VignetteIndexEntry{Differential expression analysis of count data using XBSeq package}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8} 
output: 
  BiocStyle::html_document:
    toc: true
---

```{r style, echo = FALSE, results = 'asis'}
BiocStyle::markdown()
```

## Introduction 

XBSeq is a novel algorithm for testing RNA-seq differential expression (DE), where a statistical model was established based on the assumption that observed signals are the convolution of true expression signals and sequencing noises. The mapped reads in non-exonic regions are considered as sequencing noises, which follows a Poisson distribution. Given measurable observed signal and background noise from RNA-seq data, true expression signals, assuming governed by the negative binomial distribution, can be delineated and thus the accurate detection of differential expressed genes

## Installation 

XBSeq can be installed from Bioconductor by 
```{r,eval=FALSE}
source('http://www.bioconductor.org/biocLite.R')
biocLite("XBSeq")
```
```{r,message = FALSE, warning=FALSE}
library("XBSeq")
```
If you would like to install the development version of XBSeq, it is recommended that you refer to the github page of XBSeq. 

## Use XBSeq for testing differential expression 

#### HTSeq procedure

In order to use XBSeq for testing DE, after sequence alignment, we need to run HTSeq twice to measure the reads mapped to exonic regions (observed signal) and non-exonic regions (background noise). Generally speaking, you will need to run the following code to generate observed signal and background noise. 

```{r,engine='python',eval=FALSE}
htseq-count [options] <alignment_file> <gtf_file> > Observed_count.txt
htseq-count [options] <alignment_file> <gtf_file_bg> > background_count.txt
```

Details regarding how HTSeq works can be found here: http://www-huber.embl.de/HTSeq/doc/count.html

The gtf file used to measure observed signal can be downloaded from UCSC database: http://genome.ucsc.edu. The gtf file used to measure background noise can be downloaded in the gtf folder from github: https://github.com/Liuy12/XBSeq_files. If you would like to construct the gtf file by yourself, we also have deposited the perl script we used to construct the gtf file in github. Details regarding the procedure we used to construct the background gtf file can be found in the Details section in the vignette.

#### XBSeq testing for DE 

After HTSeq procedure, then we will have two measurements for each gene, the observed signal and background noise. Here we will use a mouse RNA-seq dataset, which contains 3 replicates of wild type mouse liver tissues (WT) and 3 replicates of Myc transgenic mouse liver tissues (MYC). The dataset is obtained from Gene Expression Omnibus [(GSE61875)](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE61875). 

As a preliminary step, we have already carried out HTSeq procedure mentioned above to generate observed signal and background noise for each gene. The two datasets can be loaded into user's working space by 

```{r}
data(ExampleData)
```

We can first take a look at the two datasets:

```{r}
head(Observed)
head(Background)
```

Rows represent reads mapped to each gene or corresponding background region. Column represent samples. And differential expression analysis will be carried out as follows:

Firstly, we need to construct a XBSeqDataSet object. Conditions are the design matrix for the experiment. Observe and background are the output matrix from HTSeq (Remeber to remove the bottom few lines of summary statistics of the output matrix).  

```{r,tidy=TRUE}
conditions <- factor(c(rep('C',3), rep('T', 3)))
XB <- XBSeqDataSet(Observed, Background, conditions)
```

It is always recommended that you examine the distribution of observed signal and background noise beforehand. We provide function 'r XBplot' to achieve this. We recommended to examine the distribution in log2 reads per kilobase per million (RPKM) unit by setting argument unit equals to "RPKM". Genelength information is loaded via "ExampleData". Ideally, library size should also be provided. By default, the sum of all the reads that mapped to exonic regions are used.

```{r,tidy=TRUE,fig.width=5,fig.height=4}
XBplot(XB, Samplenum = 1, unit = "LogRPKM", Genelength = genelength[,2])
```

Then estimate the preliminary underlying signal followed by normalizing factor and dispersion estimates 

```{r, tidy=TRUE}
XB <- estimateRealCount(XB)
XB <- estimateSizeFactors(XB)
XB <- estimateSCV( XB, method='pooled', sharingMode='maximum', fitType='local' )
```

Take a look at the scv fitting information 

```{r,fig.width=3,fig.height=3}
plotSCVEsts(XB)
```

Carry out the DE test by using function XBSeqTest

```{r}
Teststas <- XBSeqTest( XB, levels(conditions)[1L], levels(conditions)[2L] )
```

Plot Maplot based on test statistics

```{r,fig.width=3,fig.height=3}
MAplot(Teststas, padj = FALSE, pcuff = 0.01, lfccuff = 1)
```

```{r,eval=FALSE,tidy=TRUE}
# Alternatively, all the codes above can be done with a wrapper function XBSeq
Teststats <- XBSeq( Observed, Background, conditions, method='pooled', sharingMode='maximum',
  fitType='local', pvals_only=FALSE )
```

#### Compare the results with DESeq

Now we will carry out DE analysis on the same dataset by using DESeq and then compare the results obtained by these two methods

If you have not installed DESeq before, DESeq is also available from Bioconductor

```{r,eval=FALSE}
biocLite("DESeq")
```

Then DE analysis for DESeq can be carried out by:

```{r,message=FALSE}
library('DESeq')
library('ggplot2')
de <- newCountDataSet(Observed, conditions)
de <- estimateSizeFactors(de)
de <- estimateDispersions(de, method = "pooled", fitType="local")
res <- nbinomTest(de, levels(conditions)[1], levels(conditions)[2])
```

Then we can compare the results from XBSeq and DESeq

```{r,warning=FALSE,message=FALSE,tidy=TRUE, fig.width=3,fig.height=3}
DE_index_DESeq <- with(res, which(pval<0.01 & abs(log2FoldChange)>1))
DE_index_XBSeq <- with(Teststas, which(pval<0.01 & abs(log2FoldChange)>1))
DE_index_inters <- intersect(DE_index_DESeq, DE_index_XBSeq)
DE_index_DESeq_uniq <- setdiff(DE_index_DESeq, DE_index_XBSeq)
DE_plot <- MAplot(Teststas, padj = FALSE, pcuff = 0.01, lfccuff = 1, shape=16)
DE_plot + geom_point( data=Teststas[DE_index_inters,], aes(x=baseMean, y=log2FoldChange),
                      color= 'green', shape=16 ) + 
  geom_point( data=Teststas[DE_index_DESeq_uniq,], aes( x=baseMean, y=log2FoldChange ),
              color= 'blue', shape=16 )
```

The red dots indicate DE genes identified only by XBSeq. Then green dots are the shared results of XBSeq and DESeq. The blue dots are DE genes identified only by DESeq. 

## Details 

#### Construction of gtf file for background region

* Exonic region annotation is obtained from UCSC database. 

* Non-exonic regions are constructed by following several criteria: 
    1. Download refFlat table from UCSC database 
    2. Several functional elements (mRNA, peudo genes, etc.) of the genome are excluded from the whole genome annotation. 
    3. Construct the background region for each gene by making the region to have the same sturcture or length as the exonic region of the gene. 

More details regarding how do we construct the background region annotation file of an real example can be found in manual page of ExampleData and also our publication of XBSeq.  

## Bug reports
Report bugs as issues on our [GitHub repository](https://github.com/Liuy12/XBSeq/issues) or you can report directly to my email: liuy12@uthscsa.edu.

## Session information 
```{r}
sessionInfo()
```

## Acknowledgements 
XBSeq is implemented in R based on the source code from DESeq and DESeq2. 

## References
H. I. Chen, Y. Liu, Y. Zou, Z. Lai, D. Sarkar, Y. Huang, et al., "Differential expression analysis of RNA sequencing data by incorporating non-exonic mapped reads," BMC Genomics, vol. 16 Suppl 7, p. S14, Jun 11 2015.