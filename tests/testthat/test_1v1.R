set.seed(1990)
Observed <- matrix(rnbinom(5000*2, 10, 0.5), 5000, 2)
Background <- matrix(rpois(5000*2, 3), 5000, 2)
rownames(Observed) <- paste('G', 1:5000, sep = '')
rownames(Background) <- paste('G', 1:5000, sep='')
conditions <- factor(c('C1', 'C2'))
XB <- XBSeqDataSet(Observed, Background, conditions)
XB <- estimateRealCount(XB)
XB <- estimateSizeFactors(XB)
expect_error(estimateSCV(XB, method='pooled', sharingMode='maximum', fitType='local'))



