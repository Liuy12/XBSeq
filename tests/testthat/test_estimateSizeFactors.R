data("ExampleData")
conditions <- factor(c(rep('C1',3), rep('C2', 3)))
XB <- XBSeqDataSet(Observed, Background, conditions)
XB <- estimateRealCount(XB)
expect_error(estimateSCV(XB, method='pooled', sharingMode='maximum', fitType='local'))
expect_error(counts(XB, normalized = TRUE))




