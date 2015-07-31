data("ExampleData")
conditions <- factor(c(rep('C',3), rep('T', 3)))
XB <- XBSeqDataSet(Observed, Background, conditions)
expect_error(XBplot(XB))

expect_error(XBplot(XB, Samplenum = 'Sample_54_WT'))

expect_error(XBplot(XB, Samplenum = 7))

expect_warning(XBplot(XB, Samplenum = 1, unit = 'LogRPKM', Genelength = genelength[,2]))

expect_error(XBplot(XB, Samplenum = 1, unit = 'LogRPKM'))

expect_error(XBplot(XB, Samplenum = 1, unit = 'LogRPKM', Genelength = genelength))

expect_error(XBplot(XB, Samplenum = 1, unit = 'LogRPKM', Genelength = genelength[-1,2]))