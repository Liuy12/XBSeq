data("ExampleData")
conditions <- factor(c(rep('C1',3), rep('C2', 3)))
XB <- XBSeqDataSet(Observed, Observed, conditions)
expect_error(counts(XB, 3))