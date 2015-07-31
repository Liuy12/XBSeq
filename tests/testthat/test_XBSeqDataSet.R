set.seed(1990)
Observed <- matrix(rnbinom(2000, 10,0.5) + 0.1, 500, 4)
Background <- matrix(rpois(2000,2) + 0.1, 500, 4)
conditions <- factor(c(rep('C1',2), rep('C2', 2)))
expect_error(XBSeqDataSet(Observed, Background, conditions))

Observed <- matrix(rnbinom(2000, 10,0.5), 500, 4)
Background <- matrix(rpois(1996,2), 499, 4)
conditions <- factor(c(rep('C1',2), rep('C2', 2)))
expect_error(XBSeqDataSet(Observed, Background, conditions))

rownames(Observed) <- 1:500
rownames(Background) <- 1:499
expect_error(XBSeqDataSet(Observed, Background, conditions))