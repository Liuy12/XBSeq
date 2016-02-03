XBplot <- function(XB, Samplenum = NULL, unit = c('counts', 'LogTPM'), Libsize = NULL, Genelength = NULL, xlab = 'log2 TPM', ylab = 'Frequencies', col = c('blue', 'red'), alpha =c(1, 0.6)){
  if(is.null(Samplenum))
    stop("You need to provide the column number of the sample you want to examine")
  Observed <- as.data.frame(counts(XB, slot = 1)) %>% select(Samplenum) %>% mutate(Group = 'Observed')
  Background <- as.data.frame(counts(XB, slot = 2)) %>% select(Samplenum) %>% mutate(Group = 'Background')
  colnames(Observed) <- c('Sample', 'Group')
  colnames(Background) <- c('Sample', 'Group')
  unit <- match.arg(unit, c('counts', 'LogTPM'))
  if(unit == 'counts'){
    xlab <- 'Counts'
    xlim <- c(0, median(Observed$Sample))
    binwidth <- 1
    Combined <- bind_rows(Observed, Background)
    }
  else{
    xlab <- 'Log2 TPM'
    if(is.null(Libsize)){
      warning("Libsize is not provided, the sum of all the read counts that mapped to exonic
regions in each sample is used as the total library size for that sample")
      Libsize <- sum(Observed$Sample)
    }
    if(is.null(Genelength))
      stop("Please provide the gene length information if you choose 'unit' equals to 'LogTPM'")
    if(!is.numeric(Genelength))
      stop("Please make sure that 'genelength' is a numeric vector of the same order and length of your XB object")
    if(nrow(Observed) != length(Genelength))
      stop("Please make sure 'genelength' information is of the same length with your XB object")
    Genelength <- as.numeric(Genelength)
    Observed$Sample <- Observed$Sample*10^9/(Genelength*Libsize)
    Background$Sample <- Background$Sample*10^9/(Genelength*Libsize)
    Observed$Sample <- Observed$Sample*10^6/sum(Observed$Sample, Background$Sample)
    Background$Sample <- Background$Sample*10^6/sum(Observed$Sample, Background$Sample)
    Combined <- bind_rows(Observed, Background)
    Combined$Sample <- log2(Combined$Sample)
    xlim <- range(Combined$Sample[!is.infinite(Combined$Sample)])
    binwidth <- 0.1
  }
  ggplot(data = Combined) + 
    geom_histogram(aes(x = Sample, fill=Group, y=..count.., alpha=Group), 
                   binwidth=binwidth, position='identity') +
    scale_fill_manual(values = col) + scale_alpha_manual(values = alpha, guide = FALSE) +
    guides(fill=guide_legend('')) +
    labs(x=xlab, y=ylab) + xlim(xlim)
}


MAplot <- function(stats, ylim, padj=TRUE, pcuff=0.1, lfccuff=1, linecol='red3',
                   xlab='mean of normalized counts', ylab=expression(log[2]~fold~change), shape)
{
  if(!(is.data.frame(stats) && all(c("baseMean", "log2FoldChange") %in% colnames(stats))))
    stop("'stats' must be a data frame with columns named 'baseMean', 'log2FoldChange'.")
  if(padj)
    col = ifelse(stats$padj>=pcuff, "gray32", "red")
  else
    col = ifelse(stats$pval>=pcuff, "gray32", "red")
  col = col[ stats$baseMean != 0 ]
  y = stats$log2FoldChange[ stats$baseMean != 0 ]
  stats = subset(stats, baseMean!=0)
  if(missing(ylim))
    ylim = c(-1,1) * quantile(abs(y[is.finite(y)]), probs=0.99) * 1.1
  if (missing(shape))
    shape = ifelse(y<ylim[1], 6, ifelse(y>ylim[2], 2, 16) )
  stats$log2FoldChange = pmax( ylim[1], pmin(ylim[2], y) )
  ggplot() + geom_point( data=stats,aes( x=baseMean, y=log2FoldChange ), color=col, shape=shape ) +
    ylim(ylim) + geom_hline(yintercept=0,colour=linecol,size=1)  + scale_x_log10() +
    labs( x=xlab, y=ylab )
}


plotSCVEsts = function( XB, name=NULL, ymin, linecol='red3',
                        xlab = "mean of normalized counts", ylab = "SCV")
{
  px = rowMeans(counts(XB, normalized=TRUE))
  sel = (px>0)
  px = px[sel]
  py = fitInfo(XB, name=name)$perGeneSCVEsts[sel]
  if(missing(ymin))
    ymin = 10^floor(log10(min(py[py>0], na.rm=TRUE))-0.1)
  py = pmax(py, ymin)
  shape = ifelse(py<ymin, 6, 16)
  sel1 = complete.cases(px)&complete.cases(py)
  px = px[sel1]
  py = py[sel1]
  shape = shape[sel1]
  fitd = data.frame(px=px,py=py)
  fx = 10^seq( -.5, 5, length.out=100 )
  fy = fitInfo(XB, name=name)$SCVFunc(fx)
  fitl = data.frame(fx=fx, fy=fy)
  ggplot() + geom_point( data=fitd, aes( x=px, y=py), shape=shape) +
    geom_line( data=fitl, aes ( x=fx, y=fy), col=linecol, size=1.5, alpha=0.6) + scale_x_log10() +
    scale_y_log10() + labs(x=xlab, y=ylab)
}