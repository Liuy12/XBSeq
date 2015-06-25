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
  px = rowMeans( counts( XB, normalized=TRUE ) )
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
    geom_line( data=fitl, aes ( x=fx, y=fy), shape=shape, col=linecol, size=1.5, alpha=0.6) + scale_x_log10() +
    scale_y_log10() + labs(x=xlab, y=ylab)
}