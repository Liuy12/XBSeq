getCountParams <- function( counts, sizeFactors ) {

  # Divides the counts by sizeFactors and calculates the estimates for
  # base means and variances for each gene

  data.frame(
    baseMean = rowMeans( t( t(counts) / sizeFactors ) ),
    baseVar = rowVars( t( t(counts) / sizeFactors ) ) )
}


getCountParamsPooled <- function( counts, sizeFactors, conditions ) {

  basecounts <- t( t(counts) / sizeFactors )
  replicated_sample <- conditions %in% names(which(table(conditions)>1))
  df <- sum(replicated_sample) - length( unique( conditions[ replicated_sample ] ) )

  data.frame(
    baseMean = rowMeans( basecounts ),
    baseVar =
      rowSums(
        sapply(
          tapply(
            ( seq_len(ncol(counts)) )[ replicated_sample ],
            factor( conditions[ replicated_sample ] ),
            function(cols)
              rowSums( ( basecounts[,cols] - rowMeans(basecounts[,cols]) )^2 ) ),
          identity ) ) / df )
}

getSCV <- function( means,
                    variances, sizeFactors, fitType = c( "parametric", "local" ),
                    locfit_extra_args=list(), lp_extra_args=list(), adjustForBias=TRUE ) {

  fitType <- match.arg( fitType )

  xim <- mean( 1/sizeFactors )
  SCVAll <- ( variances - xim * means ) / means^2

  variances <- variances[ means > 0 ]
  SCV <- SCVAll[ means > 0 ]
  means <- means[ means > 0 ]

  if( adjustForBias )
    SCV <- adjustScv( SCV, length( sizeFactors ) )

  if( fitType == "local" ) {

    fit <- do.call( "locfit", c(
      list(
        variances ~ do.call( "lp", c( list( log(means) ), lp_extra_args ) ),
        family = "gamma" ),
      locfit_extra_args ) )

    rm( means )
    rm( variances )

    if( adjustForBias )
      ans <- function( q )
        adjustScv(
          pmax( ( predict_helper( fit, log(q) ) - xim * q ) / q^2, 1e-8 ),
          length(sizeFactors) )
    else
      ans <- function( q )
        pmax( ( predict_helper( fit, log(q) ) - xim * q ) / q^2, 1e-8 )

    # Note: The 'pmax' construct above serves to limit the overdispersion to a minimum
    # of 10^-8, which should be indistinguishable from 0 but ensures numerical stability.

  } else if( fitType == "parametric" ) {

    ans <- parametricscvFit( means, SCV )

  } else
    stop( "Unknown fitType." )

  attr( ans, "fitType" ) <- fitType
  list( SCV=SCVAll, SCVfunc=ans )
}


parametricscvFit <- function( means, disps )
{
  coefs <- c( .1, 1 )
  iter <- 0
  while(TRUE) {
    residuals <- disps / ( coefs[1] + coefs[2] / means )
    good <- which( (residuals > 1e-4) & (residuals < 15) )
    fit <- glm( disps[good] ~ I(1/means[good]),
                family=Gamma(link="identity"), start=coefs )
    oldcoefs <- coefs
    coefs <- coefficients(fit)
    if( !all( coefs > 0 ) )
      stop( "Parametric dispersion fit failed. Try a local fit and/or a pooled estimation. (See '?estimateSCV')" )
    if( sum( log( coefs / oldcoefs )^2 ) < 1e-6 )
      break
    iter <- iter + 1
    if( iter > 10 ) {
      warning( "Dispersion fit did not converge." )
      break }
  }

  names( coefs ) <- c( "asymptDisp", "extraPois" )
  ans <- function( q )
    coefs[1] + coefs[2] / q
  attr( ans, "coefficients" ) <- coefs
  ans
}


predict_helper <- function( fit, x )
{
  # A wrapper around predict to avoid the issue that predict.locfit cannot
  # propagate NAs and NaNs properly.

  res <- rep.int( NA_real_, length(x) )
  res[ is.finite(x) ] <- predict( fit, x[is.finite(x)] )
  res
}


XBSeqTestForMatrices <- function( countsA, countsB, sizeFactorsA, sizeFactorsB,
                                  SCVA, SCVB )
{
  kAs <- apply( countsA,1,sum )
  kBs <- apply( countsB ,1,sum)

  mus <- rowMeans( cbind(
    t( t( countsA ) / sizeFactorsA ),
    t( t( countsB ) / sizeFactorsB ) ) )

  signalmuA <- mus*sum(sizeFactorsA)
  signalmuB <- mus*sum(sizeFactorsB)

  signalVarsA <- pmax( mus * sum( sizeFactorsA ) + SCVA * mus^2 * sum(sizeFactorsA^2),
                     mus * sum( sizeFactorsA ) * (1+1e-8) )
  signalVarsB <- pmax( mus * sum( sizeFactorsB ) + SCVB * mus^2 * sum(sizeFactorsB^2),
                     mus * sum( sizeFactorsB ) * (1+1e-8) )

  sapply( seq(along=kAs), function(i) {

    if( kAs[i] == 0 & kBs[i] == 0 )
      return( NA )

    # probability of all possible counts sums with the same total count:
    ks <- 0 : ( kAs[i] + kBs[i] )
    ps <- dnbinom(                   ks, mu = signalmuA[i], size = signalmuA[i]^2/(signalVarsA[i]-signalmuA[i]) ) *
      dnbinom( kAs[i] + kBs[i] - ks, mu = signalmuB[i], size = signalmuB[i]^2/(signalVarsB[i]-signalmuB[i]) )

    # probability of observed count sums:
    pobs <- dnbinom( kAs[i], mu = signalmuA[i], size = signalmuA[i]^2/(signalVarsA[i]-signalmuA[i]) ) *
      dnbinom( kBs[i], mu = signalmuB[i], size = signalmuB[i]^2/(signalVarsB[i]-signalmuB[i]) )

    #stopifnot( na.omit(pobs == ps[ kAs[i]+1 ]) )
    if( kAs[i] * sum( sizeFactorsB ) < kBs[i] * sum( sizeFactorsA ) )
      numer <- ps[ 1 : (kAs[i]+1) ]
    else
      numer <- ps[ (kAs[i]+1) : length(ps) ]
    min( 1, 2 * sum(numer) / sum(ps) )
  } )
}


insertRow <- function(existingDF, newrow, r) {
  for(i in 1:length(r))
    existingDF[seq(r[i]+1,nrow(existingDF)+1),] <- existingDF[seq(r[i],nrow(existingDF)),]
  existingDF[r[i],] <- newrow[i]
  existingDF
}


getSignalVars<-function( counts, bgcounts){
  if(!is.matrix(counts))
    counts <- as.matrix(counts)
  if(!is.matrix(bgcounts))
    bgcounts <- as.matrix(bgcounts)
  bgsizeFactors <- estimateSizeFactorsForMatrix(bgcounts)
  lambda <- rowMeans( sweep(bgcounts, 2, bgsizeFactors, '/' ) )
  observe_param <- getCountParams(counts, estimateSizeFactorsForMatrix(counts))
  observe_sf <- estimateSizeFactorsForMatrix(counts)
  tho <- c()
  for(i in 1:nrow(counts) )
  {
    if(lambda[i]==0){
      temp <- 0
      if ( length( which(counts[i,]==0) ) >= ncol(counts)/2 )
        temp <- NA
    }
    else {
      if(length( which(counts[i,]==0) ) >= ncol(counts)/2 ){
        temp <- NA
      }
      else {
        temp2 <- cor( counts[i,]*observe_sf, bgcounts[i,]*bgsizeFactors )
        if( is.na(temp2))
          temp <- 0
        else
          temp <-temp2
      }
    }
    tho <- c(tho, temp)
  }
  fullvar <- observe_param$baseVar
   signalvar <- fullvar + lambda - 2*tho*sqrt(fullvar)*sqrt(lambda)
  as.matrix(signalvar)
}


prepareScvBiasCorrectionFits <- function( maxnrepl=15, mu=100000, ngenes=10000,
                                          true_raw_scv = c( seq( 0, 2, length.out=100 )[-1], seq( 2, 10, length.out=20 )[-1] ) )
  lapply( 2:maxnrepl, function( m ) {
    est_raw_scv <- sapply( true_raw_scv, function( alpha ) {
      k <- matrix( rnbinom( ngenes*m, mu=mu, size=1/alpha ), ncol=m )
      k <- k[ rowSums(k)>0, ]
      mean( rowVars(k) / rowMeans(k)^2 ) } )
    locfit( true_raw_scv ~ lp( est_raw_scv, nn=.2 ) ) } )


adjustScv <- function( scv, nsamples ) {
  stopifnot( nsamples > 1 )
  if(!exists(("scvBiasCorrectionFits")))
    data( "scvBiasCorrectionFits" )
  if( nsamples - 1 > length( scvBiasCorrectionFits ) )
    scv
  else
    ifelse( scv > .02,
            pmax( predict_helper( scvBiasCorrectionFits[[ nsamples-1 ]], scv ), 1e-8 * scv ),
            scv )   # For scv < .02, our fit is too coarse, but no correction seems necessary anyway
}