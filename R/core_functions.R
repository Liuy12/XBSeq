setGeneric('estimateRealCount', function(object) standardGeneric('estimateRealCount'))


setMethod('estimateRealCount', signature(object = 'XBSeqDataSet'),
          function(object){
            signal <- assay(object, 1) - assay(object, 2)
            signal[signal < 0] <- 0
            assay(object,3) <- signal
            object
          })


setMethod('counts', signature(object = 'XBSeqDataSet'),
          function(object, slot = 3, normalized = FALSE){
            if(length(assays(object)) == 2 & slot == 3)
              stop('Only two assays exist. Call "estimateRealCount" first')
            if (!normalized) {
              return(assay(object, slot))
            } 
            else if (is.null(sizeFactors(object)) | any(is.na(sizeFactors(object)))) {
                stop("first calculate size factors, add normalizationFactors, or set normalized=FALSE")
              } else {
                return(t(t(assay(object,slot)) / sizeFactors(object)))
              }
})

setGeneric("estimateSCV",function(object, ...) standardGeneric('estimateSCV'))


setMethod("estimateSCV", signature(object="XBSeqDataSet"),
          function( object, method = c( "pooled", "per-condition", "blind" ),
                    sharingMode = c( "maximum", "fit-only", "gene-est-only" ),
                    fitType = c("local", "parametric"),
                    locfit_extra_args=list(), lp_extra_args=list(), ... )
          {
            stopifnot( is( object, "XBSeqDataSet" ) )
            if(any(c(is.null(sizeFactors(object)))))
              stop( "NAs found in size factors. Have you called already 'estimateSizeFactors'?" )
            method <- match.arg( method )
            sharingMode <- match.arg( sharingMode )
            fitType <- match.arg( fitType )
            if( length(list(...)) != 0 )
              warning( "in estimateSCV: Ignoring extra argument(s)." )
            if( sharingMode == "gene-est-only" && length(conditions(object))/length(levels(conditions(object))) <= 2 )
              warning( "in estimateSCV: sharingMode=='gene-est-only' will cause inflated numbers of false positives unless you have many replicates." )
            #Remove results from previous fits
            index <- which(! colnames(dispEst(object)) %in% paste( "disp", object@dispTable, sep="_" ))
            if(length(index)){
              dispEst(object) <- dispEst(object, index)
              object@dispTable <- character()
              object@fitInfo = new.env( hash=TRUE )
            }
            if( method == "blind" ) {
              data <- getCountParams(counts(object), sizeFactors(object))
              data_var <- getSignalVars(counts(object, 1), counts(object, 2))
              SCVf <- getSCV(data$baseMean,
                            data_var, sizeFactors(object), fitType, locfit_extra_args, lp_extra_args )
              object@fitInfo[[ "blind" ]] <- list(
                perGeneSCVEsts = SCVf$SCV,
                SCVFunc = SCVf$SCVfunc,
                fittedSCVEsts = SCVf$SCVfunc( data$baseMean ),
                df = ncol(counts(object)) - 1,
                sharingMode = sharingMode )
              a <- rep( "blind", length( levels( conditions(object) ) ) )
              names(a) <- levels( conditions(object) )
              object@dispTable <- a
            } else if( method == "per-condition" ) {
              replicated <- names( which( tapply( conditions(object), conditions(object), length ) > 1 ) )
              if( length( replicated ) < 1 )
                stop( "None of your conditions is replicated. Use method='blind' to estimate across conditions, if you have crossed factors." )
              nonreplicated <- names( which( tapply( conditions(object), conditions(object), length ) == 1 ) )
              overall_basemeans <- rowMeans( counts( object, normalized=TRUE ) )
              for( cond in replicated ) {
                cols <- conditions(object)==cond
                data <- getCountParams(counts(object)[ , cols ], sizeFactors(object)[ cols ] )
                data_var <- getSignalVars(counts(object, 1)[, cols], counts(object, 2)[, cols])
                SCVf <- getSCV( data$baseMean,
                                data_var, sizeFactors(object)[cols], fitType, locfit_extra_args, lp_extra_args )
                object@fitInfo[[ cond ]] <- list(
                  perGeneSCVEsts = SCVf$SCV,
                  SCVFunc = SCVf$SCVfunc,
                  fittedSCVEsts = SCVf$SCVfunc( overall_basemeans ),     # Note that we do not use bmv$baseMean here
                  df = sum(cols) - 1,
                  sharingMode = sharingMode ) }
              
              object@dispTable <- sapply( levels(conditions(object)), function( cond )
                ifelse( cond %in% replicated, cond, "max" ) )
              
            } else if( method == "pooled" ) {
              conds <- conditions(object)
              if( !any( duplicated( conds ) ) )
                stop( "None of your conditions is replicated. Use method='blind' to estimate across conditions, or 'pooled-CR', if you have crossed factors." )
              data <- getCountParamsPooled( counts(object), sizeFactors(object), conds )
              baseMeans <- data$baseMean
              data_var <- getSignalVars(counts(object, 1), counts(object, 2))
              SCVf <- getSCV(data$baseMean,
                            data$baseVar, sizeFactors(object), fitType, locfit_extra_args, lp_extra_args )
              df <- ncol(counts(object)) - length(unique(conds))
              object@fitInfo[[ "pooled" ]] <- list(
                perGeneSCVEsts = SCVf$SCV,
                SCVFunc = SCVf$SCVfunc,
                fittedSCVEsts = SCVf$SCVfunc( baseMeans ),
                df = df,
                sharingMode = sharingMode )
              a <- rep( "pooled", length( levels( conditions(object) ) ) )
              names(a) <- levels( conditions(object) )
              dispTable(object) = a
            } else
              stop(sprintf("Invalid method '%s'.", method))
            for( n in ls(object@fitInfo) )
              dispEst(object, paste( "disp", n, sep="_" ) ) <-
              switch( sharingMode,
                      `fit-only`      = object@fitInfo[[ n ]]$fittedSCVEsts,
                      `gene-est-only` = {
                        a <- object@fitInfo[[ n ]]$perGeneSCVEsts
                        a[ is.nan(a) ] <- 0
                        pmax( a, 1e-8 ) },
                      `maximum`       = pmax( object@fitInfo[[ n ]]$fittedSCVEsts, object@fitInfo[[ n ]]$perGeneSCVEsts, na.rm=TRUE ),
                      stop(sprintf("Invalid sharingMode '%s'.", sharingMode))
              ) ## switch
            
            if( "max" %in% object@dispTable )
              dispEst(object, "disp_max") <- do.call( pmax,
                                                c( dispEst(object, which(colnames(dispEst(object)) %in% paste( "disp", object@dispTable, sep="_" )) ), na.rm=TRUE ) )
            
            validObject( object )
            object
          })


XBSeqTest <- function( XB, condA, condB, pvals_only=FALSE )
{
  stopifnot( is( XB, "XBSeqDataSet" ) )
  if( all( is.na( dispTable(XB) ) ) )
    stop( "Call 'estimateSCV' first." )
  if( dispTable(XB)[condA] == "blind") {
    if( fitInfo( XB, "blind" )$sharingMode != "fit-only" )
      warning( 'You have used \'method="blind"\' in estimateSCV without also setting \'sharingMode="fit-only"\'. This will not yield useful results.' )
  }
  stopifnot( condA %in% levels(conditions(XB)) )
  stopifnot( condB %in% levels(conditions(XB)) )
  colA <- conditions(XB)==condA
  colB <- conditions(XB)==condB
  
  rawScvA <- dispEst(XB, paste( "disp", dispTable(XB)[condA], sep="_" ))
  rawScvB <- dispEst(XB, paste( "disp", dispTable(XB)[condB], sep="_" ))
  
  pval <- XBSeqTestForMatrices(
    counts(XB)[,colA],
    counts(XB)[,colB],
    sizeFactors(XB)[colA],
    sizeFactors(XB)[colB],
    rawScvA,
    rawScvB )
  
  if( pvals_only )
    pval
  else {
    data <- getCountParams( counts(XB), sizeFactors(XB)[colA|colB] )
    dataA <- getCountParams( counts(XB)[,colA], sizeFactors(XB)[colA] )
    dataB <- getCountParams( counts(XB)[,colB], sizeFactors(XB)[colB] )
    data.frame(
      id    = rownames( counts(XB) ),
      baseMean  = data$baseMean,
      baseMeanA = dataA$baseMean,
      baseMeanB = dataB$baseMean,
      foldChange = dataB$baseMean / dataA$baseMean,
      log2FoldChange = log2( dataB$baseMean / dataA$baseMean ),
      pval = pval,
      padj = p.adjust( pval, method="BH" ),
      stringsAsFactors = FALSE ) }
}


# a wrapper function once and for all
XBSeq <- function(counts, bgcounts, conditions, method='pooled', sharingMode='maximum', fitType='local', pvals_only=FALSE ){
  if(!is.factor(conditions))
    conditions <- as.factor(conditions)
  XB <- XBSeqDataSet(counts, bgcounts, conditions)
  XB <- estimateRealCount(XB)
  XB <- estimateSizeFactors(XB)
  XB <- estimateSCV(XB, method=method, sharingMode=sharingMode, fitType=fitType)
  Teststas <- XBSeqTest(XB, levels(conditions)[1L], levels(conditions)[2L], pvals_only=pvals_only)
  Teststas
}