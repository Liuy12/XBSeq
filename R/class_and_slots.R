setClass( "XBSeqDataSet", 
          contains = "DESeqDataSet",
          slots = c( 
            fitInfo = "environment",
            dispTable = "character",
            conditions = 'factor',
            dispEst = 'list')
)


setValidity( "XBSeqDataSet", function( object ) {
  if(!is.factor(object@conditions))
    return("conditions have to be factors")
  TRUE
} )


setGeneric("fitInfo", function(object, name=NULL) standardGeneric("fitInfo"))


setMethod('fitInfo', signature(object = "XBSeqDataSet"), 
          function( object, name){
            if( length( ls( object@fitInfo ) ) == 0 )
              stop( "No fits available. Call 'estimateSCV' first." )
            if( length( ls( object@fitInfo ) ) > 1 && is.null(name) )
              stop( "More than one fitInfo object available. Specify by name. (See 'ls(XB@fitInfo)' for a list.)" )
            if( length( ls( object@fitInfo ) ) == 1 && is.null(name) )
              name = ls( object@fitInfo )[ 1 ]
            object@fitInfo[[ name]]
          }
)


setMethod("conditions", signature(object="XBSeqDataSet"),
          function( object, ... ) {
            if(length(list(...))!=0)
              warning("in conditions: Ignoring second and/or further arguments.")
            conds <- object@conditions
            names( conds ) <- colnames( counts(object,1) )
            conds
          })   


setReplaceMethod("conditions", signature(object="XBSeqDataSet"),
                 function( object, value ) {
                   object@conditions <- factor( value )
                   validObject( object )
                   object
                 })


setMethod("dispTable", signature(object="XBSeqDataSet"),
          function( object ) {
            object@dispTable
          })   


setReplaceMethod("dispTable", signature(object="XBSeqDataSet"),
                 function( object, value ) {
                   object@dispTable <- value
                   validObject( object )
                   object
                 })   


setGeneric("dispEst", function(object, varname = NA) standardGeneric("dispEst"))


setGeneric("dispEst<-", function(object, varname = NA, value) standardGeneric("dispEst<-"))


setMethod("dispEst", signature(object="XBSeqDataSet"),
          function(object, varname) {
            if(!is.na(varname))
              object@dispEst[[varname]]
            else
              object@dispEst
          })


setReplaceMethod("dispEst", signature(object="XBSeqDataSet"),
                 function(object, varname,  value) {
                   if(!is.na(varname))
                     object@dispEst[[varname]] <- value
                   else
                     object@dispEst <- value
                   validObject( object )
                   object
                 })


XBSeqDataSet <- function(counts, bgcounts, conditions, sizeFactors=NULL, ...)
{
  counts <- as.matrix(counts)
  bgcounts <-  as.matrix(bgcounts)
  conditions <- as.factor(conditions)
  if(any(round(counts) != counts) | any(round(bgcounts) != bgcounts))
    stop("The input data have to be integer!")
  mode( counts ) <- "integer"
  mode( bgcounts ) <- "integer"
  if( nrow( counts )!= nrow( bgcounts ) ){
    MissedRecord <- which( rownames( counts) %in% setdiff( rownames( counts ),rownames( bgcounts ) ) )
    bgcounts <- insertRow( bgcounts, repmat(apply(bgcounts,2,mean),length(MissedRecord),1) ,MissedRecord)
  }
  if( is.null( sizeFactors ) ) {
    sizeFactors <- rep( NA_real_, ncol(counts) )
  }
  assays <- list(counts=counts, bgcounts=bgcounts)
  colData <- data.frame(conditions = conditions)
  rownames(colData) <- colnames(counts)
  colData <- DataFrame(colData, row.names=rownames(colData))
  se <- SummarizedExperiment(assays, colData = colData)
  stopifnot(length(conditions) == ncol(counts))
  rvft <- rep(NA_character_, length(levels(conditions)))
  XB <- DESeqDataSet(se, formula(~conditions))
  XB <- new("XBSeqDataSet",
            XB, 
            fitInfo = new.env(hash=TRUE),
            dispTable = rvft,
            conditions = conditions)
  return(XB)
}
