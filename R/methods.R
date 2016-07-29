getCountParams <- function(counts, sizeFactors) {
  # Divides the counts by sizeFactors and calculates the estimates for
  # base means and variances for each gene
  
  data.frame(baseMean = rowMeans(t(t(counts) / sizeFactors)),
             baseVar = rowVars(t(t(counts) / sizeFactors)))
}


getCountParamsPooled <-
  function(counts, sizeFactors, conditions) {
    basecounts <- t(t(counts) / sizeFactors)
    replicated_sample <-
      conditions %in% names(which(table(conditions) > 1))
    df <-
      sum(replicated_sample) - length(unique(conditions[replicated_sample]))
    
    data.frame(baseMean = rowMeans(basecounts),
               baseVar =
                 rowSums(sapply(
                   tapply((seq_len(ncol(
                     counts
                   )))[replicated_sample],
                   factor(conditions[replicated_sample]),
                   function(cols)
                     rowSums((
                       basecounts[,cols] - rowMeans(basecounts[,cols])
                     ) ^ 2)),
                   identity
                 )) / df)
  }

getSCV <- function(means,
                   variances, sizeFactors, fitType = c("parametric", "local"),
                   locfit_extra_args = list(), lp_extra_args = list(), adjustForBias =
                     TRUE) {
  fitType <- match.arg(fitType)
  
  xim <- mean(1 / sizeFactors)
  SCVAll <- (variances - xim * means) / means ^ 2
  
  variances <- variances[means > 0]
  SCV <- SCVAll[means > 0]
  means <- means[means > 0]
  
  if (adjustForBias)
    SCV <- adjustScv(SCV, length(sizeFactors))
  
  if (fitType == "local") {
    fit <- do.call("locfit", c(
      list(variances ~ do.call("lp", c(
        list(log(means)), lp_extra_args
      )),
      family = "gamma"),
      locfit_extra_args
    ))
    
    rm(means)
    rm(variances)
    
    if (adjustForBias)
      ans <- function(q)
        adjustScv(pmax((predict_helper(fit, log(
          q
        )) - xim * q) / q ^ 2, 1e-8),
        length(sizeFactors))
    else
      ans <- function(q)
        pmax((predict_helper(fit, log(q)) - xim * q) / q ^ 2, 1e-8)
    
    # Note: The 'pmax' construct above serves to limit the overdispersion to a minimum
    # of 10^-8, which should be indistinguishable from 0 but ensures numerical stability.
    
  } else if (fitType == "parametric") {
    ans <- parametricscvFit(means, SCV)
    
  } else
    stop("Unknown fitType.")
  
  attr(ans, "fitType") <- fitType
  list(SCV = SCVAll, SCVfunc = ans)
}


parametricscvFit <- function(means, disps)
{
  coefs <- c(.1, 1)
  iter <- 0
  while (TRUE) {
    residuals <- disps / (coefs[1] + coefs[2] / means)
    good <- which((residuals > 1e-4) & (residuals < 15))
    fit <- glm(disps[good] ~ I(1 / means[good]),
               family = Gamma(link = "identity"), start = coefs)
    oldcoefs <- coefs
    coefs <- coefficients(fit)
    if (!all(coefs > 0))
      stop(
        "Parametric dispersion fit failed. Try a local fit and/or a pooled estimation. (See '?estimateSCV')"
      )
    if (sum(log(coefs / oldcoefs) ^ 2) < 1e-6)
      break
    iter <- iter + 1
    if (iter > 10) {
      warning("Dispersion fit did not converge.")
      break
    }
  }
  
  names(coefs) <- c("asymptDisp", "extraPois")
  ans <- function(q)
    coefs[1] + coefs[2] / q
  attr(ans, "coefficients") <- coefs
  ans
}


predict_helper <- function(fit, x)
{
  # A wrapper around predict to avoid the issue that predict.locfit cannot
  # propagate NAs and NaNs properly.
  
  res <- rep.int(NA_real_, length(x))
  res[is.finite(x)] <- predict(fit, x[is.finite(x)])
  res
}


XBSeqTestForMatrices <-
  function(countsA, countsB, bgcountsA, bgcountsB, sizeFactorsA, sizeFactorsB,
           SCVA, SCVB , method = c('NP', 'MLE'), big_count)
  {
    method <- match.arg(method, c('NP', 'MLE'))
    if (ncol(countsA) < 5 & method == 'MLE')
      warning(
        'Non-parametric estimation method is recommended for experiments with replicates smaller than 5'
      )
    
    kAs <- apply(countsA, 1, sum)
    kBs <- apply(countsB, 1, sum)
    
    mus <- rowMeans(cbind(t(t(countsA) / sizeFactorsA),
                          t(t(countsB) / sizeFactorsB)))
    
    if (method == 'NP') {
      signalmuA <- mus * sum(sizeFactorsA)
      signalmuB <- mus * sum(sizeFactorsB)
      
      signalVarsA <-
        pmax(
          mus * sum(sizeFactorsA) + SCVA * mus ^ 2 * sum(sizeFactorsA ^ 2),
          mus * sum(sizeFactorsA) * (1 + 1e-8)
        )
      signalVarsB <-
        pmax(
          mus * sum(sizeFactorsB) + SCVB * mus ^ 2 * sum(sizeFactorsB ^ 2),
          mus * sum(sizeFactorsB) * (1 + 1e-8)
        )
      
      sizeA <- signalmuA ^ 2 / (signalVarsA - signalmuA)
      sizeB <- signalmuB ^ 2 / (signalVarsB - signalmuB)
    }
    else{
      musA <- mus * mean(sizeFactorsA)
      musB <- mus * mean(sizeFactorsB)
      
      VarsA <-
        pmax(
          mus * mean(sizeFactorsA) + SCVA * mus ^ 2 * mean(sizeFactorsA ^ 2),
          mus * mean(sizeFactorsA) * (1 + 1e-8)
        )
      VarsB <-
        pmax(
          mus * mean(sizeFactorsB) + SCVB * mus ^ 2 * mean(sizeFactorsB ^ 2),
          mus * mean(sizeFactorsB) * (1 + 1e-8)
        )
      
      lambda <- rowMeans(cbind(t(t(bgcountsA) / sizeFactorsA),
                               t(t(bgcountsB) / sizeFactorsB)))
      
      lambdaA <- lambda * mean(sizeFactorsA)
      lambdaB <- lambda * mean(sizeFactorsB)
      
      
      ParamsA <-
        sapply(1:nrow(countsA), function(i)
          estimation_param_PoissonNB_MLE(
            countsA[i,] + bgcounts[i,],
            bgcountsA[i,],
            musA[i] ^
              2 / (VarsA[i] - musA[i]),
            (VarsA[i] -
               musA[i]) / musA[i],
            lambdaA[i]
          ))
      ParamsA <- matrix(unlist(ParamsA), ncol = 3, byrow = TRUE)
      ParamsB <-
        sapply(1:nrow(countsB), function(i)
          estimation_param_PoissonNB_MLE(
            countsB[i,] + bgcounts[i,],
            bgcountsB[i,],
            musB[i] ^
              2 / (VarsB[i] - musB[i]),
            (VarsB[i] -
               musB[i]) / musB[i],
            lambdaB[i]
          ))
      ParamsB <- matrix(unlist(ParamsB), ncol = 3, byrow = TRUE)
      
      sizeA <- ParamsA[,1] * sum(sizeFactorsA)
      sizeB <- ParamsB[,1] * sum(sizeFactorsB)
      
      signalmuA <- ParamsA[,1] * ParamsB[,2] * sum(sizeFactorsA)
      signalmuB <- ParamsB[,1] * ParamsB[,2] * sum(sizeFactorsB)
    }
    
    pvals <- rep(1,nrow(countsA))
    names(pvals) <- rownames(countsA)
    
    big <- kAs > big_count & kBs > big_count
    if (any(big))
      pvals[big] <-
      exactTestBetaApprox(
        sweep(countsA[big,,drop = FALSE], 2, sizeFactorsA, '/'), sweep(countsB[big,,drop =
                                                                                 FALSE], 2, sizeFactorsB, '/'), sizeA[big], sizeB[big]
      )
    if (any(!big))
      pvals[!big] <- sapply(seq(along = kAs[!big]), function(i) {
        if (kAs[!big][i] == 0 & kBs[!big][i] == 0)
          return(NA)
        
        # probability of all possible counts sums with the same total count:
        ks <- 0:(kAs[!big][i] + kBs[!big][i])
        ps <-
          dnbinom(ks, mu = signalmuA[!big][i], size = sizeA[!big][i]) *
          dnbinom(kAs[!big][i] + kBs[!big][i] - ks, mu = signalmuB[!big][i], size = sizeB[!big][i])
        
        # probabilit y of observed count sums:
        pobs <-
          dnbinom(kAs[!big][i], mu = signalmuA[!big][i], size = sizeA[!big][i]) *
          dnbinom(kBs[!big][i], mu = signalmuB[!big][i], size = sizeB[!big][i])
        
        if (kAs[!big][i] * sum(sizeFactorsB) < kBs[!big][i] * sum(sizeFactorsA))
          numer <- ps[1:(kAs[!big][i] + 1)]
        else
          numer <- ps[(kAs[!big][i] + 1):length(ps)]
        min(1, 2 * sum(numer) / sum(ps))
      })
    return(pvals)
  }


exactTestBetaApprox <- function(y1,y2,size1, size2)
  #	Test for differences in means between two negative binomial
  #	or Poisson random variables, or between two groups of variables,
  #	using a beta distribution approximation.
  #	Test is naturally conditional on total sum.
  #	Left and right rejection regions have equal probability.
  
  # adopted and revised from source code of edgeR
{
  #	Convert matrices to vectors
  ntags <- NROW(y1)
  n1 <- NCOL(y1)
  n2 <- NCOL(y2)
  if(n1>1) y1 <- rowSums(y1)
  if(n2>1) y2 <- rowSums(y2)
  
  #	Null fitted values
  y <- y1+y2
  mu <- y/(n1+n2)
  
  #	Compute p-values
  pvals <- rep(1,ntags)
  all.zero <- y<=0
  alpha1 <- n1*mu/(1+n1/size1*mu)
  alpha2 <- n2*mu/(1+n2/size2*mu)
  med <- rep(0,ntags)
  med[!all.zero] <- qbeta(0.5,alpha1[!all.zero],alpha2[!all.zero])
  left <- (y1+0.5)/y<med & !all.zero
  if(any(left)) {
    pvals[left] <- 2*pbeta((y1[left]+0.5)/y[left],alpha1[left],alpha2[left])
  }
  right <- (y1-0.5)/y>med & !all.zero
  if(any(right)) {
    pvals[right] <- 2*pbeta((y1[right]-0.5)/y[right],alpha1[right],alpha2[right],lower.tail=FALSE)
  }
  names(pvals) <- names(y1)
  pvals
}


insertRow <- function(existingDF, newrow, r) {
  existingDF <- as.data.frame(existingDF)
  for (i in 1:length(r)){
    if(r[i] < nrow(existingDF) | r[i] == nrow(existingDF))
      existingDF[seq(r[i] + 1,nrow(existingDF) + 1),] <-
        existingDF[seq(r[i],nrow(existingDF)),]
    existingDF[r[i],] <- newrow[i,]
  }
  existingDF <- as.matrix(existingDF)
  storage.mode(existingDF) <- 'integer'
  existingDF
}


getSignalVars <- function(counts, bgcounts) {
  if (!is.matrix(counts))
    counts <- as.matrix(counts)
  if (!is.matrix(bgcounts))
    bgcounts <- as.matrix(bgcounts)
  bgsizeFactors <- estimateSizeFactorsForMatrix(bgcounts)
  lambda <- rowMeans(sweep(bgcounts, 2, bgsizeFactors, '/'))
  observe_param <-
    getCountParams(counts, estimateSizeFactorsForMatrix(counts))
  observe_sf <- estimateSizeFactorsForMatrix(counts)
  tho <- c()
  for (i in 1:nrow(counts))
  {
    if (lambda[i] == 0) {
      temp <- 0
      if (length(which(counts[i,] == 0)) >= ncol(counts) / 2)
        temp <- NA
    }
    else {
      if (length(which(counts[i,] == 0)) >= ncol(counts) / 2) {
        temp <- NA
      }
      else {
        temp2 <- cor(counts[i,] * observe_sf, bgcounts[i,] * bgsizeFactors)
        if (is.na(temp2))
          temp <- 0
        else
          temp <- temp2
      }
    }
    tho <- c(tho, temp)
  }
  fullvar <- observe_param$baseVar
  signalvar <- fullvar + lambda - 2 * tho * sqrt(fullvar) * sqrt(lambda)
  as.matrix(signalvar)
}


prepareScvBiasCorrectionFits <-
  function(maxnrepl = 15, mu = 100000, ngenes = 10000,
           true_raw_scv = c(seq(0, 2, length.out =
                                  100)[-1], seq(2, 10, length.out = 20)[-1]))
    lapply(2:maxnrepl, function(m) {
      est_raw_scv <- sapply(true_raw_scv, function(alpha) {
        k <- matrix(rnbinom(ngenes * m, mu = mu, size = 1 / alpha), ncol = m)
        k <- k[rowSums(k) > 0,]
        mean(rowVars(k) / rowMeans(k) ^ 2)
      })
      locfit(true_raw_scv ~ lp(est_raw_scv, nn = .2))
    })


adjustScv <- function(scv, nsamples) {
  stopifnot(nsamples > 1)
  if (!exists(("scvBiasCorrectionFits")))
    data("scvBiasCorrectionFits")
  if (nsamples - 1 > length(scvBiasCorrectionFits))
    scv
  else
    ifelse(scv > .02,
           pmax(predict_helper(scvBiasCorrectionFits[[nsamples - 1]], scv), 1e-8 * scv),
           scv)   # For scv < .02, our fit is too coarse, but no correction seems necessary anyway
}


Loglikhood <- function(counts, bgcounts) {
  function(para) {
    alpha <- para[1]
    beta <- para[2]
    lambda <- para[3]
    - sum(
      ddelap(
        counts, alpha = alpha, beta = beta, lambda = lambda, log = TRUE
      ), dpois(bgcounts, lambda = lambda, log = TRUE)
    )
  }
}


estimation_param_PoissonNB_MLE <-
  function(counts, bgcounts, alpha, beta, lambda) {
    if (any(is.na(c(alpha, beta, lambda))))
      list(alpha = alpha, beta = beta, lambda = lambda)
    else{
      mle <-
        try(optim(
          c(alpha, beta, lambda), Loglikhood(counts, bgcounts),
          method = 'L-BFGS-B', lower = c(0,0,0)
        ), silent = TRUE)
      if (class(mle) == 'try-error')
        list(alpha = alpha, beta = beta, lambda = lambda)
      else{
        mle <- optim(
          c(alpha, beta, lambda), Loglikhood(counts, bgcounts),
          method = 'L-BFGS-B', lower = c(0,0,0)
        )
        if (mle$convergence > 0 | any(is.na(mle$par)))
          list(alpha = alpha, beta = beta, lambda = lambda)
        else
          list(
            alpha = mle$par[1], beta = mle$par[2], lambda = mle$par[3]
          )
      }
    }
  }