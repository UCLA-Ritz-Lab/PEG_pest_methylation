makeContrasts <- function (..., contrasts = NULL, levels) 
{
  e <- substitute(list(...))
  if (is.factor(levels)) 
    levels <- levels(levels)
  if (!is.character(levels)) 
    levels <- colnames(levels)
  if (levels[1] == "(Intercept)") {
    levels[1] <- "Intercept"
    warning("Renaming (Intercept) to Intercept")
  }
  notvalid <- (levels != make.names(levels))
  if (any(notvalid)) 
    stop("The levels must by syntactically valid names in R, see help(make.names).  Non-valid names: ", 
         paste(levels[notvalid], collapse = ","))
  n <- length(levels)
  if (n < 1) 
    stop("No levels to construct contrasts from")
  indicator <- function(i, n) {
    out <- rep(0, n)
    out[i] <- 1
    out
  }
  levelsenv <- new.env()
  for (i in 1:n) assign(levels[i], indicator(i, n), pos = levelsenv)
  if (!is.null(contrasts)) {
    if (length(e) > 1) 
      stop("Can't specify both ... and contrasts")
    e <- as.character(contrasts)
    ne <- length(e)
    cm <- matrix(0, n, ne, dimnames = list(Levels = levels, 
                                           Contrasts = e))
    if (ne == 0) 
      return(cm)
    for (j in 1:ne) {
      ej <- parse(text = e[j])
      cm[, j] <- eval(ej, envir = levelsenv)
    }
    return(cm)
  }
  ne <- length(e)
  enames <- names(e)[2:ne]
  easchar <- as.character(e)[2:ne]
  if (is.null(enames)) 
    cn <- easchar
  else cn <- ifelse(enames == "", easchar, enames)
  cm <- matrix(0, n, ne - 1, dimnames = list(Levels = levels, 
                                             Contrasts = cn))
  if (ne < 2) 
    return(cm)
  for (j in 1:(ne - 1)) {
    ej <- e[[j + 1]]
    if (is.character(ej)) 
      ej <- parse(text = ej)
    ej <- eval(ej, envir = levelsenv)
    if (!is.numeric(ej)) {
      colnames(cm)[j] <- as.character(ej)[1]
      if (is.character(ej)) 
        ej <- parse(text = ej)
      ej <- eval(ej, envir = levelsenv)
    }
    cm[, j] <- ej
  }
  cm
}


lmFit <- function (object, design = NULL, ndups = NULL, spacing = NULL, 
          block = NULL, correlation, weights = NULL, method = "ls", 
          ...) 
{
  if (inherits(object, "data.frame")) {
    ColumnIsNumeric <- vapply(object, is.numeric, FUN.VALUE = TRUE)
    if (all(ColumnIsNumeric)) {
      y <- list(exprs = as.matrix(object))
      y$Amean <- rowMeans(y$exprs, na.rm = TRUE)
    }
    else {
      WhichNotNumeric <- which(!ColumnIsNumeric)
      if (identical(sum(WhichNotNumeric), 1L) && length(ColumnIsNumeric) > 
          1L) {
        y <- list()
        y$exprs <- as.matrix(object[, -1, drop = FALSE])
        y$probes <- object[, 1, drop = FALSE]
        y$Amean <- rowMeans(y$exprs, na.rm = TRUE)
        message("Converting data.frame to matrix, treating first column as gene IDs.")
      }
      else {
        stop("Expression object should be numeric, instead it is a data.frame with ", 
             length(WhichNotNumeric), " non-numeric columns")
      }
    }
  }
  else {
    y <- getEAWP(object)
  }
  if (!nrow(y$exprs)) 
    stop("expression matrix has zero rows")
  if (is.null(design)) 
    design <- y$design
  if (is.null(design)) 
    design <- matrix(1, ncol(y$exprs), 1)
  else {
    design <- as.matrix(design)
    if (!identical(mode(design), "numeric")) 
      stop("design must be a numeric matrix")
    if (!identical(nrow(design), ncol(y$exprs))) 
      stop("row dimension of design doesn't match column dimension of data object")
    if (anyNA(design)) 
      stop("NAs not allowed in design matrix")
  }
  ne <- nonEstimable(design)
  if (!is.null(ne)) 
    cat("Coefficients not estimable:", paste(ne, collapse = " "), 
        "\n")
  if (is.null(ndups)) 
    ndups <- y$printer$ndups
  if (is.null(ndups)) 
    ndups <- 1
  if (is.null(spacing)) 
    spacing <- y$printer$spacing
  if (is.null(spacing)) 
    spacing <- 1
  if (is.null(weights)) 
    weights <- y$weights
  method <- match.arg(method, c("ls", "robust"))
  if (ndups > 1) {
    if (!is.null(y$probes)) 
      y$probes <- uniquegenelist(y$probes, ndups = ndups, 
                                 spacing = spacing)
    if (!is.null(y$Amean)) 
      y$Amean <- rowMeans(unwrapdups(as.matrix(y$Amean), 
                                     ndups = ndups, spacing = spacing), na.rm = TRUE)
  }
  if (method == "robust") 
    fit <- mrlm(y$exprs, design = design, ndups = ndups, 
                spacing = spacing, weights = weights, ...)
  else if (ndups < 2 && is.null(block)) 
    fit <- lm.series(y$exprs, design = design, ndups = ndups, 
                     spacing = spacing, weights = weights)
  else {
    if (missing(correlation)) 
      stop("the correlation must be set, see duplicateCorrelation")
    fit <- gls.series(y$exprs, design = design, ndups = ndups, 
                      spacing = spacing, block = block, correlation = correlation, 
                      weights = weights, ...)
  }
  if (NCOL(fit$coefficients) > 1) {
    n <- rowSums(is.na(fit$coefficients))
    n <- sum(n > 0 & n < NCOL(fit$coefficients))
    if (n > 0) 
      warning("Partial NA coefficients for ", n, " probe(s)", 
              call. = FALSE)
  }
  fit$genes <- y$probes
  fit$Amean <- y$Amean
  fit$method <- method
  fit$design <- design
  new("MArrayLM", fit)
}

nonEstimable <- function (x) 
{
  x <- as.matrix(x)
  p <- ncol(x)
  QR <- qr(x)
  if (QR$rank < p) {
    n <- colnames(x)
    if (is.null(n)) 
      n <- as.character(1:p)
    notest <- n[QR$pivot[(QR$rank + 1):p]]
    blank <- notest == ""
    if (any(blank)) 
      notest[blank] <- as.character(((QR$rank + 1):p)[blank])
    return(notest)
  }
  else {
    return(NULL)
  }
}

lm.series <- function (M, design = NULL, ndups = 1, spacing = 1, weights = NULL) 
{
  M <- as.matrix(M)
  narrays <- ncol(M)
  if (is.null(design)) 
    design <- matrix(1, narrays, 1)
  else design <- as.matrix(design)
  nbeta <- ncol(design)
  coef.names <- colnames(design)
  if (is.null(coef.names)) 
    coef.names <- paste("x", 1:nbeta, sep = "")
  if (!is.null(weights)) {
    weights <- asMatrixWeights(weights, dim(M))
    weights[weights <= 0] <- NA
    M[!is.finite(weights)] <- NA
  }
  if (ndups > 1) {
    M <- unwrapdups(M, ndups = ndups, spacing = spacing)
    design <- design %x% rep_len(1, ndups)
    if (!is.null(weights)) 
      weights <- unwrapdups(weights, ndups = ndups, spacing = spacing)
  }
  ngenes <- nrow(M)
  stdev.unscaled <- beta <- matrix(NA, ngenes, nbeta, dimnames = list(rownames(M), 
                                                                      coef.names))
  NoProbeWts <- all(is.finite(M)) && (is.null(weights) || 
                                        !is.null(attr(weights, "arrayweights")))
  if (NoProbeWts) {
    if (is.null(weights)) 
      fit <- lm.fit(design, t(M))
    else {
      fit <- lm.wfit(design, t(M), weights[1, ])
      fit$weights <- NULL
    }
    if (fit$df.residual > 0) {
      if (is.matrix(fit$effects)) 
        fit$sigma <- sqrt(colMeans(fit$effects[(fit$rank + 
                                                  1):narrays, , drop = FALSE]^2))
      else fit$sigma <- sqrt(mean(fit$effects[(fit$rank + 
                                                 1):narrays]^2))
    }
    else fit$sigma <- rep_len(NA_real_, ngenes)
    fit$fitted.values <- fit$residuals <- fit$effects <- NULL
    fit$coefficients <- t(fit$coefficients)
    fit$cov.coefficients <- chol2inv(fit$qr$qr, size = fit$qr$rank)
    est <- fit$qr$pivot[1:fit$qr$rank]
    dimnames(fit$cov.coefficients) <- list(coef.names[est], 
                                           coef.names[est])
    stdev.unscaled[, est] <- matrix(sqrt(diag(fit$cov.coefficients)), 
                                    ngenes, fit$qr$rank, byrow = TRUE)
    fit$stdev.unscaled <- stdev.unscaled
    fit$df.residual <- rep_len(fit$df.residual, length.out = ngenes)
    dimnames(fit$stdev.unscaled) <- dimnames(fit$stdev.unscaled) <- dimnames(fit$coefficients)
    fit$pivot <- fit$qr$pivot
    return(fit)
  }
  beta <- stdev.unscaled
  sigma <- rep_len(NA_real_, ngenes)
  df.residual <- rep_len(0, ngenes)
  for (i in 1:ngenes) {
    y <- as.vector(M[i, ])
    obs <- is.finite(y)
    if (sum(obs) > 0) {
      X <- design[obs, , drop = FALSE]
      y <- y[obs]
      if (is.null(weights)) 
        out <- lm.fit(X, y)
      else {
        w <- as.vector(weights[i, obs])
        out <- lm.wfit(X, y, w)
      }
      est <- !is.na(out$coefficients)
      beta[i, ] <- out$coefficients
      stdev.unscaled[i, est] <- sqrt(diag(chol2inv(out$qr$qr, 
                                                   size = out$rank)))
      df.residual[i] <- out$df.residual
      if (df.residual[i] > 0) 
        sigma[i] <- sqrt(mean(out$effects[-(1:out$rank)]^2))
    }
  }
  QR <- qr(design)
  cov.coef <- chol2inv(QR$qr, size = QR$rank)
  est <- QR$pivot[1:QR$rank]
  dimnames(cov.coef) <- list(coef.names[est], coef.names[est])
  list(coefficients = beta, stdev.unscaled = stdev.unscaled, 
       sigma = sigma, df.residual = df.residual, cov.coefficients = cov.coef, 
       pivot = QR$pivot, rank = QR$rank)
}


contrasts.fit <- function (fit, contrasts = NULL, coefficients = NULL) 
{
  if (identical(is.null(contrasts), is.null(coefficients))) 
    stop("Must specify exactly one of contrasts or coefficients")
  if (!is.null(coefficients)) 
    return(fit[, coefficients])
  if (is.null(fit$coefficients)) 
    stop("fit must contain coefficients component")
  if (is.null(fit$stdev.unscaled)) 
    stop("fit must contain stdev.unscaled component")
  fit$t <- NULL
  fit$p.value <- NULL
  fit$lods <- NULL
  fit$F <- NULL
  fit$F.p.value <- NULL
  ncoef <- ncol(fit$coefficients)
  if (!is.numeric(contrasts)) 
    stop("contrasts must be a numeric matrix")
  if (anyNA(contrasts)) 
    stop("NAs not allowed in contrasts")
  contrasts <- as.matrix(contrasts)
  if (!identical(nrow(contrasts), ncoef)) 
    stop("Number of rows of contrast matrix must match number of coefficients in fit")
  rn <- rownames(contrasts)
  cn <- colnames(fit$coefficients)
  if (!is.null(rn) && !is.null(cn) && !identical(rn, cn)) 
    warning("row names of contrasts don't match col names of coefficients")
  fit$contrasts <- contrasts
  if (!ncol(contrasts)) 
    return(fit[, 0])
  if (is.null(fit$cov.coefficients)) {
    warning("cov.coefficients not found in fit - assuming coefficients are orthogonal", 
            call. = FALSE)
    var.coef <- colMeans(fit$stdev.unscaled^2)
    fit$cov.coefficients <- diag(var.coef, nrow = ncoef)
    cormatrix <- diag(nrow = ncoef)
    orthog <- TRUE
  }
  else {
    cormatrix <- cov2cor(fit$cov.coefficients)
    if (length(cormatrix) < 2) {
      orthog <- TRUE
    }
    else {
      orthog <- sum(abs(cormatrix[lower.tri(cormatrix)])) < 
        1e-12
    }
  }
  r <- nrow(cormatrix)
  if (r < ncoef) {
    if (is.null(fit$pivot)) 
      stop("cor.coef not full rank but pivot column not found in fit")
    est <- fit$pivot[1:r]
    if (any(contrasts[-est, ] != 0)) 
      stop("trying to take contrast of non-estimable coefficient")
    contrasts <- contrasts[est, , drop = FALSE]
    fit$coefficients <- fit$coefficients[, est, drop = FALSE]
    fit$stdev.unscaled <- fit$stdev.unscaled[, est, drop = FALSE]
    ncoef <- r
  }
  ContrastsAllZero <- which(rowSums(abs(contrasts)) == 0)
  if (length(ContrastsAllZero)) {
    contrasts <- contrasts[-ContrastsAllZero, , drop = FALSE]
    fit$coefficients <- fit$coefficients[, -ContrastsAllZero, 
                                         drop = FALSE]
    fit$stdev.unscaled <- fit$stdev.unscaled[, -ContrastsAllZero, 
                                             drop = FALSE]
    fit$cov.coefficients <- fit$cov.coefficients[-ContrastsAllZero, 
                                                 -ContrastsAllZero, drop = FALSE]
    cormatrix <- cormatrix[-ContrastsAllZero, -ContrastsAllZero, 
                           drop = FALSE]
    ncoef <- ncol(fit$coefficients)
  }
  NACoef <- anyNA(fit$coefficients)
  if (NACoef) {
    i <- is.na(fit$coefficients)
    fit$coefficients[i] <- 0
    fit$stdev.unscaled[i] <- 1e+30
  }
  fit$coefficients <- fit$coefficients %*% contrasts
  if (length(cormatrix) < 2) {
    orthog <- TRUE
  }
  else {
    orthog <- all(abs(cormatrix[lower.tri(cormatrix)]) < 
                    1e-14)
  }
  R <- chol(fit$cov.coefficients)
  fit$cov.coefficients <- crossprod(R %*% contrasts)
  if (orthog) 
    fit$stdev.unscaled <- sqrt(fit$stdev.unscaled^2 %*% 
                                 contrasts^2)
  else {
    R <- chol(cormatrix)
    ngenes <- NROW(fit$stdev.unscaled)
    ncont <- NCOL(contrasts)
    U <- matrix(1, ngenes, ncont, dimnames = list(rownames(fit$stdev.unscaled), 
                                                  colnames(contrasts)))
    o <- array(1, c(1, ncoef))
    for (i in 1:ngenes) {
      RUC <- R %*% .vecmat(fit$stdev.unscaled[i, ], contrasts)
      U[i, ] <- sqrt(o %*% RUC^2)
    }
    fit$stdev.unscaled <- U
  }
  if (NACoef) {
    i <- (fit$stdev.unscaled > 1e+20)
    fit$coefficients[i] <- NA
    fit$stdev.unscaled[i] <- NA
  }
  fit
}

.ebayes <- function (fit, proportion = 0.01, stdev.coef.lim = c(0.1, 4), 
                     trend = FALSE, robust = FALSE, winsor.tail.p = c(0.05, 0.1)) 
{
  coefficients <- fit$coefficients
  stdev.unscaled <- fit$stdev.unscaled
  sigma <- fit$sigma
  df.residual <- fit$df.residual
  if (is.null(coefficients) || is.null(stdev.unscaled) || 
      is.null(sigma) || is.null(df.residual)) 
    stop("No data, or argument is not a valid lmFit object")
  if (max(df.residual) == 0) 
    stop("No residual degrees of freedom in linear model fits")
  if (!any(is.finite(sigma))) 
    stop("No finite residual standard deviations")
  if (is.logical(trend)) {
    if (trend) {
      covariate <- fit$Amean
      if (is.null(covariate)) 
        stop("Need Amean component in fit to estimate trend")
    }
    else {
      covariate <- NULL
    }
  }
  else {
    if (!is.numeric(trend)) 
      stop("trend should be either a logical scale or a numeric vector")
    if (!identical(length(sigma), length(trend))) 
      stop("If trend is numeric then it should have length equal to the number of genes")
    covariate <- trend
  }
  out <- squeezeVar(sigma^2, df.residual, covariate = covariate, 
                    robust = robust, winsor.tail.p = winsor.tail.p)
  out$s2.prior <- out$var.prior
  out$s2.post <- out$var.post
  out$var.prior <- out$var.post <- NULL
  out$t <- coefficients/stdev.unscaled/sqrt(out$s2.post)
  df.total <- df.residual + out$df.prior
  df.pooled <- sum(df.residual, na.rm = TRUE)
  df.total <- pmin(df.total, df.pooled)
  out$df.total <- df.total
  out$p.value <- 2 * pt(-abs(out$t), df = df.total)
  var.prior.lim <- stdev.coef.lim^2/median(out$s2.prior)
  out$var.prior <- tmixture.matrix(out$t, stdev.unscaled, 
                                   df.total, proportion, var.prior.lim)
  if (anyNA(out$var.prior)) {
    out$var.prior[is.na(out$var.prior)] <- 1/out$s2.prior
    warning("Estimation of var.prior failed - set to default value")
  }
  r <- rep(1, NROW(out$t)) %o% out$var.prior
  r <- (stdev.unscaled^2 + r)/stdev.unscaled^2
  t2 <- out$t^2
  Infdf <- out$df.prior > 10^6
  if (any(Infdf)) {
    kernel <- t2 * (1 - 1/r)/2
    if (any(!Infdf)) {
      t2.f <- t2[!Infdf]
      r.f <- r[!Infdf]
      df.total.f <- df.total[!Infdf]
      kernel[!Infdf] <- (1 + df.total.f)/2 * log((t2.f + 
                                                    df.total.f)/(t2.f/r.f + df.total.f))
    }
  }
  else kernel <- (1 + df.total)/2 * log((t2 + df.total)/(t2/r + 
                                                           df.total))
  out$lods <- log(proportion/(1 - proportion)) - log(r)/2 + 
    kernel
  out
}


eBayes <- function (fit, proportion = 0.01, stdev.coef.lim = c(0.1, 4), 
          trend = FALSE, robust = FALSE, winsor.tail.p = c(0.05, 0.1)) 
{
  if (is.logical(trend) && trend && is.null(fit$Amean)) 
    stop("Need Amean component in fit to estimate trend")
  eb <- .ebayes(fit = fit, proportion = proportion, stdev.coef.lim = stdev.coef.lim, 
                trend = trend, robust = robust, winsor.tail.p = winsor.tail.p)
  fit$df.prior <- eb$df.prior
  fit$s2.prior <- eb$s2.prior
  fit$var.prior <- eb$var.prior
  fit$proportion <- proportion
  fit$s2.post <- eb$s2.post
  fit$t <- eb$t
  fit$df.total <- eb$df.total
  fit$p.value <- eb$p.value
  fit$lods <- eb$lods
  if (!is.null(fit$design) && is.fullrank(fit$design)) {
    F.stat <- classifyTestsF(fit, fstat.only = TRUE)
    fit$F <- as.vector(F.stat)
    df1 <- attr(F.stat, "df1")
    df2 <- attr(F.stat, "df2")
    if (df2[1] > 1e+06) 
      fit$F.p.value <- pchisq(df1 * fit$F, df1, lower.tail = FALSE)
    else fit$F.p.value <- pf(fit$F, df1, df2, lower.tail = FALSE)
  }
  fit
}



champ.DMP_modify <- function (beta = myNorm, pheno = myLoad$pd$Sample_Group, compare.group = NULL, 
          adjPVal = 0.05, adjust.method = "BH", arraytype = "450K") 
{
  message("[===========================]")
  message("[<<<<< ChAMP.DMP START >>>>>]")
  message("-----------------------------")
  CalculateDMP <- function(beta, pheno, tmp_compare, adjPVal = adjPVal, 
                           adjust.method = adjust.method) {
    message("  -----------------------------")
    message("  Start to Compare : ", tmp_compare[1], ", ", 
            tmp_compare[2])
    p <- pheno[which(pheno %in% tmp_compare)]
    tmpbeta <- beta[, which(pheno %in% tmp_compare)]
    design <- model.matrix(~0 + p)
    contrast.matrix <- makeContrasts(contrasts = paste(colnames(design)[2:1], 
                                                       collapse = "-"), levels = colnames(design))
    message("  Contrast Matrix")
    print(contrast.matrix)
    fit <- lmFit(tmpbeta, design)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    tryCatch(fit3 <- eBayes(fit2), warning = function(w) {
      stop("limma failed, No sample variance.\n")
    })
    DMP <- topTable(fit3, coef = 1, number = nrow(tmpbeta), 
                    adjust.method = adjust.method, p.value = adjPVal)
    message("  You have found ", sum(DMP$adj.P.Val <= adjPVal), 
            " significant MVPs with a ", adjust.method, " adjusted P-value below ", 
            adjPVal, ".")
    message("  Calculate DMP for ", tmp_compare[1], " and ", 
            tmp_compare[2], " done.")
    return(DMP)
  }
  message("!!! Important !!! New Modification has been made on champ.DMP(): \n")
  message("    (1): In this version champ.DMP() if your pheno parameter contains more than two groups of phenotypes, champ.DMP() would do pairewise differential methylated analysis between each pair of them. But you can also specify compare.group to only do comparasion between any two of them.\n")
  message("    (2): champ.DMP() now support numeric as pheno, and will do linear regression on them. So covariates like age could be inputted in this function. You need to make sure your inputted \"pheno\" parameter is \"numeric\" type.\n")
  Compare <- NULL
  message("--------------------------------")
  if (is.null(pheno) | length(unique(pheno)) <= 1) {
    stop("pheno parameter is invalid. Please check the input, pheno MUST contain at least two phenotypes.")
  }
  else {
    if (length(pheno) != ncol(beta)) 
      stop("Your Pheno's length is not in accord with your beta value's ncol.")
    message("\n[ Section 1:  Check Input Pheno Start ]\n")
    if (class(pheno) == "numeric") {
      message("  You pheno is numeric type.")
      message("    pheno parameter contains :", length(pheno), 
              " values.")
      message("    pheno parameter ranges from ", min(pheno), 
              " to ", max(pheno))
    }
    else {
      message("  You pheno is ", class(pheno), " type.")
      message("    Your pheno information contains following groups. >>")
      sapply(unique(pheno), function(x) message("    <", 
                                                x, ">:", sum(pheno == x), " samples."))
      message("    [The power of statistics analysis on groups contain very few samples may not strong.]")
      if (length(unique(pheno)) == 2) {
        message("    pheno contains only 2 phenotypes")
        if (is.null(compare.group)) {
          message("    compare.group parameter is NULL, two pheno types will be added into Compare List.")
          Compare <- list(x1 = unique(pheno))
        }
        else if (sum(compare.group %in% unique(pheno)) == 
                 2) {
          message("    Your compare.group parameter is in accord with your pheno. Two pheno types has been added into Compare List.")
          Compare <- list(x1 = unique(pheno))
        }
        else {
          stop(" You have specified compare.group, but it's not in accord with your pheno parameter. Please recheck your compare.group or pheno.")
        }
      }
      else if (length(unique(pheno)) > 2) {
        message("    pheno contains ", length(unique(pheno)), 
                " phenotypes")
        if (is.null(compare.group)) {
          message("    compare.group parameter is NULL, EACH PAIR of phenotypes will be added into Compare List.")
          Compare <- as.list(data.frame(combn(unique(pheno), 
                                              2)))
        }
        else if (sum(compare.group %in% unique(pheno)) == 
                 2) {
          message("    Your compare.group parameter is in accord with your pheno. Two pheno types has been added into Compare List.")
          Compare <- list(x1 = sort(compare.group))
        }
        else {
          stop("    Your pheno parameter contains multiple phenotypes, but values in your compare.group parameter are not all find in them. Please recheck your compare.group or pheno.")
        }
      }
      else {
        stop("    !!! Something wrong with your pheno. Please check your input.")
      }
      tmpnamelist <- vector()
      for (i in 1:length(Compare)) {
        tmpname <- paste(Compare[[i]][1], Compare[[i]][2], 
                         sep = "_to_")
        message("    ", tmpname, " compare group : ", 
                Compare[[i]][1], ", ", Compare[[i]][2])
        tmpnamelist <- c(tmpnamelist, tmpname)
      }
      names(Compare) <- tmpnamelist
    }
    message("\n[ Section 1:  Check Input Pheno Done ]\n")
  }
  DMPs <- list()
  if (is.null(Compare)) {
    message("\n[ Section 2:  Find Numeric Covariates Linear Regression CpGs Start ]\n")
    df <- data.frame(pheno = pheno)
    model.matrix <- model.matrix(~pheno, data = df)
    fit1 <- lmFit(beta, model.matrix)
    fit2 <- eBayes(fit1)
    DMP <- topTable(fit2, coef = 2, number = nrow(beta), 
                    adjust.method = adjust.method, p.value = adjPVal)
    message("  You have found ", sum(DMP$adj.P.Val <= adjPVal), 
            " significant MVPs with a ", adjust.method, " adjusted P-value below ", 
            adjPVal, ".")
    if (sum(DMP$adj.P.Val <= adjPVal) != 0) 
      DMPs[["NumericVariable"]] <- DMP
  }
  else {
    message("\n[ Section 2:  Find Differential Methylated CpGs Start ]\n")
    for (i in names(Compare)) {
      DMP <- CalculateDMP(beta, pheno, Compare[[i]], adjPVal, 
                          adjust.method)
      if (sum(DMP$adj.P.Val <= adjPVal) != 0) 
        DMPs[[i]] <- DMP
    }
  }
  message("\n[ Section 2:  Find Numeric Vector Related CpGs Done ]\n")
  # if (length(DMPs) == 0) 
  #   stop("ChAMP.DMP Have not detected even one significant CpGs. You may try other threshold.")
  message("\n[ Section 3:  Match Annotation Start ]\n")
  if (arraytype == "EPIC") 
    data(probe.features.epic)
  else data(probe.features)
  for (i in names(DMPs)) {
    com.idx <- intersect(rownames(DMPs[[i]]), rownames(probe.features))
    if (!is.null(Compare)) {
      avg <- cbind(rowMeans(beta[com.idx, which(pheno == 
                                                  Compare[[i]][1])]), rowMeans(beta[com.idx, which(pheno == 
                                                                                                     Compare[[i]][2])]))
      avg <- cbind(avg, avg[, 2] - avg[, 1])
      colnames(avg) <- c(paste(Compare[[i]], "AVG", sep = "_"), 
                         "deltaBeta")
      DMPs[[i]] <- data.frame(DMPs[[i]][com.idx, ], avg, 
                              probe.features[com.idx, ])
    }
    else {
      DMPs[[i]] <- data.frame(DMPs[[i]][com.idx, ], probe.features[com.idx, 
      ])
    }
  }
  message("\n[ Section 3:  Match Annotation Done ]\n")
  message("[<<<<<< ChAMP.DMP END >>>>>>]")
  message("[===========================]")
  message("[You may want to process DMP.GUI() or champ.GSEA() next.]\n")
  return(DMPs)
}
