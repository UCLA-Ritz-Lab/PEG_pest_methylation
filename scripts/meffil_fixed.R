meffil.ewas.summary_fix <- function (ewas.object, beta, selected.cpg.sites = character(0), 
          parameters = meffil.ewas.parameters(), verbose = T) 
{
  if (is.data.frame(beta)) 
    beta <- as.matrix(beta)
  stopifnot(is.matrix(beta))
  stopifnot(parameters$model %in% colnames(ewas.object$p.value))
  stopifnot(parameters$max.plots < nrow(ewas.object$p.value))
  stopifnot(all(selected.cpg.sites %in% rownames(ewas.object$p.value)))
  p.values <- ewas.object$p.value[, parameters$model]
  parameters$practical.threshold <- p.values[order(p.values)[parameters$max.plots + 
                                                               1]]
  if (is.na(parameters$sig.threshold)) 
    parameters$sig.threshold <- 0.05/nrow(ewas.object$p.value)
  sig.idx <- unique(which(ewas.object$p.value < parameters$sig.threshold, 
                          arr.ind = T)[, "row"])
  practical.idx <- which(p.values < parameters$practical.threshold)
  selected.idx <- match(selected.cpg.sites, rownames(ewas.object$p.value))
  significant.sites <- rownames(ewas.object$p.value)[sig.idx]
  selected.sites <- rownames(ewas.object$p.value)[selected.idx]
  cpg.idx <- union(sig.idx, union(practical.idx, selected.idx))
  cpg.sites <- rownames(ewas.object$p.value)[cpg.idx]
  cpg.stats <- data.frame(ewas.object$analyses[[1]]$table[cpg.sites, 
                                                          c("chromosome", "position")], p.value = ewas.object$p.value[cpg.idx, 
                                                          ], coefficient = ewas.object$coefficient[cpg.idx, ])
  msg("QQ plots", verbose = T)
  qq.plots <- meffil.ewas.qq.plot(ewas.object, sig.threshold = parameters$sig.threshold, 
                                  lambda.method = parameters$qq.inflation.method)
  msg("Manhattan plots", verbose = T)
  manhattan.plots <- meffil.ewas.manhattan.plot(ewas.object, 
                                                sig.threshold = parameters$sig.threshold)
  plot.sites <- rownames(ewas.object$p.value)[union(practical.idx, 
                                                    selected.idx)]
  msg("CpG site plots:", length(plot.sites), verbose = T)
  if (is.character(beta)) 
    beta <- retrieve.gds.methylation(beta, sites = plot.sites, 
                                     samples = NULL)
  cpg.plots <- sapply(plot.sites, function(cpg) {
    msg("Plotting", cpg, verbose = T)
    meffil.ewas.cpg.plot(ewas.object, cpg = cpg, beta = beta)
  }, simplify = F)
  msg("Sample characteristics", verbose = T)
  sample.characteristics <- meffil.ewas.sample.characteristics_fix(ewas.object)
  covariate.associations <- meffil.ewas.covariate.associations(ewas.object)
  parameters$winsorize.pct <- ewas.object$winsorize.pct
  parameters$outlier.iqr.factor <- ewas.object$outlier.iqr.factor
  parameters$rlm <- ewas.object$rlm
  parameters$most.variable <- ewas.object$most.variable
  parameters$random.seed <- ewas.object$random.seed
  parameters$sample.size <- length(ewas.object$samples)
  sort.by.p <- function(x) {
    sites <- x
    if (!is.character(x)) 
      sites <- names(x)
    idx <- match(sites, rownames(ewas.object$p.value))
    p <- ewas.object$p.value[idx, parameters$model]
    x[order(p)]
  }
  list(parameters = parameters, qq.plots = qq.plots, manhattan.plots = manhattan.plots, 
       cpg.stats = cpg.stats, cpg.plots = sort.by.p(cpg.plots), 
       significant.sites = sort.by.p(significant.sites), selected.sites = sort.by.p(selected.sites), 
       sample.characteristics = sample.characteristics, covariate.associations = covariate.associations)
}

msg <- function (..., verbose = T) 
{
  if (verbose) {
    x <- paste(list(...))
    name <- sys.call(sys.parent(1))[[1]]
    cat(paste("[", name, "]", sep = ""), date(), x, "\n")
  }
}

meffil.ewas.sample.characteristics_fix <- function (ewas.object) 
{
  stopifnot(is.ewas.object(ewas.object))
  msg("summarizing variables")
  summarize.variable <- function(name, variable) {
    msg(name)
    if (is.character(variable)) 
      variable <- as.factor(variable)
    if (is.factor(variable)) {
      n <- sapply(levels(variable), function(level) sum(variable == 
                                                          level, na.rm = T))
      data.frame(variable = name, value = names(n), mean = format(n), 
                 var = format(n/length(variable) * 100, digits = 3, 
                              nsmall = 1), row.names = NULL)
    }
    else {
      data.frame(variable = name, value = "", mean = format(mean(variable, na.rm = T)), 
                 var = format(sd(variable, na.rm = T)), 
                 row.names = NULL)
    }
  }
  var.summary <- summarize.variable("variable of interest", 
                                    ewas.object$variable)
  cov.summary <- NULL
  if (!is.null(ewas.object$covariates)) {
    cov.summary <- lapply(1:ncol(ewas.object$covariates), 
                          function(i) {
                            summarize.variable(colnames(ewas.object$covariates)[i], 
                                               ewas.object$covariates[, i]%>% as_vector)
                          })
    cov.summary <- do.call(rbind, cov.summary)
  }
  rbind(var.summary, cov.summary)
}


is.ewas.object <- function (object) 
  is.list(object) && all(c("class", "version", "samples", "variable", 
                           "covariates", "winsorize.pct", "outlier.iqr.factor", "analyses", 
                           "random.seed") %in% names(object))


meffil.ewas.parameters <- function (sig.threshold = NA, max.plots = 10, model = "none", 
          qq.inflation.method = "median") 
{
  list(sig.threshold = sig.threshold, max.plots = max.plots, 
       model = model, qq.inflation.method = qq.inflation.method)
}

meffil.ewas.covariate.associations <- function (ewas.object) 
{
  stopifnot(is.ewas.object(ewas.object))
  if (is.null(ewas.object$covariates)) 
    return(NULL)
  msg("covariate associations")
  ret <- lapply(1:ncol(ewas.object$covariates), function(i) {
    msg(colnames(ewas.object$covariates)[i])
    vars <- data.frame(variable = ewas.object$variable, 
                       covariate = ewas.object$covariates[, i])
    for (j in colnames(vars)) {
      if (!is.numeric(vars[[j]]) && length(unique(na.omit(vars[[j]]))) >= 
          2) 
        vars[[j]] <- as.factor(vars[[j]])
    }
    colnames(vars)[2] <- colnames(ewas.object$covariates)[i]
    meffil.summarize.relationship(vars)
  })
  names(ret) <- colnames(ewas.object$covariates)
  ret
}


meffil.ewas.qq.plot <- function (ewas.object, sig.threshold = 1e-07, sig.color = "red", 
                                 title = "QQ plot", xlab = bquote(-log[10]("expected p-values")), 
                                 ylab = bquote(-log[10]("observed p-values")), lambda.method = "median") 
{
  stopifnot(is.ewas.object(ewas.object))
  sapply(names(ewas.object$analyses), function(name) {
    p.values <- sort(ewas.object$analyses[[name]]$table$p.value, 
                     decreasing = T)
    p.values[which(p.values < .Machine$double.xmin)] <- .Machine$double.xmin
    stats <- data.frame(is.sig = p.values < sig.threshold, 
                        expected = -log(sort(ppoints(p.values), decreasing = T), 
                                        10), observed = -log(p.values, 10))
    lambda <- qq.lambda(p.values[which(p.values > sig.threshold)], 
                        method = lambda.method)
    label.x <- min(stats$expected) + diff(range(stats$expected)) * 
      0.1
    label.y <- min(stats$expected) + diff(range(stats$observed)) * 
      0.9
    lambda.label <- paste("lambda == ", format(lambda$estimate, 
                                               digits = 3), "%+-%", format(lambda$se, digits = 3), 
                          "~(", lambda.method, ")", sep = "")
    selection.idx <- scatter.thinning(stats$observed, stats$expected, 
                                      resolution = 100, max.per.cell = 100)
    lim <- range(c(0, stats$expected, stats$observed))
    sig.threshold <- format(sig.threshold, digits = 3)
    (ggplot(stats[selection.idx, ], aes(x = expected, y = observed)) + 
        geom_abline(intercept = 0, slope = 1, colour = "black") + 
        geom_point(aes(colour = factor(sign(is.sig)))) + 
        scale_colour_manual(values = c("black", "red"), 
                            name = "Significant", breaks = c("0", "1"), 
                            labels = c(paste("p-value >", sig.threshold), 
                                       paste("p-value <", sig.threshold))) + annotate(geom = "text", 
                                                                                      x = label.x, y = label.y, hjust = 0, label = lambda.label, 
                                                                                      parse = T) + xlim(lim) + ylim(lim) + xlab(xlab) + 
        ylab(ylab) + coord_fixed() + ggtitle(paste(title, 
                                                   ": ", name, sep = "")))
  }, simplify = F)
}

  
meffil.ewas.manhattan.plot <- function (ewas.object, sig.threshold = 1e-07, title = "Manhattan plot") 
{
  stopifnot(is.ewas.object(ewas.object))
  chromosomes <- paste("chr", c(1:22, "X", "Y"), sep = "")
  sapply(names(ewas.object$analyses), function(name) {
    chromosomes <- intersect(chromosomes, ewas.object$analyses[[name]]$table$chromosome)
    stats <- ewas.object$analyses[[name]]$table
    stats$chromosome <- factor(as.character(stats$chromosome), 
                               levels = chromosomes)
    stats$chr.colour <- 0
    stats$chr.colour[stats$chromosome %in% chromosomes[seq(1, 
                                                           length(chromosomes), 2)]] <- 1
    p.values <- stats$p.value
    p.values[which(p.values < .Machine$double.xmin)] <- .Machine$double.xmin
    stats$stat <- -log(p.values, 10) * sign(stats$coefficient)
    stats <- stats[order(stats$stat, decreasing = T), ]
    chromosome.lengths <- sapply(chromosomes, function(chromosome) max(stats$position[which(stats$chromosome == 
                                                                                              chromosome)]))
    chromosome.lengths <- as.numeric(chromosome.lengths)
    chromosome.starts <- c(1, cumsum(chromosome.lengths) + 
                             1)
    names(chromosome.starts) <- c(chromosomes, "NA")
    stats$global <- stats$position + chromosome.starts[stats$chromosome] - 
      1
    selection.idx <- scatter.thinning(stats$global, stats$stat, 
                                      resolution = 100, max.per.cell = 100)
    (ggplot(stats[selection.idx, ], aes(x = position, y = stat)) + 
        geom_point(aes(colour = chr.colour)) + facet_grid(. ~ 
                                                            chromosome, space = "free_x", scales = "free_x") + 
        theme(strip.text.x = element_text(angle = 90)) + 
        guides(colour = FALSE) + labs(x = "Position", y = bquote(-log[10]("p-value") * 
                                                                   sign(beta))) + geom_hline(yintercept = log(sig.threshold, 
                                                                                                              10), colour = "red") + geom_hline(yintercept = -log(sig.threshold, 
                                                                                                                                                                  10), colour = "red") + theme(axis.text.x = element_blank(), 
                                                                                                                                                                                               axis.ticks.x = element_blank()) + ggtitle(paste(title, 
                                                                                                                                                                                                                                               ": ", name, sep = "")))
  }, simplify = F)
}


retrieve.gds.methylation <- function (gds.filename, sites, samples) 
{
  retrieve.gds.matrix(gds.filename, sites, samples)
}


retrieve.gds.matrix <- function (gds.filename, sites, samples) 
{
  stopifnot(file.exists(gds.filename))
  gds.file <- openfn.gds(gds.filename)
  on.exit(closefn.gds(gds.file))
  all.sites <- read.gdsn(index.gdsn(gds.file, "row.names"))
  if (is.null(sites)) 
    sites <- all.sites
  else stopifnot(all(sites %in% all.sites))
  all.samples <- read.gdsn(index.gdsn(gds.file, "col.names"))
  if (is.null(samples)) 
    samples <- all.samples
  else stopifnot(all(samples %in% all.samples))
  matrix.node <- index.gdsn(gds.file, "matrix")
  mat <- readex.gdsn(matrix.node, sel = list(all.sites %in% 
                                               sites, all.samples %in% samples), simplify = "none")
  rownames(mat) <- all.sites[which(all.sites %in% sites)]
  colnames(mat) <- all.samples[which(all.samples %in% samples)]
  mat[sites, samples, drop = F]
}

qq.lambda <- function(p.values, method="median", B=100) {
  stopifnot(method %in% c("median","regression","robust"))
  p.values <- na.omit(p.values)
  observed <- qchisq(p.values, df=1, lower.tail = FALSE)
  observed <- sort(observed)
  expected <- qchisq(ppoints(length(observed)), df=1, lower.tail=FALSE)
  expected <- sort(expected)
  
  lambda <- se <- NA
  if (method == "median")  {
    lambda <- median(observed)/qchisq(0.5, df=1)
    boot.medians <- sapply(1:B, function(i) median(sample(observed, replace=T)))
    se <- sd(boot.medians/qchisq(0.5,df=1))
  } else if (method %in% c("regression","robust")) {
    if (method == "regression")
      coef.table <- summary(lm(observed ~ 0 + expected))$coeff
    else
      coef.table <- summary(rlm(observed ~ 0 + expected))$coef
    lambda <- coef.table["expected",1]
    se <- coef.table["expected", "Std. Error"]
  }
  list(method=method, estimate=lambda, se=se)
}

scatter.thinning <- function(x,y,resolution=100,max.per.cell=100) {
  x.cell <- floor((resolution-1)*(x - min(x,na.rm=T))/diff(range(x,na.rm=T))) + 1
  y.cell <- floor((resolution-1)*(y - min(y,na.rm=T))/diff(range(y,na.rm=T))) + 1
  z.cell <- x.cell * resolution + y.cell
  frequency.table <- table(z.cell)
  frequency <- rep(0,max(z.cell, na.rm=T))
  frequency[as.integer(names(frequency.table))] <- frequency.table
  f.cell <- frequency[z.cell]
  
  big.cells <- length(which(frequency > max.per.cell))
  sort(c(which(f.cell <= max.per.cell),
         sample(which(f.cell > max.per.cell),
                size=big.cells * max.per.cell, replace=F)),
       decreasing=F)
}
