IsingFit <-
  function(x, family = "binomial", AND = TRUE, gamma = 0.25, plot = TRUE, progressbar = TRUE, min_sum = -Inf, lowerbound.lambda = NA, ...) {
    t0 <- Sys.time()
    xx <- x
    if (family != "binomial") {
      stop("This procedure is currently only supported for binary (family='binomial') data")
    }

    ## Check to prevent error of lognet() in package glmnet
    # checklognet <- function(y){
    #   res <- c() # 0: too little variance, 1: good to go
    #   y=as.factor(y)
    #   ntab=table(y)
    #   minclass=min(ntab)
    #   if(minclass<=1) res=0 else res=1
    #   return(res)
    # }

    allowedNodes <- function(nodeValues) {
      nodeValues <- as.factor(nodeValues)
      valuesFrequency <- table(nodeValues)
      minFrequency <- min(valuesFrequency)
      maxFrequency <- max(valuesFrequency)
      if (minFrequency <= 1 || maxFrequency >= length(nodeValues) - 1) {
        return(0)
      } else {
        return(1)
      }
    }

    # NodesToAnalyze <- apply(x,2,checklognet) !=0
    NodesToAnalyze <- apply(x, 2, allowedNodes) != 0
    names(NodesToAnalyze) <- colnames(x)
    if (!any(NodesToAnalyze)) stop("No variance in dataset")
    if (any(!NodesToAnalyze)) {
      warning(paste("Nodes with too little variance (not allowed):", paste(colnames(x)[!NodesToAnalyze], collapse = ", ")))
    }
    ##

    x <- as.matrix(x)
    if (!all( x == 1 | x== 0)) {warning("IsingFit only supports data that is encoded as (0,1)")}
    allthemeans <- colMeans(x)
    x <- x[, NodesToAnalyze, drop = FALSE]
    nvar <- ncol(x)
    p <- nvar - 1
    intercepts <- betas <- lambdas <- list(vector, nvar)
    nlambdas <- rep(0, nvar)
    N <- vector()
    for (i in 1:nvar) {
      subData <- x[rowSums(replace(x[, -i, drop = FALSE], x[, -i, drop = FALSE] < 0, 0)) != (min_sum - 1), ]
      a <- glmnet(subData[, -i], subData[, i], family = "binomial")
      intercepts[[i]] <- a$a0
      betas[[i]] <- a$beta
      lambdas[[i]] <- a$lambda
      nlambdas[i] <- length(lambdas[[i]])
      N[i] <- nrow(subData)
    }

    if (progressbar == TRUE) pb <- txtProgressBar(max = nvar, style = 3)
    # penalty now different for each variable so for EBIC and penalty also empty matrix
    P <- logl <- sumlogl <- J <- EBIC <- penalty <- matrix(0, max(nlambdas), nvar)
    for (i in 1:nvar)
    {
      J[1:ncol(betas[[i]]), i] <- colSums(betas[[i]] != 0)
    }
    logl_M <- P_M <- list() # I could not get an array with differ nrows so now a list
    for (i in 1:length(N)) {
      logl_M[[i]] <- P_M[[i]] <- array(0, dim = c(N[i], max(nlambdas))) # sample size
    }

    for (i in 1:nvar) { # i <- 1
      
      subData <- x[rowSums(replace(x[, -i, drop = FALSE], x[, -i, drop = FALSE] < 0, 0)) != (min_sum - 1), ]
      sample_size <- N[i]
      betas.ii <- as.matrix(betas[[i]])
      int.ii <- intercepts[[i]]
      y <- matrix(0, nrow = sample_size, ncol = ncol(betas.ii))
      xi <- subData[, -i]
      NB <- nrow(betas.ii) # number of rows in beta
      for (bb in 1:NB) { # bb <- 1
        y <- y + betas.ii[rep(bb, sample_size), ] * xi[, bb]
      }
      y <- matrix(int.ii, nrow = sample_size, ncol = ncol(y), byrow = TRUE) + y
      # number of NAs
      n_NA <- max(nlambdas) - ncol(y)
      if (n_NA > 0) {
        for (vv in 1:n_NA) {
          y <- cbind(y, NA)
        }
      }
      # calculate P matrix
      P_M[[i]] <- exp(y * subData[, i]) / (1 + exp(y))
      logl_M[[i]] <- log(P_M[[i]])
      if (progressbar==TRUE) setTxtProgressBar(pb, i)
    }

    for (i in 1:nvar) {
      sumlogl[, i] <- colSums(logl_M[[i]])
    }
    if (progressbar==TRUE) close(pb)
    sumlogl[sumlogl == 0] <- NA

    for (i in 1:nvar) {
      penalty[, i] <- J[, i] * log(N[i]) + 2 * gamma * J[, i] * log(p)
      EBIC[, i] <- -2 * sumlogl[, i] + penalty[, i]
    }



    lambda.mat <- matrix(NA, nrow(EBIC), ncol(EBIC))
    for (i in 1:nvar) {
      lambda.mat[, i] <- c(lambdas[[i]], rep(NA, nrow(EBIC) - length(lambdas[[i]])))
    }

    if (!is.na(lowerbound.lambda)) {
      EBIC <- EBIC / (lambda.mat >= lowerbound.lambda) * 1
    }

    lambda.opt <- apply(EBIC, 2, which.min)
    lambda.val <- rep(NA, nvar)
    thresholds <- 0
    for (i in 1:length(lambda.opt)) {
      lambda.val[i] <- lambda.mat[lambda.opt[i], i]
      thresholds[i] <- intercepts[[i]][lambda.opt[i]]
    }
    weights.opt <- matrix(, nvar, nvar)
    for (i in 1:nvar) {
      weights.opt[i, -i] <- betas[[i]][, lambda.opt[i]]
    }
    asymm.weights <- weights.opt
    diag(asymm.weights) <- 0
    if (AND == TRUE) {
      adj <- weights.opt
      adj <- (adj != 0) * 1
      EN.weights <- adj * t(adj)
      EN.weights <- EN.weights * weights.opt
      meanweights.opt <- (EN.weights + t(EN.weights)) / 2
      diag(meanweights.opt) <- 0
    } else {
      meanweights.opt <- (weights.opt + t(weights.opt)) / 2
      diag(meanweights.opt) <- 0
    }
    graphNew <- matrix(0, length(NodesToAnalyze), length(NodesToAnalyze))
    graphNew[NodesToAnalyze, NodesToAnalyze] <- meanweights.opt
    colnames(graphNew) <- rownames(graphNew) <- colnames(xx)
    threshNew <- ifelse(allthemeans > 0.5, -Inf, Inf)
    threshNew[NodesToAnalyze] <- thresholds
    if (plot == TRUE) notplot <- FALSE else notplot <- TRUE
    q <- qgraph(graphNew, layout = "spring", labels = names(NodesToAnalyze), DoNotPlot = notplot, ...)
    Res <- list(
      weiadj = graphNew, thresholds = threshNew, q = q, gamma = gamma,
      AND = AND, time = Sys.time() - t0, asymm.weights = asymm.weights,
      lambda.values = lambda.val
    )
    class(Res) <- "IsingFit"
    return(Res)
  }

plot.IsingFit <- function(object,...) qgraph(object$q,DoNotPlot = FALSE, ...)

print.IsingFit <- function(x)
{
  cat("Estimated network:\n")
  
  print(round(x$weiadj,2))
  
  cat("\n\nEstimated Thresholds:\n")
  
  print(x$thresholds)  
}

summary.IsingFit <- function(object)
{
  cat("\tNetwork Density:\t\t", round(mean(object$weiadj[upper.tri(object$weiadj)]!=0),2),"\n",
      "Gamma:\t\t\t",round(object$gamma,2),"\n",
      "Rule used:\t\t",ifelse(object$AND,"And-rule","Or-rule"),"\n",
      "Analysis took:\t\t",format(object$time,format="%s"),"\n"
  )
}

