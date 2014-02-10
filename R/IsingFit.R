IsingFit <-
function(x, family=c('binomial', 'mgaussian'), AND = TRUE, gamma = 0.25, plot = TRUE, progressbar = TRUE, ...){
  t0 <- Sys.time()
  NodesToAnalyze <- apply(x,2,sd, na.rm=TRUE) != 0
  names(NodesToAnalyze) <- colnames(x)
  if (!any(NodesToAnalyze)) stop("No variance in dataset")
  if (any(!NodesToAnalyze))
  {
    warning(paste("Nodes without variance:",paste(colnames(x)[!NodesToAnalyze],collapse = ", ")))
  }
  x <- as.matrix(x)
  allthemeans <- colMeans(x)
  x <- x[,NodesToAnalyze,drop=FALSE]
  nvar <- ncol(x)
  p <- nvar - 1
  intercepts <- betas <- lambdas <- list(vector,nvar)
  nlambdas <- rep(0,nvar)
  for (i in 1: nvar){
    a <- glmnet(x[,-i], x[,i], family = family)
    intercepts[[i]] <- a$a0
    betas[[i]] <- a$beta
    lambdas[[i]] <- a$lambda
    nlambdas[i] <- length(lambdas[[i]])
  }
  if (progressbar==TRUE) pb <- txtProgressBar(max=nrow(x), style = 3)
  P <- logl <- sumlogl <- J <- matrix(0, max(nlambdas), nvar)
  for (i in 1:nvar)
  {
    J[1:ncol(betas[[i]]),i] <- colSums(betas[[i]]!=0)
  }
  for (n in 1:nrow(x)){
    for (i in 1: nvar){
      y <- intercepts[[i]] + colSums(betas[[i]]*x[n,-i])
      y <- c(y,rep(NA,max(nlambdas)-length(y)))
      P[,i] <- exp(y*x[n,i])/(1+exp(y))
      logl[,i] <- log(P[,i])
    }
    sumlogl <- sumlogl + logl
    if (progressbar==TRUE) setTxtProgressBar(pb, n)
  }
  if (progressbar==TRUE) close(pb)
  sumlogl[sumlogl==0]=NA
  penalty <- J * log(nrow(x)) + 2 * gamma * J * log(p)
  EBIC <- -2 * sumlogl + penalty
  lambda.opt <- apply(EBIC,2,which.min)
  thresholds <- 0
  for(i in 1:length(lambda.opt))
    thresholds[i] <- intercepts[[i]][lambda.opt[i]]
  weights.opt <- matrix(,nvar,nvar)
  for (i in 1:nvar){
    weights.opt[i,-i] <- betas[[i]][,lambda.opt[i]]
  }
  if (AND==TRUE) {
    adj <- weights.opt
#     adj[weights.opt!=0] <- 1
#     adj[weights.opt=0] <- 0
    adj <- (adj!=0)*1
    EN.weights <- adj * t(adj)
    EN.weights <- EN.weights * weights.opt
    meanweights.opt <- (EN.weights+t(EN.weights))/2
    diag(meanweights.opt) <- 0 
  } else {
    meanweights.opt <- (weights.opt+t(weights.opt))/2
    diag(meanweights.opt) <- 0
  }
  graphNew <- matrix(0,length(NodesToAnalyze),length(NodesToAnalyze))
  graphNew[NodesToAnalyze,NodesToAnalyze] <- meanweights.opt
  colnames(graphNew) <- rownames(graphNew) <- colnames(x)
  threshNew <- ifelse(allthemeans > 0.5, -Inf, Inf)
  threshNew[NodesToAnalyze] <- thresholds
  if (plot==TRUE) notplot=FALSE else notplot=TRUE
  q <- qgraph(graphNew,layout='spring',labels=names(NodesToAnalyze),DoNotPlot=notplot,...)
  Res <- list(weiadj = graphNew, thresholds = threshNew, q = q, gamma = gamma, AND = AND, time = Sys.time() - t0)
  class(Res) <- "IsingFit"
  return(Res)
}

## Methods:
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

exportNetLogo.IsingFit <- function(x)
{
  if (is.character(x))
  {
    x[] <- paste0('"',x,'"')
  }
  if (is.vector(x))
  {
    return(paste("\n",paste(paste0("[",paste(x, collapse = " "),"]"),collapse="\n"),"\n\n"))
  } else if (is.matrix(x))
  {
    return(paste("[\n",paste(apply(x,1,function(s)paste0("[",paste(s, collapse = " "),"]")),collapse="\n"),"\n]")  )
  } else stop("Object not supported")
  
}

# 
# print(fit)
# plot(fit)
# summary(fit)
# export(fit)

