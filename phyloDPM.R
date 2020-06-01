# ***************************************************#
#
# SCRIPT: phyloDPM.R
# AUTHOR: Denis Awany; Emile R. Chimusa
# DATE: 30/03/2020
# REV: 1.0
# VERSION: v1.0
# PLATFORM: Not Platform Dependent
# 
# PURPOSE: This script is used to analyse micobiome data. 
#          It based on the logistic Dirichlet process random effects model.
#          Produces Bayesian posterior estimates for OTU effects.


# ===================================================
#                    BEGINNING OF MAIN
# ===================================================
phyloDPM.main <- function(Y, #case-control phenotype [vector]
                          Z, # microbiome abundance data [matrix or dataframe object]
                          mTree, #phylogenetic tree[matrix or phylo class object]
                          Nsim, #number of MCMC iterations [integer]
                          thin, # thinning value for posterior samples [integer]
                          burnin, #number of burn-in samples [integer]
                          post.plots, # whether or not to display plots of posterior samples
                          mu, # tuning parameter for phylogenetic and RBF kernels.
                          alpha.par,#first hyper-parameter for Gamma prior of precision parameter
                          sig2.bet.par,#variance parameter on beta
                          tau2.par, # variance parameter on the base parameter of the DP.
                          rbf.par, # kernel coefficient for the RBF kernel.
                          saveSim, # save raw simulated samples to a file in the current directory.
                          evol.rate, # evolutionary rate parameter for the phylogenetic kernel.
                          ...){

start.time = proc.time()[3] # time the duration of the code execution.

welcome=function()
{
  message("\n");
  message("******************************************");
  message("*              phyloDPM v1.0             *");
  message("*                                        *");
  message("*      Awany Denis; Emile R. Chimusa     *");
  message("*   Human Genetics, Dept. of Pathology   *");
  message("*         University of Cape Town        *");
  message("******************************************");
  message("\n");
}

welcome()
  

# ----------------------------------------------------
#                   Check user inputs
# ----------------------------------------------------

if(missing(Y) | missing (Z)) {
  stop("Response, and OTU data are required!")}
if(!is.vector(Y)){
  stop('Response must be a vector.')}
if(class(Z)!='data.frame' & class(Z)!='matrix'){
  stop('OTU data must be given as a dataframe or matrix.')}
if(length(Y) != nrow(Z)) {
  stop("Response vector must have the same length as the microbiome sample size!")}
if(!is.finite(mu) | !is.numeric(mu) | mu<0 | mu>1){stop("invalid entry: mu must be a numeric between 0 and 1.")}
if(is.null(mTree)){
  message("Phylogenetic tree data is missing; taxa similarity will be based on the radial basis function.")
  if(mu >0){message("Warning: mu ignored, as mTree is missing.")}
}else{
    if(class(mTree)!="phylo" & class(mTree)!="matrix"){
      stop("mTree must be a distance matrix or an object of the class 'phylo' ")}
    if(mu==1){message("Warning: radial basis function information will not be used, as 'mu' is set to 1.")}
    if(mu==0){message("Warning: phylogenetic tree information will not be used, as 'mu' is set to 0.")}
}
if (thin >= Nsim ) {stop("Error: thin cannot be equal to or greater than the number of posterior draws (Nsim).")}
if (thin < 1) {stop("Error: thin should be greater than or equal to 1.")}
if(is.null(burnin)){
    burnin = round(0.1 * Nsim) # if burn-in not supplied, set it to 10% of number of posterior draws (Nsim).
}else{
    if(burnin < 0){stop("Error! burnin should be a positive number.") }
    if(burnin >= Nsim){stop("Error! burnin should be less than the number of posterior draws (Nsim).") }
}
if(!is.finite(evol.rate)| !is.numeric(evol.rate) | evol.rate < 0){stop("Error: 'evol.rate' parameter must be a numberic, positive value.")}

message("\nInput data succesfully verified ... \n")


# Load required packages
if(post.plots==TRUE) {
  message("The 'coda' package is required for some posterior plots. It will be installed and/or loaded.")
  if (!require("coda")) { # Needed for Bayesian posterior plots
  install.packages("coda")
  library(coda)
  }
} 

  # ===================================================
  #                DEFINE FUNCTIONS HERE
  # ====================================================

# thinning Bayesian posterior samples. From boom 20/07/2019.
thin.samples <- function (x, thin) {
  ## Thins the data in x, so that every thin'th observation is
  ## returned.  This is useful for obtaining aprximately uncorrelated
  ## posterior samples from Bayesian posterior.
  ##
  ## Args:
  ##   x: A numeric vector, matrix, or array to be thinned.  If a
  ##     matrix or array then the first dimension corresponds to the
  ##     observation number, and will be used for thinning.
  ##   thin: The frequency with which to thin.  E.g. if thin == 10
  ##     then every 10th observation will be returned.
  ##
  ## Returns:
  ##    The thinned subset of x.
  stopifnot(is.numeric(thin) && length(thin) == 1)
  if (thin <= 1)
    return(x)
  if (is.array(x)) {
    nr <- dim(x)[1]
  } else if (is.numeric(x)) {
    nr <- length(x)
  } else stop("x must be a numeric type in thin()")

  top <- floor(nr/thin)
  indx <- (1:top) * thin

  if (is.matrix(x)) {
    return(x[indx, ])
  } else if (is.numeric(x)) {
    return(x[indx])
  } else if (is.array(x)) {
    stop("Error! you cannot drop all samples.")
  }
}

  # rdirichlet from LearnBayes (24/03/2019)
  rdirichlet <- function(n, alpha)
  ## sampling from Dirichlet process using the Polya Urn scheme.
  ##
  ## Args:
  ##   n: number of samples; likened to number of balls in the Urn
  ##   alpha: concentration (precision) parameter.
  ##
  ## Returns:
  ##    samples distributed as a DP.
  {
      k = length(alpha)
      z = array(0, dim = c(n, k))
      s = array(0, dim = c(n, 1))
      for (i in 1:k) {
          z[, i] = rgamma(n, shape=alpha[i])
          s = s + z[, i]
      }
      for (i in 1:k) {
          z[, i] = z[, i]/s
      }
      return(z)
  }
  
  #Posterior predictive p-value, adapted from bayespval func in 'contig' package by Antony O.
  # 21/07/2019.
  bayespval <- function(Y, Z, Beta, statistic="X2"){
    MU   <- exp(as.matrix(Z )%*% as.vector(Beta))
    Pred <- matrix(rbinom(n=prod(dim(MU)), size=1, prob=MU/(1+MU)),ncol=dim(MU)[2])
    Obs <- Y

    Tpred<-apply(((Pred-MU)^2)/MU,1,sum)
    Tobs<-apply(((Obs-MU)^2)/MU,1,sum)
     
    pval<-mean(as.numeric(Tpred>Tobs))
    pval
  }

  HPDinterval.Bayesian <- function(obj, prob = 0.95, ...) {
    obj <- as.matrix(obj)
    vals <- apply(obj, 2, sort)
    if (!is.matrix(vals)) stop("obj must have nsamp > 1")
    nsamp <- nrow(vals)
    npar <- ncol(vals)
    gap <- max(1, min(nsamp - 1, round(nsamp * prob)))
    init <- 1:(nsamp - gap)
    inds <- apply(vals[init + gap, ,drop=FALSE] - vals[init, ,drop=FALSE],
                  2, which.min)
    ans <- cbind(vals[cbind(inds, 1:npar)],
                 vals[cbind(inds + gap, 1:npar)])
    dimnames(ans) <- list(colnames(obj), c("lower", "upper"))
    attr(ans, "Probability") <- gap/nsamp
    ans 
  }

  # cophenetic.phylo from 'ape' (06/04/2019)
  dist.nodes <- function(x)
  {
      x <- reorder(x) # required for the C code
      n <- Ntip(x)
      m <- x$Nnode
      nm <- n + m
      d <- .C("dist_nodes", as.integer(n), as.integer(m), # "dist_nodes", from a C file named 'dist_nodes.c'.
              as.integer(x$edge[, 1] - 1L), as.integer(x$edge[, 2] - 1L),
              as.double(x$edge.length), as.integer(Nedge(x)),
              double(nm * nm), NAOK = TRUE)[[7]]
      dim(d) <- c(nm, nm)
      dimnames(d) <- list(1:nm, 1:nm)
      d
  }
  
  cophenetic.phylo <- function(x)
  {
      n <- length(x$tip.label)
      ans <- dist.nodes(x)[1:n, 1:n]
      dimnames(ans)[1:2] <- list(x$tip.label)
      ans
  }
  
  # rbfkernel from rdetools (02/04/2019)
  rbfkernel <- function(X, sigma = 1, Y = NULL)
  {
      # test if X is a matrix
      if(!is.matrix(X))
      {
          print("X must be a matrix containing samples in its rows")
          return()
      }
      # test if sigma is a number and > 0
      if(length(sigma) != 1 || sigma <= 0)
      {
          print("sigma must be a number > 0 specifying the rbf-kernel width")
          return()
      }
      if(!is.null(Y))
      {
          # test if Y is a matrix
          if(!is.matrix(Y))
          {
              print("Y must be a matrix containing samples in its rows or NULL if it should not be used")
              return()
          }
          # test if vectors in X and Y have same dimension
          if(ncol(X) != ncol(Y))
          {
              print("The samples in the rows of X and Y must be of same dimension")
              return()
          }
      }
  
      n <- nrow(X) # number of samples in X
  
      if(is.null(Y))
      {
          # calculate distance matrix
          XtX <- tcrossprod(X)
          XX <- matrix(1, n) %*% diag(XtX)
          D <- XX - 2*XtX + t(XX) # distance matrix
      }
      else
      {
          m <- nrow(Y) # number of samples in Y
          # calculate distance matrix (between vectors of X and Y)
          XX <- matrix(apply(X ^ 2, 1, sum), n, m)
          YY <- matrix(apply(Y ^ 2, 1, sum), n, m, byrow = TRUE)
          XY <- tcrossprod(X, Y)
          D <- XX - 2*XY + YY
      }
  
      # calculate rbf-kernel matrix
      K <- exp(-D/(2*sigma^2))
  
      return(K)
  }
  
  # rdiscrete from e1071 (24/03/2019)
  rdiscrete <- function (n, probs, values = 1:length(probs), ...)
  {
      sample(values, size = n, replace = TRUE, prob = probs)
  }
  

  # rinvgamma from MCMCpack (24/03/2019)
  rinvgamma <- function (n, shape, scale = 1)
  {
      return(1/rgamma(n = n, shape = shape, rate = scale))
  }
  
  # rmvnorm from mvtnorm (24/03/2019)
  rmvnorm <- function (n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)),
      method = c("eigen", "svd", "chol"))
  {
      if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps))) {
          stop("sigma must be a symmetric matrix")
      }
      if (length(mean) != nrow(sigma)) {
          stop("mean and sigma have non-conforming size")
      }
      sigma1 <- sigma
      dimnames(sigma1) <- NULL
      if (!isTRUE(all.equal(sigma1, t(sigma1)))) {
          warning("sigma is numerically not symmetric")
      }
      method <- match.arg(method)
      if (method == "eigen") {
          ev <- eigen(sigma, symmetric = TRUE)
          if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
              warning("sigma is numerically not positive definite")
          }
          retval <- ev$vectors %*% diag(sqrt(ev$values), length(ev$values)) %*%
              t(ev$vectors)
      }
      else if (method == "svd") {
          sigsvd <- svd(sigma)
          if (!all(sigsvd$d >= -sqrt(.Machine$double.eps) * abs(sigsvd$d[1]))) {
              warning("sigma is numerically not positive definite")
          }
          retval <- t(sigsvd$v %*% (t(sigsvd$u) * sqrt(sigsvd$d)))
      }
      else if (method == "chol") {
          retval <- chol(sigma, pivot = TRUE)
          o <- order(attr(retval, "pivot"))
          retval <- retval[, o]
      }
      retval <- matrix(rnorm(n * ncol(sigma)), nrow = n) %*% retval
      retval <- sweep(retval, 2, mean, "+")
      colnames(retval) <- names(mean)
      retval
  }
  
  
  # rtnorm from msm (06/04/2019); truncated normal sampler
  rtnorm <- function (n, mean = 0, sd = 1, lower = -Inf, upper = Inf)
  {
      if (length(n) > 1)
          n <- length(n)
      mean <- rep(mean, length = n)
      sd <- rep(sd, length = n)
      lower <- rep(lower, length = n)
      upper <- rep(upper, length = n)
      lower <- (lower - mean)/sd
      upper <- (upper - mean)/sd
      ind <- seq(length = n)
      ret <- numeric(n)
      alg <- ifelse(lower > upper, -1, ifelse(((lower < 0 & upper == Inf) | (lower == 
                   -Inf & upper > 0) | (is.finite(lower) & is.finite(upper) & (lower < 0) & 
                    (upper > 0) & (upper - lower > sqrt(2 * pi)))), 0, ifelse((lower >= 0 & 
                    (upper > lower + 2 * sqrt(exp(1))/(lower + sqrt(lower^2 + 4)) * 
                    exp((lower * 2 - lower * sqrt(lower^2 + 4))/4))), 1, ifelse(upper <= 0 & 
                    (-lower > -upper + 2 * sqrt(exp(1))/(-upper + sqrt(upper^2 + 4)) * 
                       exp((upper * 2 - -upper * sqrt(upper^2 + 4))/4)), 2, 3))))
      ind.nan <- ind[alg == -1]
      ind.no <- ind[alg == 0]
      ind.expl <- ind[alg == 1]
      ind.expu <- ind[alg == 2]
      ind.u <- ind[alg == 3]
      ret[ind.nan] <- NaN
      while (length(ind.no) > 0) {
          y <- rnorm(length(ind.no))
          done <- which(y >= lower[ind.no] & y <= upper[ind.no])
          ret[ind.no[done]] <- y[done]
          ind.no <- setdiff(ind.no, ind.no[done])
      }
      stopifnot(length(ind.no) == 0)
      while (length(ind.expl) > 0) {
          a <- (lower[ind.expl] + sqrt(lower[ind.expl]^2 + 4))/2
          z <- rexp(length(ind.expl), a) + lower[ind.expl]
          u <- runif(length(ind.expl))
          done <- which((u <= exp(-(z - a)^2/2)) & (z <= upper[ind.expl]))
          ret[ind.expl[done]] <- z[done]
          ind.expl <- setdiff(ind.expl, ind.expl[done])
      }
      stopifnot(length(ind.expl) == 0)
      while (length(ind.expu) > 0) {
          a <- (-upper[ind.expu] + sqrt(upper[ind.expu]^2 + 4))/2
          z <- rexp(length(ind.expu), a) - upper[ind.expu]
          u <- runif(length(ind.expu))
          done <- which((u <= exp(-(z - a)^2/2)) & (z <= -lower[ind.expu]))
          ret[ind.expu[done]] <- -z[done]
          ind.expu <- setdiff(ind.expu, ind.expu[done])
      }
      stopifnot(length(ind.expu) == 0)
      while (length(ind.u) > 0) {
          z <- runif(length(ind.u), lower[ind.u], upper[ind.u])
          rho <- ifelse(lower[ind.u] > 0, exp((lower[ind.u]^2 -
              z^2)/2), ifelse(upper[ind.u] < 0, exp((upper[ind.u]^2 -
              z^2)/2), exp(-z^2/2)))
          u <- runif(length(ind.u))
          done <- which(u <= rho)
          ret[ind.u[done]] <- z[done]
          ind.u <- setdiff(ind.u, ind.u[done])
      }
      stopifnot(length(ind.u) == 0)
      ret * sd + mean
  }
  

  # 'is.empty' from sjmisc: to check whether string, list or vector is empty (09/04/2019)
  is.empty <- function(x, first.only = FALSE) {
    # do we have a valid vector?
    if (!is.null(x)) {
      # if it's a character, check if we have only one element in that vector
      if (is.character(x)) {
        # characters may also be of length 0
        if (length(x) == 0) return(TRUE)
        # else, check all elements of x
        zero_len <- nchar(x) == 0
        # return result for multiple elements of character vector
        if (first.only) {
          zero_len <- .is_true(zero_len[1])
          if (length(x) > 0) x <- x[!is.na(x)][1]
        } else {
          return(unname(zero_len))
        }
        # we have a non-character vector here. check for length
      } else if (is.list(x)) {
        x <- purrr::compact(x)
        zero_len <- length(x) == 0
      } else {
        zero_len <- length(x) == 0
      }
    }
  
    any(is.null(x) || zero_len || all(is.na(x)))
  
  }
  
  # 'is.empty' from rapporttools: to check whether, vector value is empty (12/04/2019)
  is_empty <- function(x, trim = TRUE, ...) {
      if (length(x) <= 1) {
          if (is.null(x))
              return (TRUE)
          if (length(x) == 0)
              return (TRUE)
          if (is.na(x) || is.nan(x))
              return (TRUE)
          if (is.character(x) && nchar(ifelse(trim, trim.space(x), x)) == 0)
              return (TRUE)
          if (is.logical(x) && !isTRUE(x))
              return (TRUE)
          if (is.numeric(x) && x == 0)
              return (TRUE)
          return (FALSE)
      } else
          sapply(x, is_empty, trim = trim, ...)
  }
  
  # con_par; concentration parameter of DP (29/03/2019)
  con_par <- function(initial_m_value, n, K, prior_m_params)
  {
    x <- rbeta(1, initial_m_value + 1, n)
    pi1 <- prior_m_params[1] + K - 1
    pi2 <- n * (prior_m_params[2] - log(x))
    pi1 <- pi1/(pi1 + pi2)
    if (runif(1) < pi1) {
      g1 <- rgamma(1, prior_m_params[1] + K, prior_m_params[2] - log(x))
      new_m_value <- g1
    } else {
      g2 <- rgamma(1, prior_m_params[1] + K - 1, prior_m_params[2] - log(x))
      new_m_value <- g2
    }
    #new_concentration parameter
    return(new_m_value)
  }

 hpd_interval <- function(obj, prob = 0.95, ...) {
   obj <- as.matrix(obj)
   vals <- apply(obj, 2, sort)
   if (!is.matrix(vals)) stop("obj must have nsamp > 1")
   nsamp <- nrow(vals)
   npar <- ncol(vals)
   gap <- max(1, min(nsamp - 1, round(nsamp * prob)))
   init <- 1:(nsamp - gap)
   inds <- apply(vals[init + gap, ,drop=FALSE] - vals[init, ,drop=FALSE],
                 2, which.min)
   ans <- cbind(vals[cbind(inds, 1:npar)],
                vals[cbind(inds + gap, 1:npar)])
   dimnames(ans) <- list(colnames(obj), c('lower', 'upper'))
   attr(ans, "Probability") <- gap/nsamp
   round(ans,3) 
 }


 is.invertible <- function(m) class(try(solve(m),silent=T))=="matrix" #Test invertibility of matrix.
  
  # ---------------------------------------------
  # Check user input and get basic info from data
  # ---------------------------------------------
  
  # ----------------------------------------------
  # Get basic info from data/ convert to convenient forms for computations.
  # ----------------------------------------------
  Y       <- as.vector(Y)
  Z       <- as.matrix(Z)
  n       <- num.sam <- length(Y)
  num.otu <- ncol(Z)

  # -----------------------------------------------
  # Initialize the model parameters
  # -----------------------------------------------
  precis_param  <- alpha.par #Gamma (a.b) parameters for the precision parameter of the DP.
  if(is.null(precis_param)) {precis_param = c(exp(-0.033*num.sam), exp(-0.033*num.sam))} # these are the parameters of Gamma(a,b) prior on m.
  a.tau         <- tau2.par[1]  # gamma distr parameter for tau - the parameter for base measure of DP.
  b.tau         <- tau2.par[2]  
  a.beta        <- sig2.bet.par[1]  # inv-gamma parameter for sig.beta2 - the parameter for variance component of Beta.
  b.beta        <- sig2.bet.par[2]  #

  tau.sq        <- rinvgamma(1,shape=a.tau,scale=b.tau)
  if(!is.finite(tau.sq)) {tau.sq <- 1.0e+50} # If Infinity, use 1.0e+50 to represent the large enough infinite value!!
  tau           <- sqrt(tau.sq)
  sig.beta      <- 1
  
  # Also initialize the DP parameters.
  # Algorithm to simulate DP parameters adapted from Kyung et al (2010). 
  q             <- c(rdirichlet(1, rep(1, num.sam)))      # PROBABILITY OF ASSIGNMENT
  A.n           <- sample(1:num.sam,num.sam,replace=TRUE,prob=q)      # CREATES THE ASSIGNMENTS.
  A.K           <- table(A.n)                             # CREATES THE LIST OF OCCUPIED
  A.K.labels    <- as.numeric(names(A.K))                 # LOCATIONS OF OCCUPIED
  K             <- length(A.K)                            # NUMBER OF OCCUPIED TABLES
  B             <- rep(0,length=num.sam)                        # LENGTH num.sam FOR ALL CASES
  B[A.K.labels] <- 1                                      # 1 WHERE OCCUPIED >0
  B             <- B*rnorm(num.sam,0,tau)                       # ASSIGN psi VALUES TO OCCUPIED
  nk            <- rep(0,num.sam); for (i in 1:num.sam) nk[A.n[i]] <- nk[A.n[i]] + 1  # COUNTS OCCUPANTS AT TABLES
  psi           <- B[A.n]                               # MAP TABLE psi'S TO CASES 
  eta           <- as.vector(rmvnorm(1, mean=rep(0,K), diag(rep(tau.sq), K)))
  

  ## ------------------------------------------------
  # Compute OTU similarity measures

  # Create 'RBF kernel
  #X must be a matrix containing samples in its rows
  rbf.kernel = as.matrix(rbfkernel(X=t(Z[,]), sigma=rbf.par)) #transpose to apply across taxa

  # Create Phylogenetic tree kernel
  if(is.null(mTree)){
  	JointKernel = rbf.kernel
  }else{
  	if(class(mTree)=='matrix'){
      # if mTree supplied as a similarity matrix, convert to dissimilarity matrix
      # as exponentiation converts back to similarity matrix. 
      if(mTree[1,1]==1){ # testing if a similarity matrix.
        mTree.dissim.matrix <- (1-mTree) 
      }else if(mTree[1,1]==0){ # testing if a dissimilarity matrix
        mTree.dissim.matrix <- mTree # if mTree given as dissimilarity matrix
      }else if (class(mTree)=='matrix' & mTree[1,1]!=0 & mTree[1,1]!=1 ) {
        stop("Error: the mTree matrix object supplied is not a similarity measure. Please check!")
      }
    }
    if(class(mTree)=='phylo'){
      mTree.dissim.matrix <- as.matrix(cophenetic.phylo(mTree))
	  }

    phylo.kernel <- exp(-2.0*evol.rate*mTree.dissim.matrix^2.0) #
    JointKernel  <- mu*phylo.kernel + (1-mu)*rbf.kernel #Composite kernel
  }

  joint.cov <- JointKernel 


  #scaling data, avoid need for intercept#
  mean.z <- drop(rep(1, num.sam) %*% Z)/num.sam
  Z      <- scale(Z, center=mean.z, scale=FALSE)

  Beta1  <- rmvnorm (1, mean = rep(0, ncol(Z)), sigma = diag(ncol(Z)))


  # -------------------------------------------------
  # Declare storage matrix for the sampled parameters.
  # --------------------------------------------------
  Beta     <- rbind(Beta1, matrix(NA, nrow=Nsim, ncol=ncol(Z)))
  tau      <- c(tau, rep(NA,Nsim))
  sig.beta <- c(sig.beta, rep(NA,Nsim))
  
  #name the matrices and vectors [variables]:
  colnames(Beta)  <- paste0("OTU",1:num.otu)
  names(tau)      <- "tau"
  names(sig.beta) <- "sig.beta"

  # ----------------------------------------------
  # calculate the likelihood.
  pr <- 1 - 1/(1 + exp(Z%*%t(Beta1) + psi))  #
  like.K <- 0 #
  for (i in 1:n){
    like.K <- like.K + (dbinom(Y[i], size=1, prob=pr[i], log=TRUE)) + (dnorm(psi[i], mean=0, sd=tau[1],log=TRUE))
    like.K <- exp(like.K + sum(lgamma(A.K)))  #
  }

  # For the prior of m, we use Gamma(a,b)
  m        = n/log(n) # precision parameter of the DP part. Initialize by n/log(n)
  m        = c(m, rep(NA,Nsim))
  names(m) = "m"

  Mu   <- rnorm(ncol(Z),0,1)
  Gama <- rnorm(ncol(Z),0,1)
  tau2.gam <- 1
  
  cat('\nSampling using slice-Gibbs algorithm ... \n\n')
  
  # -------------------------------------------------------
  # Perform MCMC (Gibbs-slice) sampling for Nsim iterations.
  # -------------------------------------------------------
  for (iter in 2:(Nsim+1)){  # simulate samples (particles).
  
    ## print iteration to screen, if print.iter==TRUE
       if ( iter%%8==0 | iter%%20==0 ){
         cat('Completed', format(round((iter/(Nsim+1))*100,1),nsmall=1),"%."," Wait...", "\r")
       }

    # -------------------------------------------------------
    # Performing M-H algorithm to generate 'A' matrix.
    # -------------------------------------------------------
    p.A.old  <- ((m[iter-1])^K)
    f.y.old  <- like.K 
    mult.old <- dmultinom(x=A.K, prob=q[A.K.labels])
    
    # Create new 'A' matrix; 'can' stands for candidate
    pq                <- nk +1
    new.q             <- rdirichlet(1, pq)
    A.n.can           <- sample(1:num.sam,num.sam,replace=TRUE,new.q)
    A.K.can           <- table(A.n.can)
    A.K.labels.can    <- as.numeric(names(A.K.can))
    K.can             <- length(A.K.can)
    B                 <- rep(0,length=num.sam)
    B[A.K.labels.can] <- 1
    B                 <- B*rnorm(num.sam,0, tau[iter-1])
    nk.can            <- rep(0,num.sam); for (i in 1:num.sam) nk.can[A.n.can[i]] <- nk.can[A.n.can[i]] + 1
    psi.can           <- B[A.n.can]
    
    pr <- 1 - 1/(1 + exp(Z%*%Beta[iter-1,] + psi))
    like.K.can        <- 0
    for(i in 1:n){
      like.K.can <- like.K.can + (dbinom(Y[i], size=1, prob=pr[i], log=TRUE)) + (dnorm(psi[i], mean=0, sd=tau[iter-1], log=TRUE)) #
      like.K.can <- exp(like.K.can + sum(lgamma(A.K.can))) #
    }
    
    p.A.can  <- (m[iter-1]^K.can)
    f.y.can  <- like.K.can
    mult.can <- dmultinom(x=A.K.can, prob=new.q[A.K.labels.can])

    # Define the M-H acceptance probability to ACCEPT or REJECT simulated candidates.
    p.ratio <- p.A.can/p.A.old; f.ratio <- f.y.can/f.y.old; mult.ratio <- mult.can/mult.old
    if (!is.finite(p.ratio)) {p.ratio       <- 1} # If value NOT finite, set it to 1.
    if (!is.finite(f.ratio)) {f.ratio       <- 1}
    if (!is.finite(mult.ratio)) {mult.ratio <- 1}
    rho <- p.ratio * f.ratio * mult.ratio # Acceptance probability.
    if (rho>runif(1)) {
        A.n        <- A.n.can
        A.K        <- A.K.can
        A.K.labels <- A.K.labels.can
        K          <- K.can
        nk         <- nk.can
        psi        <- psi.can
        like.K     <- like.K.can
    }

    ## ==================
    # Simulate uniform random variates; delta and lamb.
    func     <- Z%*%Beta[iter-1,] + psi
    prob.fun <- exp(func)/(1+exp(func)) #  
    delta    <- runif(num.sam, rep(0,num.sam), prob.fun)

    #And for lamb, we have:
    lamb <- runif(num.sam, rep(0,num.sam), 1-prob.fun)

    # Simulate Beta 
    for (j in 1:ncol(Z)){
        Gama[j] <- rnorm(n=1, mean=(1/(as.numeric(solve(t(joint.cov[,j])%*%joint.cov[,j])) + 1/tau2.gam)) * 
        	(as.numeric(solve(t(joint.cov[,j])%*%joint.cov[,j]))*Mu[j]) , 
        	sd=as.numeric(solve(t(joint.cov[,j])%*%joint.cov[,j])) + 1/tau2.gam)

        Mu[j] <- rnorm(n=1, mean=1/(1/sig.beta[iter-1]^2 + as.numeric(solve(t(joint.cov[,j])%*%joint.cov[,j]))) * 
        	((1/sig.beta[iter-1]^2)*Beta[iter-1,j] + as.numeric(solve(t(joint.cov[,j])%*%joint.cov[,j]))%*%Gama[j]), 
        	sd=1/sig.beta[iter-1]^2 + as.numeric(solve(t(joint.cov[,j])%*%joint.cov[,j])) )          

        # Define the truncation limits (intervals) to simulate 'Beta'. 
        S1 <- (1/Z[Y==1,j]) * (log(delta[Y==1]/(1-delta[Y==1])) - Z[Y==1,-j] %*% Beta[iter-1,-j] - psi[Y==1])
        S2 <- (1/Z[Y==0,j]) * (log((1-lamb[Y==0])/lamb[Y==0]) - Z[Y==0,-j] %*% Beta[iter-1,-j] - psi[Y==0])   

        #left truncation limit
        indx1 = Z[Y==1,j] > 0 
        A = indx1*S1 
        A[A==0] = NaN #
        left1 = suppressWarnings(max(A, na.rm=TRUE))

        indx2 = Z[Y==0,j] < 0 
        B = indx2*S2
        B[B==0] = NaN
        left2 = suppressWarnings(max(B, na.rm=TRUE))
        beta.left = max(left1, left2)

        #right truncation limit
        indx3 = Z[Y==1,j] < 0 
        C = indx3*S1
        C[C==0] = NaN
        right1 = suppressWarnings(min(C, na.rm=TRUE))

        indx4 = Z[Y==0,j] >0 
        D = indx4*S2
        D[D==0] = NaN
        right2  = suppressWarnings(min(D, na.rm=TRUE))
        beta.right = min(right1, right2)

        # simulating from truncated normal distribution.
        #Either:
        unif.random1 <- runif(1,0,1)
        a1 <- pnorm(beta.left, Mu[j], sig.beta[iter-1])   #
        a2 <- pnorm(beta.right, Mu[j], sig.beta[iter-1])
        pa <- a1 + unif.random1*(a2-a1)
        if(pa > 0.99999) {pa=0.99999} # Avoid numerical issues, since qnorm(1,mean,sd)==Inf
        current.beta.j <- qnorm(pa, Mu[j], sig.beta[iter-1])
        #Or: 
        #current.beta.j <- rtnorm(n=1, mean = Mu[j], sd = sig.beta[iter-1], lower = beta.left, upper = beta.right)

        Beta[iter,j] <- current.beta.j
       
    }

    ## ====================
    # Simulate 'eta'.
    A.m1   = matrix(0, nrow=num.sam, ncol=num.sam)
    for (i in 1:n){ A.m1[i,A.n[i]] <- 1 }
    A.K.m  = A.m1[,which(apply(A.m1,2,sum)>0)]

    for (k in 1:K){
      # Lower limit for eta
       eta.lim1 = log(delta[Y==1]/(1-delta[Y==1])) - Z[Y==1,] %*% Beta[iter,] - A.K.m[Y==1,-k] %*% eta[1:K][-k] #
       eta.lim2 = log((1-lamb[Y==0])/lamb[Y==0]) - Z[Y==0,] %*% Beta[iter,] - A.K.m[Y==0,-k] %*% eta[1:K][-k] #

       eta.left  = suppressWarnings(max(eta.lim1, na.rm=TRUE)) # If set empty, set limit to -Inf
       eta.right = suppressWarnings(min(eta.lim2, na.rm=TRUE)) # If set empty, set limit to  Inf

       #simulating from truncated normal distribution.
       #Either:
       unif.random2 <- runif(1,0,1)
       b1 <- pnorm(eta.left, 0, tau[iter-1])
       b2 <- pnorm(eta.right, 0, tau[iter-1])
       pb <- b1 + unif.random2*(b2-b1)
       if(pb > 0.99999) {pb=0.99999} #Avoid numerical issues, since qnorm(1,mean,sd)==Inf
       current.eta.k <- qnorm(pb, 0,  tau[iter-1])
       #Or:
       #current.eta.k <- rtnorm(n=1, mean=0, sd=tau[iter-1], lower=eta.left, upper=eta.right)
       eta[k] <- current.eta.k 
    }

    ## =====================
    # Use simulated eta to generate psi.
    psi   <- A.K.m %*% eta[1:K]

    ## =====================
    # Perform conditional update on tau2, sig.beta2.
    # simulate tau2
    tau.sq <- rinvgamma(1,shape=(0.5*K)+a.tau, scale=(0.5*sum(eta[1:K]^2))+b.tau )
    
    if (is.finite(tau.sq)){
      tau[iter] <- sqrt(tau.sq) 
    }else{  #since tau.sq can go to infinity.
    	tau.sq <- 1#1.7e100 # If infinite, 1.7e100 is big enough to denote infinity!!
      tau[iter] <- sqrt(tau.sq) #mean(tau[1:iter-1])
      cat("Numerical error encountered in simulating 'tau.square'; infinite value during iteration ",iter,".")
    } #


    ## ====================
    # simulate sig.beta.sq
    sig2.beta <- rinvgamma(1, shape=(0.5*ncol(Z)+a.beta), scale=(0.5*sum((Beta[iter,]-Mu)^2)+b.beta))
    if (is.finite(sig2.beta)){ 
      sig.beta[iter] <- sqrt(sig2.beta)
    }else{ #since sig2.beta can go to infinity.
      sig2.beta <- 1 #
      sig.beta[iter] <- sqrt(sig2.beta) # mean(sig.beta[1:iter-1])
      cat("Numerical error encountered in simulating 'sigma2.beta'; infinite value during iteration ",iter,".")
    } #


    ## ===================
    # Update concentration paremeter,'m':
    # inputs: [previous m, sample size (n), cluster size (K), Gamma prior on m].
    m.cur  <- con_par(m[iter-1], n, K, precis_param)  # 
    m[iter]<- m.cur

    if(iter==(Nsim+1)){print('Sampling completed.', quote = FALSE)}

  } # End of Nsim.


    # Write (if save is TRUE) the raw simulation output to a file.
   if (saveSim){
      cat("\nWriting raw simulation data to a file in the present directory ...\n")
      out.df <- cbind(Beta, sig.beta, m, tau)
      column.names <- c(paste0("OTU", 1:num.otu), "sigma.beta", "m", "tau")
      #In writing the output to the file, exclude first row (initialization values for the simulations)
      write.table(out.df[-1,], file="simParams.txt", quote=FALSE, row.names = FALSE, col.names = column.names)
    }

  cat('\nPerforming posterior inference ...\n')

  # Estimate of Parameters
  # Discard burn-in samples and do thinning (or not, when thin =1) to posterior samples before conducting 
  # inference on the parameters.

  # Burn-in
  Beta.Burn         <- Beta[(burnin+1):Nsim, ] 
  tau.Burn          <- tau[(burnin+1):Nsim]
  sig.beta.Burn     <- sig.beta[(burnin+1):Nsim]
  m.Burn            <- m[(burnin+1):Nsim]

  # thinning
  Beta.BurnThin     <- thin.samples(as.matrix(Beta.Burn), thin) 
  tau.BurnThin      <- thin.samples(as.vector(tau.Burn),thin)
  sig.beta.BurnThin <- thin.samples(as.vector(sig.beta.Burn),thin)
  m.BurnThin        <- thin.samples(as.vector(m.Burn),thin)


  cat(paste0("\nNumber of effective posterior samples ", nrow(Beta.BurnThin),".\n"))


# =================================================
#              Posterior summary
# =================================================
posterior.summary <- function(object,...){
#
# S3 method to compute and print posterior summaries for a matrix of draws
# Modified from package bayeslm, author Hahn et al.(2017).
#
    numEff= 
    function(x,m=as.integer(min(length(x),(100/sqrt(5000))*sqrt(length(x)))))
    {
    # purpose:
    #  compute N-W std error and relative numerical efficiency
    # 
    # Arguments:
    #  x is vector of draws
    #  m is number of lags to truncate acf
    #  def is such that m=100 if length(x)= 5000 and grows with sqrt(length)
    # 
    # Output:
    #  list with numerical std error and variance multiple (f)
    #
    wgt=as.vector(seq(m,1,-1))/(m+1)
    z=acf(x,lag.max=m,plot=FALSE)
    f=1+2*wgt %*% as.vector(z$acf[-1])
    stderr=sqrt(var(x)*f/length(x))
    list(stderr=stderr,f=f,m=m)
    }
  X = object
  if(mode(X) == "list") stop("list entered \n Please extract sampling matrix from the list \n")
  if(mode(X) !="numeric") stop("Requires numeric argument \n")
  if(is.null(attributes(X)$dim)) X=as.matrix(X)
  nx=ncol(X)
  
  if(is.null(colnames(X))){
      names=as.character(1:nx)
  }else{
      names = colnames(X)
  }

  
  post.mean.beta=matrix(format(round(apply(X,2,mean), digits=3), nsmall=3), nrow=1)

  sds.beta=format(round(apply(X,2,sd),3), nsmall=3) # posterior stand devs.
  quantile95.beta=format(round(apply(X, 2, function(x){quantile(x,c(0.025, 0.975), na.rm=TRUE)}),3), nsmall=3) # quantiles (Credible Interval).
  hpd=format(round(hpd_interval(object, prob = 0.95),3), nsmall=3) # Credible Interval

  num_se=double(nx); rel_eff=double(nx); eff_s_size=double(nx)
  for(i in 1:nx)  # posterior numerical standard errors.
     {out=numEff(X[,i])
      if(is.nan(out$stderr)) 
          {num_se[i]=-9999; rel_eff[i]=-9999; eff_s_size[i]=-9999} 
      else
          {num_se[i]=out$stderr; rel_eff[i]=out$f; eff_s_size[i]=nrow(X)/ceiling(out$f)}
     }
  num_se = format(round(num_se,3), nsmall=3); rel_eff= format(round(rel_eff,3), nsmall=3)

  mat = rbind(post.mean.beta, sds.beta, num_se, rel_eff, quantile95.beta )
  mat = t(mat)
  
  mat = data.frame(mat) 

  #names of variables
  rownames(mat)=c(paste0("OTU",1:num.otu), "sig2.beta", "m")
  colnames(mat)[1]="mean" #posterior mean
  colnames(mat)[2]="std.dev" #posterior stan devs
  colnames(mat)[3]="num.se" # numerical standard errors
  colnames(mat)[4]="rel.eff" # relative efficiency
  colnames(mat)[5]="95% CI(lo)" # Credible Interval (lower)
  colnames(mat)[6]="CI(up)" # (upper)
  
  # format the output results
  mat=format(mat, width = 4, justify = "centre", digits=3, scientific=FALSE)  

  cat("\nEstimates of the parameters:\n")
  print((mat),digits=2)

  #return(invisible(1))
  return(mat)
  # invisible(t(mat))
} # End of 'posterior.summary' function 


# ============================================
#  Posterior Plots (traceplot, ACF plot) on the full simulated data
# ============================================
graphs.function <- function(X, tvalues,...){
    
  nx=ncol(X)

  if(nx==1) par(mfrow=c(1,1)) 
  if(nx==2) par(mfrow=c(2,1))
  if(nx==3) par(mfrow=c(3,1))
  if(nx==4) par(mfrow=c(2,2))
  if(nx>=5) par(mfrow=c(5,4))

  if(post.plots){
     if(nx==1) par(mfrow=c(1,4)) 
     if(nx==2) par(mfrow=c(2,4))
     if(nx>=3) par(mfrow=c(5,4))
     for(index in 1:nx){
        plot(as.vector(X[,index]),xlab="Iteration",ylab="",main=paste0("Trace plot of ",colnames(X)[index]),type="l",col="red") #traceplot
        plot(density(X[,index]), xlab="Value",ylab="Density",main=paste0("Density plot of ",colnames(X)[index]), type="l",col="blue") #density plot
        cumuplot(mcmc(X[,index]),probs=c(0.025,0.5,0.975),ylab="",xlab="Iteration",main=paste0("Cumulative quantile plot of ",colnames(X)[index]),
          lty=c(2,1),lwd=c(1,2),type="l",ask=FALSE,auto.layout=FALSE, col=1) #cumulative quantile plot
        if(!missing(tvalues)) abline(h=tvalues[index],lwd=2,col="blue")
        if(var(X[,index])>1.0e-20) {acf(as.vector(X[,index]),xlab="Lag",ylab="Autocorrelation",main=paste0("ACF plot of ",colnames(X)[index]))} #ACF plot
        else {plot.default(X[,index],xlab="",ylab="",type="n",main="No ACF Produced")}
     }
   }
  invisible()
} # end of graphs.function


# =============================================
#            Posterior Inference
# =============================================
# # Do posterior inference on 'effective' data (after thinning and burn-in).
post.est <- posterior.summary(cbind(Beta.BurnThin, sig.beta.BurnThin,m.BurnThin))
            

# # Make posterior plots using ALL the simulated data (no thin and no burn-in)
if(post.plots==TRUE){
  graphs.function(cbind(Beta)) #posterior plots for OTU effects..
}


# Stop the clock
end.time <- proc.time()[3] # '3' in seconds.
time.taken <- end.time - start.time # calculate difference
#In realistic use (with large Nsim), 'time.taken' will be several minutes to hours.
if(time.taken < 60){
  cat("\nTotal time elaspsed:", round(time.taken,digits=0), "seconds.\n")
}else if(time.taken >=60 & time.taken < 3600){
    cat("\nTotal time elaspsed:", round(time.taken/60, digits=0), 
      ifelse(round(time.taken/60, digits=0)==1,"minute\n", "minutes.\n"))
  }else{ # for time in hours (>3600 seconds)
    cat("\nTotal time elaspsed:", (time.taken %/%3600), 
      ifelse((time.taken %/%3600) < 2,"hour","hours"), 
      round((time.taken %% 3600)/60,digits=0),
      ifelse(round((time.taken %% 3600)/60,digits=0) < 2,"minute.\n","minutes.\n"))
  }

cat("\nDone! Thank you.\n")
cat("----------------\n")


#return posterior means of parameters
return(list(all.post = data.frame(post.est), 
	          beta = round(apply(Beta.BurnThin,2, mean),3),
            beta.95q = round(apply(Beta.BurnThin, 2, function(x){quantile(x,c(0.025, 0.975))}),3),
            beta.sd = round(apply(Beta.BurnThin, 2, sd),3),
            beta.se = round(sqrt(diag(var(Beta.BurnThin))),3),
            precis.param = round(mean(m.BurnThin),3),
            tau = round(mean(tau.BurnThin),3),
            sigma.beta = round(mean(sig.beta.BurnThin),3),
            raw.data.beta = as.matrix(Beta),
            raw.data.sig.beta = as.vector(sig.beta),
            raw.data.tau = as.vector(tau), 
            raw.data.m = as.vector(m) ))
} # End of 'phyloDPM.main' function.

#================================================
# We create a generic user-input function called "PhyloDPM" (18/07/2019)
#================================================
phyloDPM <- function(pheno, otu.data, mTree, Nsim, thin=1, burnin=NULL, post.plots=FALSE, mu=0, 
                      alpha.par=NULL, tau2.par=c(20,20),sig2.bet.par=c(20,15), rbf.par=0.001, saveSim=FALSE,evol.rate=10^3,...) {

  est <- phyloDPM.main(Y = pheno, Z = otu.data, mTree = mTree, tau2.par=tau2.par, sig2.bet.par=sig2.bet.par,
  			alpha.par=alpha.par, Nsim=Nsim, thin=thin, burnin=burnin, 
  			post.plots=post.plots, mu=mu, rbf.par=rbf.par,saveSim=saveSim,evol.rate=evol.rate,...)
  #class(est) <- "phyloDPM"
  est	
}

# ================================================
#             END OF phyloDPM.R SCRIPT
# ================================================
