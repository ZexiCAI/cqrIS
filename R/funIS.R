#' A Functional Censored Quantile Regression Model (The Interface)
#'
#' An interface to fit a censored quantile regression model to the provided dataset where
#' the covariates include both non-functional and functional types. The parameters
#' are obtained by minimizing an equivalent objective $l_1$ function of an estimating
#' equation. It also provides an interface for variance estimation via multiplier bootstrap.
#'
#' @param Z Non-functional covariates
#' @param W Transformed functional covariates
#' @param X Survival time, not log-transformed
#' @param cen Censoring indicator, 1 for failure, 0 for censored
#' @param xi Weights for samples, default: a 1-vector
#' @param tau.max Upper bound of estimated quantile
#' @param grid The grid size
#' @param diverge Norm of the estimate for the update to be considered divergent
#'
#' @import quantreg
#'
#' @return This function returns a list of lists with each list containing three elements:
#' \itemize{
#'   \item beta.est, the estimates for the non-functional coefficients
#'   \item gamma.est, the estimates for the functional coefficients
#'   \item idenLim, the identifiability limit for the data
#' }
#'
#' @references Jiang, F., Cheng, Q., Yin, G. and Shen, H. (2020),
#' "Generalizing Quantile Regression for Counting Processes with Applications to Recurrent Events,"
#' \emph{Journal of the American Statistical Association}, \bold{115}, 931-944.
#' @seealso Cai, Z. and Sit, T. (2020+),
#' "Censored Quantile Regression with Induced Smoothing,"
#' \emph{Working Paper}.

estjiang_ <- function(Z, W, X, cen, xi=rep(1,length(X)), tau.max, grid, diverge){
  R <- 1e10
  qtile.list <- seq(0, tau.max, by = grid)
  n <- dim(Z)[1]
  p <- dim(Z)[2]
  d <- dim(W)[2]
  H.qtile <- -log(1-qtile.list)
  diff.H <- diff(H.qtile)
  identLimit <- max(qtile.list)

  resp <- c(xi*cen*log(X), R, R)
  jiang.est <- array(NA, c(length(qtile.list), p+d), dimnames = list(qtile.list, c(paste0("beta",1:p),paste0("gamma",1:d))))
  jiang.est[1,] <- rep(NA, p+d)
  minus <- (log(X) >= -Inf) * diff.H[1]

  covr <- array(NA, dim = c(n+2,p+d))
  bZ <- cbind(Z, W)
  covr[1:n,] <- cen*(xi*bZ)
  covr[dim(covr)[1]-1,] <- crossprod(xi*bZ,-cen)
  # estimation
  for (i in 1:(length(qtile.list)-1)){
    qtile <- qtile.list[i+1]

    covr[dim(covr)[1],] <- crossprod(xi*bZ,2*minus)
    jiang.est[i+1,] <- rq(resp~0+covr,tau=0.5)$coef

    if (anyNA(jiang.est[i+1,]) | any(is.infinite(jiang.est[i+1,])) | t(jiang.est[i+1,])%*%jiang.est[i+1,]>diverge){
      jiang.est[i+1,] <- jiang.est[i, ]
      identLimit <- min(identLimit, qtile)
      warning(paste0("jiang procedure does not converge at tau = ", qtile))
    }
    minus = minus + (log(X) >= bZ %*% jiang.est[i+1,]) * diff.H[i+1]
  }
  return(list(beta.est=jiang.est[,1:p], gamma.est=jiang.est[,-(1:p)], identLimit))
}

#' Adaptation of cqrIS Method to a Functional Censored Quantile Regression Model (The Interface)
#'
#' An interface to fit a censored quantile regression model to the provided dataset where
#' the covariates include both non-functional and functional types. The parameters
#' are obtained by solving an induced-smoothed version of the estimating equation used in
#' Jiang et al. (2020). It also provides an interface for variance estimation via multiplier
#' boostrap approach.
#'
#' @param Z Non-functional covariates
#' @param W Transformed functional covariates
#' @param X Survival time, not log-transformed
#' @param cen Censoring indicator, 1 for failure, 0 for censored
#' @param xi Weights for samples, default: a 1-vector
#' @param tau.max Upper bound of estimated quantile
#' @param grid The grid size
#' @param tol Norm of the estimate for the update to be considered convergent
#' @param diverge Norm of the estimate for the update to be considered divergent
#' @param maxiter Maximum number of iteration
#'
#' @import MASS
#' @import quantreg
#' @import stats
#'
#' @return This function returns a list of lists with each list containing three elements:
#' \itemize{
#'   \item beta.ISest, the estimates for the non-functional coefficients
#'   \item gamma.ISest, the estimates for the functional coefficients
#'   \item idenLim, the identifiability limit for the data
#' }
#'
#' @references Cai, Z. and Sit, T. (2020+),
#' "Censored Quantile Regression with Induced Smoothing,"
#' \emph{Working Paper}.

funIS_ <- function(Z, W, X, cen, xi=rep(1,length(X)), tau.max, grid, tol, diverge, maxiter){
  R <- 1e10
  qtile.list <- seq(0, tau.max, by = grid)
  n <- dim(Z)[1]
  p <- dim(Z)[2]
  d <- dim(W)[2]
  H.qtile <- -log(1-qtile.list)
  diff.H <- diff(H.qtile)
  identLimit <- max(qtile.list)

  jiang.res <- estjiang_(Z, W, X, cen, tau.max=tau.max, grid=grid, diverge=diverge)
  jiang.est <- cbind(jiang.res[[1]],jiang.res[[2]])

  # funIS estimation
  funIS.est <- array(NA, c(length(qtile.list), p+d), dimnames = list(qtile.list, c(paste0("beta",1:p),paste0("gamma",1:d))))
  funIS.est[1,] <- jiang.est[1,]
  Sigma <- diag(c(rep(n^(-1),p),rep((n/(d-2))^(-1),d)))#n^(-1) * diag(p)
  minus <- (log(X) >= -Inf) * diff.H[1]
  bZ <- cbind(Z, W)

  for (i in 1:(length(qtile.list)-1)){
    qtile <- qtile.list[i+1]

    funIS.est.tmp <- jiang.est[i+1,]

    loop <- 1
    iter <- 0
    alter <- 0
    while (loop){
      iter <- iter + 1
      funIS.est.now <- funIS.est.tmp

      wAn <- dnorm((bZ%*%funIS.est.now-log(X))/sqrt(diag(bZ%*%Sigma%*%t(bZ))))*(xi*cen)/sqrt(diag(bZ%*%Sigma%*%t(bZ)))
      An <- crossprod(bZ, bZ*as.numeric(wAn))/n

      Sn <- n^(-1)*t(xi*bZ)%*%(cen*pnorm((bZ%*%funIS.est.now-log(X))/sqrt(diag(bZ%*%Sigma%*%t(bZ))))-minus)

      # update beta and Sigma
      funIS.est.tmp <- funIS.est.now - MASS::ginv(An)%*%Sn
      if (anyNA(funIS.est.tmp) | any(is.infinite(funIS.est.tmp)) | t(funIS.est.tmp)%*%funIS.est.tmp>diverge){
        funIS.est.tmp <- funIS.est[i, ]
        loop <- 0
        warning(paste0("funfunIS procedure does not converge at q = ", qtile))
        idenLim <- min(idenLim, qtile)
      } else if (iter > maxiter | norm(funIS.est.now-funIS.est.tmp,type = "2")<tol){
        loop <- 0
      }
    }
    funIS.est[i+1,] <- funIS.est.tmp
    minus = minus + (log(X) >= bZ %*% funIS.est[i+1,]) * diff.H[i+1]
  }
  return(list(beta.est=funIS.est[,1:p], gamma.est=funIS.est[,-(1:p)], identLimit))
}




#' A Functional Censored Quantile Regression Model
#'
#' This function fits a censored quantile regression model to the provided dataset where
#' the covariates include both non-functional and functional types. The parameters
#' are obtained by minimizing an equivalent objective $l_1$ function of an estimating
#' equation.
#'
#' @param Z A numeric matrix. Non-functional covariates
#' @param W A numeric matrix. Transformed functional covariates
#' @param X A numeric vector. Survival time, not log-transformed
#' @param cen A numeric vector. Censoring indicator, 1 for failure, 0 for censored
#' @param tau.max A numeric number. Upper bound of estimated quantile, default: 0.8
#' @param grid A numeric number. The grid size, default: 0.01
#' @param diverge A numeric number. Norm of the estimate for the update to be considered divergent, default: 1e4
#' @param boot Logical. Whether to perform variance estimation, default: FALSE
#' @param nboot A numeric number. The number of bootstrap samples if boot=TRUE, default: 250
#'
#' @import quantreg
#' @export
#'
#' @examples
#' \donttest{
#' dat <- funct.sam200.v10.nb5.cen20
#' res.jiang <- estjiang(Z=dat[,c(3:4)], W=dat[,c(5:9)], X=dat[,1], cen=dat[,2])
#' }
#'
#' @return This function returns a list of lists with each list containing three elements:
#' \itemize{
#'   \item beta.est, a matrix giving the estimates for the non-functional coefficients. If boot=FALSE, it will
#'   contain p columns giving the coefficient estimates where p is the dimension of the non-functional covariate;
#'   if boot=TRUE, it will contain 2*p columns giving the coefficient estimates in the first p columns and the
#'   estimated standard errors in the remaining columns.
#'   \item gamma.est, a matrix giving the estimates for the functional coefficients. If boot=FALSE, it will
#'   contain d columns giving the coefficient estimates where d is the dimension of the basis-transformed functional
#'   covariate; if boot=TRUE, it will contain 2*d columns giving the coefficient estimates in the first d columns
#'   and the estimated standard errors in the remaining columns.
#'   \item idenLim, a numeric giving the identifiability limit for the data
#' }
#'
#' @references Jiang, F., Cheng, Q., Yin, G. and Shen, H. (2020),
#' "Generalizing Quantile Regression for Counting Processes with Applications to Recurrent Events,"
#' \emph{Journal of the American Statistical Association}, \bold{115}, 931-944.
#' @seealso Cai, Z. and Sit, T. (2020+),
#' "Censored Quantile Regression with Induced Smoothing,"
#' \emph{Working Paper}.

estjiang <- function(Z, W, X, cen, tau.max=0.8, grid=0.01, diverge=1e4, boot=FALSE, nboot=250){
  if (!is.numeric(Z) |
      !is.numeric(W) |
      !is.numeric(X) |
      !is.numeric(cen) |
      !is.numeric(tau.max) |
      !is.numeric(grid) |
      !is.numeric(diverge) |
      !is.logical(boot) |
      !is.numeric(nboot)){
    stop("Incorrect variable type(s). Try \"?estjiang\" for more information.")
  }

  if (!is.matrix(Z) |
      !is.matrix(W)){
    stop("Z or W is of incorrect type. Make sure it is a matrix. Try \"?estjiang\" for more information.")
  }

  if (!is.vector(X) |
      !is.vector(cen) |
      !is.vector(tau.max) |
      !is.vector(grid) |
      !is.vector(diverge) |
      !is.vector(boot) |
      !is.vector(nboot)){
    stop("Any of the variables X, cen, tau.max, grid, diverge, boot, nboot is of incorrect type. Make sure they are vectors. Try \"?estjiang\" for more information.")
  } else {
    if(dim(Z)[1] != dim(W)[1]){
      stop("Z or W is of incorrect size. Make sure the length of its first dimension is equal to the sample size.")
    }

    if(length(X) != dim(Z)[1] |
       length(cen) != dim(Z)[1]){
      stop("X or cen is of incorrect length. Make sure it is equal to the sample size.")
    }

    if(length(tau.max)!= 1 |
       length(grid)!= 1 |
       length(diverge)!= 1 |
       length(boot)!= 1 |
       length(nboot)!= 1){
      stop("Any of the variables tau.max, grid, diverge, boot, nboot is of incorrect length. Make sure they are of length 1.")
    }
  }

  if (tau.max < 0.1 | tau.max > 1){
    stop("The specified upper limit is either too low or is larger than 1. Make sure it is in the range (0.1, 1].")
  }

  if (1/grid<20){
    warning("The number of grids might be too small to produce stable results.")
  }

  if (boot ==TRUE & nboot < 30){
    warning("Bootstrap times might be too small to produce reliable standard error estimates.")
  }



  n <- dim(Z)[1]

  jiang.res <- estjiang_(Z, W, X, cen, tau.max=tau.max, grid=grid, diverge=diverge)
  jiang.beta.res <- jiang.res[[1]]
  jiang.gamma.res <- jiang.res[[2]]
  if (jiang.res[[3]]!=tau.max){
    warning(paste0("The identifiability limit might not reach the provided upper bound. Suggested limit: ", jiang.res[[3]]))
  }

  if (boot==FALSE){
    jiang.ret <- list(beta.est=jiang.beta.res, gamma.est=jiang.gamma.res, idenLim=jiang.res[[3]])
  } else {
    jiang.beta.boot.res <- array(NA,dim=c(dim(jiang.beta.res),nboot),dimnames=list(dimnames(jiang.beta.res),1:nboot))
    jiang.gamma.boot.res <- array(NA,dim=c(dim(jiang.gamma.res),nboot),dimnames=list(dimnames(jiang.gamma.res),1:nboot))
    for (iboot in 1:nboot){
      jiang.boot.res <- estjiang_(Z, W, X, cen, xi=rexp(n), tau.max=tau.max, grid=grid, diverge=diverge)
      jiang.beta.boot.res[,,iboot] <- jiang.boot.res[[1]]
      jiang.gamma.boot.res[,,iboot] <- jiang.boot.res[[1]]
    }
    jiang.beta.sd <- apply(jiang.beta.boot.res,c(1,2),sd,na.rm=TRUE)
    jiang.gamma.sd <- apply(jiang.gamma.boot.res,c(1,2),sd,na.rm=TRUE)
    jiang.ret <- list(beta.est=cbind(jiang.beta.res,jiang.beta.sd), gamma.est=cbind(jiang.gamma.res,jiang.gamma.boot.res), idenLim=jiang.res[[3]])
  }

  return(jiang.ret)
}

#' Adaptation of cqrIS Method to a Functional Censored Quantile Regression Model
#'
#' This function fits a censored quantile regression model to the provided dataset where
#' the covariates include both non-functional and functional types. The parameters
#' are obtained by solving an induced-smoothed version of the estimating equation used in
#' Jiang et al. (2020). It also provides an interface for variance estimation via multiplier
#' boostrap approach.
#'
#' @param Z A numeric matrix. Non-functional covariates
#' @param W A numeric matrix. Transformed functional covariates
#' @param X A numeric vector. Survival time, not log-transformed
#' @param cen A numeric vector. Censoring indicator, 1 for failure, 0 for censored
#' @param tau.max A numeric number. Upper bound of estimated quantile, default: 0.8
#' @param grid A numeric number. The grid size, default: 0.01
#' @param tol A numeric number. Norm of the estimate for the update to be considered convergent, default: 1e-4
#' @param diverge A numeric number. Norm of the estimate for the update to be considered divergent, default: 1e4
#' @param maxiter A numeric number. Maximum number of iteration, default: 200
#' @param boot Logical. Whether to perform variance estimation, default: FALSE
#' @param nboot A numeric number. The number of bootstrap samples if boot=TRUE, default: 250
#'
#' @import MASS
#' @import quantreg
#' @import stats
#' @export
#'
#' @examples
#' \donttest{
#' dat <- funct.sam200.v10.nb5.cen20
#' res.funIS <- estjiang(Z=dat[,c(3:4)], W=dat[,c(5:9)], X=dat[,1], cen=dat[,2])
#' }
#'
#' @return This function returns a list of lists with each list containing three elements:
#' \itemize{
#'   \item beta.funISest, a matrix giving the estimates for the non-functional coefficients. If boot=FALSE, it will
#'   contain p columns giving the coefficient estimates where p is the dimension of the non-functional covariate;
#'   if boot=TRUE, it will contain 2*p columns giving the coefficient estimates in the first p columns and the
#'   estimated standard errors in the remaining columns.
#'   \item gamma.funISest, a matrix giving the estimates for the functional coefficients. If boot=FALSE, it will
#'   contain d columns giving the coefficient estimates where d is the dimension of the basis-transformed functional
#'   covariate; if boot=TRUE, it will contain 2*d columns giving the coefficient estimates in the first d columns
#'   and the estimated standard errors in the remaining columns.
#'   \item idenLim, a numeric giving the identifiability limit for the data
#' }
#'
#' @references Cai, Z. and Sit, T. (2020+),
#' "Censored Quantile Regression with Induced Smoothing,"
#' \emph{Working Paper}.

funIS <- function(Z, W, X, cen, tau.max=0.8, grid=0.01, tol=1e-4, diverge=1e4, maxiter=200, boot=FALSE, nboot=250){
  if (!is.numeric(Z) |
      !is.numeric(W) |
      !is.numeric(X) |
      !is.numeric(cen) |
      !is.numeric(tau.max) |
      !is.numeric(grid) |
      !is.numeric(tol) |
      !is.numeric(diverge) |
      !is.numeric(maxiter) |
      !is.logical(boot) |
      !is.numeric(nboot)){
    stop("Incorrect variable type(s). Try \"?funIS\" for more information.")
  }

  if (!is.matrix(Z) |
      !is.matrix(W)){
    stop("Z or W is of incorrect type. Make sure it is a matrix. Try \"?funIS\" for more information.")
  }

  if (!is.vector(X) |
      !is.vector(cen) |
      !is.vector(tau.max) |
      !is.vector(grid) |
      !is.vector(tol) |
      !is.vector(diverge) |
      !is.vector(maxiter) |
      !is.vector(boot) |
      !is.vector(nboot)){
    stop("Any of the variables X, cen, tau.max, grid, tol, diverge, maxiter, boot, nboot is of incorrect type. Make sure they are vectors. Try \"?funIS\" for more information.")
  } else {
    if(dim(Z)[1] != dim(W)[1]){
      stop("Z or W is of incorrect size. Make sure the length of its first dimension is equal to the sample size.")
    }

    if(length(X) != dim(Z)[1] |
       length(cen) != dim(Z)[1]){
      stop("X or cen is of incorrect length. Make sure it is equal to the sample size.")
    }

    if(length(tau.max)!= 1 |
       length(grid)!= 1 |
       length(tol)!= 1 |
       length(diverge)!= 1 |
       length(maxiter)!= 1 |
       length(boot)!= 1 |
       length(nboot)!= 1){
      stop("Any of the variables tau.max, grid, tol, diverge, maxiter, boot, nboot is of incorrect length. Make sure they are of length 1.")
    }
  }

  if (tau.max < 0.1 | tau.max > 1){
    stop("The specified upper limit is either too low or is larger than 1. Make sure it is in the range (0.1, 1].")
  }

  if (1/grid<20){
    warning("The number of grids might be too small to produce stable results.")
  }

  if (boot ==TRUE & nboot < 30){
    warning("Bootstrap times might be too small to produce reliable standard error estimates.")
  }



  n <- dim(Z)[1]

  funIS.res <- funIS_(Z, W, X, cen, tau.max=tau.max, grid=grid, tol=tol, diverge=diverge, maxiter=maxiter)
  funIS.beta.res <- funIS.res[[1]]
  funIS.gamma.res <- funIS.res[[2]]
  if (funIS.res[[3]]!=tau.max){
    warning(paste0("The identifiability limit might not reach the provided upper bound. Suggested limit: ", funIS.res[[3]]))
  }

  if (boot==FALSE){
    funIS.ret <- list(beta.funISest=funIS.beta.res, gamma.funISest=funIS.gamma.res, idenLim=funIS.res[[3]])
  } else {
    funIS.beta.boot.res <- array(NA,dim=c(dim(funIS.beta.res),nboot),dimnames=list(dimnames(funIS.beta.res),1:nboot))
    funIS.gamma.boot.res <- array(NA,dim=c(dim(funIS.gamma.res),nboot),dimnames=list(dimnames(funIS.gamma.res),1:nboot))
    for (iboot in 1:nboot){
      funIS.boot.res <- funIS_(Z, W, X, cen, xi=rexp(n), tau.max=tau.max, grid=grid, tol=tol, diverge=diverge, maxiter=maxiter)
      funIS.beta.boot.res[,,iboot] <- funIS.boot.res[[1]]
      funIS.gamma.boot.res[,,iboot] <- funIS.boot.res[[1]]
    }
    funIS.beta.sd <- apply(funIS.beta.boot.res,c(1,2),sd,na.rm=TRUE)
    funIS.gamma.sd <- apply(funIS.gamma.boot.res,c(1,2),sd,na.rm=TRUE)
    funIS.ret <- list(beta.funISest=cbind(funIS.beta.res,funIS.beta.sd), gamma.funISest=cbind(funIS.gamma.res,funIS.gamma.boot.res), idenLim=funIS.res[[3]])
  }

  return(funIS.ret)
}
