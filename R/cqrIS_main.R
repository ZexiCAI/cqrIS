#' Censored Quantile Regression (The Interface)
#'
#' This function fits a censored quantile regression model to the provided survival dataset.
#' The parameters are obtained by minimizing an equivalent objective $l_1$ function of an
#' estimating equation. It also provides an interface for variance estimation via multiplier
#' bootstrap.
#'
#' @param Z Covariates
#' @param X Survival time, not log-transformed
#' @param cen Censoring indicator, 1 for failure, 0 for censored
#' @param xi Weights for samples, default: a 1-vector
#' @param tau.max Upper bound of estimated quantile
#' @param grid The grid size
#' @param diverge Norm of the estimate for the update to be considered divergent
#'
#' @import quantreg
#'
#' @return This function returns a list of lists with each list containing two elements:
#' \itemize{
#'   \item beta.est, the estimates for the parameters
#'   \item idenLim, the identifiability limit for the data
#' }
#'
#' @references Peng, L. and Huang, Y. (2008),
#' "Survival Analysis with Quantile Regression Models,"
#' \emph{Journal of the American Statistical Association}, \bold{103}, 637-649.
#' @seealso Cai, Z. and Sit, T. (2020+),
#' "Censored Quantile Regression with Induced Smoothing,"
#' \emph{Working Paper}.

estPH_ <- function(Z, X, cen, xi=rep(1,length(X)), tau.max, grid, diverge){
  R <- 1e10
  qtile.list <- seq(0, tau.max, by = grid)
  n <- dim(Z)[1]
  p <- dim(Z)[2]
  H.qtile <- -log(1-qtile.list)
  diff.H <- diff(H.qtile)
  identLimit <- max(qtile.list)

  resp <- c(xi*cen*log(X), R, R)
  beta.est <- array(NA, c(length(qtile.list), p), dimnames = list(qtile.list, paste0("beta",1:p)))
  beta.est[1,] <- c(-R,rep(0,p-1))
  minus <- (log(X) >= Z %*% beta.est[1,]) * diff.H[1]

  covr <- array(NA, dim = c(n+2,p))
  covr[1:n,] <- cen*(xi*Z)
  covr[dim(covr)[1]-1,] <- crossprod(xi*Z,-cen)
  # estimation
  for (i in 1:(length(qtile.list)-1)){
    qtile <- qtile.list[i+1]

    covr[dim(covr)[1],] <- crossprod(xi*Z,2*minus)
    beta.est[i+1,] <- rq(resp~0+covr,tau=0.5)$coef

    if (anyNA(beta.est[i+1,]) | any(is.infinite(beta.est[i+1,])) | t(beta.est[i+1,])%*%beta.est[i+1,]>diverge){
      beta.est[i+1,] <- beta.est[i, ]
      identLimit <- min(identLimit, qtile)
      warning(paste0("PH procedure does not converge at tau = ", qtile))
    }

    minus = minus + (log(X) >= Z %*% beta.est[i+1,]) * diff.H[i+1]
  }
  return(list(beta.est, identLimit))
}

#' Censored Quantile Regression with Induced Smoothing (The Interface)
#'
#' This function fits a censored quantile regression model to the provided survival dataset.
#' The parameters are obtained by solving an induced-smoothed version of the estimating
#' equation used in Peng and Huang (2008). It also provides an interface for variance estimation
#' via multiplier boostrap approach.
#'
#' @param Z Covariates
#' @param X Survival time, not log-transformed
#' @param cen Censoring indicator, 1 for failure, 0 for censored
#' @param xi Weights for samples, default: a 1-vector
#' @param init Whether Peng and Huang (2008) estimates are used as initial guess
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
#' @return This function returns a list of lists with each list containing two elements:
#' \itemize{
#'   \item beta.ISest, the estimates for the parameters
#'   \item idenLim, the identifiability limit for the data
#' }
#'
#' @references Cai, Z. and Sit, T. (2020+),
#' "Censored Quantile Regression with Induced Smoothing,"
#' \emph{Working Paper}.

cqrIS_ <- function(Z, X, cen, xi=rep(1,length(X)), init, tau.max, grid, tol, diverge, maxiter){
  R <- 1e10
  qtile.list <- seq(0, tau.max, by = grid)
  H.qtile <- -log(1-qtile.list)
  n <- dim(Z)[1]
  p <- dim(Z)[2]
  diff.H <- diff(H.qtile)
  identLimit <- max(qtile.list)
  trun <- 5

  resp <- c(xi*cen*log(X), R, R)
  beta.est <- array(NA, c(length(qtile.list), p), dimnames = list(qtile.list, paste0("beta",1:p)))
  beta.est[1,] <- c(-R,rep(0,p-1))
  minus_PH <- (log(X) >= Z %*% beta.est[1,]) * diff.H[1]

  covr <- array(NA, dim = c(n+2,p))
  covr[1:n,] <- cen*(xi*Z)
  covr[dim(covr)[1]-1,] <- crossprod(xi*Z,-cen)

  # PH initial values
  for (i in 1:(length(qtile.list)-1)){
    qtile <- qtile.list[i+1]

    covr[dim(covr)[1],] <- crossprod(xi*Z,2*minus_PH)
    beta.est[i+1,] <- rq(resp~0+covr,tau=0.5)$coef

    minus_PH = minus_PH + (log(X) >= Z %*% beta.est[i+1,]) * diff.H[i+1]
  }

  # IS estimation
  beta.ISest <- array(NA, c(length(qtile.list), p), dimnames = list(qtile.list, paste0("beta",1:p)))
  beta.ISest[1:trun,] <- beta.est[1:trun,]
  Sigma <- n^(-1) * diag(p)

  minus <- rowSums((matrix(rep(log(X),trun),ncol=trun)>=Z%*%t(array(beta.ISest[1:trun,],dim=c(trun,p))))*matrix(rep(diff.H[1:trun],n),nrow=n,byrow=T))

  for (i in trun:(length(qtile.list)-1)){
    qtile <- qtile.list[i+1]

    if (init==TRUE){
      beta.ISest.tmp <- beta.est[i+1,]
    } else{
      beta.ISest.tmp <- beta.ISest[i,]
    }

    loop <- 1
    iter <- 0
    alter <- 0
    while (loop){
      iter <- iter + 1
      beta.ISest.now <- beta.ISest.tmp

      wAn <- dnorm((Z%*%beta.ISest.now-log(X))/sqrt(diag(Z%*%Sigma%*%t(Z))))*(xi*cen)/sqrt(diag(Z%*%Sigma%*%t(Z)))
      An <- crossprod(Z, Z*as.numeric(wAn))/n

      Sn <- n^(-1)*t(xi*Z)%*%(cen*pnorm((Z%*%beta.ISest.now-log(X))/sqrt(diag(Z%*%Sigma%*%t(Z))))-minus)

      # update beta and Sigma
      beta.ISest.tmp <- beta.ISest.now - MASS::ginv(An)%*%Sn
      if ((anyNA(beta.ISest.tmp) | any(is.infinite(beta.ISest.tmp)) | t(beta.ISest.tmp)%*%beta.ISest.tmp>diverge) & alter == 0){
        warning(paste0("IS procedure does not converge for the current initial estimate at tau = ", qtile, ". Change to another choice."))
        beta.ISest.tmp <- beta.est[i+1,]*(init==FALSE) + beta.ISest[i,]*(init==TRUE)
        iter <- 0
        alter <- 1
      } else if (iter > maxiter | norm(beta.ISest.now-beta.ISest.tmp,type = "2")<tol){
        loop <- 0
      } else if (alter == 1) {
        beta.ISest.tmp <- beta.ISest[i, ]
        identLimit <- min(identLimit, qtile)
        loop <- 0
        warning(paste0("IS procedure does not converge at tau = ", qtile))
      }
    }
    beta.ISest[i+1,] <- beta.ISest.tmp
    minus = minus + (log(X) >= Z %*% beta.ISest[i+1,]) * diff.H[i+1]
  }
  return(list(beta.ISest, identLimit))
}

#' Censored Quantile Regression with Induced Smoothing, the Second Formulation (The Interface)
#'
#' This function fits a censored quantile regression model to the provided survival dataset.
#' The parameters are obtained by solving an induced-smoothed version of the estimating
#' equation used in Peng and Huang (2008). It also provides an interface for variance estimation
#' via multiplier boostrap approach.
#'
#' @param Z Covariates
#' @param X Survival time, not log-transformed
#' @param cen Censoring indicator, 1 for failure, 0 for censored
#' @param xi Weights for samples, default: a 1-vector
#' @param init Whether Peng and Huang (2008) estimates are used as initial guess
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
#' @return This function returns a list of lists with each list containing two elements:
#' \itemize{
#'   \item beta.sISest, the estimates for the parameters
#'   \item idenLim, the identifiability limit for the data
#' }
#'
#' @references Cai, Z. and Sit, T. (2020+),
#' "Censored Quantile Regression with Induced Smoothing,"
#' \emph{Working Paper}.

cqrsIS_ <- function(Z, X, cen, xi=rep(1,length(X)), init, tau.max, grid, tol, diverge, maxiter){
  R <- 1e10
  qtile.list <- seq(0, tau.max, by = grid)
  n <- dim(Z)[1]
  p <- dim(Z)[2]
  H.qtile <- -log(1-qtile.list)
  diff.H <- diff(H.qtile)
  identLimit <- max(qtile.list)
  trun <- 5

  resp <- c(cen*log(X), R, R)
  beta.est <- array(NA, c(length(qtile.list), p), dimnames = list(qtile.list, paste0("beta",1:p)))
  beta.est[1,] <- c(-R,rep(0,p-1))
  minus_PH <- (log(X) >= Z %*% beta.est[1,]) * diff.H[1]

  covr <- array(NA, dim = c(n+2,p))
  covr[1:n,] <- cen*Z
  covr[dim(covr)[1]-1,] <- crossprod(Z,-cen)

  # PH initial values
  for (i in 1:(length(qtile.list)-1)){
    qtile <- qtile.list[i+1]

    covr[dim(covr)[1],] <- crossprod(Z,2*minus_PH)
    beta.est[i+1,] <- rq(resp~0+covr,tau=0.5)$coef

    minus_PH = minus_PH + (log(X) >= Z %*% beta.est[i+1,]) * diff.H[i+1]
  }

  # sIS estimation
  beta.sISest <- array(NA, c(length(qtile.list), p), dimnames = list(qtile.list, paste0("beta",1:p)))
  beta.sISest[1:trun,] <- beta.est[1:trun,]
  Sigma <- n^(-1) * diag(p)

  minus <- rep(0, n)
  plus <- rep(0, n)

  for (i in 1:(trun-1)){
    minus <- minus + pnorm((-Z%*%beta.sISest[i,]+log(X))/sqrt(diag(Z%*%Sigma%*%t(Z))))*diff.H[i]
    plus <- plus + dnorm((-Z%*%beta.sISest[i,]+log(X))/sqrt(diag(Z%*%Sigma%*%t(Z))))/sqrt(diag(Z%*%Sigma%*%t(Z)))*diff.H[i]
  }

  for (i in trun:(length(qtile.list)-1)){
    qtile <- qtile.list[i+1]

    if (init==TRUE){
      beta.sISest.tmp <- beta.est[i+1,]
    } else{
      beta.sISest.tmp <- beta.sISest[i,]
    }

    minus <- minus + pnorm((-Z%*%beta.sISest[i,]+log(X))/sqrt(diag(Z%*%Sigma%*%t(Z))))*diff.H[i]
    plus <- plus + dnorm((-Z%*%beta.sISest[i,]+log(X))/sqrt(diag(Z%*%Sigma%*%t(Z))))/sqrt(diag(Z%*%Sigma%*%t(Z)))*diff.H[i]

    loop <- 1
    iter <- 0
    alter <- 0
    while (loop){
      iter <- iter + 1
      beta.sISest.now <- beta.sISest.tmp

      wAn <- dnorm((Z%*%beta.sISest.now-log(X))/sqrt(diag(Z%*%Sigma%*%t(Z))))*(xi*cen)/sqrt(diag(Z%*%Sigma%*%t(Z)))+plus
      An <- crossprod(Z, Z*as.numeric(wAn))/n

      Sn <- n^(-1)*t(xi*Z)%*%(cen*pnorm((Z%*%beta.sISest.now-log(X))/sqrt(diag(Z%*%Sigma%*%t(Z))))-minus)

      # update beta and Sigma
      beta.sISest.tmp <- beta.sISest.now - MASS::ginv(An)%*%Sn
      if ((anyNA(beta.sISest.tmp) | any(is.infinite(beta.sISest.tmp)) | t(beta.sISest.tmp)%*%beta.sISest.tmp>diverge) & alter == 0){
        warning(paste0("sIS procedure does not converge for the current initial estimate at tau = ", qtile, ". Change to another choice."))
        beta.sISest.tmp <- beta.est[i+1,]*(init==FALSE) + beta.sISest[i,]*(init==TRUE)
        iter <- 0
        alter <- 1
      } else if (iter > maxiter | norm(beta.sISest.now-beta.sISest.tmp,type = "2")<tol){
        loop <- 0
      } else if (alter == 1) {
        beta.sISest.tmp <- beta.sISest[i, ]
        identLimit <- min(identLimit, qtile)
        loop <- 0
        warning(paste0("sIS procedure does not converge at tau = ", qtile))
      }
    }
    beta.sISest[i+1,] <- beta.sISest.tmp
  }

  return(list(beta.sISest, identLimit))
}



#' Censored Quantile Regression
#'
#' This function fits a censored quantile regression model to the provided survival dataset.
#' The parameters are obtained by minimizing an equivalent objective $l_1$ function of an
#' estimating equation.
#'
#' @param Z A numeric matrix. Covariates
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
#' dat <- ordin.sam200.cen25.homo
#' res.PH <- estPH(Z=dat[,-c(1:2)], X=dat[,1], cen=dat[,2])
#' }
#'
#' @return This function returns a list of lists with each list containing two elements:
#' \itemize{
#'   \item beta.est, a matrix giving the estimates for the parameters. If boot=FALSE, it will contain p columns
#'   giving the coefficient estimates where p is the dimension of the covariate; if boot=TRUE, it will contain
#'   2*p columns giving the coefficient estimates in the first p columns and the estimated standard errors in
#'   the remaining columns.
#'   \item idenLim, a numeric giving the identifiability limit for the data
#' }
#'
#' @references Peng, L. and Huang, Y. (2008),
#' "Survival Analysis with Quantile Regression Models,"
#' \emph{Journal of the American Statistical Association}, \bold{103}, 637-649.
#' @seealso Cai, Z. and Sit, T. (2020+),
#' "Censored Quantile Regression with Induced Smoothing,"
#' \emph{Working Paper}.

estPH <- function(Z, X, cen, tau.max=0.8, grid=0.01, diverge=1e4, boot=FALSE, nboot=250){
  if (!is.numeric(Z) |
      !is.numeric(X) |
      !is.numeric(cen) |
      !is.numeric(tau.max) |
      !is.numeric(grid) |
      !is.numeric(diverge) |
      !is.logical(boot) |
      !is.numeric(nboot)){
    stop("Incorrect variable type(s). Try \"?estPH\" for more information.")
  }

  if (!is.matrix(Z)){
    stop("Z is of incorrect type. Make sure it is a matrix. Try \"?estPH\" for more information.")
  }

  if (!is.vector(X) |
      !is.vector(cen) |
      !is.vector(tau.max) |
      !is.vector(grid) |
      !is.vector(diverge) |
      !is.vector(boot) |
      !is.vector(nboot)){
    stop("Any of the variables X, cen, tau.max, grid, diverge, boot, nboot is of incorrect type. Make sure they are vectors. Try \"?estPH\" for more information.")
  } else {
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

  PH.res <- estPH_(Z, X, cen, tau.max=tau.max, grid=grid, diverge=diverge)
  if (PH.res[[2]]!=tau.max){
    warning(paste0("The identifiability limit might not reach the provided upper bound. Suggested limit: ", PH.res[[2]]))
  }

  if (boot==FALSE){
    PH.ret <- list(beta.est=PH.res[[1]], idenLim=PH.res[[2]])
  } else {
    PH.boot.res <- array(NA,dim=c(dim(PH.res[[1]]),nboot),dimnames=list(dimnames(PH.res[[1]]),1:nboot))
    for (iboot in 1:nboot){
      PH.boot.res[,,iboot] <- estPH_(Z, X, cen, xi=rexp(n), tau.max=tau.max, grid=grid, diverge=diverge)[[1]]
    }
    PH.res.sd <- apply(PH.boot.res,c(1,2),sd,na.rm=TRUE)
    PH.ret <- list(beta.est=cbind(PH.res[[1]],PH.res.sd), idenLim=PH.res[[2]])
  }

  return(PH.ret)
}

#' Censored Quantile Regression with Induced Smoothing
#'
#' This function fits a censored quantile regression model to the provided survival dataset.
#' The parameters are obtained by solving an induced-smoothed version of the estimating
#' equation used in Peng and Huang (2008).
#'
#' @param Z A numeric matrix. Covariates
#' @param X A numeric vector. Survival time, not log-transformed
#' @param cen A numeric vector. Censoring indicator, 1 for failure, 0 for censored
#' @param init Logical. Whether Peng and Huang (2008) estimates are used as initial guess, default: TRUE for yes
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
#' dat <- ordin.sam200.cen25.homo
#' res.cqrIS <- cqrIS(Z=dat[,-c(1:2)], X=dat[,1], cen=dat[,2])
#' }
#'
#' @return This function returns a list of lists with each list containing two elements:
#' \itemize{
#'   \item beta.ISest, a matrix giving the estimates for the parameters. If boot=FALSE, it will contain p columns
#'   giving the coefficient estimates where p is the dimension of the covariate; if boot=TRUE, it will contain
#'   2*p columns giving the coefficient estimates in the first p columns and the estimated standard errors in
#'   the remaining columns.
#'   \item idenLim, a numeric giving the identifiability limit for the data
#' }
#'
#' @references Cai, Z. and Sit, T. (2020+),
#' "Censored Quantile Regression with Induced Smoothing,"
#' \emph{Working Paper}.

cqrIS <- function(Z, X, cen, init=TRUE, tau.max=0.8, grid=0.01, tol=1e-4, diverge=1e4, maxiter=200, boot=FALSE, nboot=250){
  if (!is.numeric(Z) |
      !is.numeric(X) |
      !is.numeric(cen) |
      !is.logical(init) |
      !is.numeric(tau.max) |
      !is.numeric(grid) |
      !is.numeric(tol) |
      !is.numeric(diverge) |
      !is.numeric(maxiter) |
      !is.logical(boot) |
      !is.numeric(nboot)){
    stop("Incorrect variable type(s). Try \"?cqrIS\" for more information.")
  }

  if (!is.matrix(Z)){
    stop("Z is of incorrect type. Make sure it is a matrix. Try \"?cqrIS\" for more information.")
  }

  if (!is.vector(X) |
      !is.vector(cen) |
      !is.vector(init) |
      !is.vector(tau.max) |
      !is.vector(grid) |
      !is.vector(tol) |
      !is.vector(diverge) |
      !is.vector(maxiter) |
      !is.vector(boot) |
      !is.vector(nboot)){
    stop("Any of the variables X, cen, init, tau.max, grid, tol, diverge, maxiter, boot, nboot is of incorrect type. Make sure they are vectors. Try \"?cqrIS\" for more information.")
  } else {
    if(length(X) != dim(Z)[1] |
       length(cen) != dim(Z)[1]){
      stop("X or cen is of incorrect length. Make sure it is equal to the sample size.")
    }

    if(length(init)!= 1 |
       length(tau.max)!= 1 |
       length(grid)!= 1 |
       length(tol)!= 1 |
       length(diverge)!= 1 |
       length(maxiter)!= 1 |
       length(boot)!= 1 |
       length(nboot)!= 1){
      stop("Any of the variables init, tau.max, grid, tol, diverge, maxiter, boot, nboot is of incorrect length. Make sure they are of length 1.")
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

  cqrIS.res <- cqrIS_(Z, X, cen, init=init, tau.max=tau.max, grid=grid, tol=tol, diverge=diverge, maxiter=maxiter)
  if (cqrIS.res[[2]]!=tau.max){
    warning(paste0("The identifiability limit might not reach the provided upper bound. Suggested limit: ", cqrIS.res[[2]]))
  }

  if (boot==FALSE){
    cqrIS.ret <- list(beta.ISest=cqrIS.res[[1]], idenLim=cqrIS.res[[2]])
  } else {
    cqrIS.boot.res <- array(NA,dim=c(dim(cqrIS.res[[1]]),nboot),dimnames=list(dimnames(cqrIS.res[[1]]),1:nboot))
    for (iboot in 1:nboot){
      cqrIS.boot.res[,,iboot] <- cqrIS_(Z, X, cen, xi=rexp(n), init=init, tau.max=tau.max, grid=grid, tol=tol, diverge=diverge, maxiter=maxiter)[[1]]
    }
    cqrIS.res.sd <- apply(cqrIS.boot.res,c(1,2),sd,na.rm=TRUE)
    cqrIS.ret <- list(beta.ISest=cbind(cqrIS.res[[1]],cqrIS.res.sd), idenLim=cqrIS.res[[2]])
  }

  return(cqrIS.ret)
}

#' Censored Quantile Regression with Induced Smoothing, the Second Formulation
#'
#' This function fits a censored quantile regression model to the provided survival dataset.
#' The parameters are obtained by solving an induced-smoothed version of the estimating
#' equation used in Peng and Huang (2008).
#'
#' @param Z A numeric matrix. Covariates
#' @param X A numeric vector. Survival time, not log-transformed
#' @param cen A numeric vector. Censoring indicator, 1 for failure, 0 for censored
#' @param init Logical. Whether Peng and Huang (2008) estimates are used as initial guess, default: TRUE for yes
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
#' dat <- ordin.sam200.cen25.homo
#' res.cqrsIS <- cqrsIS(Z=dat[,-c(1:2)], X=dat[,1], cen=dat[,2])
#' }
#'
#' @return This function returns a list of lists with each list containing two elements:
#' \itemize{
#'   \item beta.sISest, a matrix giving the estimates for the parameters. If boot=FALSE, it will contain p columns
#'   giving the coefficient estimates where p is the dimension of the covariate; if boot=TRUE, it will contain
#'   2*p columns giving the coefficient estimates in the first p columns and the estimated standard errors in
#'   the remaining columns.
#'   \item idenLim, a numeric giving the identifiability limit for the data
#' }
#'
#' @references Cai, Z. and Sit, T. (2020+),
#' "Censored Quantile Regression with Induced Smoothing,"
#' \emph{Working Paper}.

cqrsIS <- function(Z, X, cen, init=TRUE, tau.max=0.8, grid=0.01, tol=1e-4, diverge=1e4, maxiter=200, boot=FALSE, nboot=250){
  if (!is.numeric(Z) |
      !is.numeric(X) |
      !is.numeric(cen) |
      !is.logical(init) |
      !is.numeric(tau.max) |
      !is.numeric(grid) |
      !is.numeric(tol) |
      !is.numeric(diverge) |
      !is.numeric(maxiter) |
      !is.logical(boot) |
      !is.numeric(nboot)){
    stop("Incorrect variable type(s). Try \"?cqrsIS\" for more information.")
  }

  if (!is.matrix(Z)){
    stop("Z is of incorrect type. Make sure it is a matrix. Try \"?cqrsIS\" for more information.")
  }

  if (!is.vector(X) |
      !is.vector(cen) |
      !is.vector(init) |
      !is.vector(tau.max) |
      !is.vector(grid) |
      !is.vector(tol) |
      !is.vector(diverge) |
      !is.vector(maxiter) |
      !is.vector(boot) |
      !is.vector(nboot)){
    stop("Any of the variables X, cen, init, tau.max, grid, tol, diverge, maxiter, boot, nboot is of incorrect type. Make sure they are vectors. Try \"?cqrsIS\" for more information.")
  } else {
    if(length(X) != dim(Z)[1] |
       length(cen) != dim(Z)[1]){
      stop("X or cen is of incorrect length. Make sure it is equal to the sample size.")
    }

    if(length(init)!= 1 |
       length(tau.max)!= 1 |
       length(grid)!= 1 |
       length(tol)!= 1 |
       length(diverge)!= 1 |
       length(maxiter)!= 1 |
       length(boot)!= 1 |
       length(nboot)!= 1){
      stop("Any of the variables init, tau.max, grid, tol, diverge, maxiter, boot, nboot is of incorrect length. Make sure they are of length 1.")
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

  cqrsIS.res <- cqrsIS_(Z, X, cen, init=init, tau.max=tau.max, grid=grid, tol=tol, diverge=diverge, maxiter=maxiter)
  if (cqrsIS.res[[2]]!=tau.max){
    warning(paste0("The identifiability limit might not reach the provided upper bound. Suggested limit: ", cqrsIS.res[[2]]))
  }

  if (boot==FALSE){
    cqrsIS.ret <- list(beta.sISest=cqrsIS.res[[1]], idenLim=cqrsIS.res[[2]])
  } else {
    cqrsIS.boot.res <- array(NA,dim=c(dim(cqrsIS.res[[1]]),nboot),dimnames=list(dimnames(cqrsIS.res[[1]]),1:nboot))
    for (iboot in 1:nboot){
      cqrsIS.boot.res[,,iboot] <- cqrsIS_(Z, X, cen, xi=rexp(n), init=init, tau.max=tau.max, grid=grid, tol=tol, diverge=diverge, maxiter=maxiter)[[1]]
    }
    cqrsIS.res.sd <- apply(cqrsIS.boot.res,c(1,2),sd,na.rm=TRUE)
    cqrsIS.ret <- list(beta.sISest=cbind(cqrsIS.res[[1]],cqrsIS.res.sd), idenLim=cqrsIS.res[[2]])
  }

  return(cqrsIS.ret)
}
