#' Recurrent Event Modelling via Generalizing Quantile Regression for Counting Processes (The Interface)
#'
#' An interface to fit a recurrent event model to the provided dataset where the covariates
#' are related to the "expected frequency" via a generalized linear model. The parameters
#' are obtained by minimizing an equivalent objective $l_1$ function of an estimating
#' equation. It also provides an interface for variance estimation via multiplier bootstrap.
#' The codes are adapted from Prof. Limin Peng's website (with modifications).
#'
#' @param T Recurrent event time
#' @param x Covariates for samples
#' @param ID ID of those who experience recurrent events
#' @param grids Number of grids
#' @param U Upper bound of estimation
#' @param L Left censoring variable
#' @param R Right censoring variable
#' @param ZETA Weights for those who experience recurrent events, default: a 1-vector
#' @param zeta Weights for samples, default: a 1-vector
#' @param link Link function, default: "1"
#'
#' @import quantreg
#'
#' @return This function returns a list of lists with each list containing two elements:
#' \itemize{
#'   \item beta.est, the estimates for the parameters
#'   \item idenLim, the identifiability limit for the data
#' }
#'
#' @references Sun, X., Peng, L., Huang, Y. and Lai, H.J. (2016),
#' "Generalizing Quantile Regression for Counting Processes with Applications to Recurrent Events,"
#' \emph{Journal of the American Statistical Association}, \bold{111}, 145-156.
#' @seealso Cai, Z. and Sit, T. (2020+),
#' "Censored Quantile Regression with Induced Smoothing,"
#' \emph{Working Paper}.

estsun_ <- function(T, x, ID, grids, U, L, R, ZETA=rep(1,dim(X)[1]), zeta=rep(1,dim(x)[1]), link){
  X <- x[ID, ]
  N = dim(X)[1]
  p = dim(X)[2]
  n = dim(x)[1]
  idenLim <- U

  rqy = matrix(NA, N+2, 1)
  rqx =  matrix(NA, N+2, p)
  rqw =  matrix(NA, N+2, 1)

  rqy[1:N] = log(T)
  rqy[N+1] = 1e10
  rqy[N+2] = 1e10

  rqw[1:N] = ZETA
  rqw[(N+1):(N+2)] =  1

  rqx[1:N, ] = X
  rqx[N+1, ] = -apply(X*matrix(rep(rqw[1:N, ], times=p), ncol=p), 2, sum)

  beta.est =  matrix(NA, grids, p)

  lc = matrix(0, n, 1)

  for (k in 1:grids) {
    if (k == 1) {
      lc = matrix(0, n, 1)
      if (link == "1") {
        lc[which(L==0)] =  U/grids
      }
      else if (link == "u") {
        lc[which(L==0)] =  (U/grids)^2/2
      }
      else if (link == "exp") {
        lc[which(L==0)] = exp(U/grids)-1
      }
    }
    if (k > 1) {
      tmp = which(L <= exp(x%*%beta.est[k - 1, ]) & exp(x%*%beta.est[k - 1, ])<= R)
      if (link == "1") {
        lc[tmp] = lc[tmp] + U/grids
      }
      else if (link == "u") {
        lc[tmp] = lc[tmp] + (U/grids*k)^2/2 - (U/grids*(k-1))^2/2
      }
      else if (link == "exp") {
        lc[tmp] = lc[tmp] + exp(U/grids * k) - exp(U/grids*(k - 1))
      }
    }
    rqx[N + 2, ] = 2*apply(x*matrix(rep(lc*zeta, times=p), ncol=p), 2, sum)

    est <- rq(rqy ~ rqx - 1, tau = 0.5, weights = rqw)$coefficients
    if (max(abs(est)) < 10^5) {
      beta.est[k, ] = est
    } else {
      if (k!=1){
        beta.est[k, ] = beta.est[k-1, ]
      } else {
        beta.est[k, ] = rep(0,p)
      }
      warning(paste0("Sun procedure does not converge at u = ", U/grids*k))
      idenLim <- min(idenLim, U/grids*k)
    }
  }
  beta.est <- rbind(c(-1e10,rep(0,p-1)),beta.est)
  return(list(beta.est,idenLim))
}

#' Adaptation of cqrIS Method to a Recurrent Event Model (The Interface)
#'
#' An interface to fit a recurrent event model to the provided dataset where the covariates
#' are related to the "expected frequency" via a generalized linear model. The parameters
#' are obtained by solving an induced-smoothed version of the estimating equation used in
#' Sun et al. (2016). It also provides an interface for variance estimation via multiplier
#' boostrap approach.
#'
#' @param T Recurrent event time
#' @param x Covariates for samples
#' @param ID ID of those who experience recurrent events
#' @param grids Number of grids
#' @param U Upper bound of estimation
#' @param L Left censoring variable
#' @param R Right censoring variable
#' @param zeta Weights for samples, default: a 1-vector
#' @param link Link function, default: "1"
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

recurIS_ <- function(T, x, ID, grids, U, L, R, zeta=rep(1,dim(x)[1]), link, tol, diverge, maxiter){
  X <- x[ID, ]
  u.list <- seq(0, U, by=U/grids)
  G.utile <- u.list
  diff.G <- diff(G.utile)

  N = dim(X)[1]
  p = dim(X)[2]
  n = dim(x)[1]
  idenLim <- U

  IDtran <- matrix(0,nrow=n,ncol=N)
  j<-1
  for (iii in 1:n){
    if (j>length(ID)) break;
    while (ID[j]==iii){
      IDtran[iii,j] <- 1
      j <- j+1
      if (j>length(ID)) break;
    }
  }

  beta.est <- estsun_(T=T, x=x, ID=ID, grids=grids, U=U, L=L, R=R, link=link)[[1]]
  # recurIS estimation
  beta.ISest <- array(NA, c(length(u.list), p), dimnames = list(u.list, paste0("beta",1:p)))
  beta.ISest[1,] <- beta.est[1,]
  Sigma <- n^(-1) * diag(p)

  minus <- (log(R) >= x %*% beta.est[1,]) *(x %*% beta.est[1,] >= log(L)) * diff.G[1]
  for (i in 1:(length(u.list)-1)){
    utile <- u.list[i+1]

    beta.ISest.tmp <- beta.est[i+1,]

    loop <- 1
    iter <- 0
    while (loop){
      iter <- iter + 1
      beta.ISest.now <- beta.ISest.tmp

      wAn <- zeta*(IDtran%*%(dnorm((X%*%beta.ISest.now-log(T))/sqrt(diag(X%*%Sigma%*%t(X))))/sqrt(diag(X%*%Sigma%*%t(X)))))
      An <- crossprod(x, x*as.numeric(wAn))/n

      Sn <- n^(-1)*t(zeta*x)%*%(IDtran%*%pnorm((X%*%beta.ISest.now-log(T))/sqrt(diag(X%*%Sigma%*%t(X))))-minus)

      # update beta and Sigma
      beta.ISest.tmp <- beta.ISest.now - MASS::ginv(An)%*%Sn
      if (anyNA(beta.ISest.tmp) | any(is.infinite(beta.ISest.tmp)) | t(beta.ISest.tmp)%*%beta.ISest.tmp>diverge){
        beta.ISest.tmp <- beta.ISest[i, ]
        loop <- 0
        warning(paste0("recurIS procedure does not converge at u = ", utile))
        idenLim <- min(idenLim, utile)
      } else if (iter > maxiter | norm(beta.ISest.now-beta.ISest.tmp,type = "2")<tol){
        loop <- 0
      }
    }
    beta.ISest[i+1,] <- beta.ISest.tmp
    minus <- minus + (log(R) >= x %*% beta.ISest[i+1,]) *(x %*% beta.ISest[i+1,] >= log(L)) * diff.G[i+1]
  }
  return(list(beta.ISest,idenLim))
}





#' Recurrent Event Modelling via Generalizing Quantile Regression for Counting Processes
#'
#' This function fits a recurrent event model to the provided dataset where the covariates
#' are related to the "expected frequency" via a generalized linear model. The parameters
#' are obtained by minimizing an equivalent objective $l_1$ function of an estimating
#' equation.
#'
#' @param T A numeric vector. Recurrent event time
#' @param x A numeric matrix. Covariates for samples
#' @param ID A numeric vector. ID of those who experience recurrent events
#' @param grids A numeric number. Number of grids
#' @param U A numeric number. Upper bound of estimation
#' @param L A numeric vector. Left censoring variable
#' @param R A numeric vector. Right censoring variable
#' @param link A character string. Link function, default: "1"
#' @param boot Logical. Whether to perform variance estimation, default: FALSE
#' @param nboot A numeric number. The number of bootstrap samples if boot=TRUE, default: 250
#'
#' @import quantreg
#' @export
#'
#' @examples
#' \donttest{
#' evt <- recur[[1]]
#' sam <- recur[[2]]
#' res.estsun <- estsun(T=evt[,1], x=sam[,-c(1:2)], ID=evt[,2], grids=100, U=4, L=sam[,1], R=sam[,2])
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
#' @references Sun, X., Peng, L., Huang, Y. and Lai, H.J. (2016),
#' "Generalizing Quantile Regression for Counting Processes with Applications to Recurrent Events,"
#' \emph{Journal of the American Statistical Association}, \bold{111}, 145-156.
#' @seealso Cai, Z. and Sit, T. (2020+),
#' "Censored Quantile Regression with Induced Smoothing,"
#' \emph{Working Paper}.

estsun <- function(T, x, ID, grids, U, L, R, link="1", boot=FALSE, nboot=250){
  if (!is.numeric(T) |
      !is.numeric(x) |
      !is.numeric(ID) |
      !is.numeric(grids) |
      !is.numeric(U) |
      !is.numeric(L) |
      !is.numeric(R) |
      !is.character(link) |
      !is.logical(boot) |
      !is.numeric(nboot)){
    stop("Incorrect variable type(s). Try \"?estsun\" for more information.")
  }

  if (!is.matrix(x)){
    stop("x is of incorrect type. Make sure it is a matrix. Try \"?estsun\" for more information.")
  }

  if (!is.vector(T) |
      !is.vector(ID) |
      !is.vector(grids) |
      !is.vector(U) |
      !is.vector(L) |
      !is.vector(R) |
      !is.vector(link) |
      !is.vector(boot) |
      !is.vector(nboot)){
    stop("Any of the variables T, ID, grids, U, L, R, link, boot, nboot is of incorrect type. Make sure they are vectors. Try \"?estsun\" for more information.")
  } else {
    if(length(T) != length(ID)){
      stop("T or ID is of incorrect length. Make sure their lengths are equal.")
    }

    if (length(L) != dim(x)[1] |
        length(R) != dim(x)[1] ){
      stop("L or R is of incorrect length. Make sure it is equal to the sample size.")
    }

    if(length(grids)!= 1 |
       length(U)!= 1 |
       length(link)!= 1 |
       length(boot)!= 1 |
       length(nboot)!= 1){
      stop("Any of the variables grids, U, link, boot, nboot is of incorrect length. Make sure they are of length 1.")
    }
  }

  if (U < 0.1){
    warning("The specified upper limit might be too low.")
  }

  if (U/grids>0.1){
    warning("The number of grids might be too small to produce stable results.")
  }

  if (boot ==TRUE & nboot < 30){
    warning("Bootstrap times might be too small to produce reliable standard error estimates.")
  }



  X <- x[ID, ]
  n <- dim(x)[1]

  sun.res <- estsun_(T=T, x=x, ID=ID, grids=grids, U=U, L=L, R=R, link=link)
  if (sun.res[[2]]!=U){
    warning(paste0("The identifiability limit might not reach the provided upper bound. Suggested limit: ", sun.res[[2]]))
  }
  if (boot==FALSE){
    sun.ret <- list(beta.est=sun.res[[1]], idenLim=sun.res[[2]])
  } else {
    sun.boot.res <- array(NA,dim=c(dim(sun.res[[1]]),nboot),dimnames=list(dimnames(sun.res[[1]]),1:nboot))
    for (iboot in 1:nboot){
      zeta <- rexp(n)
      ZETA <- zeta[ID]
      sun.boot.res[,,iboot] <- estsun_(T=T, x=x, ID=ID, grids=grids, U=U, L=L, R=R, ZETA=ZETA, zeta=zeta, link=link)[[1]]
    }
    sun.res.sd <- apply(sun.boot.res,c(1,2),sd,na.rm=TRUE)
    sun.ret <- list(beta.est=cbind(sun.res[[1]],sun.res.sd), idenLim=sun.res[[2]])
  }

  return(sun.ret)
}

#' Adaptation of cqrIS Method to a Recurrent Event Model
#'
#' This function fits a recurrent event model to the provided dataset where the covariates
#' are related to the "expected frequency" via a generalized linear model. The parameters
#' are obtained by solving an induced-smoothed version of the estimating equation used in
#' Sun et al. (2016).
#'
#' @param T A numeric vector. Recurrent event time
#' @param x A numeric matrix. Covariates for samples
#' @param ID A numeric vector. ID of those who experience recurrent events
#' @param grids A numeric number. Number of grids
#' @param U A numeric number. Upper bound of estimation
#' @param L A numeric vector. Left censoring variable
#' @param R A numeric vector. Right censoring variable
#' @param link A character string. Link function, default: "1"
#' @param tol A numeric number. Norm of the estimate for the update to be considered convergent, default: 1e-5
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
#' evt <- recur[[1]]
#' sam <- recur[[2]]
#' res.recurIS <- recurIS(T=evt[,1], x=sam[,-c(1:2)], ID=evt[,2], grids=100, U=4, L=sam[,1], R=sam[,2])
#' }
#'
#' @return This function returns a list of lists with each list containing two elements:
#' \itemize{
#'   \item beta.recurISest, a matrix giving the estimates for the parameters. If boot=FALSE, it will contain p columns
#'   giving the coefficient estimates where p is the dimension of the covariate; if boot=TRUE, it will contain
#'   2*p columns giving the coefficient estimates in the first p columns and the estimated standard errors in
#'   the remaining columns.
#'   \item idenLim, a numeric giving the identifiability limit for the data
#' }
#'
#' @references Cai, Z. and Sit, T. (2020+),
#' "Censored Quantile Regression with Induced Smoothing,"
#' \emph{Working Paper}.

recurIS <- function(T, x, ID, grids, U, L, R, link="1", boot=FALSE, nboot=250, tol=1e-5, diverge=1e4, maxiter=200){
  if (!is.numeric(T) |
      !is.numeric(x) |
      !is.numeric(ID) |
      !is.numeric(grids) |
      !is.numeric(U) |
      !is.numeric(L) |
      !is.numeric(R) |
      !is.character(link) |
      !is.logical(boot) |
      !is.numeric(nboot) |
      !is.numeric(tol) |
      !is.numeric(diverge) |
      !is.numeric(maxiter)){
    stop("Incorrect variable type(s). Try \"?recurIS\" for more information.")
  }

  if (!is.matrix(x)){
    stop("x is of incorrect type. Make sure it is a matrix. Try \"?recurIS\" for more information.")
  }

  if (!is.vector(T) |
      !is.vector(ID) |
      !is.vector(grids) |
      !is.vector(U) |
      !is.vector(L) |
      !is.vector(R) |
      !is.vector(link) |
      !is.vector(boot) |
      !is.vector(nboot) |
      !is.vector(tol) |
      !is.vector(diverge) |
      !is.vector(maxiter)){
    stop("Any of the variables T, ID, grids, U, L, R, link, boot, nboot, tol, diverge, maxiter is of incorrect type. Make sure they are vectors. Try \"?recurIS\" for more information.")
  } else {
    if(length(T) != length(ID)){
      stop("T or ID is of incorrect length. Make sure their lengths are equal.")
    }

    if (length(L) != dim(x)[1] |
        length(R) != dim(x)[1] ){
      stop("L or R is of incorrect length. Make sure it is equal to the sample size.")
    }

    if(length(grids)!= 1 |
       length(U)!= 1 |
       length(link)!= 1 |
       length(boot)!= 1 |
       length(nboot)!= 1 |
       length(tol)!= 1 |
       length(diverge)!= 1 |
       length(maxiter)!= 1){
      stop("Any of the variables grids, U, link, boot, nboot, tol, diverge, maxiter is of incorrect length. Make sure they are of length 1.")
    }
  }

  if (U < 0.1){
    warning("The specified upper limit might be too low.")
  }

  if (U/grids>0.1){
    warning("The number of grids might be too small to produce stable results.")
  }

  if (boot ==TRUE & nboot < 30){
    warning("Bootstrap times might be too small to produce reliable standard error estimates.")
  }



  X <- x[ID, ]
  n <- dim(x)[1]

  recurIS.res <- recurIS_(T, x, ID, grids, U, L, R, link=link, tol=tol, diverge=diverge, maxiter=maxiter)
  if (recurIS.res[[2]]!=U){
    warning(paste0("The identifiability limit might not reach the provided upper bound. Suggested limit: ", recurIS.res[[2]]))
  }
  if (boot==FALSE){
    recurIS.ret <- list(beta.recurISest=recurIS.res[[1]], idenLim=recurIS.res[[2]])
  } else {
    recurIS.boot.res <- array(NA,dim=c(dim(recurIS.res[[1]]),nboot),dimnames=list(dimnames(recurIS.res[[1]]),1:nboot))
    for (iboot in 1:nboot){
      recurIS.boot.res[,,iboot] <- recurIS_(T, x, ID, grids, U, L, R, zeta=rexp(n), link=link, tol=tol, diverge=diverge, maxiter=maxiter)[[1]]
    }
    recurIS.res.sd <- apply(recurIS.boot.res,c(1,2),sd,na.rm=TRUE)
    recurIS.ret <- list(beta.recurISest=cbind(recurIS.res[[1]],recurIS.res.sd), idenLim=recurIS.res[[2]])
  }

  return(recurIS.ret)
}
