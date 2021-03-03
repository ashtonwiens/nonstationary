
#' Logit function
#'
#' @param x a probability between 0 and 1
#'
#' @return the log odds
#' @export
#'
#' @examples
logit <- function(x) {
  log(x / (1 - x))
}


#' Inverse logit / expit function
#'
#' @param x any number (log odds)
#'
#' @return the probability between 0 and 1
#' @export
#'
#' @examples
ilogit <- expit <- function(x) {
  exp(x) / (1 + exp(x))
}


# library(fields)
# x <- seq(0,1,,35)
# grd <- expand.grid(x, x)
# N <- 30
# y <- matrix(NA, nr=dim(grd)[1], nc=N)
# angl <- pi/6
# U <- matrix(c(cos(angl), -sin(angl),
#               sin(angl), cos(angl)), nrow=2, ncol=2, byrow=TRUE)
# D <- matrix(c(8, 0, 0, 1), nrow=2, ncol=2, byrow=TRUE)
# S <- U %*% D %*% t(U)
# Linv <- t(chol(solve(S)))
# grd2 <- t(Linv %*% t(grd) )
# #grd2 <- t( D %*% U %*% t(grd) )
# Sig <- 1*Matern(rdist(grd2), range=1, smoothness=1)
# Sc <- chol(Sig)
# for(i in 1:N){
#   y[,i] <- t(Sc) %*% rnorm(dim(grd)[1]) + rnorm(dim(grd)[1], sd=0.01)
# }
#
# library(ellipse)
# par(mfrow=c(2,2))
# E2 <- ellipse( S, level=0.05)  # / kappa[indx[1], indx[2]]^2
# plot(E2)
# image.plot(x, x, matrix(y[,1], nrow=length(x)))
# image.plot(x, x, matrix(y[,2], nrow=length(x)))
#

#' Wrapper function for fitting an anisotropic Matern GP
#'
#' @param grd a two columns matrix with the spatial locations of the data
#' @param y the value of the response variable at the locations in grd
#' @param NU a number > 0, the smoothness parameter
#'
#' @return an optim object with the parameter estimates of the anisotropic Matern
#' @export
#'
#' @examples
fit.Matern.aniso <- function(grd, y, NU = 1) {
  theta <- 1
  sigma2 <- 0.1
  tau2 <- 0.01
  lambda_y <- 1
  angle <- 0 # logit(0)
  stats::optim(
    par = c(logtheta = log(theta), logsigma2 = log(sigma2), logtau2 = log(tau2), loglambda2 = log(lambda_y), logitangle = angle),
    fn = Matern.logLik.aniso, nu = NU, grd = grd, y = y,
    hessian = TRUE,
    # control=list(trace=1),
    method = "L-BFGS-B", lower = c(log(0.1), -Inf, -Inf, log(0.1), -pi / 4 + 0.1), upper = c(log(150), Inf, Inf, log(150), pi / 4)
  )
}


#' Log likelihood function for anisotropic Matern GP
#'
#' @param pars covariance parameters: range in x, Matern variance, white noise
#' variance, range in y, and rotation angle
#' @param grd a two columns matrix with the spatial locations of the data
#' @param y the value of the response variable at the locations in grd
#' @param nu a number > 0, the smoothness parameter
#'
#' @return the negative log likelihood
#' @export
#'
#' @examples
Matern.logLik.aniso <- function(pars = c(log(1.5), log(0.1), log(0.01), log(0.5), 0),
                                grd, y, nu = 1.0) {
  nobs <- dim(grd)[1]
  nreps <- dim(y)[2]
  # print("theta sigma2 tau2")
  # print( round( pars, 3) )
  angl <- pars[5] # (pi/2)*ilogit(pars[5]) - (pi/4)
  U <- matrix(c(
    cos(angl), -sin(angl),
    sin(angl), cos(angl)
  ), nrow = 2, ncol = 2, byrow = TRUE)
  D <- matrix(c(exp(pars[1]), 0, 0, exp(pars[4])), nrow = 2, ncol = 2, byrow = TRUE)
  S <- U %*% D %*% t(U)
  Linv <- t(chol(solve(S)))
  grd2 <- t(Linv %*% t(grd))
  d <- fields::rdist(grd2)
  Sigma <- exp(pars[2]) * fields::Matern(d, range = 1, smoothness = nu) + exp(pars[3]) * diag(nobs)
  Sigma.c <- t(chol(Sigma))
  out <- forwardsolve(Sigma.c, y)
  quad.form <- sum(out^2)
  det.part <- nreps * 2 * sum(log(diag(Sigma.c)))
  nLL <- 0.5 * det.part + 0.5 * quad.form
  # print('log det :: quadratic form')
  # print(c(logDetCov, quad.form) )
  # cat('\n')
  print(nLL)
  return(nLL)
}

# myLik <- fit.Matern.aniso(grd, y, NU=1)
#
# print(exp(myLik$par[c(1:4)]))
# print(myLik$par[5])
# #print((pi/2)*ilogit(myLik$par[5]) - (pi/4) )
#
# #angl <- (pi/2)*ilogit(myLik$par[5]) - (pi/4)
# angl <- myLik$par[5]
# U <- matrix(c(cos(angl), -sin(angl),
#               sin(angl), cos(angl)), nrow=2, ncol=2, byrow=TRUE)
# D <- matrix(c(exp(myLik$par[1]), 0, 0, exp(myLik$par[4])), nrow=2, ncol=2, byrow=TRUE)
# S <- U %*% D %*% t(U)
#
# E2 <- ellipse( S, level=0.05)  # / kappa[indx[1], indx[2]]^2
# plot(E2)


#' Log likelihood function for isotropic Matern GP
#'
#' @param pars covariance parameters: range parameter, Matern variance, and
#' white noise variance
#' @param grd a two columns matrix with the spatial locations of the data
#' @param y the value of the response variable at the locations in grd
#' @param nu a number > 0, the smoothness parameter
#'
#' @return the negative log likelihood
#' @export
#'
#' @examples
Matern.logLik <- function(pars = c(log(1.5), log(3), log(0.3)),
                          grd, y, nu = 1.0) {
  nobs <- dim(grd)[1]
  nreps <- dim(y)[2]
  # print("theta sigma2 tau2")
  # print( round( pars, 3) )
  d <- fields::rdist(grd)
  Sigma <- exp(pars[2]) * fields::Matern(d, range = exp(pars[1]), smoothness = nu) +
    exp(pars[3]) * diag(nobs)
  Q <- solve(Sigma)
  logDetCov <- nreps * sum(log(eigen(Sigma, only.values = TRUE)$values))
  quad.form <- sum(apply(y, 2, function(z) {
    t(z) %*% Q %*% z
  }))
  nLL <- logDetCov + quad.form
  # print('log det :: quadratic form')
  # print(c(logDetCov, quad.form) )
  # cat('\n')
  return(nLL)
}

#' Wrapper function for fitting an isotropic Matern GP
#'
#' @param grd a two columns matrix with the spatial locations of the data
#' @param y the value of the response variable at the locations in grd
#' @param NU a number > 0, the smoothness parameter
#'
#' @return an optim object with the parameter estimates of the isotropic Matern
#' @export
#'
#' @examples
fit.Matern <- function(grd, y, NU = 1) {
  theta <- 1
  sigma2 <- 0.1
  tau2 <- 0.01
  stats::optim(
    par = c(logtheta = log(theta), logsigma2 = log(sigma2), logtau2 = log(tau2)),
    fn = Matern.logLik, nu = NU, grd = grd, y = y,
    hessian = TRUE,
    # control=list(trace=1),
    method = "L-BFGS-B", lower = c(log(0.1), -Inf, -Inf), upper = c(log(150), Inf, Inf)
  )
}







#' Residual log likelihood function for isotropic Matern GP
#'
#' @param pars covariance parameters: range parameter, Matern variance, and
#' white noise variance
#' @param grd a two columns matrix with the spatial locations of the data
#' @param y the value of the response variable at the locations in grd
#' @param Z covariates as columns in a matrix
#' @param nu a number > 0, the smoothness parameter
#'
#' @return the negative log likelihood
#' @export
#'
#' @examples
Matern.reml.loglik <- function(pars = c(log(1.5), log(3), log(0.3)),
                               grd, y, Z, nu = 1.0) {
  nobs <- dim(grd)[1]
  nreps <- ifelse(is.null(dim(y)[2]), 1, dim(y)[2])
  # print("theta sigma2 tau2")
  print(round(exp(pars), 3))
  d <- fields::rdist(grd)
  Sigma <- exp(pars[2]) * fields::Matern(d, range = exp(pars[1]), smoothness = nu) +
    exp(pars[3]) * diag(nobs)
  Q <- solve(Sigma)
  logDetCov <- nreps * sum(log(eigen(Sigma, only.values = TRUE)$values))
  logDetCov2 <- log(det(t(Z) %*% Q %*% Z))
  # quad.form <- sum(apply(y, 2, function(z){t(z) %*% Q %*% z}))
  quad.form <- t(y) %*% (Q - (Q %*% Z %*% solve(t(Z) %*% Q %*% Z) %*% t(Z) %*% Q)) %*% y
  nLL <- logDetCov + logDetCov2 + quad.form
  # print('log det :: quadratic form')
  # print(c(logDetCov, quad.form) )
  # cat('\n')
  return(nLL)
}

#' Wrapper function for fitting an isotropic Matern GP with REML
#'
#' @param grd a two columns matrix with the spatial locations of the data
#' @param y the value of the response variable at the locations in grd
#' @param Z covariates as columns in a matrix
#' @param NU a number > 0, the smoothness parameter
#'
#' @return an optim object with the parameter estimates of the isotropic Matern
#' @export
#'
#' @examples
fit.Matern.reml <- function(grd, y, Z, NU = 1) {
  theta <- 1
  sigma2 <- 0.1
  tau2 <- 0.01
  stats::optim(
    par = c(logtheta = log(theta), logsigma2 = log(sigma2), logtau2 = log(tau2)),
    fn = Matern.reml.loglik, nu = NU, grd = grd, y = y, Z = Z,
    hessian = TRUE,
    # control=list(trace=1),
    method = "L-BFGS-B", lower = c(log(0.1), -Inf, -Inf), upper = c(log(150), Inf, Inf)
  )
}





#' Residual log likelihood function for anisotropic Matern GP
#'
#' @param pars covariance parameters: range in x, Matern variance, white noise
#' variance, range in y, and rotation angle
#' @param grd a two columns matrix with the spatial locations of the data
#' @param y the value of the response variable at the locations in grd
#' @param Z covariates as columns in a matrix
#' @param nu a number > 0, the smoothness parameter
#'
#' @return the negative log likelihood
#' @export
#'
#' @examples
Matern.reml.aniso.loglik <- function(pars = c(log(1.5), log(3), log(0.3), log(1.5), 0),
                                     grd, y, Z, nu = 1.0) {
  nobs <- dim(grd)[1]
  nreps <- ifelse(is.null(dim(y)[2]), 1, dim(y)[2])
  # print("theta sigma2 tau2")
  print(round(exp(pars), 3))
  angl <- pars[5] # (pi/2)*ilogit(pars[5]) - (pi/4)
  U <- matrix(c(
    cos(angl), -sin(angl),
    sin(angl), cos(angl)
  ), nrow = 2, ncol = 2, byrow = TRUE)
  D <- matrix(c(exp(pars[1]), 0, 0, exp(pars[4])), nrow = 2, ncol = 2, byrow = TRUE)
  S <- U %*% D %*% t(U)
  Linv <- t(chol(solve(S)))
  grd2 <- t(Linv %*% t(grd))
  d <- fields::rdist(grd2)
  Sigma <- exp(pars[2]) * fields::Matern(d, range = 1, smoothness = nu) + exp(pars[3]) * diag(nobs)
  # Sigma.c <- t(chol(Sigma))
  # out <- forwardsolve(Sigma.c, y)
  # quad.form <- sum(out^2)
  # det.part <- nreps * 2*sum(log(diag(Sigma.c)))
  # nLL <- 0.5*det.part + 0.5*quad.form
  Q <- solve(Sigma)
  logDetCov <- nreps * sum(log(eigen(Sigma, only.values = TRUE)$values))
  logDetCov2 <- log(det(t(Z) %*% Q %*% Z))
  # quad.form <- sum(apply(y, 2, function(z){t(z) %*% Q %*% z}))
  quad.form <- t(y) %*% (Q - (Q %*% Z %*% solve(t(Z) %*% Q %*% Z) %*% t(Z) %*% Q)) %*% y
  nLL <- logDetCov + logDetCov2 + quad.form
  # print('log det :: quadratic form')
  # print(c(logDetCov, quad.form) )
  # cat('\n')
  return(nLL)
}

#' Wrapper function for fitting an anisotropic Matern GP with REML
#'
#' @param grd a two columns matrix with the spatial locations of the data
#' @param y the value of the response variable at the locations in grd
#' @param Z covariates as columns in a matrix
#' @param NU a number > 0, the smoothness parameter
#'
#' @return an optim object with the parameter estimates of the anistropic Matern
#' @export
#'
#' @examples
fit.Matern.aniso.reml <- function(grd, y, Z, NU = 1) {
  theta <- 1
  sigma2 <- 0.1
  tau2 <- 0.01
  lambda_y <- 1
  angle <- 0
  stats::optim(
    par = c(logtheta = log(theta), logsigma2 = log(sigma2), logtau2 = log(tau2), loglambda2 = log(lambda_y), logitangle = angle),
    fn = Matern.reml.loglik, nu = NU, grd = grd, y = y, Z = Z,
    hessian = TRUE,
    # control=list(trace=1),
    method = "L-BFGS-B",
    lower = c(log(0.1), -Inf, -Inf, log(0.1), -pi / 4 + 0.1),
    upper = c(log(150), Inf, Inf, log(150), pi / 4)
  )
}
