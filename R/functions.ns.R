
#' Encode nonstationary GMRF precision matrix
#'
#' @param x x-coordinates of grid
#' @param y y-coordinates of grid
#' @param grd two column location of grid
#' @param symm if TRUE, defines H using the symmetric square root, otherwise uses the Cholesky
#'
#' @return a precision matrix
#' @export
#'
#' @examples
fivepointB.J <- function(x, y, grd, symm = FALSE) {
  ### this function was used to encode the SAR precision in the nonstationary paper
  # grd columns 1 and 2 are expand.grid(x,y)
  # grd columns 3,4,5,6 are kappa, J1, J2, J3
  m <- length(x)
  n <- length(y)
  G <- expand.grid(1:m, 1:n)
  N <- dim(grd)[1]
  out <- matrix(0, N, N)
  # specification of kappa, J1, J2, J3 as columns 3:6 of grd
  # print('Using 8 columns specification with kappa, J1, J2, and J3')
  for (i in 1:N) {
    O <- matrix(0, m, n)
    I <- matrix(unlist(G[i, ]), ncol = 2)
    hx <- diff(x)[1]
    hy <- diff(y)[1]
    kappa <- grd[i, 3] # 0 to Inf
    J1 <- grd[i, 4] # -Inf to Inf
    J2 <- grd[i, 5] # -Inf to Inf
    J3 <- grd[i, 6] # -Inf to Inf
    if (symm == TRUE) {
      J <- matrix(c(J1, J2, J2, J3), nrow = 2, byrow = TRUE)
    } else {
      J <- matrix(c(J1, J2, 0, J3), nrow = 2, byrow = TRUE)
    }
    H <- J %*% t(J)
    H1 <- H[1, 1]
    H2 <- H[1, 2]
    H3 <- H[2, 2]
    Ins <- rbind(
      I + cbind(0, 1),
      I - cbind(0, 1)
    )
    Iew <- rbind(
      I + cbind(1, 0),
      I - cbind(1, 0)
    )
    I_ne_sw <- rbind(
      I + cbind(1, 1),
      I - cbind(1, 1)
    )
    I_se_nw <- rbind(
      I + cbind(1, -1),
      I - cbind(1, -1)
    )
    Ins <- matrix(Ins[Ins[, 2] <= n & Ins[, 2] > 0, ], ncol = 2, byrow = FALSE)
    Iew <- matrix(Iew[Iew[, 1] <= m & Iew[, 1] > 0, ], ncol = 2, byrow = FALSE)
    I_ne_sw <- matrix(I_ne_sw[I_ne_sw[, 1] <= m &
      I_ne_sw[, 1] > 0 &
      I_ne_sw[, 2] <= n & I_ne_sw[, 2] > 0, ],
    ncol = 2, byrow = FALSE
    )
    I_se_nw <- matrix(I_se_nw[I_se_nw[, 1] <= m &
      I_se_nw[, 1] > 0 &
      I_se_nw[, 2] <= n & I_se_nw[, 2] > 0, ],
    ncol = 2, byrow = FALSE
    )
    O[Iew] <- -H1 / (hx^2)
    O[Ins] <- -H3 / (hy^2)
    O[I_ne_sw] <- -(2 * H2) / (4 * hx * hy)
    O[I_se_nw] <- (2 * H2) / (4 * hx * hy)
    middle.term <- kappa^2 + (nrow(Iew) * H1 / (hx^2)) + (nrow(Ins) * H3 / (hy^2))
    if (!((nrow(I_se_nw) + nrow(I_ne_sw)) %% 2 == 0)) {
      print("Corner adjustment")
      if (nrow(I_se_nw) == 0 & nrow(I_ne_sw) == 1) {
        middle.term <- middle.term + ((2 * H2) / (4 * hx * hy))
      } else {
        middle.term <- middle.term - ((2 * H2) / (4 * hx * hy))
      }
    }
    O[I] <- middle.term
    out[i, ] <- c(O)
  }
  return(out)
}

#' Log likelihood function using a GMRF precision matrix
#'
#' @param pars parameters defining the precision matrix; the range
#' @param x x-coordinates of grid
#' @param y y-coordinates of grid
#' @param Z response variable of data at locations
#' @param symm if TRUE, defines H using the symmetric square root, otherwise uses the Cholesky
#'
#' @return the negative log likelihood
#' @export
#'
#' @examples
fivepoint.logLik <- function(
                             pars = c(
                               log(0.1), 1, 0,
                               # kappa, J1, J2, J3,
                               log(0.1), log(1)
                             ),
                             # sigma2=1, tau2=0.1,
                             x, y, Z, symm) {
  # pars = c( kappa, J1, J2, J3, tau2)
  nobs <- dim(Z)[1]
  nreps <- dim(Z)[2]
  grd <- expand.grid(x, y)
  grd$kappa <- exp(pars[1])
  grd$J1 <- pars[2]
  grd$J2 <- pars[3]
  grd$J3 <- 1
  B <- spam::as.spam(fivepointB.J(x, y, grd, symm = symm))
  Q <- t(B) %*% B
  S <- spam::chol2inv(spam::chol.spam(Q))
  C <- exp(pars[5]) * S + exp(pars[4]) * diag(nobs)
  C.c <- t(chol(C))
  out <- forwardsolve(C.c, Z)
  quad.form <- sum(out^2)
  log.det <- nreps * 2 * sum(log(diag(C.c)))
  nLL <- 0.5 * log.det + 0.5 * quad.form
  cat(
    "kappa: ", exp(pars[1]), ",  J1: ", pars[2], ", J2: ", pars[3],
    #  ', J3: ', pars[4],
    ", tau^2: ", exp(pars[4]),
    ", sigma^2: ", exp(pars[5]),
    ", nLL: ", nLL, "\n"
  )
  return(nLL)
}

#' Fit nonstationary GMRF precision matrix from data
#'
#' @param xx x-coordinates of data
#' @param yy y-coordinates of data
#' @param ZZ response variable of data at locations
#' @param symm if TRUE, defines H using the symmetric square root, otherwise uses the Cholesky
#' @param init initial values in optimization
#'
#' @return results from optim
#' @export
#'
#' @examples
fit.fivepoint.J.s <- function(xx, yy, ZZ, symm, init = NULL) {
  ## xx, yy are vectors of coordinates, ZZ is a matrix with rows = obs and cols = replicates
  if (is.null(init)) {
    init <- c(log(0.15), 1, 0, log(0.01), log(1))
  }
  stats::optim(
    par = init,
    fn = fivepoint.logLik,
    x = xx, y = yy, Z = ZZ, symm = symm,
    method = "L-BFGS-B", hessian = TRUE
  ) # ,
  # lower=c(-10,-10,-10,-10, -10), upper=c(10,10,10,10,10) )
}
