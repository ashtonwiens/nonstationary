
#' Define a nonstationary covariance matrix based on a generalization of the Matern with geometric anisotropy
#'
#' Encode (local) spatially varying parameters (ranges Lx and Ly, smoothness nu,
#' rotation angle, variance sigma) into a global nonstationarity covariance with
#' local geometric anisotropy. See Paciorek (2003) and Stein (2005).
#'
#' @param x x-coordinate of locations on a regular rectangular grid
#' @param y y-coordinate of locations on a regular rectangular grid
#' @param Lx range in x direction (geometric anistotropy)
#' @param Ly range in x direction (geometric anistotropy)
#' @param ang angle of rotation (geometric anistotropy)
#' @param nu smoothness parameter
#' @param sigma standard deviation parameter
#' @param verbose whether to print progress or not
#'
#' @return a covariance matrix
#' @export
#'
#' @examples
Anisotropic.Matern.Nonstationary.Cov <- function(x, y, Lx, Ly, ang, nu, sigma=1, verbose=FALSE){
  m <- length(x)
  n <- length(y)
  stopifnot(dim(Lx)==c(m,n) & dim(nu)==c(m,n))
  par.grd <- expand.grid(x, y)
  par.grd$Lx <- c(Lx)
  par.grd$Ly <- c(Ly)
  par.grd$ang <- c(ang)
  par.grd$nu <- c(nu)
  par.grd$sigma <- c(sigma)
  cov.NS <- matrix(NA, nrow=m*n, ncol=m*n)

  #chunk.ind <- seq(1, 12601, 200)
  #chunk.ind2 <- chunk.ind+199
  #chunks <- map2(chunk.ind, chunk.ind2, ~ seq(.x, .y))
  #chunks[[44]] <- 12901:13056

  for(count in 1){ #rev( seq_along(chunks)) ){
      for(i in 1:(m*n) ){#chunks[[count]]){
        for(j in 1:(m*n)){
          if(i >= j){
            nu.avg <- (par.grd$nu[i] + par.grd$nu[j])/2 # U+03BD
            const <- 1/(gamma(nu.avg)*2^(nu.avg-1))
            #const <- 1 # for nu==1
            U <- matrix(c(cos(par.grd$ang[i]), -sin(par.grd$ang[i]),
                          sin(par.grd$ang[i]), cos(par.grd$ang[i])), nrow=2, ncol=2, byrow=TRUE)
            D <- matrix(c(par.grd$Lx[i]^2, 0, 0, par.grd$Ly[i]^2), nrow=2, ncol=2, byrow=TRUE)
            sig.i <- U %*% D %*% t(U)
            U <- matrix(c(cos(par.grd$ang[j]), -sin(par.grd$ang[j]),
                          sin(par.grd$ang[j]), cos(par.grd$ang[j])), nrow=2, ncol=2, byrow=TRUE)
            D <- matrix(c(par.grd$Lx[j]^2, 0, 0, par.grd$Ly[j]^2), nrow=2, ncol=2, byrow=TRUE)
            sig.j <- U %*% D %*% t(U)
            #sig.i <- matrix(c(par.grd$Lx[i]^2,0,0,par.grd$Lx[i]^2), nrow=2, ncol=2)
            #sig.j <- matrix(c(par.grd$Lx[j]^2,0,0,par.grd$Lx[j]^2), nrow=2, ncol=2)
            d.ij <- as.matrix( abs(par.grd[i,1:2]-par.grd[j,1:2] ) )
            sig.av <- (1/2)*(sig.i+sig.j)
            sp <- eigen(sig.av)
            U <- -sp$vectors
            D <- diag(sp$values)
            sig.av <- solve(sig.av)
            L <- sqrt(solve(D))%*%t(U)
            #sig.av <- matrix( c(2/(par.grd$Lx[i] + par.grd$Lx[j]),0,0,2/(par.grd$Lx[i] + par.grd$Lx[j])),
            #                  ncol=2) # for Sigma diagonal matrices
            Q <- sum( (t(L%*%t(par.grd[i,1:2]-par.grd[j,1:2])) )^2) #d.ij %*% (sig.av) %*% t(d.ij)
            if(dplyr::near(Q, 0, tol=1e-6)){
              cov.NS[i,j] <- par.grd$sigma[i] * par.grd$sigma[j]
            }else{
              sig.dets <- det(sig.i)^(1/4) * det(sig.j)^(1/4) * det(sig.av)^(1/2)
              #sig.dets <- sqrt(par.grd$Lx[i]) * sqrt(par.grd$Lx[j]) * 2/(par.grd$Lx[i] + par.grd$Lx[j])
              cov.NS[i,j] <-  par.grd$sigma[i] * par.grd$sigma[j] * const * sig.dets *
                (sqrt(Q))^nu.avg * besselK(sqrt(Q),nu.avg)
              #             (2*sqrt(nu.avg*Q))^nu.avg * besselK(2*sqrt(nu.avg*Q),nu.avg)
            }
          }

        }
        if(verbose)print(i/(m*n))
      }
  }
  for(i in 1:(m*n)){
    for(j in 1:(m*n)){
      if( i < j){
        cov.NS[i,j] <- cov.NS[j,i] #C[i,j] <- C[j,i]
      }
    }
  }
  return(cov.NS)
}



#' Define a nonstationary covariance matrix based on a generalization of the
#' Matern with geometric anisotropy
#'
#' Encode (local) spatially varying parameters (range theta, smoothness nu,
#' variance sigma) into a global nonstationarity covariance with local isotropy.
#' See Paciorek (2003) and Stein (2005).
#'
#' @param x x-coordinate of locations on a regular rectangular grid
#' @param y y-coordinate of locations on a regular rectangular grid
#' @param theta correlation range parameter
#' @param nu smoothness parameter
#' @param sigma standard deviation parameter
#' @param family Matern, exp, or RQ (rational quadratic)
#' @param verbose whether to print progress or not
#'
#' @return a covariance matrix
#' @export
#'
#' @examples
Isotropic.Matern.Exp.RQ.Nonstationary.Covs <- function(x, y, theta, nu, sigma=1,
                                                       family='matern',
                                                       verbose=FALSE){
  m <- length(x)
  n <- length(y)
  stopifnot(dim(theta)==c(m,n) & dim(nu)==c(m,n))
  par.grd <- expand.grid(x, y)
  par.grd$theta <- c(theta)
  par.grd$nu <- c(nu)
  par.grd$sigma <- c(sigma)
  if(plot==TRUE){
    graphics::par(mfrow=c(1,2))
    # plots input parameters
    fields::quilt.plot(par.grd[,1:2], par.grd$theta, nx=m, ny=n, main="theta")
    fields::quilt.plot(par.grd[,1:2], par.grd$nu, nx=m, ny=n, main="nu")
  }
  cov.NS <- matrix(NA, nrow=m*n, ncol=m*n)

  #chunk.ind <- seq(1, 12601, 200)
  #chunk.ind2 <- chunk.ind+199
  #chunks <- map2(chunk.ind, chunk.ind2, ~ seq(.x, .y))
  #chunks[[44]] <- 12901:13056

  for(count in 1){ #rev( seq_along(chunks)) ){
    if(family=='exp'){
      for(i in 1:(m*n)){
        for(j in 1:(m*n)){
          if(i >= j){
            nu.avg <- (par.grd$nu[i] + par.grd$nu[j])/2 # U+03BD
            const <- 1/(gamma(nu.avg)*2^(nu.avg-1))
            sig.i <- matrix(c(par.grd$theta[i],0,0,par.grd$theta[i]), nrow=2, ncol=2)
            sig.j <- matrix(c(par.grd$theta[j],0,0,par.grd$theta[j]), nrow=2, ncol=2)
            print(par.grd[i,1:2])
            print(par.grd[j,1:2])
            d.ij <- as.matrix( abs(par.grd[i,1:2]-par.grd[j,1:2] ))
            sig.av <- solve((1/2)*(sig.i+sig.j))
            Q <- d.ij %*% sig.av %*% t(d.ij)
            if(Q == 0){
              Q <- 1e-10
            }
            sig.dets <- det(sig.i)^(1/4) * det(sig.j)^(1/4) * det(sig.av)^(1/2)
            cov.NS[i,j] <-  sigma[i] * sigma[j] * const * sig.dets * exp(-Q)
          }
        }
        if(verbose)print(i/(m*n))
      }
    }
    if(family=='matern'){
      print(family)
      for(i in 1:(m*n) ){#chunks[[count]]){
        for(j in 1:(m*n)){
          if(i >= j){
            nu.avg <- (par.grd$nu[i] + par.grd$nu[j])/2 # U+03BD
            const <- 1/(gamma(nu.avg)*2^(nu.avg-1))
            #const <- 1 # for nu==1
            sig.i <- matrix(c(par.grd$theta[i]^2,0,0,par.grd$theta[i]^2), nrow=2, ncol=2)
            sig.j <- matrix(c(par.grd$theta[j]^2,0,0,par.grd$theta[j]^2), nrow=2, ncol=2)
            d.ij <- as.matrix( abs(par.grd[i,1:2]-par.grd[j,1:2] ) )
            sig.av <- solve((1/2)*(sig.i+sig.j))
            #sig.av <- matrix( c(2/(par.grd$theta[i] + par.grd$theta[j]),0,0,2/(par.grd$theta[i] + par.grd$theta[j])),
            #                  ncol=2) # for Sigma diagonal matrices
            Q <- d.ij %*% (sig.av) %*% t(d.ij)
            if(dplyr::near(Q, 0, tol=1e-6)){
              cov.NS[i,j] <- par.grd$sigma[i] * par.grd$sigma[j]
            }else{
              sig.dets <- det(sig.i)^(1/4) * det(sig.j)^(1/4) * det(sig.av)^(1/2)
              #sig.dets <- sqrt(par.grd$theta[i]) * sqrt(par.grd$theta[j]) * 2/(par.grd$theta[i] + par.grd$theta[j])
              cov.NS[i,j] <-  par.grd$sigma[i] * par.grd$sigma[j] * const * sig.dets *
                (sqrt(Q))^nu.avg * besselK(sqrt(Q),nu.avg)
              #             (2*sqrt(nu.avg*Q))^nu.avg * besselK(2*sqrt(nu.avg*Q),nu.avg)
            }
          }

        }
        if(verbose)print(i/(m*n))
      }
    }
    if(family=='rq'){
      for(i in 1:(m*n)){
        for(j in 1:(m*n)){
          if(i >= j){
            nu.avg <- (par.grd$nu[i] + par.grd$nu[j])/2 # U+03BD
            sig.i <- matrix(c(par.grd$theta[i],0,0,par.grd$theta[i]), nrow=2, ncol=2)
            sig.j <- matrix(c(par.grd$theta[j],0,0,par.grd$theta[j]), nrow=2, ncol=2)
            d.ij <- as.matrix( abs(par.grd[i,1:2]-par.grd[j,1:2] ) )
            sig.av <- solve((1/2)*(sig.i+sig.j))
            Q <- d.ij %*% sig.av %*% t(d.ij)
            if(Q == 0){
              Q <- 1e-10
            }
            sig.dets <- det(sig.i)^(1/4) * det(sig.j)^(1/4) * det(sig.av)^(1/2)
            cov.NS[i,j] <- sigma[i] * sigma[j] * sig.dets * (Q+1)^(-nu.avg)
          }

        }
        if(verbose)print(i/(m*n))
      }
    }

  }
  for(i in 1:(m*n)){
    for(j in 1:(m*n)){
      if( i < j){
        cov.NS[i,j] <- cov.NS[j,i] #C[i,j] <- C[j,i]
      }
    }
  }
  # if(plot==TRUE){
  #   if(is.null(plotsize)){
  #     plotsize <- dim(Q)[1]
  #   }
  #   ## compare "distance" matrices
  #   dd <- fields::rdist(par.grd)
  #   fields::image.plot(rotate(dd[1:plotsize,1:plotsize]), main='euclidean distance matrix of spatial locations')
  #   #image.plot(rotate(Q[1:plotsize,1:plotsize]), main='"Q" distance')
  #   fields::image.plot( rotate(cov.NS[1:plotsize,1:plotsize]), main='Nonstationary Cov Matrix' )
  # }
  # if(!is.null(savefile)){
  #   save( x, y, par.grd, C, file=savefile )
  # }
  return(cov.NS)
}


