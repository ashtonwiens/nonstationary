## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(nonstationary)

## ---- fig.width=8-------------------------------------------------------------
#load("../nsMaternEstimates.rda")
data("nsMaternEstimates")
grd <- nsMaternEstimates$grd
lonNS <- nsMaternEstimates$lonNS
latNS <- nsMaternEstimates$latNS
thx <- nsMaternEstimates$thx
thy <- nsMaternEstimates$thy
ang <- nsMaternEstimates$ang
s2 <- nsMaternEstimates$s2
t2 <- nsMaternEstimates$t2

# getting cholesky decomposition from eigen decomposition estimates
J1 <- J2 <- J3 <- lambda <- ecc <- c()

for(i in 1:nrow(grd)){
  U <- matrix(c(cos(ang[i]), -sin(ang[i]),
                sin(ang[i]), cos(ang[i])), nrow=2, ncol=2, byrow=TRUE)
  lambda[i] <- sqrt((thx[i] * thy[i]) )
  D <- matrix(c(thx[i]/lambda[i], 0, 0, thy[i]/lambda[i]),
              nrow=2, ncol=2, byrow=TRUE)
  #      D <- matrix(c(theta[i], 0, 0, Ly[i]), nrow=2, ncol=2, byrow=TRUE)
  H <- U %*% D^2 %*% t(U) 
  ecc[i] <- (thx[i] / thy[i]) 
  
  J <- chol(H)
  J1[i] <- J[1,1]
  J2[i] <- J[1,2]
  J3[i] <- J[2,2]
}


par(mfrow=c(1,4), mar=c(1,0.3,2,0.7))

#fields::image.plot(lonNS, latNS, theta, main='theta'); fields::world(add=TRUE) #, zlim=c(0,range.thresh)
#fields::image.plot(lonNS, latNS, Ly, main='Ly'); fields::world(add=TRUE) #, zlim=c(0,range.thresh)
#fields::image.plot(lonNS, latNS, ang, main='angle'); fields::world(add=TRUE) #, zlim=c(0,range.thresh)

fields::image.plot(lonNS, latNS, matrix(s2, nrow=length(lonNS), ncol=length(latNS)), main=expression(paste('Sill ', sigma^2 )),
           xaxt='n', yaxt='n'); fields::world(add=TRUE, col='grey') #, zlim=c(0,range.thresh)
fields::image.plot(lonNS, latNS, matrix(t2, nrow=length(lonNS), ncol=length(latNS)), main=expression(paste('Nugget ', tau^2 )),
           xaxt='n', yaxt='n'); fields::world(add=TRUE, col='grey') #, zlim=c(0,range.thresh)
fields::image.plot(lonNS, latNS, matrix(lambda, nrow=length(lonNS), ncol=length(latNS)), main=expression(paste('Geom Avg Range ',sqrt(paste(lambda[x],lambda[y]) ) )),
           xaxt='n', yaxt='n'); fields::world(add=TRUE, col='grey') #, zlim=c(0,range.thresh)

K <- matrix(s2, nrow=length(lonNS), ncol=length(latNS)); K[] <- NA
par(mar=c(2.2,0.5,4.3,0.5))
par(mar=c(1.8,0.7,3.4,0.7))

image(lonNS, latNS, K, main=expression(paste('Anisotropy ', H )),
      xaxt='n', yaxt='n',
      col=fields::two.colors(100, start = 'white', end='black'), zlim=c(0,1)); fields::world(add=TRUE, col='grey')
#fields::image.plot(lonNS, latNS, kappa, main='kappa and H', col=fields::tim.colors()); fields::world(add=TRUE)

mn <- 15
xx <- seq(2, 102,,mn) # for making grid of ellipses
yy <- seq(3, 123,,mn)
Fg <- expand.grid(xx,yy)
G <- expand.grid(lonNS[xx], latNS[yy])

for(i in 1:dim(Fg)[1]){
  indx <- unlist(Fg[i,])
  
  ang.i <- matrix(ang, nrow=length(lonNS), ncol=length(latNS))[indx[1], indx[2]]
  thx.i <- matrix(thx, nrow=length(lonNS), ncol=length(latNS))[indx[1], indx[2]]
  thy.i <- matrix(thy, nrow=length(lonNS), ncol=length(latNS))[indx[1], indx[2]]
  
  U <- matrix(c(cos(ang.i), -sin(ang.i),
                sin(ang.i), cos(ang.i)), nrow=2, ncol=2, byrow=TRUE)
  D <- matrix(c(thx.i, 0, 0, thy.i), nrow=2, ncol=2, byrow=TRUE)# * (1/lambda[indx[1], indx[2]])
  H <- U %*% D^2 %*% t(U) 
  
  J <- matrix(NA, 2,2)
  J[1,1] <- matrix(J1, nrow=length(lonNS), ncol=length(latNS))[indx[1], indx[2]]
  J[1,2] <- J[2,1] <- matrix(J2, nrow=length(lonNS), ncol=length(latNS))[indx[1], indx[2]]
  J[2,2] <- matrix(J3, nrow=length(lonNS), ncol=length(latNS))[indx[1], indx[2]]
  H2 <- t(J)%*%J #* lambda[indx[1], indx[2]]
  
  #E2 <- ellipse( H2 , level=0.95)  # / kappa[indx[1], indx[2]]^2
  E2 <- ellipse::ellipse( H , level=0.001)#0.015  # / kappa[indx[1], indx[2]]^2
  E2 <- E2 + matrix( rep(unlist(G[i,]), nrow(E2)), nc=2, byrow = TRUE)
  lines(E2, col='black')
} 
# par(mar=c(1,0.5,2,0.5))
# fields::image.plot(lonNS, latNS, lambda, main=expression(paste('Geom Avg Range ',sqrt(paste(lambda[x],lambda[y]) ) )),
#            xaxt='n', yaxt='n'); fields::world(add=TRUE) #, zlim=c(0,range.thresh)
# fields::image.plot(lonNS, latNS, log(ecc), main=expression(paste('Eccentricity ',log( lambda[x]/lambda[y]) )),
#            xaxt='n', yaxt='n'); fields::world(add=TRUE) #, zlim=c(0,range.thresh)
# fields::image.plot(lonNS, latNS, nu)

## ---- eval=FALSE--------------------------------------------------------------
#  
#  ############### ENCODE SAR
#  #######################
#  grd <- expand.grid(lonNS, latNS)
#  grd$kappa <- 1/sqrt(c(lambda)) #c(1/theta)
#  adjst <- function(x){
#    x^1.158 # empirical fit based on simulation study
#  }
#  grd$kappa <- 1/adjst ( sqrt(c(lambda)) ) #c(1/theta)
#  grd$J1 <- c(J1)
#  grd$J2 <- c(J2)
#  grd$J3 <- c(J3)
#  grd$nu <- 2
#  
#  system.time({
#    Bs <- spam::as.spam(fivepointB.J(lonNS, latNS, grd)) #.order2
#    #fields::image.plot(matrix(diag(Bs), nrow=length(lonNS), ncol=length(latNS)))
#    Qs1 <- spam::crossprod(Bs)
#    Qs <- spam::crossprod(Qs1) #t(Bs) %d*% Bs %d*% t(Bs) %d*% Bs
#    Sigma.s <- spam::chol2inv( spam::chol.spam(Qs))
#    s1 <- diag(Sigma.s)
#    var.adjust <- diag( c(sqrt(s1)/c(sqrt(s2))) )
#    vas <- spam::as.spam(var.adjust)
#    Q1n <- Qs1 %d*% vas
#    Qn <- crossprod(Q1n)
#    Sigma.G <- spam::chol2inv( spam::chol.spam(Qn))
#    Sigma <- Sigma.G + diag(c(t2))
#    Sigma.c <- chol(Sigma)
#  })

## ---- eval=FALSE, fig.width=8-------------------------------------------------
#  S <- Sigma.c
#  set.seed(94)
#  YM <- t(S) %*% rnorm(dim(S)[1])
#  YM <- matrix(YM, nr=length(lonNS))
#  f <- 1 #sqrt(sdRaw)
#  image(lonNS, latNS, YM*f, zlim=c(-3.6, 3.6),xlab='', ylab='', xaxt='n', yaxt='n', col=fields::tim.colors()); fields::world(add=TRUE, col='gray')

## ---- eval=FALSE, fig.width=8-------------------------------------------------
#  a <- c(3010, 3025, 3030, 3040)#, 7928)
#  
#  par(mfrow=c(1,4), mar=c(1,1,1,1))
#  S.n <- cov2cor(Sigma)
#  rw <- a[1]
#  cond.ind <- matrix( S.n[rw,], nr=length(lonNS), nc=length(latNS)) #Re(gB) #emp.cor
#  pt <- grd[rw, ]
#  image(lonNS, latNS, cond.ind, ylim=c(-58, -0), xaxt='n', yaxt='n', col=fields::tim.colors()) #c(-50,0) c(-20, 60)
#  fields::world(add=TRUE, col='gray70')
#  points(pt, col='white', pch=1, cex=0.5)
#  
#  rw <- a[2]
#  cond.ind <- matrix( S.n[rw,], nr=length(lonNS), nc=length(latNS)) #Re(gB) #emp.cor
#  pt <- grd[rw, ]
#  image(lonNS, latNS, cond.ind, ylim=c(-58, -0), xaxt='n', yaxt='n', col=fields::tim.colors()) #c(-50,0) c(-20, 60)
#  fields::world(add=TRUE, col='gray70')
#  points(pt, col='white', pch=1, cex=0.5)
#  
#  rw <- a[3]
#  cond.ind <- matrix( S.n[rw,], nr=length(lonNS), nc=length(latNS)) #Re(gB) #emp.cor
#  pt <- grd[rw, ]
#  image(lonNS, latNS, cond.ind, ylim=c(-58, -0), xaxt='n', yaxt='n', col=fields::tim.colors()) #c(-50,0) c(-20, 60)
#  fields::world(add=TRUE, col='gray70')
#  points(pt, col='white', pch=1, cex=0.5)
#  
#  rw <- a[4]
#  cond.ind <- matrix( S.n[rw,], nr=length(lonNS), nc=length(latNS)) #Re(gB) #emp.cor
#  pt <- grd[rw, ]
#  image(lonNS, latNS, cond.ind, ylim=c(-58, -0), xaxt='n', yaxt='n', col=fields::tim.colors()) #c(-50,0) c(-20, 60)
#  fields::world(add=TRUE, col='gray70')
#  points(pt, col='white', pch=1, cex=0.5)
#  
#  

