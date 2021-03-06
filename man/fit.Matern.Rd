% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Matern.logLik.aniso.R
\name{fit.Matern}
\alias{fit.Matern}
\title{Wrapper function for fitting an isotropic Matern GP}
\usage{
fit.Matern(grd, y, NU = 1)
}
\arguments{
\item{grd}{a two columns matrix with the spatial locations of the data}

\item{y}{the value of the response variable at the locations in grd}

\item{NU}{a number > 0, the smoothness parameter}
}
\value{
an optim object with the parameter estimates of the isotropic Matern
}
\description{
Wrapper function for fitting an isotropic Matern GP
}
