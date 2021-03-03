# #' Climate model ensemble
# #'
# #' A dataset containing 30 simulations of a climate model ensemble
# #'
# #' @format An array Ufield with the simulations, their mean Umean, defined on a
# #' grid expand.grid(Longitude, Latitude) (which are vectors), and covariate
# #' landMask
# "BRACEUfields.rda"


#' Local Matern covariance parameter estimates
#'
#' A dataset containing 30 simulations of a climate model ensemble
#'
#' @format Local parameter estimates: thx, thy, and ang defining the local
#' geometric anisotropy, s2 and t2 defining the process and white noise
#' variances. Estimated on a grd (two columns matrix of locations)
"nsMaternEstimates"
