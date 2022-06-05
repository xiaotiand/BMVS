#' Estimate the marginal posterior probability
#'
#' P(Z = 1)
#'
#' @param iter_res output from BMVS()
#' @param burnin burn-in iterations
#' @return estimated marginal posterior probability
#' @examples
#' estProb(BMVS_res, burnin = 5000)
estProb = function(iter_res, burnin) {
  z_prob_vec = iter_res$z_prob_vec[(burnin + 1):length(iter_res$z_prob_vec)]
  z_prob = Reduce('+', z_prob_vec) / length(z_prob_vec)
  return(z_prob)
}


#' Estimate the regression coefficients
#'
#' Regression coefficients
#'
#' @param iter_res output from BMVS()
#' @param burnin burn-in iterations
#' @return estimate regression coefficients
#' @examples
#' estBeta(BMVS_res, burnin = 5000)
estBeta = function(iter_res, burnin) {
  beta_vec = iter_res$beta_vec[(burnin + 1):length(iter_res$beta_vec)]
  beta = Reduce('+', beta_vec) / length(beta_vec)
  return(beta)
}


#' Select an optimal model
#'
#' Select an optimal model using BIC.
#' The inputs are the output from BMVS() function call.
#'
#' @param iter_res output from BMVS()
#' @param Y input response variables
#' @param X input predictors
#' @param burnin burn-in iterations
#' @return indices of selected predictors
#' @examples
#' BMVS::select(BMVS_res, Y, X, burnin = 5000)
select = function(iter_res, Y, X, burnin) {
  n = dim(Y)[1]
  q = dim(Y)[2]
  z_prob = estProb(iter_res, burnin)
  beta = estBeta(iter_res, burnin)
  tbeta = matrix(0, dim(beta)[1], dim(beta)[2])
  rank = unlist(lapply(z_prob, function(x) min(which(x == sort(z_prob, decreasing = TRUE)))))
  nmodel = dim(X)[2]
  AIC = rep(0, nmodel)
  for (k in 1:nmodel) {
    tbeta[rank == k, ] = beta[rank == k, ]
    AIC[k] = n * log(det(t(Y) %*% (Y - X %*% tbeta) / n)) + q * (n + k) * n / (n - (k + q + 1))
  }
  best = which(AIC == min(AIC))
  return(which(rank <= best))
}
