library(mvtnorm)
library(LearnBayes)
library(MCMCpack)

#' Run one iteration of the MCMC chain
#'
#' This is an internal function called by BMVS()
#'
#' @param Y input response variables
#' @param X input predictors
#' @param n sample size
#' @param p the number of predictors
#' @param q the number of response variables
#' @return The results of one MCMC iteration
iter = function(Y, X, n, p, q, par_list, tao0, tao1, qn,
                # vectors of parameter values
                beta_vec, z_vec, z_prob_vec, sigma2_vec, invY_vec,
                # values of prior distribution parameters
                alpha1, alpha2, vk, deltak, cur) {
  beta = par_list$beta
  z = par_list$z
  z_prob = par_list$z_prob
  sigma2 = par_list$sigma2
  invY = par_list$invY

  t0 = proc.time()
  mu_m = X %*% beta
  for (i in 1:p) {
    ## beta
    Dk = matrix(0, q, q)
    if (z[i] == 1) {
      diagn = rep(1 / tao1, q)
    } else {
      diagn = rep(1 / tao0, q)
    }
    diag(Dk) = diagn / sigma2
    mu_m = mu_m - sweep(matrix(rep(beta[i, ], n), n, q, byrow = TRUE), MARGIN = 1, X[, i], "*")
    mu_i = solve(Dk + as.numeric(t(X[, i]) %*% X[, i]) * invY) %*%
      (apply(sweep(Y - mu_m, MARGIN = 1, X[, i], "*") %*% invY, 2, sum))
    Sigma_i = solve(Dk + as.numeric(t(X[, i]) %*% X[, i]) * invY)
    #Sigma_ai = Sigma_ai + diag(ncol(Sigma_ai)) * 0.1
    beta[i, ] = mvrnorm(1, mu_i, Sigma_i)
    mu_m = mu_m + sweep(matrix(rep(beta[i, ], n), n, q, byrow = TRUE), MARGIN = 1, X[, i], "*")

    ## z & z_prob
    z_prob[i] = (qn * dmvnorm(beta[i, ], rep(0, q), diag(q) * sigma2 * tao1)) /
      (qn * dmvnorm(beta[i, ], rep(0, q), diag(q) * sigma2 * tao1) +
         (1 - qn) * dmvnorm(beta[i, ], rep(0, q), diag(q) * sigma2 * tao0))
    #za_prob[is.na(za_prob)] = 1
    z[i] = rbinom(1, 1, z_prob[i])
  }

  ## sigma2
  Dk = matrix(0, p * q, p * q)
  diagn = rep(1 / tao0, p * q)
  diagn[as.vector(t(z)) == 1] = 1 / tao1
  betasum = as.numeric(as.vector(t(beta)) %*% Dk %*% as.vector(t(beta)))
  #meansum = sum(diag((Y - mu_m) %*% invY %*% t(Y - mu_m)))
  sigma2 = rigamma(1, alpha1 + p * q, alpha2 + betasum / 2)

  ## invY
  Z = Y - X %*% beta
  #lambdak
  #Smat = (L %*% t(Z)) %*% t(L %*% t(Z))
  #for (k in 1:q) {
  #  D[k, k] = rigamma(1, (vk[k] + n) / 2, (Smat[k, k] / sigma2 + deltak[k]) / 2)
  #}
  #ak
  #for (k in 2:q) {
  #  ak = rep(0, k - 1)
  #  invAk = matrix(0, k - 1, k - 1)
  #  diag(invAk) = 1 / D[k, k]
  #  sigma_ak = solve(invAk + (1 / D[k, k]) * t(Z[, 1:(k - 1)]) %*% Z[, 1:(k - 1)] / sigma2)
  #  #sigma_ak = solve(sigma_ak + diag(ncol(sigma_ak)) * 0.1)
  #  mu_ak = sigma_ak %*% (invAk %*% ak + (1 / D[k, k]) * t(Z[, 1:(k - 1)]) %*% Z[, k] / sigma2)
  #  L[k, 1:(k - 1)] = mvrnorm(1, mu_ak, sigma_ak)
  #}
  #invY = t(L) %*% solve(D) %*% L
  Sigma_invY = t(Z) %*% Z
  diag(Sigma_invY) = diag(Sigma_invY) + 1
  invY = solve(riwish(n + q + 1, Sigma_invY))
  proc.time() - t0

  beta_vec[[cur]] = beta
  z_vec[[cur]] = z
  z_prob_vec[[cur]] = z_prob
  sigma2_vec[cur] = sigma2
  invY_vec[[cur]] = invY
  cur = cur + 1

  par_list$beta = beta
  par_list$z = z
  par_list$z_prob = z_prob
  par_list$sigma2 = sigma2
  par_list$invY = invY

  result = list(cur = cur,
                par_list = par_list,
                beta_vec = beta_vec,
                z_vec = z_vec,
                z_prob_vec = z_prob_vec,
                sigma2_vec = sigma2_vec,
                invY_vec = invY_vec)
  return(result)
}


#' BMVS main function
#'
#' This is the main function of BMVS algorithm
#'
#' @param Y input response variables
#' @param X input predictors
#' @param iteration  the number of iterations of the MCMC chain
#' @return The MCMC simulation resutls. The results can be analyzed by estProb(), estBeta() and select().
#' @examples
#' data(tobacco)
#' Y = as.matrix(tobacco[, 1:3])
#' X = as.matrix(tobacco[, 4:9])
#' BMVS_res = BMVS(Y, X, iteration = 10000)
BMVS = function(Y, X, iteration = 3000) {
  n = dim(X)[1]
  p = dim(X)[2]
  q = dim(Y)[2]
  # overall mean
  #mu = apply(Y, 2, mean)
  #Yt = Y - matrix(rep(mu, n), n, q, byrow = TRUE)
  # tao & qn
  tao0 = sd(Y) ^ 2 / (10 * n)
  tao1 = sd(Y) ^ 2 * ((p * q) ^ 2.1) / (100 * n)
  qn = 0.1 #* (p / 500)

  # initialize parameters
  alpha1 = 0.1
  alpha2 = 0.1
  vk = 10 + 1:q
  deltak = rep(1, q)
  invY0 = cov(Y)
  par_list = list("beta" = matrix(0, p, q),
                  "z" = rep(0, p), "z_prob" = rep(0, p),
                  "sigma2" = 1, "invY" = invY0)
  beta_vec = list()
  z_vec = list()
  z_prob_vec = list()
  sigma2_vec = c()
  invY_vec = list()
  beta_vec[[1]] = par_list$beta
  z_vec[[1]] = par_list$z
  z_prob_vec[[1]] = par_list$z_prob
  sigma2_vec[1] = par_list$sigma2
  invY_vec[[1]] = par_list$invY
  cur = 1

  #t0 = proc.time()
  #iteration = 100
  for (index in 1:iteration) {
    iter_res = iter(Y, X, n, p, q, par_list, tao0, tao1, qn,
                    beta_vec, z_vec, z_prob_vec, sigma2_vec, invY_vec,
                    alpha1, alpha2, vk, deltak, cur)
    cur = iter_res$cur
    par_list = iter_res$par_list
    beta_vec = iter_res$beta_vec
    z_vec = iter_res$z_vec
    z_prob_vec = iter_res$z_prob_vec
    sigma2_vec = iter_res$sigma2_vec
    invY_vec = iter_res$invY_vec
  }
  #proc.time() - t0

  return(iter_res)
}
