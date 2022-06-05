source("functions.R")
n = 200
p = 5000

DCmat = matrix(NA, 500, 10)
MMmat = matrix(NA, 500, 10)
HUmat = matrix(NA, 500, 10)
CAmat = matrix(NA, 500, 10)

for (set in 1:500) {
  pr = runif(1, 0.05, 0.95)
  Y = sample(0:1, n, replace = TRUE, prob = c(pr, 1 - pr))
  prob = matrix(c(0.3, 0.8, 0.6, 0.7, 0.2, 0.6, 0.4, 0.6, 0.4, 0.7,
                  0.6, 0.4, 0.2, 0.4, 0.7, 0.1, 0.7, 0.2, 0.7, 0.4), nrow = 2, byrow = TRUE)
  X = matrix(NA, n, p)
  for (i in 1:10) {
    X[Y == 0, i] = rbinom(length(which(Y == 0)), 2, prob[1, i])
    X[Y == 1, i] = rbinom(length(which(Y == 1)), 2, prob[2, i])
  }
  for (i in 11:p) {
    pr = runif(1, 0.05, 0.95)
    X[, i] = rbinom(n, 2, pr)
  }
  data = cbind(Y, X)
  
  # compare results
  DCres = DCSIS(data)
  DCmat[set, ] = (p - rank(DCres, ties.method = "first") + 1)[1:10]
  MMres = MMLE(data)
  MMmat[set, ] = (p - rank(MMres, ties.method = "first") + 1)[1:10]
  HUres = Huang(data)
  HUmat[set, ] = (p - rank(HUres, ties.method = "first") + 1)[1:10]
  CAres = CATrend(data)
  CAmat[set, ] = (p - rank(CAres, ties.method = "first") + 1)[1:10]
  
  write.csv(data, paste("Simu_data", set, ".csv", sep = ""), row.names = FALSE)
  write.csv(DCmat, "DCmat", row.names = FALSE)
  write.csv(MMmat, "MMmat", row.names = FALSE)
  write.csv(HUmat, "HUmat", row.names = FALSE)
  write.csv(CAmat, "CAmat", row.names = FALSE)
}
