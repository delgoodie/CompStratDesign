library(plotrix)
library(matlib)
library(MASS)

mu <- c(4.23, 3.01, 2.91)
mu2 <- rep(mu, 10)


stddev <- c(1.23, 0.92, 1.32)

corMat <- matrix(c(1, 0.78, 0.23, 0.78, 1, 0.27, 0.23, 0.27, 1), ncol = 3)

corMat2 <- matrix(0, nrow=30, ncol=30)


covMat <- stddev %*% t(stddev) * corMat

covMat2 <- matrix(0, nrow=30, ncol=30)


for(k in 1:10) {
    covMat2[((k-1)*3+1):(k*3),((k-1)*3+1):(k*3)] <- covMat
}

test_data <- mvrnorm(n=9, mu = mu2, Sigma = covMat2, empirical = FALSE)

fDist <- function(m) { return(sqrt(sum(m ** 2))) }

x <- 0:100 / 50- 1
y <- c()

for(lambda in x) {
    y <- c(y, fDist(EstimateCov(test_data, lambda) - covMat2))
}

plot(x, y)
