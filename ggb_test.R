library(MASS)


generate_gb_covariance <- function(g, b, sigma = 0.01, cor = TRUE) {
    p <- length(b)
    stopifnot(igraph::vcount(g) == p)
    stopifnot(b >= 0, b == round(b))
    D <- igraph::shortest.paths(g)
    b <- pmin(b, apply(D, 1, max))
    A <- D <= b # A_jk indicates whether k is in N_j
    B <- A | t(A)
    Sig <- B / D
    diag(Sig) <- 0
    # add to diag so min eval=sig^2
    diag(Sig) <- max(-min(eigen(Sig)$val), 0) + sigma^2
    if (cor) {
        sig <- sqrt(diag(Sig))
        D <- diag(1 / sig)
        Sig <- D %*% Sig %*% D
    }
    return(Sig)
}


# Generate Covariance Matrix and Graph
g <- igraph::graph.lattice(c(3, 2))
b <- rep(1, 3 * 2)
Sig <- generate_gb_covariance(g, b)
p <- nrow(Sig)
adj_mat <- as.matrix(igraph::as_adjacency_matrix(g))


# create observations
n <- 30
# optional fixed seed for testing
# set.seed(123)
eig <- eigen(Sig)
A <- diag(sqrt(eig$values)) %*% t(eig$vectors)
x <- matrix(rnorm(n * p), n, p) %*% A

# Sample Covariance Matrix
S <- cov(x)

# Error of Sample Covariance Matrix
covNorm <- norm(S - Sig, 'F')

# our ggb covariance matrix (call to C++ function in rcpparma.cpp)
lambda <- .05
Sighat <- CompStratDesign::bcd_method(S, adj_mat, lambda)

# ggb package covariance matrix
packageSighat <- as.matrix(ggb::ggb(S, g, 'local', lambda)$Sig[[1]])

print('Dif of our ggb vs ggb package')
print(norm(packageSighat - Sighat, type='F'))


# ggb Frobenius Norm
ggbNorm <- norm(Sighat - Sig, type='F')


print('')
print('')
print('Dif of our ggb norm vs sample cov norm (higher number is better)')
print(covNorm - ggbNorm)


#starting S for iter
#inverse of diagonal of S
