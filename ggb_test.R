library(MASS)
library(matlib)

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




#starting S for iter
#inverse of diagonal of S


iter <- function(S, adj_mat, t,  tol, maxiter, lambda, ggb_maxiter=500, ggb_tol=1e-4, verbose=1) {
    theta_old <- diag(diag(inv(S)))
    theta <- CompStratDesign::bcd_method(theta_old - t * (-inv(theta_old) + S), adj_mat, lambda, ggb_maxiter, ggb_tol, verbose)
    
    i <- 0
    
    while (norm(theta - theta_old, 'F') > tol && i < maxiter || i < 1) {
        i <- i + 1
        
        theta_old <- theta
        
        new_S <- theta_old - t * (-inv(theta_old) + S)
        
        theta <- CompStratDesign::bcd_method(new_S, adj_mat, lambda * t, ggb_maxiter, ggb_tol, verbose)
    }
    
    print(paste('R Iter Method converged after ', i, ' iterations\n'))
    
    return(inv(theta))
}













    # Error of Sample Covariance Matrix
covNorm <- norm(S - Sig, 'F')

# our ggb covariance matrix (call to C++ function in rcpparma.cpp)
t <- .2
lambda <- .05
tol <- 1e-4
maxiter <- 500


Sighat <- CompStratDesign::iter_method(S, adj_mat, t, tol, maxiter, lambda)

print(norm(Sighat -  iter(S, adj_mat, t, tol, maxiter, lambda), type='F'))
proxNorm <- norm(Sighat - Sig, type='F')


# ggb package covariance matrix
packageSighat <- as.matrix(ggb::ggb(S, g, 'local', lambda)$Sig[[1]])


# ggb Frobenius Norm
ggbNorm <- norm(packageSighat - Sig, type='F')


print('')
print('')
print('Dif of our prox norm vs sample cov norm (higher number is better)')
print(covNorm - proxNorm)


print('')
print('')
print('Dif of ggb package norm vs sample cov norm (higher number is better)')
print(covNorm - ggbNorm)















