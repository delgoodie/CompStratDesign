library(MASS)
# library(ggb)
# library(igraph)


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
    Sig
}

ggb_method <- function(S, g, lambda, max_depths = NULL, maxiter = 500, tol = 1e-4, verbose = 1) {
    
    adj <- as.matrix(igraph::as_adjacency_matrix(g))
    
    sumvv <- CompStratDesign::bcd_method(S, adj, lambda, maxiter, tol, verbose)

    sumV <- Matrix::Matrix(0, p, p, sparse = TRUE)
    sumV[lower.tri(sumV)] <- sumvv
    sumV <- Matrix::forceSymmetric(sumV, uplo = "L")

    diag(sumV) <- diag(S)

    return(as.matrix(sumV))
}


test <- function() {
    # generate Sig
    g <- igraph::graph.lattice(c(3, 2))
    b <- rep(1, 3 * 2)
    Sig <- generate_gb_covariance(g, b)
    p <- nrow(Sig)
    n <- 30
    
    # create observations
    # set.seed(123)
    eig <- eigen(Sig)
    A <- diag(sqrt(eig$values)) %*% t(eig$vectors)
    x <- matrix(rnorm(n * p), n, p) %*% A
    
    # sample covariance matrix
    S <- cov(x)
    
    # sample Frobenius Norm
    covNorm <- norm(S - Sig, 'F')
    
    # ggb package covariance matrix
    Sighat <- ggb_method(S, g, .05)
    
    actualSighat <- as.matrix(ggb::ggb(S, g, 'local', .05)$Sig[[1]])
    
    print('Dif of our ggb vs ggb package')
    print(norm(actualSighat - Sighat))
    
    
    # ggb Frobenius Norm
    ggbNorm <- norm(Sighat - Sig, type='F')
    
    
    print('')
    print('')
    print('Dif of our ggb norm vs sample cov norm (higher number is better)')
    print(covNorm - ggbNorm)
}


test()




#inverse of diagonal of S
