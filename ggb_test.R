library(MASS)
library(matlib)

generate_gb_covariance <- function(g, b, sigma = 0.01, cor = TRUE) {
    set.seed(0)
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

iter <- function(S, adj_mat, t,  tol, maxiter, lambda, ggb_maxiter=500, ggb_tol=1e-4, verbose=1) {
    theta_old <- diag(diag(inv(S)))
    
    penalty <- .1
    
    
    theta <- CompStratDesign::bcd_method(theta_old - t * (-inv(theta_old) + S), adj_mat, lambda, penalty, ggb_maxiter, ggb_tol, verbose)
    
    i <- 0
    
    while (norm(theta - theta_old, 'F') > tol && i < maxiter || i < 1) {
        i <- i + 1
        
        theta_old <- theta
        
        new_S <- theta_old - t * (-inv(theta_old) + S)
        
        theta <- CompStratDesign::bcd_method(new_S, adj_mat, lambda * t, penalty, ggb_maxiter, ggb_tol, verbose)
        
        print(penalty)
    }
    
    print(paste('R Iter Method converged after ', i, ' iterations\n'))
    
    return(inv(theta))
}



# penalty <- 0
# 
# package_ggb <- ggb::ggb(S, g, "local", lambda)
# 
# our_ggb <- CompStratDesign::bcd_method(S, adj_mat, lambda, penalty)
# 
# dif <- norm(as.matrix(package_ggb$Sig[[1]]) - our_ggb, "F")
# 
# print("Difference:")
# print(dif)


# Generate Covariance Matrix and Graph
g <- igraph::graph.lattice(c(3, 2))
b <- rep(1, 3 * 2)
Sig <- generate_gb_covariance(g, b)
p <- nrow(Sig)
print('Sigma')
print(Sig)
print('')
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

# objective <- rep(0, maxiter)
# 
# Theta_hat <- CompStratDesign::iter_method(S, adj_mat, t, tol, B, maxiter, lambda, objective)
# 
# print(norm(inv(Theta_hat) - S, 'F'))
# 
# plot(log(objective))
# 
# #print(objective)
# 
# err <-norm(Theta_hat -  inv(Sig), type='F')
# rel_err <- err / norm(inv(Sig), type='F')
# 
# print(rel_err)




# iterative method
t <- 1
lambda <- .1
tol <- 1e-4
maxiter <- 16
B <- .9



objective <- rep(0, maxiter)

new_S <- CompStratDesign::iter_method(S, adj_mat, t, tol, B, maxiter, lambda, objective)

package_ggb <- ggb::ggb(new_S, g, "local", lambda, verbose=2)

print(package_ggb$obj)


print(norm(inv(Theta_hat) - S, 'F'))

plot(log(objective))

#print(objective)

err <-norm(Theta_hat -  inv(Sig), type='F')
rel_err <- err / norm(inv(Sig), type='F')

print(rel_err)# ggb package covariance matrix
# packageSighat <- as.matrix(ggb::ggb(S, g, 'local', lambda)$Sig[[1]])


# ggb Frobenius Norm
# ggbNorm <- norm(packageSighat - Sig, type='F')


# print('')
# print('')
# print('Dif of our prox norm vs sample cov norm (higher number is better)')
# print(covNorm - proxNorm)


# print('')
# print('')
# print('Dif of ggb package norm vs sample cov norm (higher number is better)')
# print(covNorm - ggbNorm)




#shrink t if obj does not decrease - watch carnegie mellon video on backtracking - email Tianxi if doesn't make sense

#check if last sum of v_jb matricies == ouptut matrix (Covariance Matrix)

#collect and return vector of objective values for each iteration in arma::vec (also used by backtracking)










