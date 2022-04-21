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

ggb_mm_R <- function(S, sig_til, adj_mat, t,  tol, maxiter, lambda, prox_maxiter=500, prox_tol=1e-4, verbose=1) {

    # initialize variables
    penalty <- 0 # this is a pass by reference value that the C function prox() will set
    sig_k <- S # sigma of current iteration
    sig_k_1 <- S #sigma of last iteration
    first_iter <- TRUE
    init_t <- t
    obj <- c() # objective list
    it <- 0 # iteration number
    
    # Helper functions
    G_t <- function(E_k_1, E_k, t) { (E_k_1 - E_k) / t }
    
    calc_obj <- function(sig) { CompStratDesign::Majorizor_linear(sig_til, sig, S) + lambda * penalty }
    

    
    while((first_iter || (norm(sig_k - sig_k_1, 'F') / norm(sig_k_1, 'F') > tol)) && it < maxiter) {
        first_iter<-FALSE
        it <- it + 1
        t <- init_t
    
        grad <- CompStratDesign::Majorizor_linear_grad(sig_til, sig_k, S)
        
        #print(norm(grad - t(grad), 'F'))
        
        mm_sig_k_1 <- CompStratDesign::Majorizor_linear(sig_til, sig_k, S)
        
        
        new_S <- sig_k_1 - t *grad
        #print(norm(new_S - t(new_S), 'F'))
        
        sig_bar <- CompStratDesign::ggb_psd(new_S, adj_mat, t * lambda, penalty, prox_maxiter, prox_tol, verbose)
        
        
        # ===== BACKTRACK Line Search =============
        
        while(TRUE) {
            mm_sig_bar <- CompStratDesign::Majorizor_linear(sig_til, sig_bar, S)
            
            term_2 <- -t * tr(t(grad) * G_t(sig_k, sig_bar, t))
            
            term_3 <- t/2 * (norm(G_t(sig_k, sig_bar, t), 'F'))**2
            
            if (mm_sig_bar > mm_sig_k_1 + term_2 + term_3) 
            {
                new_S <- sig_k - t * grad
                sig_bar = CompStratDesign::ggb_psd(new_S, adj_mat, t * lambda, penalty, prox_maxiter, prox_tol, verbose)
                # print("=============================Shrinking t=============================")
                # print(t)
                t <- t * B
                
                # print(norm(sig_bar - t(sig_bar), 'F'))
                
            } else {
                break
            }
        }
        
        # ============ set sig for next iteration ===============

        sig_k_1 <- sig_k
        sig_k <- sig_bar

        
        obj <- append(obj, calc_obj(sig_bar))    
    }
    
    plot(obj)
    
    return(sig_k)
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
n <- 300
# optional fixed seed for testing
# set.seed(123)
eig <- eigen(Sig)
A <- diag(sqrt(eig$values)) %*% t(eig$vectors)
x <- matrix(rnorm(n * p), n, p) %*% A

# Sample Covariance Matrix
S <- cor(x)

#s <- diag(Sig)
#Sig <- t(t(Sig / sqrt(s)) / sqrt(s))


#ggb_out <- ggb::ggb(Sig, g, "local", 1)

#eigen(ggb_out$Sig[[1]])$values

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
t <- .5
lambda <- 1e-4
tol <- 1e-9
maxiter <- 100
B <- .2
delta <- .01

objective <- rep(0, maxiter)

#theta_hat <- CompStratDesign::iter_method(Sig, adj_mat, t, tol, B, maxiter, lambda, objective)


penalty <- 0

sig_hat_R <- ggb_mm_R(S, S, adj_mat, t, tol, maxiter, lambda)

# sigma_hat <- CompStratDesign::ggb_mm(S, S, adj_mat, t, tol, B, maxiter, lambda)



print(norm(sigma_hat - sig_hat_R, 'F'))

# package_ggb <- ggb::ggb(new_S, g, "local", lambda, verbose=2)
# 
# print(package_ggb$obj)


# print(norm(inv(Theta_hat) - S, 'F'))
# 
# plot(log(objective))

#print(objective)

# err <-norm(Theta_hat -  inv(Sig), type='F')
# rel_err <- err / norm(inv(Sig), type='F')

#print(rel_err)# ggb package covariance matrix
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










