library(MASS)
library(matlib)

generate_gb_covariance <- function(g, b, sigma = 0.01, cor = TRUE) {
    # set.seed(0)
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

ggb_mm_R <- function(S, sig_til, adj_mat, t,  tol, maxiter, lambda, tau, prox_maxiter=500, prox_tol=1e-4, verbose=1) {

    # initialize variables
    penalty <- 0 # this is a pass by reference value that the C function prox() will set
    sig_k <- S # diag(nrow(S)) # sigma of current iteration
    sig_k_1 <- S #sigma of last iteration
    first_iter <- TRUE
    init_t <- t
    #obj <- c() # objective list
    it <- 0 # iteration number
    
    # Helper functions
    G_t <- function(E_k_1, E_k, t) { (E_k_1 - E_k) / t }
    
    calc_obj <- function(pen, sig) { CompStratDesign::Majorizor_linear(sig_til, sig, S) + lambda * pen }
    
    objective <- c()
    
    while((first_iter || (norm(sig_k - sig_k_1, 'F') / norm(sig_k_1, 'F') > tol)) && it < maxiter) {
        
        #print(norm(sig_k - sig_k_1, 'F') / norm(sig_k_1, 'F'))
        
        first_iter<-FALSE
        it <- it + 1
        #t <- init_t
        

    
        grad <- CompStratDesign::Majorizor_linear_grad(sig_til, sig_k, S)
        
        #print(norm(grad - t(grad), 'F'))
        
        mm_sig_k_1 <- CompStratDesign::Majorizor_linear(sig_til, sig_k, S)
        
        
        new_S <- sig_k - t * grad
        #print(norm(new_S - t(new_S), 'F'))
        
        
        ggb_out <- ggb::ggb(new_S, adj_mat, "local", c(lambda * t), tau, out_iter = maxiter, in_iter = maxiter, out_tol = tol, in_tol = tol, verbose=0)
        
        sig_bar <- as.matrix(ggb_out$Sig[[1]])
       #print(norm(sig_bar - sig_k_1, 'F') / norm(sig_k_1, 'F'))
        
        
        penalty <- (ggb_out$obj - 1/2 * norm(new_S - sig_bar, 'F') ** 2) / (lambda * t)
        
        objective <- c(objective, calc_obj(penalty, sig_bar))
        
        #print(calc_obj(penalty, sig_bar))
        
        
            
        #CompStratDesign::ggb_psd(new_S, adj_mat, t * lambda, penalty, prox_maxiter, prox_tol, tau, verbose)

        
                
        # ===== BACKTRACK Line Search =============
        
        t_iter <- 0
        while(t_iter < 5) {
            t_iter <- t_iter + 1
            mm_sig_bar <- CompStratDesign::Majorizor_linear(sig_til, sig_bar, S)

            term_2 <- -t * tr(t(grad) %*% G_t(sig_k, sig_bar, t))

            term_3 <- t/2 * (norm(G_t(sig_k, sig_bar, t), 'F'))**2

            if (mm_sig_bar > mm_sig_k_1 + term_2 + term_3)
            {
                t <- t * B

                new_S <- sig_k - t * grad

                ggb_out <- ggb::ggb(new_S, adj_mat, "local", c(lambda * t), tau, out_iter = maxiter, in_iter = maxiter, out_tol = tol, in_tol = tol, verbose=0)

                sig_bar <- as.matrix(ggb_out$Sig[[1]])

                penalty <- (ggb_out$obj - 1/2 * norm(new_S - sig_bar, 'F') ** 2) / (lambda * t)

                append(objective, calc_obj(penalty, sig_bar))

                print("=============================Shrinking t=============================")
                print(t)

                # print(norm(sig_bar - t(sig_bar), 'F'))

            } else {
                break
            }
        }
        
        # ============ set sig for next iteration ===============
    
        sig_k_1 <- sig_k
        sig_k <- sig_bar

        
        #obj <- append(obj, calc_obj(sig_bar))    
    }
    
    print(objective)
    plot(objective)
    
    return(sig_k)
}





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


















# iterative method
t <- .05
lambda <- 0.001
tol <- 1e-14
maxiter <- 40
B <- .2
delta <- 1e-1
tau <- .3

objective <- rep(0, maxiter)

#theta_hat <- CompStratDesign::iter_method(Sig, adj_mat, t, tol, B, maxiter, lambda, objective)


penalty <- 0

sig_hat_R <- ggb_mm_R(S, Sig, g, t, tol, maxiter, lambda, tau)


#sigma_hat <- CompStratDesign::ggb_mm(S, S, adj_mat, t, tol, B, maxiter, lambda)



# our_prox_psd <- CompStratDesign::ggb_psd(S, adj_mat, lambda, penalty, maxiter, tol, delta)
# 
# ggb_prox_psd <- as.matrix(ggb::ggb(S, g, "local", c(lambda), delta, out_iter = maxiter, in_iter = maxiter, out_tol = tol, in_tol = tol, verbose=1)$Sig[[1]])
# 
# print(norm(our_prox_psd - ggb_prox_psd, 'F'))




