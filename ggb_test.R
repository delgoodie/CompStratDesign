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
    #add to diag so min eval=sig^2
    diag(Sig) <- max(-min(eigen(Sig)$val), 0) + sigma^2
    if (cor) {
        sig <- sqrt(diag(Sig))
        D <- diag(1 / sig)
        Sig <- D %*% Sig %*% D
    }
    return(Sig)
}


ggb_mm_linear_R <- function(S, sig_til, adj_mat, lambda, tau, t, B, maxiter, tol, verbose) {

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
        
        
        new_obj <- CompStratDesign::Majorizor_linear(sig_til, sig_bar, S) + lambda * (ggb_out$obj - 1/2 * norm(new_S - sig_bar, 'F') ** 2) /  t
        
        objective <- c(objective, new_obj)
        
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

                
                new_obj <- CompStratDesign::Majorizor_linear(sig_til, sig_bar, S) + lambda * (ggb_out$obj - 1/2 * norm(new_S - sig_bar, 'F') ** 2) /  t
                objective <- c(objective, new_obj)
                
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
    
    #print(objective)
    plot(objective)
    
    print('R obj')
    print(objective[length(objective)])
    
    return(sig_k)
}


#========================== Generate Covariance Matrix and Graph ==========================
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







#========================== Initialize Parameters ==========================

S <- cor(x)
sig_til <- Sig
lambda <- 0.001
tau <- .3
t <- .05
B <- .2

tol <- 1e-4
maxiter <- 100

prox_tol <- 1e-2
prox_maxiter <- 500

verbose <- 0





#========================== Test ggb_mm_linear method C vs R ==========================

print(system.time(sig_hat_R <- ggb_mm_linear_R(S, sig_til, g, lambda, tau, t, B, maxiter, tol, verbose)))

print(system.time(ggb_mm_C_out <- ggb_mm_linear(S, sig_til, adj_mat, lambda, tau, t, B, maxiter, tol, verbose)))

print('C objective')
print(ggb_mm_C_out$objective[ggb_mm_C_out$iterations])

sig_hat_C <- ggb_mm_C_out$sigma

plot(ggb_mm_C_out$objective)

dif <- norm(sig_hat_R - sig_hat_C, 'F')
print(dif)







#========================== Compare our ggb_psd with ggb package ==========================

# our_prox_psd <- unpack_vector(CompStratDesign::ggb_psd(S, adj_mat, lambda, delta, maxiter, tol, verbose), nrow(S))$sigma
# ggb_prox_psd <- as.matrix(ggb::ggb(S, g, "local", c(lambda), delta, out_iter = maxiter, in_iter = maxiter, out_tol = tol, in_tol = tol, verbose=1)$Sig[[1]])
# print(norm(our_prox_psd - ggb_prox_psd, 'F'))



#========================== Test Main Algo ==========================
res <- main_algo(S, adj_mat, lambda, tol=1e-6)
print(res$sigma)
print(res$iterations)
plot(res$objective)
