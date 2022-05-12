
unpack_vector <- function(packed_vec, mat_length) {
    i <- 1
    
    did_converge <- packed_vec[i]
    i <- i+1
    iterations <- packed_vec[i] 
    i <- i+1
    penalty <- packed_vec[i]
    i <- i+1
    
    sigma <- matrix(packed_vec[i:(i+mat_length**2 - 1)], nrow=mat_length, ncol=mat_length)
    
    i <- i + mat_length ** 2
    
    objective <- packed_vec[i:length(packed_vec)]
    
    unpack_list <- list(did_converge = did_converge, iterations = iterations, penalty=penalty, sigma=sigma, objective=objective)
    
    class(unpack_list) <- 'CPPResult'
    
    return(unpack_list)
}


prox <- function(S, G, lambda, maxiter, tol, verbose) 
    { unpack_vector(.Call(`_CompStratDesign_c_prox`, S, G, lambda, maxiter, tol, verbose), nrow(S)) }


ggb_psd <- function(S, adj_mat, lambda, delta, maxiter, tol, verbose) 
    { unpack_vector(.Call(`_CompStratDesign_c_ggb_psd`, S, G, lambda, delta, maxiter, tol, verbose), nrow(S)) }


ggb_mm_linear <- function(S, sig_til, adj_mat, lambda, tau, t, B, maxiter, tol, verbose) 
    { unpack_vector(.Call(`_CompStratDesign_c_ggb_mm_linear`, S, sig_til, adj_mat, lambda, tau, t, B, maxiter, tol, verbose), nrow(S)) }


ggb_mm_EM <- function(S, sig_til, adj_mat, lambda, tau, t, B, maxiter, tol, verbose) 
    { unpack_vector(.Call(`_CompStratDesign_c_ggb_mm_EM`, S, sig_til, adj_mat, lambda, tau, t, B, maxiter, tol, verbose), nrow(S)) }


main_algo <- function(S, adj_mat, lambda, tau=.01, t=.5, B=.2, maxiter=50, tol=1e-4, prox_maxiter=100, prox_tol=1e-4, verbose=0) {
    
    sig_til <- S
    
    objective <- c()
    
    penalty <- 0
    
    for (iter in 1:maxiter) {
        
        out <- ggb_mm_linear(S, sig_til, adj_mat, lambda, tau, t, B, prox_maxiter, prox_tol, verbose)
        sig_hat <- out$sigma
        objective <- c(objective, out$objective[out$iterations])
        penalty <- out$penalty
        
        if (norm(sig_til - sig_hat, 'F') / norm(sig_til, 'F') < tol) {
            unpack_list <- list(did_converge = TRUE, iterations = iter, penalty=penalty, sigma=sig_hat, objective=objective)
            class(unpack_list) <- 'CPPResult'
            return(unpack_list)
        } else {
            sig_til <- sig_hat
        }
    }
    
    unpack_list <- list(did_converge = FALSE, iterations = maxiter, penalty=penalty, sigma=sig_til, objective=objective)
    class(unpack_list) <- 'CPPResult'
    return(unpack_list)
}