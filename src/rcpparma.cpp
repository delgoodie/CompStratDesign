#include "RcppArmadillo.h"

// modified BCD Function from Jacob Bien github: https://github.com/jacobbien/ggb
extern "C" void pathwiseprox_bcd3(double *rr, int *dd, int *M, double *lambda, int *p, double *sumvv, double *obj, double* penalty, int *maxiter, double *tol, int *verbose, int* out_iter, int* did_converge);


// Helper Functions

// finding minimum distance
int min_dist(const arma::imat& dist, const arma::ivec& t_set, int r) {
    int mi = 0, mv = INT_MAX;
    
    for(int i=0;i<(int)dist.n_cols;i++) 
    {
        if(!t_set(i) && dist(r, i) <= mv)      
        {
            mv = dist(r, i);
            mi = i;
        }
    }
    return mi;
}

// adjacency matrix 
arma::imat shortest_path(arma::imat graph) {
    int rows = graph.n_rows, cols=rows;
    if (rows != cols) return arma::imat(1, 1);
    
    arma::imat dist(rows, cols);                             
    arma::ivec t_set(cols);
    
    for(int i=0;i<rows;i++)
        for(int j=0;j<cols;j++)
        dist(i, j) = INT_MAX;
    

    for (int r=0;r<rows;r++) {
        dist(r, r) = 0;
        for(int i=0;i<cols;i++) t_set(i)=0;
        
        for(int i=0;i<cols;i++)                           
        {
            int m = min_dist(dist,t_set, r);
            t_set(m) = 1;
            
            for(int k=0;k<rows;k++)                  
            {
                if(!t_set(k) && graph(m, k) && dist(r, m) != INT_MAX && dist(r, m) + graph(m, k) < dist(r, k))
                    dist(r, k) = dist(r, m) + graph(m, k);
            }
        }
    }
    
    return dist;
}

// Get lower or upper triangle of matrix in vector form
arma::vec tri(const arma::mat& matrix, char s) {
    int rows = matrix.n_rows;
    if (rows != (int)matrix.n_cols || rows < 2) return arma::vec(1);

    int size = (rows * rows - rows) / 2 + rows % 2;

    arma::vec vec(size);

    for(int i=0,c=0;i<rows-1;i++) {
        for(int j=1+i;j<rows;j++,c++) {
            if (s == 'L')
                vec(c) = matrix(j, i);
            else
                vec(c) = matrix(i, j);
        }
    }
    return vec;
}

// Get lower or upper triangle of integer matrix in integer vector form
arma::ivec tri(const arma::imat& matrix, char s) {
    int rows = matrix.n_rows;
    if (rows != (int)matrix.n_cols || rows < 2) return arma::ivec(1);

    int size = (rows * rows - rows) / 2 + rows % 2;

    arma::ivec vec(size);

    for(int i=0,c=0;i<rows-1;i++) {
        for(int j=1+i;j<rows;j++,c++) {
            if (s == 'L')
                vec(c) = matrix(j, i);
            else
                vec(c) = matrix(i, j);
        }
    }
    return vec;
}

// Create symetrical matrix from upper or lower triangle
arma::mat from_tri(const arma::vec& vec, int rows) {
    arma::mat mat(rows, rows);

    for(int i=0,c=0;i<rows;i++)
        for(int j=i+1;j<rows;j++,c++) {
            mat(i, j) = vec(c);
            mat(j, i) = vec(c);
        }
    return mat;
}

static inline float max(float a, float b) { return a > b ? a : b; }



// Majorizor Functions

// [[Rcpp::export]]
float Majorizor_EM(arma::mat sig_til, arma::mat sig, arma::mat S, float tau) {
    arma::mat W_1 = sig_til + tau * arma::mat(sig_til.n_rows, sig_til.n_cols, arma::fill::eye);
    arma::mat W_2 = sig_til;
    
    arma::mat W_1_inv = W_1.i();
    arma::mat W_2_inv = W_2.i();
    
    // are these parenthesis in correct order of operations?
    arma::mat Q = W_2_inv - W_1_inv + W_1_inv * (S * W_1_inv);
    
    return -log(arma::det(sig))
        + 1/tau * arma::trace(sig * (Q * sig))
        + arma::trace(sig * Q)
        - 2/tau * arma::trace(sig * (W_1_inv * S));
}

// [[Rcpp::export]]
arma::mat Majorizor_EM_grad(arma::mat sig_til, arma::mat sig, arma::mat S, float tau) {
    arma::mat W_1 = sig_til + tau * arma::mat(sig_til.n_rows, sig_til.n_cols, arma::fill::eye);
    arma::mat W_2 = sig_til;
    
    arma::mat W_1_inv = W_1.i();
    arma::mat W_2_inv = W_2.i();
    
    arma::mat Q = W_2_inv - W_1_inv + W_1_inv * (S * W_1_inv);
    
    return -1/tau * (S * W_1_inv + (S * W_1_inv).t() - sig * Q - (sig * Q).t()) + Q - sig.i();
}


// [[Rcpp::export]]
float Majorizor_linear(arma::mat sig_til, arma::mat sig, arma::mat S) {
    return arma::trace(sig_til.i() * sig) + arma::trace(S * sig.i());
}

// [[Rcpp::export]]
arma::mat Majorizor_linear_grad(arma::mat sig_til, arma::mat sig, arma::mat S) {
    arma::mat sig_inv = sig.i();
    arma::mat out = sig_til.i() - sig_inv * (S * sig_inv);
    
    out = (out + out.t()) / 2;
    
    return out;
}


// C++ Functions

/* Design:
 * The C implementation for each function is extracted to a C-only function (marked _implementation) which fills the CPPResult struct
 * Then a separate function 'arma::vec pack_vector(CPPResult* cppres)' converts this struct into a vector (since R can't handle a C-struct)
 * In R we convert the vector back into an R class to get the final result
 * 
 * The functions that are called by R are implemented at the end of the script and are marked with [[Rcpp::export]]
 * 
 */


// C++ Struct which stores result data for each algorithm
struct CPPResult {
    arma::mat sigma;
    arma::vec objective;
    int iterations;
    double penalty;
    bool did_converge;
};

// Converts the C++ struct into a vector by placing each variable from the struct into a vector in order
arma::vec pack_vector(CPPResult* cppres) {
    arma::vec result(3 + (int)cppres->sigma.n_elem + cppres->iterations);
    
    int i = 0;
    
    result(i++) = cppres->did_converge;
    result(i++) = cppres->iterations;
    result(i++) = cppres->penalty;
    
    for(int x = 0; x < (int)cppres->sigma.n_rows; x++) for (int y = 0; y < (int)cppres->sigma.n_cols; y++) result(i++) = cppres->sigma(x, y);
    
    for(int j = 0; j < cppres->iterations; j++) result(i++) = cppres->objective(j);
    
    return result;
}

// Solves the proximal problem (basically a direct call the the C bcd function in bcd.c written by Jacob Bien)
void prox_implementation(arma::mat S, arma::imat G, double lambda, int maxiter, double tol, int verbose, CPPResult* result) {
    int p = S.n_rows;
    
    // Finds the shortest paths between every node in the adjacency matrix of the graph (stored as an integer matrix)
    arma::imat D = shortest_path(G);
    
    arma::ivec M(D.n_cols);
    
    // Finds the max depths for each variable
    for(int i=0;i<(int)D.n_cols;i++) {
        int mi = 0;
        for(int j=0;j<(int)D.n_rows;j++)
            if (D(j, i) != INT_MAX && D(j, i) >= D(mi, i))
                mi = j;
        M(i) = D(i, mi);
    }
    
    int M_max = M.max();
    for(int i=0;i<(int)D.n_rows;i++)
        for(int j=0;j<(int)D.n_cols;j++)
            if (D(i, j) == INT_MAX)
                D(i, j) = M_max + 1;
    
    // Converts S into a lower triangle vector (only form accepted by bcd.c)
    arma::vec ss_vec = tri(S, 'L');
    
    // Converts D into a lower triangle vector (only form accepted by bcd.c)
    arma::ivec dd_vec = tri(D, 'L');
    
    // allocates sumvv memory (for use by bcd.c)
    int n_sumvv = p * (p - 1) / 2;
    double *sumvv = new double[n_sumvv];

    // allocates objective memory (for use by bcd.c)    
    double *obj = new double[maxiter];
    double penalty = 0.f;
    
    // allocates return parameters to be set by bcd.c
    int out_iter = -1;
    int did_converge = false;
    
    // call to C function in bcd.c, all arguments are passed by pointer
    pathwiseprox_bcd3(ss_vec.memptr(), dd_vec.memptr(), M.memptr(), &lambda, &p, sumvv, obj, &penalty, &maxiter, &tol, &verbose, &out_iter, &did_converge);
    
    
    // Fill CPPResult struct (passed by pointer)
    
    result->did_converge = did_converge;
    result->iterations = out_iter;
    result->penalty = penalty;
    
    result->objective = arma::vec(obj, out_iter);
    
    free(obj);
        
    arma::vec sumvv_vec(n_sumvv);
    for (int i = 0; i < n_sumvv; i++) sumvv_vec(i) = sumvv[i];
    
    free(sumvv);
    
    arma::mat sig = from_tri(sumvv_vec, S.n_rows);
    for (int i=0;i<(int)sig.n_rows;i++) sig(i, i) = S(i, i);
    
    result->sigma = sig;
}

// Solves the proximal problem w/ psd constriant (this is the same problem which is solved in the main R script ggb.R in Jacob Bien's ggb package)
void ggb_psd_implementation(arma::mat S, arma::imat G, double lambda, float delta, int maxiter, double tol, int verbose, CPPResult* result) {

    arma::mat oldSig;
    arma::mat R, C(S.n_rows, S.n_cols, arma::fill::zeros);
    arma::mat out;
    arma::mat B;
    double penalty = 0.f;
    
    bool did_converge = false;
    double *obj = new double[maxiter];
    int i = 0;
    
    for (i = 0; i < maxiter; i++) {
        // update over B: solve Prox on S + C
        R = S + C;
        
        CPPResult res;
        prox_implementation(R, G, lambda, maxiter, tol, verbose, &res);
        out = res.sigma;
        penalty = res.penalty;
        
        //Rcpp::Rcout << "|out|" << arma::norm(out - out.t(), "fro") << "\n";
        
        B = R - out;
        
        
        // checks max dif between each element to determine if we should stop iterating
        float max_dif = -1e10; 
        
        if ((int)oldSig.n_elem > 0)
            for(int j=0;j<(int)out.n_rows;j++)for(int k=0;k<(int)out.n_cols;k++) max_dif = max(max_dif, abs(out(j,k)-oldSig(j,k)));
        
        
        if (i > 1 && max_dif < tol) {
            if (verbose > 0) Rcpp::Rcout << "ggb_psd Converged after " << i << " iterations\n";
            did_converge = true;
            i--;
            break;
        }
        
        oldSig = out;
        
        // update over C: adjust eigenvalues
        arma::vec eigen_val;
        arma::mat eigen_vec;
        
        
        //Rcpp::Rcout << "|B|" << arma::norm(B - B.t(), "fro") << "\n";
        
        
        arma::eig_sym(eigen_val, eigen_vec, B - S);

        
        
        arma::mat diag_mat(eigen_val.n_elem, eigen_val.n_elem);
        for(int j = 0; j < (int)eigen_val.n_elem; j++) diag_mat(j, j) = max(eigen_val(j) + delta, 0.f);
        
        C = eigen_vec * diag_mat * eigen_vec.t();
        
        // Calculate objective for current iterations
        // obj[l] <- out$obj + (sum((S - Sig[[l]])^2) - sum((R - Sig[[l]])^2))/2
        // return penalty + pow(arma::accu(S - oldSig), 2) - sum((R - oldSig)^2))/2;
        // res.objective(res.iterations - 1) + (arma::accu(pow(S - res.sigma, 2)) - arma::accu(pow(R - res.sigma, 2))) / 2;
        obj[i] = log(arma::det(res.sigma)) + arma::trace(S * res.sigma.i()) + lambda * penalty; 
    }
    
    // Fill CPPResult struct (passed by pointer)
    result->did_converge = did_converge;
    result->iterations = i;
    result->penalty = penalty;
    result->sigma = S + C - B;
    result->objective = arma::vec(obj, i);
}

// Solve proximal decent with a linear majorizer function, iterating over the ggb_psd method 
void ggb_mm_linear_implementation(arma::mat S, arma::mat sig_til, arma::imat G, double lambda, double tau, double t, double B, int maxiter, double tol, int verbose, CPPResult* result) {
    
    // initialize variables
    arma::mat sig_k = S; // diag(nrow(S)) # sigma of current iteration
    arma::mat sig_k_1 = S; // sigma of last iteration
    bool did_converge = false;
    int it = 0;
    double* obj = new double[maxiter];
    double penalty = 0.f;
    
        
    for (it = 0; it < maxiter; it++) {
        arma::mat grad = Majorizor_linear_grad(sig_til, sig_k, S); // cache gradient
        
        float mm_sig_k_1 = Majorizor_linear(sig_til, sig_k, S); // cache majorizer
        
        arma::mat new_S = sig_k - t * grad;
        
        CPPResult res;
        ggb_psd_implementation(new_S, G, t * lambda, tau, maxiter, tol, verbose, &res);
        arma::mat sig_bar = res.sigma;
        
        penalty = res.penalty;
        obj[it] = Majorizor_linear(sig_til, sig_bar, S) + lambda * penalty; // (res.objective(res.iterations - 1) - 1/2 * pow(arma::norm(new_S - sig_bar, "fro"), 2)) /  t;
        
        // penalty <- (ggb_out$obj - 1/2 * norm(new_S - sig_bar, 'F') ** 2) / (lambda * t)
        // objective <- c(objective, calc_obj(penalty, sig_bar))
        
        //print(calc_obj(penalty, sig_bar))
        
        // ===== BACKTRACK Line Search =============
        
        for (int t_it = 0; t_it < 5; t_it++) {
            
            arma::mat G_t = (sig_k - sig_bar) / t;
            
            float mm_sig_bar = Majorizor_linear(sig_til, sig_bar, S);
            
            float term_2 = -t * arma::trace(grad.t() * G_t);
            
            float term_3 = t/2 * pow(arma::norm(G_t, "fro"), 2);
            
            if (mm_sig_bar > mm_sig_k_1 + term_2 + term_3) // backtracking condition as specified in paper
            {
                t *= B;
                
                new_S = sig_k - t * grad;
                
                ggb_psd_implementation(new_S, G, t * lambda, tau, maxiter, tol, verbose, &res);
                sig_bar = res.sigma;        
                
                penalty = res.penalty;
                obj[it] = Majorizor_linear(sig_til, sig_bar, S) + lambda * (res.objective(res.iterations - 1) - 1/2 * pow(arma::norm(new_S - sig_bar, "fro"), 2)) /  t;
            }
            else break;
        }
        
        // ============ set sig for next iteration ===============
        sig_k_1 = sig_k;
        sig_k = sig_bar;
        
        
        if (it > 0 && arma::norm(sig_k - sig_k_1, "fro") / arma::norm(sig_k_1, "fro") < tol) {
            did_converge = true;
            break;
        }
    }
    
    // Fill CPPResult struct (passed by pointer)
    result->did_converge = did_converge;
    result->iterations = it;
    result->penalty = penalty;
    result->sigma = sig_k;
    result->objective = arma::vec(obj, it);
}

// Solve proximal decent with an EM majorizer function, iterating over the ggb_psd method
void ggb_mm_EM_implementation(arma::mat S, arma::mat sig_til, arma::imat G, double lambda, double tau, double t, double B, int maxiter, double tol, int verbose, CPPResult* result) {
    
    // THIS FUNCTION SAME AS GGB_MM_LINEAR but all calls to majorizer functions are swapped from linear to EM
    
    // initialize variables
    arma::mat sig_k = S; // diag(nrow(S)) # sigma of current iteration
    arma::mat sig_k_1 = S; // sigma of last iteration
    bool did_converge = false;
    int it = 0;
    double* obj = new double[maxiter];
    double penalty = 0.f;
    
    
    for (it = 0; it < maxiter; it++) {
        arma::mat grad = Majorizor_EM_grad(sig_til, sig_k, S, tau);
        
        float mm_sig_k_1 = Majorizor_EM(sig_til, sig_k, S, tau);
        
        arma::mat new_S = sig_k - t * grad;
        
        CPPResult res;
        ggb_psd_implementation(new_S, G, t * lambda, tau, maxiter, tol, verbose, &res);
        arma::mat sig_bar = res.sigma;
        
        penalty = res.penalty;
        obj[it] = Majorizor_EM(sig_til, sig_bar, S, tau) + lambda * penalty; // (res.objective(res.iterations - 1) - 1/2 * pow(arma::norm(new_S - sig_bar, "fro"), 2)) /  t;
        
        // penalty <- (ggb_out$obj - 1/2 * norm(new_S - sig_bar, 'F') ** 2) / (lambda * t)
        // objective <- c(objective, calc_obj(penalty, sig_bar))
        
        //print(calc_obj(penalty, sig_bar))
        
        // ===== BACKTRACK Line Search =============
        
        for (int t_it = 0; t_it < 5; t_it++) {
            
            arma::mat G_t = (sig_k - sig_bar) / t;
            
            float mm_sig_bar = Majorizor_EM(sig_til, sig_bar, S, tau);
            
            float term_2 = -t * arma::trace(grad.t() * G_t);
            
            float term_3 = t/2 * pow(arma::norm(G_t, "fro"), 2);
            
            if (mm_sig_bar > mm_sig_k_1 + term_2 + term_3)
            {
                t *= B;
                
                new_S = sig_k - t * grad;
                
                ggb_psd_implementation(new_S, G, t * lambda, tau, maxiter, tol, verbose, &res);
                sig_bar = res.sigma;        
                
                penalty = res.penalty;
                obj[it] = Majorizor_EM(sig_til, sig_bar, S, tau) + lambda * (res.objective(res.iterations - 1) - 1/2 * pow(arma::norm(new_S - sig_bar, "fro"), 2)) /  t;
            }
            else break;
        }
        
        // ============ set sig for next iteration ===============
        sig_k_1 = sig_k;
        sig_k = sig_bar;
        
        
        if (it > 0 && arma::norm(sig_k - sig_k_1, "fro") / arma::norm(sig_k_1, "fro") < tol) {
            did_converge = true;
            break;
        }
    }
    
    result->did_converge = did_converge;
    result->iterations = it;
    result->penalty = penalty;
    result->sigma = sig_k;
    result->objective = arma::vec(obj, it);
}

// R Export Functions

// [[Rcpp::export]]
arma::vec c_prox(arma::mat S, arma::imat G, double lambda, int maxiter, double tol, int verbose) {
    CPPResult cppres;
    prox_implementation(S, G, lambda, maxiter, tol, verbose, &cppres); // call to prox implementation to fill cppres, passed by pointer
    return pack_vector(&cppres); // convert cppres to packed vector for R
}

// [[Rcpp::export]]
arma::vec c_ggb_psd(arma::mat S, arma::imat G, double lambda, float delta, int maxiter, double tol, int verbose) {
    CPPResult cppres;
    ggb_psd_implementation(S, G, lambda, delta, maxiter, tol, verbose, &cppres); // call to ggb psd implementation to fill cppres, passed by pointer
    return pack_vector(&cppres); // convert cppres to packed vector for R
}

// [[Rcpp::export]]
arma::vec c_ggb_mm_linear(arma::mat S, arma::mat sig_til, arma::imat G, double lambda, double tau, double t, double B, int maxiter, double tol, int verbose) {
    CPPResult cppres;
    ggb_mm_linear_implementation(S, sig_til, G, lambda, tau, t, B, maxiter, tol, verbose, &cppres); // call to ggb mm linear implementation to fill cppres, passed by pointer
    return pack_vector(&cppres); // convert cppres to packed vector for R
}

// [[Rcpp::export]]
arma::vec c_ggb_mm_EM(arma::mat S, arma::mat sig_til, arma::imat G, double lambda, double tau, double t, double B, int maxiter, double tol, int verbose) {
    CPPResult cppres;
    ggb_mm_EM_implementation(S, sig_til, G, lambda, tau, t, B, maxiter, tol, verbose, &cppres); // call to ggb mm EM implementation to fill cppres, passed by pointer
    return pack_vector(&cppres); // convert cppres to packed vector for R
}