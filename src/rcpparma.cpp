// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// include RcppArmadillo.h which pulls Rcpp.h
#include "RcppArmadillo.h"

//#include "bcd.h"

extern "C" void pathwiseprox_bcd3(double *rr, int *dd, int *M, double *lambda, int *p, double *sumvv, double *obj, float* penalty, int *maxiter, double *tol, int *verbose);
// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]
double interpolate(double a, double b, double t) {
    return a + (b-a) * t;
}

arma::mat interpolate_matrix(arma::mat m1, arma::mat m2, double t) {
    if (m1.n_rows != m2.n_rows || m1.n_cols != m2.n_cols) return m1;
    
    arma::mat res(m1.n_rows, m1.n_cols);
    
    for(int i=0; i < m1.n_rows; i++)
        for(int j=0;j<m1.n_cols;j++)
            res(i, j) = interpolate(m1(i, j), m2(i, j), t);
    
    return res;
}


// [[Rcpp::export]]
arma::mat EstimateCov(arma::mat data, double lambda) {
    arma::mat cov = arma::cov(data);
    arma::mat res = interpolate_matrix(cov, arma::mat(cov.n_rows, cov.n_cols, arma::fill::zeros), lambda);
    
    for(int i=0;i<cov.n_rows;i++) res(i, i) =  cov(i, i);
    
    return res;
}

int min_dist(const arma::imat& dist, const arma::ivec& t_set, int r) // finding minimum distance
{
    int mi = 0, mv = INT_MAX;
    
    for(int i=0;i<dist.n_cols;i++) 
    {
        if(!t_set(i) && dist(r, i) <= mv)      
        {
            mv = dist(r, i);
            mi = i;
        }
    }
    return mi;
}


arma::imat shortest_path(arma::imat graph) // adjacency matrix 
{
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


arma::vec tri(const arma::mat &matrix, char s) {
    int rows = matrix.n_rows;
    if (rows != matrix.n_cols || rows < 2) return arma::vec(1);
    
    int size = (rows * rows - rows) / 2 + rows % 2;
    
    arma::vec vec(size);
    
    for(int i=0,c=0;i<rows-1;i++)
        for(int j=1+i;j<rows;j++,c++) {
            if (s == 'L')
                vec(c) = matrix(j, i);
            else
                vec(c) = matrix(i, j);
        }
    
    return vec;
}

arma::ivec tri(const arma::imat &matrix, char s) {
    int rows = matrix.n_rows;
    if (rows != matrix.n_cols || rows < 2) return arma::ivec(1);
    
    int size = (rows * rows - rows) / 2 + rows % 2;
    
    arma::ivec vec(size);
    
    for(int i=0,c=0;i<rows-1;i++)
        for(int j=1+i;j<rows;j++,c++) {
            if (s == 'L')
                vec(c) = matrix(j, i);
            else
                vec(c) = matrix(i, j);
        }
        
        return vec;
}

arma::mat from_tri(const arma::vec& vec, int rows) {
    arma::mat mat(rows, rows);

    for(int i=0,c=0;i<rows;i++)
        for(int j=i+1;j<rows;j++,c++) {
            mat(i, j) = vec(c);
            mat(j, i) = vec(c);
        }
    return mat;
}



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

// [[Rcpp::export]]
arma::mat prox(arma::mat S, arma::imat G, double lambda, float& penalty, int maxiter=500, double tol=1e-4, int verbose=1) {
    int p = S.n_rows;
    
    arma::imat D = shortest_path(G);
    
    arma::ivec M(D.n_cols);
    
    for(int i=0;i<D.n_cols;i++) {
        int mi = 0;
        for(int j=0;j<D.n_rows;j++)
            if (D(j, i) != INT_MAX && D(j, i) >= D(mi, i))
                mi = j;
        M(i) = D(i, mi);
    }
    
    int M_max = M.max();
    for(int i=0;i<D.n_rows;i++)
        for(int j=0;j<D.n_cols;j++)
            if (D(i, j) == INT_MAX)
                D(i, j) = M_max + 1;
    
    arma::vec ss_vec = tri(S, 'L');
    
    arma::ivec dd_vec = tri(D, 'L');
    
    int n_sumvv = p * (p - 1) / 2;
    double *sumvv = new double[n_sumvv];
    
    double obj = 0;
    
    pathwiseprox_bcd3(ss_vec.memptr(), dd_vec.memptr(), M.memptr(), &lambda, &p, sumvv, &obj, &penalty, &maxiter, &tol, &verbose);
    
    arma::vec sumvv_vec(n_sumvv);
    for (int i = 0; i < n_sumvv; i++) sumvv_vec(i) = sumvv[i];
    
    free(sumvv);
    
    arma::mat res = from_tri(sumvv_vec, S.n_rows);
    for (int i=0;i<res.n_rows;i++) res(i, i) = S(i, i);
    
    return res;
}


static inline float max(float a, float b) { return a > b ? a : b; }

// [[Rcpp::export]]
arma::mat ggb_psd(arma::mat S, arma::imat G, double lambda, float& penalty, int maxiter, double tol, float delta, int verbose = 1) {

    arma::mat oldSig;
    arma::mat R, C(S.n_rows, S.n_cols, arma::fill::zeros);
    
    arma::mat out;
    
    arma::mat B;
    
    for (int i = 0; i < maxiter; i++) { 
        // update over B: solve Prox on S + C
        R = S + C;

        out = prox(R, G, lambda, penalty, maxiter, tol, verbose);
        
        //Rcpp::Rcout << "|out|" << arma::norm(out - out.t(), "fro") << "\n";
        
        B = R - out;
        
        float max_dif = -1e10; 
        
        if (oldSig.n_elem > 0)
            for(int j=0;j<out.n_rows;j++)for(int k=0;k<out.n_cols;k++) max_dif = max(max_dif, abs(out(j,k)-oldSig(j,k)));
        
        
        if (i > 1 && max_dif < tol) {
            if (verbose > 0) Rcpp::Rcout << "ggb_psd Converged after " << i << " iterations\n";
                break;
        }
        
        oldSig = out;
        
        // update over C: adjust eigenvalues
        arma::vec eigen_val;
        arma::mat eigen_vec;
        
        
        //Rcpp::Rcout << "|B|" << arma::norm(B - B.t(), "fro") << "\n";
        
        
        arma::eig_sym(eigen_val, eigen_vec, B - S);

        
        
        // am I allowed to ask for real?
        arma::mat diag_mat(eigen_val.n_elem, eigen_val.n_elem);
        for(int j = 0; j < eigen_val.n_elem; j++) diag_mat(j, j) = max(eigen_val(j) + delta, 0.f);
        
        C = eigen_vec * diag_mat * eigen_vec.t();
        
        // keep using same w without recomputing
        // w <- out$w
    }
    
    return S + C - B; // S + C - B
    // obj[l] <- out$obj + (sum((S - Sig[[l]])^2) - sum((R - Sig[[l]])^2))/2
}


static double objective(arma::mat S, arma::mat theta, double lambda, double penalty) {
    return -log(arma::det(theta)) + arma::trace(S * theta) + lambda * penalty;
}

static float g(arma::mat x, arma::mat S) {
    return -log(arma::det(x)) + arma::trace(S * x);
}

static arma::mat G_x(arma::mat x, arma::mat x_p, float t) {
    return (x - x_p) / t;
}

static arma::mat G_t(arma::mat E_k_1, arma::mat E_k, float t) {
    return (E_k_1 - E_k) / t;
}

// [[Rcpp::export]]
arma::mat ggb_mm(arma::mat S, arma::mat sig_til, arma::imat G, double t, double tol, double B, int maxiter, double lambda, int prox_maxiter=500, double prox_tol=1e-4, int verbose=1) {

    std::vector<arma::mat> sigma;
    sigma.push_back(S);
    
    float penalty;
    
    float init_t = t;
        
    while(sigma.size() < 2 || arma::norm(sigma.back() - sigma[sigma.size() - 2], "fro") / arma::norm(sigma[sigma.size() - 2]) > tol) {
        
        
        arma::mat sig_bar = prox(sigma.back() - t * Majorizor_linear_grad(sig_til, sigma.back(), S), G, lambda, penalty, prox_maxiter, prox_tol, verbose);
        
        arma::mat sig_k = sigma.back();
        
        t = init_t;
        
        while(true) {
            float mm_sig_bar = Majorizor_linear(sig_til, sig_bar, S);
            
            float mm_sig_k_1 = Majorizor_linear(sig_til, sig_k, S);
            
            float term_2 = -t * arma::trace(Majorizor_linear_grad(sig_til, sig_k, S).t() * G_t(sig_k, sig_bar, t));
            
            float term_3 = t/2 * pow(arma::norm(G_t(sig_k, sig_bar, t), "fro"), 2);
            
            if (mm_sig_bar <= mm_sig_k_1 + term_2 + term_3) 
            {
                sig_bar = prox(sig_k - t * Majorizor_linear_grad(sig_til, sig_k, S), G, lambda, penalty, prox_maxiter, prox_tol, verbose);
                Rcpp::Rcout << "=============================Shrinking t: " << t << "=============================\n";
                t *= B;
            } else {
                break;
            }
        }
        
        sigma.push_back(sig_bar);
    }
    
    return sigma.back();
}


// [[Rcpp::export]]
arma::mat iter_method(arma::mat S, arma::imat G, double t, double tol, double B, int maxiter, double lambda, arma::vec& objective_vec, int ggb_maxiter=500, double ggb_tol=1e-4, int verbose=1) {
    float orig_t = t;
    
    arma::mat theta_old = arma::diagmat(S.i());
    
    float penalty;
    
    arma::mat theta = prox(theta_old - t * (-theta_old.i() + S), G, lambda, penalty);
    
    double prev_obj = 1000;
    int i;
    for (i = 0; i < maxiter && arma::norm(theta - theta_old, "fro") > tol; ++i) {
        
        Rcpp::Rcout << "norm: " << arma::norm(theta - theta_old, "fro") << " > tol: " << tol << "\n";
        
        t = orig_t;
        theta_old = theta;
        
        // Shrink t
        // store inverse of theta for repetitive use, only write in terms of old theta
        arma::mat new_S = theta_old - t * (-theta_old.i() + S);
        
        arma::mat next_theta = prox(new_S, G, lambda * t, penalty, ggb_maxiter, ggb_tol, verbose);
        
        
        int b_maxiter = maxiter;
        int b_i = 0;
        // use theta old
        while(g(next_theta, S) > g(theta, S) - t * arma::trace((-arma::inv(theta) + S).t() * G_x(theta, next_theta, t)) + t / 2 * pow(arma::norm(G_x(theta, next_theta, t), "fro"), 2) && b_i < b_maxiter) {
            t *= B;
            new_S = theta_old - t * (-theta_old.i() + S);
            next_theta = prox(new_S, G, lambda * t, penalty, ggb_maxiter, ggb_tol, verbose);
            Rcpp::Rcout << "=============================Shrinking t: " << t << "=============================\n";
            b_i++;
        }
        
        if (b_i > b_maxiter) {
            Rcpp::Rcout << "Hit B_Maxiter";
            
        }
        
        theta = next_theta;
        
        objective_vec[i] = objective(S, theta, lambda, penalty);
        
        
        if (i > 0 && objective_vec[i] > objective_vec[i-1]) {
            Rcpp::Rcout << "Objective not decreasing on " << i << " iteration\n";
        }
    }
    
    //return theta_old - t * (-theta_old.i() + S);

        
    if (verbose)
        Rcpp::Rcout << "Iterative Method converged after " << i << " iterations\n";
    return theta;
}






