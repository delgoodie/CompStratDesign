// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// include RcppArmadillo.h which pulls Rcpp.h
#include "RcppArmadillo.h"

//#include "bcd.h"

extern "C" void pathwiseprox_bcd3(double *rr, int *dd, int *M, double *lambda, int *p, double *sumvv, double *obj, int *maxiter, double *tol, int *verbose);
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
arma::mat bcd_method(arma::mat S, arma::imat G, double lambda, int maxiter=500, double tol=1e-4, int verbose=1) {
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
    
    pathwiseprox_bcd3(ss_vec.memptr(), dd_vec.memptr(), M.memptr(), &lambda, &p, sumvv, &obj, &maxiter, &tol, &verbose);
    
    arma::vec sumvv_vec(n_sumvv);
    for (int i = 0; i < n_sumvv; i++) sumvv_vec(i) = sumvv[i];
    
    free(sumvv);
    
    arma::mat res = from_tri(sumvv_vec, S.n_rows);
    for (int i=0;i<res.n_rows;i++) res(i, i) = S(i, i);
    
    return res;
}


// [[Rcpp::export]]
arma::mat iter_method(arma::mat S, arma::imat G, double t, double tol, int maxiter, double lambda, int ggb_maxiter=500, double ggb_tol=1e-4, int verbose=1) {
    arma::mat theta_old = arma::diagmat(S.i());
    arma::mat theta = bcd_method(theta_old - t * (-theta_old.i() + S), G, lambda);
    int i=0;
    
    do {
        i++;
        theta_old = theta;
        arma::mat new_S = theta_old - t * (-theta_old.i() + S);
        theta = bcd_method(new_S, G, lambda * t, ggb_maxiter, ggb_tol, verbose);
    } while (i < maxiter && arma::norm(theta - theta_old, 2) > tol);
    
    if (verbose)
        Rcpp::Rcout << "Iterative Method converged after " << i << " iterations\n";
    return theta.i();
}


