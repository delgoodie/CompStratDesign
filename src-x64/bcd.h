

void allocate_groups(int *dd, int p, int *M, int ****elem, int ***nelem, int ***ngelem);

void free_groups(int p, int *M, int ****elem, int ***nelem, int ***ngelem);

void allocate_v(int *M, int ***elem, int **nelem, int p, double ****v);

void free_v(int *M, int ***elem, int p, double ****v);

void compute_lamlist(double *rr, int p, int *M, int ***elem, int **nelem, int **ngelem, int nlam, double flmin, double *lamlist);

void compute_sumvv(double ***v, int p, int *M, int ***elem, int **nelem, int **ngelem, double *sumvv);

double objective(double ***v, double lam, int p, int *M, int ***elem, int **nelem, int **ngelem, double *rr0);

void print_matrix(double *mat, int nr, int nc);

void print_lowertri(double *mm, int nr, int nc);

void prox_bcd3(double *rr, int *M, double lam, int p, int ***elem, int **nelem, int **ngelem, double ***v, double *rr0, int maxiter, double tol, int verbose);

void pathwiseprox_bcd3(double *rr, int *dd, int *M, double *lambda, double *flmin,int *p, double *sumvv, double *obj, int *maxiter, double *tol, int *verbose);