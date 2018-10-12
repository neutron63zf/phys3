#include "../lib/cmatrix.h"
#include "../lib/dgemm.h"
#include "../lib/mersenne_twister.h"
#include <stdio.h>
#include <math.h>

# define M_PI  3.14159265358979323846

#define vec_elem(vec, i) (vec)[i]

extern void dgetrf_(int *M, int *N, double *A, int *LDA, int*IPIV, int *INFO);
extern void dgetrs_(char *TRANS, int *N, int *NRHS, double *A, int *LDA, int *IPIV,
                    double *B, int *LDB, int *INFO);

void dgemm_(char *transa, char *transb, int *m, int *n, int *k,
double *alpha, double *A, int *ldA, double *B, int *ldB,
double *beta , double *C, int *ldC);
void dgemv_(char *trans, int *m, int *n, 
double *alpha, double *A, int *ldA, double *x, int *incx, 
double *beta , double *y, int *incy);

char *parse_filename (int argc, char** argv) {
  if (argc < 2) {
    fprintf(stderr, "Usage: %s inputfile\n", argv[0]);
    exit(1);
  }
  return argv[1];
}

int read_data (char* filename, double** x, double** y) {
    FILE *fp;
    int mdim, dim;
    double **mat;
    fp = fopen(filename, "r");
    if (fp == NULL) {
        fprintf(stderr, "Error: file can not open\n");
        exit(1);
    }
    read_dmatrix(fp, &mdim, &dim, &mat);
    *x = alloc_dvector(mdim);
    *y = alloc_dvector(mdim);
    for (int i = 0; i < mdim; ++i){
      vec_elem(*x, i) = mat_elem(mat, i, 0);
      vec_elem(*y, i) = mat_elem(mat, i, 1);
    }
    return mdim;
}

double **make_design_matrix(int dim, int n, double* x) {
  double **dm = alloc_dmatrix(n, dim + 1);
  for (int r = 0; r < n; r++) {
    for (int c = 0; c < dim + 1; c++) {
      mat_elem(dm, r, c) = pow(x[r], c);
    }
  }
  return dm;
}

double *make_right(int dim, int n, double* y, double** dm) {
  char normal = 'T';
  double alpha = 1.0;
  double beta = 0.0;
  int inc = 1;
  int d = dim + 1;

  double *ry = alloc_dvector(n);
  dgemv_(
    &normal, &n, &d,
    &alpha, mat_ptr(dm), &n, vec_ptr(y), &inc,
    &beta, vec_ptr(ry), &inc
  );
  return ry;
}

double **make_left(int dim, int n, double** dm) {
  char trans = 'T';
  char normal = 'N';
  double alpha = 1.0;
  double beta = 0.0;
  int d = dim + 1;
  double **C = alloc_dmatrix(d,d);

  dgemm_(
    &trans, &normal, &d, &d, &n,
    &alpha, mat_ptr(dm), &n, mat_ptr(dm), &n,
    &beta , mat_ptr(C), &d
  );
  return C;
}

int *lu_decomp (int dim, double*** mat) {
  char normal = 'N';
  int info;
  int *ipiv = alloc_ivector(dim);

  dgetrf_(&dim, &dim, mat_ptr(*mat), &dim, vec_ptr(ipiv), &info);
  if (info != 0) {
    fprintf(stderr, "Error: LAPACK::dgetrf failed\n");
    exit(1);
  }
  return ipiv;
}

void solve(int dim, double*** mat, double** vec, int *ipiv) {
  char normal = 'N';
  int info;
  int nrhs = 1;

  dgetrs_(&normal, &dim, &nrhs, mat_ptr(*mat), &dim, vec_ptr(ipiv), vec_ptr(*vec), &dim, &info);
  if (info != 0) {
    fprintf(stderr, "Error: LAPACK::dgetrs failed\n");
    exit(1);
  }

}

int *solve_by_lu_decomp(int dim, double*** mat, double** vec) {
    int* ipiv = lu_decomp(dim, mat);
    solve(dim, mat, vec, ipiv);
}

 
int main(int argc, char** argv) {
  int dim = 30;
  int d = dim + 1;
  char* filename = parse_filename(argc, argv);
  double *x;
  double *y;
  int n = read_data(filename, &x, &y);
  double **dm = make_design_matrix(dim, n, x);
  double *r = make_right(dim, n, y, dm);
  double **L = make_left(dim, n, dm);
  solve_by_lu_decomp(d, &L, &r);

  for (int i = 0; i < d; ++i) {
    printf("(%lf) * (x ** %d) + ", r[i], i);
  }
  printf("0\n");
  return 0;
}