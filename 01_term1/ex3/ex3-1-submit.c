#include "../lib/cmatrix.h"
#include "../lib/dgemm.h"
#include "../lib/mersenne_twister.h"
#include <stdio.h>
#include <math.h>

#define vec_elem(vec, i) (vec)[i]

/* http://www.netlib.org/lapack/explore-html/d3/d6a/dgetrf_8f.html */
extern void dgetrf_(int *M, int *N, double *A, int *LDA, int*IPIV, int *INFO);

double detA(int dim, double** mat, int *ipiv) {
  double det = 1;
  for (int i = 0; i < dim; ++i) {
    det *= mat_elem(mat, i, i);
    if (vec_elem(ipiv, i) - 1 != i) {
      det *= -1;
    }
  }
  return det;
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

double **make_vandermonde(int n, double* det) {
  double **mat = alloc_dmatrix(n, n);
  double base;
  for (int j = 0; j < n; ++j) {
    base = 5 * genrand_real3();
    for (int i = 0; i < n; ++i) {
      mat_elem(mat, i, j) = pow(base, i);
    }
  }
  *det = 1;
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < j; ++i) {
      *det *= mat_elem(mat, 1, j) - mat_elem(mat, 1, i);
    }
  }
  return mat;
}

int main() {
  int base_seed = 5153;
  int n = 5;
  double vandet;
  for(int i = 0; i < 1000; i++){
    int seed = base_seed + i;
    init_genrand(seed);
    double **van = make_vandermonde(n, &vandet);
    int *ipiv = lu_decomp(n, &van);
    double det = detA(n, van, ipiv);
    printf("error : %lf\n", (det - vandet)/vandet);
  }
}
