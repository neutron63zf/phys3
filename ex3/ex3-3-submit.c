#include "../lib/cmatrix.h"
#include "../lib/mersenne_twister.h"
#include <stdio.h>
#include <math.h>

void dgemm_(char *transa, char *transb, int *m, int *n, int *k,
double *alpha, double *A, int *ldA, double *B, int *ldB,
double *beta , double *C, int *ldC);
void dgemv_(char *trans, int *m, int *n, 
double *alpha, double *A, int *ldA, double *x, int *incx, 
double *beta , double *y, int *incy);
void daxpy_(int *n, double *alpha, double *x, int *incx, double *y, int *incy);
void dcopy_(int *n, double *x, int *incx, double *y, int *incy);
double dnrm2_(int *n, double *x, int *incx);

# define M_PI  3.14159265358979323846

#define vec_elem(vec, i) (vec)[i]

#define flatten(d, i, j) ( d - 1 ) * i + j

void calc_AB(int n, double **A, double **B, double ***C) {
  char normal = 'N';
  double alpha = 1.0;
  double beta = 0.0;

  dgemm_(
    &normal, &normal, &n, &n, &n,
    &alpha, mat_ptr(A), &n, mat_ptr(B), &n,
    &beta , mat_ptr(*C), &n
  );
}

double calc_lpA(int i, int j, int _i, int _j) {
  double m = 0;
  if (
    (_i == i) &&
    (_j == j)
  ){
    m = -4;
  }else if(
    (
      (_i == i + 1) &&
      (_j == j)
    ) ||
    (
      (_i == i - 1) &&
      (_j == j)
    ) ||
    (
      (_i == i) &&
      (_j == j + 1)
    ) ||
    (
      (_i == i) &&
      (_j == j - 1)
    )
  ){
    m = 1;
  }
  return m;
}

double calc_lpB(int d, int i, int j) {
  double v = 0;
  double y = (j + 1.0) / d;
  if (i + 1 == d - 1) {
    v = v - sin(M_PI * y );
  } else if (i == 0) {
    v = v - sin(M_PI * 2 * y);
  }
  return v;
}

double **jacobi_mtx(int d) {
  int l = d - 1;
  int n = pow(l, 2) + 1;
  double** A = alloc_dmatrix(n, n);
  double** D = alloc_dmatrix(n, n);
  for (int i = 0; i < l; ++i){
    for (int j = 0; j < l; ++j){
      int ridx = flatten(d, i, j);
      double v = calc_lpB(d, i, j);
      mat_elem(A, ridx, n - 1) = -v;
      mat_elem(D, ridx, n - 1) = 0; 
      for (int _i = 0; _i < l; ++_i){
        for (int _j = 0; _j < l; ++_j){
          int cidx = flatten(d, _i, _j);
          double m = calc_lpA(i, j, _i, _j);
          if (ridx == cidx) {
            mat_elem(A, ridx, cidx) = 0;
            mat_elem(D, ridx, cidx) = - 1.0 / m;
          } else {
            mat_elem(A, ridx, cidx) = m;
            mat_elem(D, ridx, cidx) = 0;
          }
        }
      }
    }
  }
  for (int k = 0; k < n - 1; ++k){
    mat_elem(A, n - 1, k) = 0;
    mat_elem(D, n - 1, k) = 0;
  }
  mat_elem(A, n - 1, n - 1) = 1;
  mat_elem(D, n - 1, n - 1) = 1;
  double** G = alloc_dmatrix(n, n);
  calc_AB(n, D, A, &G);
  free_dmatrix(D);
  free_dmatrix(A);
  return G;
}

void print_splot_dmatrix(int d, double* v) {
  printf("# splot data started\n");
  int k = 0;
  for (int i = 0; i < d + 1; ++i) {
    for (int j = 0; j < d + 1; ++j){
      double x = (i + 0.0) / d;
      double y = (j + 0.0) / d;
      double z = 0;
      if (j != 0 && j != d) {
        if (i == 0) {
          z = sin(2 * M_PI * j / d);
        } else if (i == d) {
          z = sin(M_PI * j / d);
        } else {
          z = vec_elem(v, k);
          k++;
        }
      }
      printf("%10.5f ", x);
      printf("%10.5f ", y);
      printf("%10.5f ", z);
      printf("\n");
    }
  }
}

double* make_randv(int d) {
  int l = d - 1;
  int n = pow(l, 2) + 1;
  double *v = alloc_dvector(n);
  for (int i = 0; i < n; ++i) {
    vec_elem(v, i) = 2.0 * genrand_real3() - 1;
  }
  vec_elem(v, n - 1) = 1;
  return v;
}

double *calc_Ax(int n, double **A, double *x) {
  char normal = 'N';
  double alpha = 1.0;
  double beta = 0.0;
  int inc = 1;

  double *y = alloc_dvector(n);
  dgemv_(
    &normal, &n, &n,
    &alpha, mat_ptr(A), &n, vec_ptr(x), &inc,
    &beta, vec_ptr(y), &inc
  );
  return y;
}

double diff_ratio(int n, double *x, double *y) {
  double *_x = alloc_dvector(n);
  int inc = 1;
  dcopy_(&n, vec_ptr(x), &inc, vec_ptr(_x), &inc);
  double alpha = -1;
  daxpy_(&n, &alpha, vec_ptr(y), &inc, vec_ptr(_x), &inc);
  double x_norm2 = dnrm2_(&n, vec_ptr(x), &inc);
  double delta_norm2 = dnrm2_(&n, vec_ptr(_x), &inc);
  return delta_norm2;
}

double *jacobi(int d, double threshold) {
  printf("# random vector generating...\n");
  double *x = make_randv(d);
  printf("# jacobi_mtx generating...\n");
  double **G = jacobi_mtx(d);
  printf("# start calculation\n");
  int n = pow(d - 1, 2) + 1;
  double limit = pow(10, 10);
  int i = 0;
  double *y;
  double diff;
  do {
    y = calc_Ax(n, G, x);
    printf("# %lf \n", diff_ratio(n, x, y));
    diff = diff_ratio(n, x, y);
    x = y;
    i++;
  } while(
    diff > threshold &&
    i < limit
  );
  if (i >= limit) {
    printf("# limit exceeded");
  }
  return x;
}

void main() {
  int d = 40;
  int seed = 1287;
  double threshold = 0.01;

  init_genrand(seed);
  double *x = jacobi(d, threshold);
  // print_splot_dmatrix(d, x);
}