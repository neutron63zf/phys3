#include "../lib/cmatrix.h"
#include "../lib/dgemm.h"
#include "../lib/mersenne_twister.h"
#include <stdio.h>
#include <math.h>
void dgemv_(char *trans, int *m, int *n, 
double *alpha, double *A, int *ldA, double *x, int *incx, 
double *beta , double *y, int *incy);
double ddot_(int *n, double *x, int *incx, double *y, int *incy);
void dscal_(int *n, double *a, double *x, int *incx);
#define vec_elem(vec, i) (vec)[i]
double **make_matrix(int n) {
  double **mat = alloc_dmatrix(n, n);
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i) {
      mat_elem(mat, i, j) = fmin(i+1,j+1);
    }
  }
  return mat;
}
double *make_rand_vector(int n) {
  int seed = 1287;
  init_genrand(seed);
  double *vec = alloc_dvector(n);
  for (int i = 0; i < n ; ++i) {
    vec_elem(vec, i) = genrand_real3();
  }
  return vec;
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
double dot(int n, double *x, double *y) {
  int inc = 1;
  return ddot_(
    &n,
    vec_ptr(x), &inc,
    vec_ptr(y), &inc);
}
void *calc_ax(int n, double a, double **x) {
  int inc = 1;
  dscal_(&n, &a, vec_ptr(*x), &inc);
}
double lambda(int n) {
  return (1.0 / 2.0)/(
    1 - cos(M_PI / (2.0 * n + 1))
  );
}
double calc_lamdba(int n, double *v_np1, double *v_n) {
  return
    dot(n, v_np1, v_np1) /
    dot(n, v_np1, v_n);
}
int main() {
  int n = 100;
  double t = pow(2, -32);
  double **A = make_matrix(n);
  double *v_n = make_rand_vector(n);
  double *v_np1;
  double vlambda;
  double climax = lambda(n);
  int i = 1;
  do {
    v_np1 = calc_Ax(n, A, v_n);
    vlambda = calc_lamdba(n, v_np1, v_n);
    printf("%d %lf\n", i, fabs(vlambda - climax));
    v_n = v_np1;
    calc_ax(n, 1.0 / vlambda, &v_n);
    i++;
  } while (
    fabs(vlambda - climax) >= t &&
    i < 30
  );
  return 0;
}