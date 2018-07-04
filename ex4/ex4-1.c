#include "../lib/cmatrix.h"
#include "../lib/dgemm.h"
#include "../lib/mersenne_twister.h"
#include <stdio.h>
#include <math.h>

/* void dgemv_(char *trans, int *m, int *n, 
double *alpha, double *A, int *ldA, double *x, int *incx, 
double *beta , double *y, int *incy); */

// mat_elem で簡単に行列要素を表現できるように、
// ベクトルについても同じように表現できるようにする
#define vec_elem(vec, i) (vec)[i]

void test(double ***target,int n) {
  double **mat_target = *target;
  mat_target = alloc_dmatrix(n, n);
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i) {
      mat_elem(mat_target, i, j) = fmin(i+1,j+1);
    }
  }
  *target = mat_target;
}

int main() {
  int n = 5;
  int seed = 114514;
  double **mat_target;
  double *vec_trial;
  int i, j;
  char trans = 'N';
  double alpha = 1.0;
  double beta = 0.0;

  /* generate target matrix */
  /* mat_target = alloc_dmatrix(n, n);
  for (j = 0; j < n; ++j) {
    for (i = 0; i < n; ++i) {
      mat_elem(mat_target, i, j) = fmin(i+1,j+1);
    }
  } */
  test(&mat_target, n);
  // output
  fprint_dmatrix(stdout, n, n, mat_target);

  // generate random vector
  vec_trial = alloc_dvector(n);
  for (i = 0; i < n; ++i) {
      printf("iii");
      vec_elem(vec_trial, i)= 1;
  }
  // output
  fprint_dvector(stdout, n, vec_trial);


  /* calculate C = A * B */
  /* matC = alloc_dmatrix(n, n);
  dgemm_(&trans, &trans, &n, &n, &n, &alpha, mat_ptr(matA), &n, mat_ptr(matB), &n,
         &beta, mat_ptr(matC), &n);

  printf("n = %d\n", n);
  printf("matC[0][0] = %15.10f\n", mat_elem(matC, 0, 0));
  printf("matC[0][n-1] = %15.10f\n", mat_elem(matC, 0, n-1));
  printf("matC[n-1][n-1] = %15.10f\n", mat_elem(matC, n-1, n-1)); */
  return 0;
}