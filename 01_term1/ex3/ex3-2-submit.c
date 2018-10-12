#include "../lib/cmatrix.h"
#include "../lib/dgemm.h"
#include "../lib/mersenne_twister.h"
#include <stdio.h>
#include <math.h>

# define M_PI  3.14159265358979323846 /* pi */
#define vec_elem(vec, i) (vec)[i]
#define flatten(d, i, j) ( d - 1 ) * i + j
extern void dgetrf_(int *M, int *N, double *A, int *LDA, int*IPIV, int *INFO);
extern void dgetrs_(char *TRANS, int *N, int *NRHS, double *A, int *LDA, int *IPIV,
                    double *B, int *LDB, int *INFO);

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

void make_laplace_eq(int d, double*** mat, double** vec) {
  int l = d - 1;
  *mat = alloc_dmatrix(pow(l,2),pow(l,2));
  *vec = alloc_dvector(pow(l,2));
  for (int i = 0; i < l; ++i){
    for (int j = 0; j < l; ++j) {
      int target_vidx = flatten(d,i,j);
      double v = 0;
      if (i + 1 == l) {
        v = v - sin(M_PI * (j + 1) / d );
      } else if (i == 0) {
        v = v - sin(M_PI * 2 * (j + 1) / d);
      }
      vec_elem(*vec, target_vidx) = v;
      int target_m_ridx = target_vidx;
      for (int _i = 0; _i < l; ++_i) {
        for (int _j = 0; _j < l; ++_j) {
          int target_m_cidx = flatten(d, _i, _j);
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
          mat_elem(*mat,target_m_ridx, target_m_cidx) = m;
        }
      }
    }
  }
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

void main() {
  int d = 100;
  int l = d - 1;
  int p = pow(l, 2);
  double **m;
  double *v;
  make_laplace_eq(d, &m, &v);
  solve_by_lu_decomp(p, &m, &v);
  // print_splot_dmatrix(d, v);
}