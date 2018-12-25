#include "cmatrix.h"
#include "dgemm.h"
#include <stdio.h>
#include <math.h>

/* v = D^{1/2} v */
void product_D(int M, double temperature, double *v) {
  int dim = 1 << M; /* dim = 2^M */
  int i, j, s;
  double elem, w[2];
  w[0] = 1.0;
  w[1] = exp(-1.0 / temperature);
  for (s = 0; s < dim; ++s) {
    elem = v[s];
    for (i = 0; i < M; ++i) {
      j = (i+1) % M;
      elem *= w[((s >> i) & 1) ^ ((s >> j) & 1)];
    }
    v[s] = elem;
  }
}

/* v = U v */
void product_U(int M, double temperature, double *v, double *vtmp) {
  int dim = 1 << M; /* dim = 2^M */
  int i, s;
  double w[2];
  w[0] = 1.0;
  w[1] = exp(-2.0 / temperature);
  for (i = 0; i < M; ++i) {
    for (s = 0; s < dim; ++s)
      vtmp[s] = w[0] * v[s] + w[1] * v[s ^ (1 << i)];
    for (s = 0; s < dim; ++s)
      v[s] = vtmp[s];
  }
}

/* v = D^{1/2} U D^{1/2} v */
void product(int M, double temperature, double *v, double *vtmp) {
  product_D(M, temperature, v);
  product_U(M, temperature, v, vtmp);
  product_D(M, temperature, v);
}

void product_L(int M, int L, double temperature, double *v, double *vtmp) {
  for (int i = 0; i < L; ++i) {
    product(M, temperature, v, vtmp);
  }
}

// void fill_T(int M, double temperature, double **T) {
//   int dim = 1 << M; /* dim = 2^M */
//   int i, j;
//   double *v0, *v1;
//   v0 = alloc_dvector(dim);
//   v1 = alloc_dvector(dim);
//   for (i = 0; i < dim; ++i) {
//     for (j = 0; j < dim; ++j) v0[j] = 0.0;
//     v0[i] = 1.0;
//     product(M, temperature, v0, v1);
//     for (j = 0; j < dim; ++j) T[i][j] = v0[j];
//   }
//   free_dvector(v0);
//   free_dvector(v1);
// }
  
int main(int argc, char** argv) {
  double temperature = 2.0;
  int L = 4;
  int M = L;
  int num_sites = L * M;
  int dim = 1 << M;
  double **T, **T0, **T1;
  int i, j, k;
  double sum;
  char trans = 'N';
  double alpha = 1.0;
  double beta = 0.0;

  /* prepare transfer matrix T */
  // T = alloc_dmatrix(dim, dim);
  // fill_T(M, temperature, T);
  
  /* calculate T^L */
  // T0 = alloc_dmatrix(dim, dim);
  // T1 = alloc_dmatrix(dim, dim);
  // /* T0 = identity matrix */
  // for (i = 0; i < dim; ++i) {
  //   for (j = 0; j < dim; ++j) mat_elem(T0, i, j) = 0.0;
  //   mat_elem(T0, i, i) = 1.0;
  // }
  // for (k = 0; k < L; ++k) {
  //   /* T1 = T * T0 */
  //   dgemm_(&trans, &trans, &dim, &dim, &dim, &alpha, mat_ptr(T),
  //          &dim, mat_ptr(T0), &dim, &beta, mat_ptr(T1), &dim);
  //   /* T0 = T1 */
  //   for (j = 0; j < dim; ++j) {
  //     for (i = 0; i < dim; ++i) {
  //       mat_elem(T0, i, j) = mat_elem(T1, i, j);
  //     }
  //   }
  // }

  /* take trace of T^L */
  sum = 0.0;
  // for (i = 0; i < dim; ++i) sum += mat_elem(T0, i, i);
  for (i = 0; i < dim; ++i) {
    double *v0, *v1;
    v0 = alloc_dvector(dim);
    v1 = alloc_dvector(dim);
    for (j = 0; j < dim; ++j) {
      v0[j] = 0.0;
      v1[j] = 0.0;
    } 
    v0[i] = 1.0;
    v1[i] = 1.0;
    product_L(M, L, temperature, v0, v1);
    sum += v0[i];
  }
  
  printf("temperature = %15.10f\n", temperature);
  printf("L = %d\n", L);
  printf("Z = %15.10f * exp(%14.10f)\n", sum, 2 * num_sites / temperature);
  printf("free energy = %15.10f\n", -2 * num_sites - temperature * log(sum));
  printf("free energy density = %15.10f\n", -2 - temperature * log(sum) / num_sites);

  // free_dmatrix(T);
  // free_dmatrix(T0);
  // free_dmatrix(T1);
}
