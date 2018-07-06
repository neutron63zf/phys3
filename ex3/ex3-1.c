#include "../lib/cmatrix.h"
#include "../lib/dgemm.h"
#include "../lib/mersenne_twister.h"
#include <stdio.h>
#include <math.h>

/* http://www.netlib.org/lapack/explore-html/d3/d6a/dgetrf_8f.html */
extern void dgetrf_(int *M, int *N, double *A, int *LDA, int*IPIV, int *INFO);

/* http://www.netlib.org/lapack/explore-html/d6/d49/dgetrs_8f.html */
extern void dgetrs_(char *TRANS, int *N, int *NRHS, double *A, int *LDA, int *IPIV,
                    double *B, int *LDB, int *INFO);

// ファイル名を受け取って、ベクトルと行列を読み取り、次元を返す
int read_data (char* filename, double*** mat, double** vec) {
    FILE *fp;
    int mdim, dim;
    fp = fopen(filename, "r");
    if (fp == NULL) {
        fprintf(stderr, "Error: file can not open\n");
        exit(1);
    }
    read_dmatrix(fp, &mdim, &dim, mat);
    if (mdim != dim) {
        fprintf(stderr, "Error: inconsistent number of equations\n");
        exit(1);
    }
    read_dvector(fp, &dim, vec);
    if (mdim != dim) {
        fprintf(stderr, "Error: inconsistent number of equations\n");
        exit(1);
    }
    return dim;
}

int main(int argc, char** argv) {
  char* filename;


  int *ipiv;
  int info;
  char trans = 'N';
  int nrhs = 1;

  if (argc < 2) {
    fprintf(stderr, "Usage: %s inputfile\n", argv[0]);
    exit(1);
  }
  filename = argv[1];

  double **mat;
  double *vec;
  int n = read_data(filename, &mat, &vec);

  printf("Matrix A:\n");
  fprint_dmatrix(stdout, n, n, mat);
  printf("Vector B (transposed):\n");
  fprint_dvector(stdout, n, vec);

  /* perform LU decomposition */
  ipiv = alloc_ivector(n);
  dgetrf_(&n, &n, mat_ptr(mat), &n, vec_ptr(ipiv), &info);
  if (info != 0) {
    fprintf(stderr, "Error: LAPACK::dgetrf failed\n");
    exit(1);
  }
  printf("Result of LU decomposition:\n");
  fprint_dmatrix(stdout, n, n, mat);
  printf("Pivot for LU decomposition:\n");
  fprint_ivector(stdout, n, ipiv);

  /* solve equations */
  dgetrs_(&trans, &n, &nrhs, mat_ptr(mat), &n, vec_ptr(ipiv), vec_ptr(vec), &n, &info);
  if (info != 0) {
    fprintf(stderr, "Error: LAPACK::dgetrs failed\n");
    exit(1);
  }
  printf("Solution X (transposed):\n");
  fprint_dvector(stdout, n, vec);
  
  free_dmatrix(mat);
  free_dvector(vec);
  free_ivector(ipiv);
}
