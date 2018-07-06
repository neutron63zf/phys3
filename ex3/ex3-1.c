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

// argvとargcからファイル名を弾き出す
char *parse_filename (int argc, char** argv) {
  if (argc < 2) {
    fprintf(stderr, "Usage: %s inputfile\n", argv[0]);
    exit(1);
  }
  return argv[1];
}

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

// LU分解をして、ピボットベクトルを返す。もとの行列は書き換えられてしまう
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

// ベクトルと行列を与え、方程式を解く。もとのベクトルは書き換えられてしまう。
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

int main(int argc, char** argv) {

  char* filename = parse_filename(argc, argv);

  // データを読み込む
  double **mat;
  double *vec;
  int n = read_data(filename, &mat, &vec);

  /*
  printf("Matrix A:\n");
  fprint_dmatrix(stdout, n, n, mat);
  printf("Vector B (transposed):\n");
  fprint_dvector(stdout, n, vec);
  */

  // LU分解を行う
  int *ipiv = lu_decomp(n, &mat);
  
  /*
  printf("Result of LU decomposition:\n");
  fprint_dmatrix(stdout, n, n, mat);
  printf("Pivot for LU decomposition:\n");
  fprint_ivector(stdout, n, ipiv);
  */

  // 方程式を解く
  solve(n, &mat, &vec, ipiv);

  /*
  printf("Solution X (transposed):\n");
  fprint_dvector(stdout, n, vec);
  */
  
  free_dmatrix(mat);
  free_dvector(vec);
  free_ivector(ipiv);
}
