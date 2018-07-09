#include "../lib/cmatrix.h"
#include "../lib/dgemm.h"
#include "../lib/mersenne_twister.h"
#include <stdio.h>
#include <math.h>

// mat_elem で簡単に行列要素を表現できるように、
// ベクトルについても同じように表現できるようにする
#define vec_elem(vec, i) (vec)[i]

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
// 与える行列は、もとの行列をLU分解したものである必要がある
// また、ピボットベクトルも、LU分解の際に生成されたものである必要がある
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

// 行列式を返す関数
// LU分解されていることが前提で、ipivはその際に生成されたピボットベクトル
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

// ファンデルモンド行列を返し、理論的な行列式も出す
double **make_vandermonde(int n, double* det) {
  int seed = 5153;
  init_genrand(seed);
  double **mat = alloc_dmatrix(n, n);
  double base;
  for (int j = 0; j < n; ++j) {
    base = genrand_real3();
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

int main(int argc, char** argv) {

  char* filename = parse_filename(argc, argv);

  // データを読み込む
  double **mat;
  double *vec;
  int n = read_data(filename, &mat, &vec);

  printf("Matrix A:\n");
  fprint_dmatrix(stdout, n, n, mat);
  printf("Vector B (transposed):\n");
  fprint_dvector(stdout, n, vec);

  // LU分解を行う
  int *ipiv = lu_decomp(n, &mat);
  
  printf("Result of LU decomposition:\n");
  fprint_dmatrix(stdout, n, n, mat);
  printf("Pivot for LU decomposition:\n");
  fprint_ivector(stdout, n, ipiv);

  // det mat
  double det = detA(n, mat, ipiv);
  printf("det A: %lf\n", det);

  // 方程式を解く
  solve(n, &mat, &vec, ipiv);
  
  printf("Solution X (transposed):\n");
  fprint_dvector(stdout, n, vec);

  double vandet;
  double **van = make_vandermonde(n, &vandet);
  fprint_dmatrix(stdout, n, n, van);
  printf("det Vand(n): %lf\n", vandet);
  
  // メモリ解放
  free_dmatrix(mat);
  free_dvector(vec);
  free_ivector(ipiv);
}
