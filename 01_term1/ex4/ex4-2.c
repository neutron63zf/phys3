#include "../lib/cmatrix.h"
#include "../lib/dgemm.h"
#include "../lib/mersenne_twister.h"
#include <stdio.h>
#include <math.h>

# define M_PI  3.14159265358979323846 /* pi */

// mat_elem で簡単に行列要素を表現できるように、
// ベクトルについても同じように表現できるようにする
#define vec_elem(vec, i) (vec)[i]

/* http://www.netlib.org/lapack/explore-html/d3/d6a/dgetrf_8f.html */
extern void dgetrf_(int *M, int *N, double *A, int *LDA, int*IPIV, int *INFO);

/* http://www.netlib.org/lapack/explore-html/d6/d49/dgetrs_8f.html */
extern void dgetrs_(char *TRANS, int *N, int *NRHS, double *A, int *LDA, int *IPIV,
                    double *B, int *LDB, int *INFO);

// BLAS宣言
void dgemm_(char *transa, char *transb, int *m, int *n, int *k,
double *alpha, double *A, int *ldA, double *B, int *ldB,
double *beta , double *C, int *ldC);
void dgemv_(char *trans, int *m, int *n, 
double *alpha, double *A, int *ldA, double *x, int *incx, 
double *beta , double *y, int *incy);

// argvとargcからファイル名を弾き出す
char *parse_filename (int argc, char** argv) {
  if (argc < 2) {
    fprintf(stderr, "Usage: %s inputfile\n", argv[0]);
    exit(1);
  }
  return argv[1];
}

// ファイル名を受け取って、xとyのベクトルを格納する
// 次元を返す
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
    // fprint_dmatrix(stdout, mdim, dim, mat);
    // fprint_dvector(stdout, mdim, *x);
    // fprint_dvector(stdout, mdim, *y);
    // printf("%d %d\n", mdim, dim);
    return mdim;
}

// dimが近似の次数。nが個数。xが入力
double **make_design_matrix(int dim, int n, double* x) {
  double **dm = alloc_dmatrix(n, dim + 1);
  for (int r = 0; r < n; r++) {
    for (int c = 0; c < dim + 1; c++) {
      mat_elem(dm, r, c) = pow(x[r], c);
    }
  }
  return dm;
}

// dimが近似の次数。nが個数。yが出力
double *make_right(int dim, int n, double* y, double** dm) {
  char normal = 'T';
  double alpha = 1.0;
  double beta = 0.0;
  int inc = 1;
  int d = dim + 1;

  double *ry = alloc_dvector(n);
  // http://azalea.s35.xrea.com/blas/gemv.html
  // 適宜改行をいれてあげると読みやすい
  dgemv_(
    &normal, &n, &d,
    &alpha, mat_ptr(dm), &n, vec_ptr(y), &inc,
    &beta, vec_ptr(ry), &inc
  );
  return ry;
}

// 係数行列を作成する。
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

/* http://www.netlib.org/lapack/explore-html/d3/d6a/dgetrf_8f.html */
extern void dgetrf_(int *M, int *N, double *A, int *LDA, int*IPIV, int *INFO);

/* http://www.netlib.org/lapack/explore-html/d6/d49/dgetrs_8f.html */
extern void dgetrs_(char *TRANS, int *N, int *NRHS, double *A, int *LDA, int *IPIV,
                    double *B, int *LDB, int *INFO);

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
// ベクトルには答えが格納される（多分）
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

// 上の２つの関数を組み合わせただけの簡単なコード
// もとの行列とベクトルは書き換えられてしまうが、ベクトルには答えが格納される
// ピボットベクトルを返す（必要ないはず）
int *solve_by_lu_decomp(int dim, double*** mat, double** vec) {
    int* ipiv = lu_decomp(dim, mat);
    solve(dim, mat, vec, ipiv);
}

 
int main(int argc, char** argv) {
  int dim = 2;
  int d = dim + 1;
  char* filename = parse_filename(argc, argv);
  double *x;
  double *y;
  int n = read_data(filename, &x, &y);
  double **dm = make_design_matrix(dim, n, x);
  double *r = make_right(dim, n, y, dm);
  double **L = make_left(dim, n, dm);
  solve_by_lu_decomp(d, &L, &r);
  // fprint_dvector(stdout, d, r);

  // pltファイル出力部
  for (int i = 0; i < d; ++i) {
    printf("(%lf) * (x ** %d) + ", r[i], i);
  }
  printf("0\n");
  return 0;
}