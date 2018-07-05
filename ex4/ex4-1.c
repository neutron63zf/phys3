#include "../lib/cmatrix.h"
#include "../lib/dgemm.h"
#include "../lib/mersenne_twister.h"
#include <stdio.h>
#include <math.h>

// BLASのライブラリの関数を使うための宣言
// （これを入れなくても動くが、警告が消えてすっきりする）
void dgemv_(char *trans, int *m, int *n, 
double *alpha, double *A, int *ldA, double *x, int *incx, 
double *beta , double *y, int *incy);
double ddot_(int *n, double *x, int *incx, double *y, int *incy);
void dscal_(int *n, double *a, double *x, int *incx);

// mat_elem で簡単に行列要素を表現できるように、
// ベクトルについても同じように表現できるようにする
#define vec_elem(vec, i) (vec)[i]

// ターゲットの行列を作成
double **make_matrix(int n) {
  double **mat = alloc_dmatrix(n, n);
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i) {
      mat_elem(mat, i, j) = fmin(i+1,j+1);
    }
  }
  return mat;
}

// べき乗法に必要な最初のベクトルを作成
double *make_rand_vector(int n) {
  // シードで乱数を初期化
  int seed = 334;
  init_genrand(seed);
  double *vec = alloc_dvector(n);
  for (int i = 0; i < n ; ++i) {
    vec_elem(vec, i) = genrand_real3();
  }
  return vec;
}

// ベクトルに行列を作用させる関数
// BLAS Level 2関数である、dgemv_を使用
double *calc_Ax(int n, double **A, double *x) {
  char normal = 'N';
  double alpha = 1.0;
  double beta = 0.0;
  int inc = 1;

  double *y = alloc_dvector(n);
  // http://azalea.s35.xrea.com/blas/gemv.html
  // 適宜改行をいれてあげると読みやすい
  dgemv_(
    &normal, &n, &n,
    &alpha, mat_ptr(A), &n, vec_ptr(x), &inc,
    &beta, vec_ptr(y), &inc
  );
  return y;
}

// 2つのベクトルの内積を計算する関数
// BLAS Level 1関数である、ddot_を使用
double dot(int n, double *x, double *y) {
  int inc = 1;
  return ddot_(
    &n,
    vec_ptr(x), &inc,
    vec_ptr(y), &inc);
}

// ベクトルを定数倍する関数
// BLAS Level 1関数である、dscal_を使用
// もとのベクトルは書き換えられてしまう
double *calc_ax(int n, double a, double *x) {
  int inc = 1;
  dscal_(&n, &a, vec_ptr(x), &inc);
  return x;
}


int main() {
  int n = 5;
  double **mat_target = make_matrix(n);
  // fprint_dmatrix(stdout, n, n, mat_target);
  double *vec_trial = make_rand_vector(n);
  fprint_dvector(stdout, n, vec_trial);
  double *y = calc_Ax(n, mat_target, vec_trial);
  fprint_dvector(stdout, n, y);
  double _dot = dot(n, vec_trial, y);
  printf("%lf\n", _dot);
  return 0;
}