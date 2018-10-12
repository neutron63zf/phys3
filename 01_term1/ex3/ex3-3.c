#include "../lib/cmatrix.h"
#include "../lib/mersenne_twister.h"
#include <stdio.h>
#include <math.h>

// かなり遅いので注意!

// BLAS宣言
void dgemm_(char *transa, char *transb, int *m, int *n, int *k,
double *alpha, double *A, int *ldA, double *B, int *ldB,
double *beta , double *C, int *ldC);
void dgemv_(char *trans, int *m, int *n, 
double *alpha, double *A, int *ldA, double *x, int *incx, 
double *beta , double *y, int *incy);
void daxpy_(int *n, double *alpha, double *x, int *incx, double *y, int *incy);
void dcopy_(int *n, double *x, int *incx, double *y, int *incy);
double dnrm2_(int *n, double *x, int *incx);

# define M_PI  3.14159265358979323846 /* pi */

// mat_elem で簡単に行列要素を表現できるように、
// ベクトルについても同じように表現できるようにする
#define vec_elem(vec, i) (vec)[i]

// ラプラス方程式を解く行列を構成するときに便利
// 行列成分を平らにしたときに
#define flatten(d, i, j) ( d - 1 ) * i + j


// 行列の積ABを計算し、Cに
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

// ラプラス方程式の係数行列の成分を計算する
// 分割数、境界条件によらない
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

// ラプラス方程式の境界とかそういうものの値の計算
// 分割数、境界条件に依存
double calc_lpB(int d, int i, int j) {
  double v = 0;
  double y = (j + 1.0) / d;
  if (i + 1 == d - 1) {
    // x正側が境界
    v = v - sin(M_PI * y );
  } else if (i == 0) {
    // x負側が境界
    v = v - sin(M_PI * 2 * y);
  }
  return v;
}

// 分割数だけからヤコビ法に用いる行列Gを構成する。
// ラプラス方程式の情報は前提としているので、引数には分割数のみが入る
double **jacobi_mtx(int d) {
  // 行列Gを構成するためには、A.Bの情報が必要なので、
  // まずは2つ用意し、同時に初期化する
  // 行の長さは (d-1)^2 + 1
  // 最後の余分な一次元を用いて、漸化式を簡略化する
  int l = d - 1;
  int n = pow(l, 2) + 1;
  // 対角成分がほぼ空っぽ
  double** A = alloc_dmatrix(n, n);
  // 対角行列
  double** D = alloc_dmatrix(n, n);
  for (int i = 0; i < l; ++i){
    for (int j = 0; j < l; ++j){
      // 現在の行インデックスを計算
      int ridx = flatten(d, i, j);
      // Aのベクトル部分（符号反転）
      double v = calc_lpB(d, i, j);
      // 符号反転
      mat_elem(A, ridx, n - 1) = -v;
      mat_elem(D, ridx, n - 1) = 0; 
      for (int _i = 0; _i < l; ++_i){
        for (int _j = 0; _j < l; ++_j){
          // 現在の列インデックスを計算
          int cidx = flatten(d, _i, _j);
          // Aの係数行列部分
          double m = calc_lpA(i, j, _i, _j);
          // この時点で成分が判明しているが...
          // 対角成分だったら別の箇所に代入（それも符号絶対値反転）
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
  // この時点で一番下の行がうめられていないので、埋める
  for (int k = 0; k < n - 1; ++k){
    mat_elem(A, n - 1, k) = 0;
    mat_elem(D, n - 1, k) = 0;
  }
  // 右下は1
  mat_elem(A, n - 1, n - 1) = 1;
  mat_elem(D, n - 1, n - 1) = 1;
  // この時点でやっと行列が2つできるので、あとはDAを計算してやるだけ
  double** G = alloc_dmatrix(n, n);
  calc_AB(n, D, A, &G);

  // fprint_dmatrix(stdout, n, n, A);
  // fprint_dmatrix(stdout, n, n, D);

  // 解放してあげる
  free_dmatrix(D);
  free_dmatrix(A);
  return G;
}

// splot用の表示
// dはベクトルの次元ではなく、分割数
// コードはcmatrixを参考に作成
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

// ランダムなベクトルXを作成する
double* make_randv(int d) {
  // 長さを計算
  int l = d - 1;
  int n = pow(l, 2) + 1;
  double *v = alloc_dvector(n);
  for (int i = 0; i < n; ++i) {
    // (-1, 1)
    vec_elem(v, i) = 2.0 * genrand_real3() - 1;
  }
  vec_elem(v, n - 1) = 1;
  return v;
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

// 2つのベクトルの差分を取り、|x-y|/|x|を計算
double diff_ratio(int n, double *x, double *y) {
  // xをコピー（daxpy_で破壊されると困るため）
  double *_x = alloc_dvector(n);
  int inc = 1;
  dcopy_(&n, vec_ptr(x), &inc, vec_ptr(_x), &inc);
  // _x - yを計算
  // _xは破壊される
  double alpha = -1;
  daxpy_(&n, &alpha, vec_ptr(y), &inc, vec_ptr(_x), &inc);
  // L2ノルムの比を計算
  double x_norm2 = dnrm2_(&n, vec_ptr(x), &inc);
  double delta_norm2 = dnrm2_(&n, vec_ptr(_x), &inc);
  return delta_norm2;
}

// jacobi法本体
double *jacobi(int d, double threshold) {
  printf("# random vector generating...\n");
  double *x = make_randv(d);
  printf("# jacobi_mtx generating...\n");
  double **G = jacobi_mtx(d);
  // double **bG = G;
  printf("# start calculation\n");
  int n = pow(d - 1, 2) + 1;
  double limit = pow(10, 10);
  int i = 0;
  double *y;
  double diff;
  // double **_G = alloc_dmatrix(n,n);
  do {
    y = calc_Ax(n, G, x);
    // 比率の表示
    printf("# %lf \n", diff_ratio(n, x, y));
    diff = diff_ratio(n, x, y);
    // y = calc_Ax(n, G, x);
    // GをG^2に変換することで収束までの速度を上げる
    /*
    calc_AB(n,G,G,&_G);
    G = _G;
    //*/
    x = y;
    i++;
  } while(
    diff > threshold &&
    // 無限ループされると困るので
    i < limit
  );
  if (i >= limit) {
    printf("# limit exceeded");
  }
  return x;
}

// コンソール確認用の表示
// 2次元上にちゃんと表示される
// メッシュ数が多いと死ぬ
// あくまでも確認で、そのままgnuplot用データにされてもいいように、#をつける
void print2d(int d, double* v) {
  int k = 0;
  for (int i = 0; i < d + 1; ++i) {
    printf("# ");
    for (int j = 0; j < d + 1; ++j){
      double p = 0;
      if (j != 0 && j != d) {
        if (i == 0) {
          p = sin(2 * M_PI * j / d);
        } else if (i == d) {
          p = sin(M_PI * j / d);
        } else {
          // printf("%d %d %d %d\n",k, i- 1,j-1,flatten(d, i - 1, j - 1));
          p = vec_elem(v, k);
          k++;
        }
      }
      printf("%10.5f ", p);
    }
    printf("\n");
  }
}

void main() {
  int d = 4;
  int seed = 1287;
  double threshold = 0.01;

  init_genrand(seed);
  // double *x = jacobi(d, threshold);
  double **A = jacobi_mtx(d);
  int n = pow(d-1,2)+1;
  fprint_dmatrix(stdout, n, n, A);
  double **B;
  double **C;
  printf("#lp\n");
  B = alloc_dmatrix(n,n);
  calc_AB(n,A,A,&B);
  A = B;
  fprint_dmatrix(stdout, n, n, A);
  C = alloc_dmatrix(n,n);
  calc_AB(n,A,A,&C);
  A = C;
  fprint_dmatrix(stdout, n, n, A);

  printf("#lp\n");
  B = alloc_dmatrix(n,n);
  calc_AB(n,A,A,&B);
  A = B;
  fprint_dmatrix(stdout, n, n, A);
  C = alloc_dmatrix(n,n);
  calc_AB(n,A,A,&C);
  A = C;
  fprint_dmatrix(stdout, n, n, A);
  // print2d(d, x);
  // print_splot_dmatrix(d, x);
}