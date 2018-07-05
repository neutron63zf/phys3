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
void *calc_ax(int n, double a, double **x) {
  int inc = 1;
  dscal_(&n, &a, vec_ptr(*x), &inc);
}

// 最大固有値の理論値を計算する関数
// ちなみに、オーダーはn^2
double lambda(int n) {
  return (1.0 / 2.0)/(
    1 - cos(M_PI / (2.0 * n + 1))
  );
}

// v_nとv_{n+1}から暫定的に固有値を計算する関数
double calc_lamdba(int n, double *v_np1, double *v_n) {
  return
    dot(n, v_np1, v_np1) /
    dot(n, v_np1, v_n);
}

 
int main() {
  // 次元の設定
  // 今の所113くらいが限界
  int n = 100;
  // 収束判定に用いる。理論値との差がここまで縮んだらループを止める
  double t = pow(2, -32);

  // 初期値
  double **A = make_matrix(n);
  double *v_n = make_rand_vector(n);
  double *v_np1;
  double vlambda;

  // 収束先の値
  double climax = lambda(n);

  // メインループ
  int i = 1;
  do {
    // v_{n+1}を計算
    v_np1 = calc_Ax(n, A, v_n);
    // lambdaの暫定値を計算
    vlambda = calc_lamdba(n, v_np1, v_n);
    // 現在のlambdaを出力
    printf("%d %lf\n", i, vlambda);

    // v_nにv_{n+1}を代入
    v_n = v_np1;
    // 計算のたびにv_nの絶対値はlambdaだけ大きくなるので、それを打ち消す
    calc_ax(n, 1.0 / vlambda, &v_n);
    // iを加算
    i++;
  } while (
    // 計算結果の差がt以上の場合は計算を続ける
    fabs(vlambda - climax) >= t
  );
  // ついでに収束先も表示してあげる
  printf("∞  %lf\n", climax);

  return 0;
}