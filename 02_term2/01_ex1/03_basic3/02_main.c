#include "../../../00_lib/cmatrix.h"
#include "../../../00_lib/mersenne_twister.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define M_PI  3.14159265358979323846 /* pi */
#define vec_elem(vec, i) (vec)[i]

// BLASのライブラリの関数を使うための宣言
// （これを入れなくても動くが、警告が消えてすっきりする）
double ddot_(int *n, double *x, int *incx, double *y, int *incy);
void dscal_(int *n, double *a, double *x, int *incx);

double potential(double x, double v, double width, double offset) {
  double value;
  if (x > (0.5 - 0.5 * width) && x < (0.5 + 0.5 * width)) {
    value = v;
  } else {
    value = 0.0;
  }
  value += offset;
  return value;
}

void extract_argv(
  int* n, double* v, double* width, double* offset, int* mode, double* tpow,
  int argc, char** argv
) {
  if (argc < 6) {
    fprintf(stderr, "Usage: %s n v width offset mode tpow\n", argv[0]);
    exit(1);
  }
  *n = atoi(argv[1]);
  *v = atof(argv[2]);
  *width = atof(argv[3]);
  *offset = atof(argv[4]);
  *mode = atoi(argv[5]);
  *tpow = atof(argv[6]);
}

double** prepare_mat(int dim, double v, double width, double offset) {
  int n = dim + 1;
  double h = 1.0 / n;
  double diag = 2.0 / (h*h);
  double offdiag = -1.0 / (h*h);
  
  /* preparation of Mat */
  double **mat = alloc_dmatrix(dim, dim);
  int i;
  double x;
  for (i = 0; i < dim; ++i) {
    x = 1.0 * (i+1) / n; /* x = 0 for i = -1 and x = 1 for i = dim */
    if (i > 0) {
      mat_elem(mat, i, i-1) = offdiag;
    }
    mat_elem(mat, i, i) = diag + potential(x, v, width, offset);
    if (i < (dim - 1)) {
      mat_elem(mat, i, i+1) = offdiag;
    }
  }

  return mat;
}

// Hに関する事前の知識を仮定する
/* x = 0 for i = -1 and x = 1 for i = dim */
// ハミルトニアン行列とベクトルの積を計算する
double* Hx(int n, double** H, double* phi) {
  int dim = n;
  double *new_phi = alloc_dvector(dim);
  double x;
  for (int i = 0; i < dim; ++i) {
    x = 1.0 * (i+1) / n;
    vec_elem(new_phi, i) = 0;
    if (i > 0) {
      vec_elem(new_phi, i) += vec_elem(phi, i - 1) * mat_elem(H, i, i - 1);
    }
    vec_elem(new_phi, i) += vec_elem(phi, i) * mat_elem(H, i, i);
    if (i < (dim - 1)) {
      vec_elem(new_phi, i) += vec_elem(phi, i + 1) * mat_elem(H, i, i + 1);
    }
  } 
  return new_phi;
}

// べき乗法に必要な最初のベクトルを作成
double *make_rand_vector(int n) {
  // シードで乱数を初期化
  int seed = 1287;
  init_genrand(seed);
  double *vec = alloc_dvector(n);
  for (int i = 0; i < n ; ++i) {
    vec_elem(vec, i) = genrand_real3();
  }
  return vec;
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

// v_nとv_{n+1}から暫定的に固有値を計算する関数
double calc_lamdba(int n, double *v_np1, double *v_n) {
  return
    dot(n, v_np1, v_np1) /
    dot(n, v_np1, v_n);
}

// べき乗法の本体
/**
 * t 閾値
 * A 行列
 * n 次元
 * eigen_vec 返り値、固有ベクトル
 * eigen_value 返り値、固有値
 */
void power_method (double t, double** A, int n, double** eigen_vec, double* eigen_value, double offset, int mode) {
  double *v_n = make_rand_vector(n);
  double *v_np1;
  double vlambda_n = 0;
  double vlambda_np1 = 0;

  // メインループ
  int i = 1;
  do {
    // \lambda_nに\lambda_{n+1}を代入
    vlambda_n = vlambda_np1;

    // v_{n+1}を計算
    v_np1 = Hx(n, A, v_n);
    // lambdaの暫定値を計算
    vlambda_np1 = calc_lamdba(n, v_np1, v_n);
    // 現在のlambdaを出力
    if (mode == 1) {
      printf("%d %lf\n", i, fabs(vlambda_np1 - vlambda_n));
    }

    // v_nにv_{n+1}を代入
    v_n = v_np1;
    // 計算のたびにv_nの絶対値はlambdaだけ大きくなるので、それを打ち消す
    calc_ax(n, 1.0 / vlambda_np1, &v_n);
    // iを加算
    i++;
  } while (
    // 計算結果の差がt以上の場合は計算を続ける
    fabs(vlambda_np1 - vlambda_n) >= t &&
    // i < 30
    1
  );
  *eigen_value = vlambda_np1 - offset;
  *eigen_vec = v_np1;
}

// ベクトルをgnuplotで表示できるような形式で出力する
void print_vector (int n, double* vec, double start, double end) {
  double x;
  for (int i = 0; i < n; ++i) {
    // i = 0     -> x = start
    // i = n - 1 -> x = end
    x = start + ( end - start ) / ( n - 1 ) * i;
    printf("%lf %lf\n", x, vec_elem(vec, i));
  }
}

// ベクトルを1に規格化
double* normalize_vector(int n, double* vec) {
  double l2_norm = sqrt(dot(n, vec, vec));
  double* new_vec = alloc_dvector(n);
  for (int i = 0; i < n; ++i) {
    vec_elem(new_vec, i) = vec_elem(vec, i) / l2_norm;
  }
  return new_vec;
}

int main(int argc, char** argv) {
  int n; /* partition number */
  double v; /* height of potential between two wells */
  double width; /* width of wall between two wells */
  double offset; /* potential offset */
  int mode; /* execution mode */
  double tpow;
  extract_argv(
    &n, &v, &width, &offset, &mode, &tpow,
    argc, argv
  );

  int dim = n - 1; /* dimension of Hamiltonian */
  
  printf("#------\n# n: %10d\n", n);
  printf("# v: %10.5f\n", v);
  printf("# width: %10.5f\n", width);
  if (mode == 1) {
    printf("# mode: 1 -> print convergence\n");
  }

  /* preparation of Mat */
  double **mat = prepare_mat(dim, v, width, offset);

  /* power method */
  double t = pow(2, tpow);  /* power method threshold */
  double *evec = alloc_dvector(dim); /* eigenvector */
  double eval; /* eigenvalue */
  power_method(t, mat, dim, &evec, &eval, offset, mode);

  /* normalize vector */
  evec = normalize_vector(dim, evec);

  /* print out eigenstate/eigenvalue */
  if (mode == 0) {
    printf("# mode: 0 -> print eigenstate/eigenvalue\n\n");
    double print_start = 1.0 / n;
    double print_end = 1.0 - print_start;
    printf("%lf %lf\n", 0.0, 0.0);
    print_vector(dim, evec, print_start, print_end);
    printf("%lf %lf\n", 1.0, 0.0);
    printf("\n# eigenvalue ->\n# %lf\n", eval);
    printf("\n# (for changing offet)\n# -E/2 ->\n# %lf\n", -eval / 2);
  }
}
