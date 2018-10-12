#include "../lib/cmatrix.h"
#include "../lib/dgemm.h"
#include "../lib/mersenne_twister.h"
#include <stdio.h>
#include <math.h>

# define M_PI  3.14159265358979323846 /* pi */

// mat_elem で簡単に行列要素を表現できるように、
// ベクトルについても同じように表現できるようにする
#define vec_elem(vec, i) (vec)[i]

// ラプラス方程式を解く行列を構成するときに便利
// 行列成分を平らにしたときに
#define flatten(d, i, j) ( d - 1 ) * i + j

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

//ラプラス方程式を差分方程式にしたものの行列とベクトルを作成する。
// dは平面の分割数
void make_laplace_eq(int d, double*** mat, double** vec) {
  // 行列、ベクトル用領域を確保
  int l = d - 1;
  *mat = alloc_dmatrix(pow(l,2),pow(l,2));
  *vec = alloc_dvector(pow(l,2));
  // mat_elem(mat, i, j) を u((i+1)/d, (j+1)/d)としている
  for (int i = 0; i < l; ++i){
    for (int j = 0; j < l; ++j) {
      // まずはベクトル成分を計算
      int target_vidx = flatten(d,i,j);
      double v = 0;
      if (i + 1 == l) {
        // x正側が境界
        v = v - sin(M_PI * (j + 1) / d );
      } else if (i == 0) {
        // x負側が境界
        v = v - sin(M_PI * 2 * (j + 1) / d);
      }
      vec_elem(*vec, target_vidx) = v;
      // 次に行列の成分を計算
      // 行一定の状態で、列を変えていく
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

// splot用の表示
// dはベクトルの次元ではなく、分割数
// コードはcmatrixを参考に作成
void print_splot_dmatrix(int d, double* v) {
  printf("# splot data started\n");
  int k = 0;
  for (int i = 0; i < d + 1; ++i) {
    // printf("# ");
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

// コンソール確認用の表示
// 2次元上にちゃんと表示される
// メッシュ数が多いと死ぬ
// あくまでも確認で、そのままgnuplot用データにされてもいいように、
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
  // 分割数
  int d = 4;
  int l = d - 1;
  int p = pow(l, 2);
  double **m;
  double *v;
  make_laplace_eq(d, &m, &v);
  // fprint_dmatrix(stdout, p, p, m);
  // fprint_dvector(stdout, p, v);
  solve_by_lu_decomp(p, &m, &v);
  // fprint_dvector(stdout, p, v);
  print2d(d, v);
  // print_splot_dmatrix(d, v);
}