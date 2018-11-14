#include "../../../00_lib/mersenne_twister.h"
#include "../../../00_lib/cmatrix.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define M_PI 3.14159265358979323846
#define vec_elem(vec, i) (vec)[i]

// 乱数を初期化
void init_rand () {
  // シードで乱数を初期化
  int seed = time(NULL);
  init_genrand(seed);
}

// (0,1)区間での一様乱数
double ran () {
  double x = genrand_real3();
  return x;
}

// 引数を抜き取る
// コマンドラインから引数を抜き取る
void extract_args(
  // mainの引数
  int argc, char** argv,
  // -> システムサイズ、総ステップ数、最初の無視する値
  int* L, int* MS, int* ignore, double* T
) {
  if (argc < 4) {
    fprintf(stderr, "Usage: %s L MS ignore T\n", argv[0]);
    exit(1);
  }
  *L=atoi(argv[1]);
  *MS=atoi(argv[2]);
  *ignore=atoi(argv[3]);
  *T=atof(argv[4]);
}

// イジングモデルのランダムな格子を作成する
int* ising_random_grid(int L) {
  int *grid = alloc_ivector(L * L);

  for (int i = 0; i < L * L; ++i) {
    int spin;
    if (ran() > 1.0 / 2.0) {
      spin = +1;
    } else {
      spin = -1;
    }
    vec_elem(grid, i) = spin;
  }
  
  return grid;
}

// イジングモデル用の結合行列
int** square_lattice(int L) {
  int **neighbor = alloc_imatrix(L * L, 4);
  /* size of neighbors should be (L*L) x 4 */
  int x, y, xn, yn, i, j;
  for (y = 0; y < L; ++y) {
    for (x = 0; x < L; ++x) {
      i = L * y + x;
      /* right */
      xn = (x + 1) % L;
      yn = y;
      j = L * yn + xn;
      mat_elem(neighbor, i, 0) = j;
      /* left */
      xn = (L + x - 1) % L;
      yn = y;
      j = L * yn + xn;
      mat_elem(neighbor, i, 1) = j;
      /* up */
      xn = x;
      yn = (y + 1) % L;
      j = L * yn + xn;
      mat_elem(neighbor, i, 2) = j;
      /* down */
      xn = x;
      yn = (L + y - 1) % L;
      j = L * yn + xn;
      mat_elem(neighbor, i, 3) = j;
    }
  }
  //fprint_imatrix(stdout, L*L, 4, neighbor);
  //printf("# %d\n", mat_elem(neighbor,L*L -1, 2));
  return neighbor;
}

// sのスピンに最近接しているj個目のsiteの番号
int neighbor(
  int s, int j, int** lattice
) {
  return mat_elem(lattice, s, j);
}

// イジングモデルの各スピンにおいて更新を試してみる
void each_montecarlo_step(
  double T, int L, int* spin,int** lattice,
  double* e, double* m
) {
  double J = -1.0;
  double beta = 1.0 / T;
  double delta = 0.0;

  int num_sites = L * L;
  int num_neighbors = 4;
  int s,j,v;
  for (s = 0; s < num_sites; ++s) {
    delta = 0.0;
    for (j = 0; j < num_neighbors; ++j) {
      v = neighbor(s, j, lattice);
      delta += 2 * J * spin[s] * spin[v];
    }
    if (ran() < exp(-beta * delta)) {
      spin[s] = -1 * spin[s];
      *m += (2.0 * spin[s]) / num_sites; 
      *e += delta;
    }
  }
  
}

// 最初のEとMを計算するのに使う
void calc_em(
  int L, int* spin, int** lattice,
  double* e, double* m
) {
  *e = 0.0;
  *m = 0.0;
  double J = -1.0;
  int num_sites = L * L;
  int num_neighbors = 4;
  int s,j,v;
  for (s = 0; s < num_sites; ++s) {
    for (j = 0; j < num_neighbors; ++j) {
      v = neighbor(s, j, lattice);
      *e += ((J / 2.0) * spin[s]) * spin[v];
    }
    *m += 1.0 * spin[s] / num_sites;
  }
}

// メイン処理
// メトロポリス法に従ってモンテカルロ
void metropolis_montecarlo(
  double T, int L, int* spin,
  int MS, int ignore
) {
  // 計算に用いるステップ数
  int sample_steps = MS - ignore;
  // エネルギー、磁場の格納用の配列
  double* E = alloc_dvector(sample_steps);
  double* M = alloc_dvector(sample_steps);

  // 上のやつの自乗
  double* E2 = alloc_dvector(sample_steps);
  double* M2 = alloc_dvector(sample_steps);

  // イジングモデル用の結合行列
  int **lattice = square_lattice(L);
  // エネルギー、磁化
  double e, m;
  // とその自乗
  double e2, m2;

  calc_em(L, spin, lattice, &e, &m);
  
  // エネルギー、エネルギーの自乗、磁化の自乗
  for (int i = 0; i < MS; ++i) {
    each_montecarlo_step(
      T, L, spin, lattice,
      &e, &m
    );
    e2 = pow(e, 2);
    m2 = pow(m, 2);
    fprint_ivector(stdout, L*L, spin);

    // 最初のignoreステップが過ぎたら記録をしていく
    if (i >= ignore) {
      int eq_index = i - ignore;
      vec_elem(E, eq_index) = e;
      vec_elem(M, eq_index) = m;
      vec_elem(E2, eq_index) = e2;
      vec_elem(M2, eq_index) = m2;
      printf("%d %lf %lf\n", eq_index, e, m2);
    }
  }
}

// 実行エントリーポイント
int main(int argc, char** argv) {
  init_rand();

  int L, MS, ignore;
  double T;
  extract_args(
    argc, argv,
    &L, &MS, &ignore, &T
  );

  int *spin = ising_random_grid(L);
  
  // fprint_ivector(stdout, L*L, spin);

  metropolis_montecarlo(
    T, L, spin,
    MS, ignore
  );
}