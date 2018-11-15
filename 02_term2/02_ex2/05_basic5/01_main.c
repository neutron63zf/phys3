#include "../../../00_lib/mersenne_twister.h"
#include "../../../00_lib/cmatrix.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define M_PI 3.14159265358979323846
#define vec_elem(vec, i) (vec)[i]

/** 入力系 **/

  // 引数を抜き取る
  // コマンドラインから引数を抜き取る
  void extract_args(
    // mainの引数
    int argc, char** argv,
    // -> システムサイズ、総ステップ数、最初の無視する値
    int* L, int* MS, int* ignore, double* T, int* mode
  ) {
    if (argc < 5) {
      fprintf(stderr, "Usage: %s L MS ignore T mode\n", argv[0]);
      exit(1);
    }
    *L=atoi(argv[1]);
    *MS=atoi(argv[2]);
    *ignore=atoi(argv[3]);
    *T=atof(argv[4]);
    *mode=atoi(argv[5]);
  }

//

/** 出力系 **/

  // スピンの±を表示
  void print_spin(int* spin, int L) {
    int num_sites = L * L;
    printf("# ");
    for (int i = 0; i < num_sites; ++i) {
      if (vec_elem(spin, i) > 0) {
        printf("+");
      } else {
        printf("-");
      }
      printf(" ");
    }
    printf("\n");
  }

  // モンテカルロステップについての出力
  void print_sequence(
    int sample_steps, double* E, double* M
  ) {
    double e, m2;
    for (int i = 0; i < sample_steps; ++i) {
      e = vec_elem(E, i);
      m2 = pow(vec_elem(M, i), 2);
      printf("%d %lf %lf\n", i, e, m2);    
    }
  }

  // 種々の物理量を表示
  void print_values (
    double T,
    double C, double C_s,
    double M2, double M2_s,
    double E, double E_s
  ) {
    printf(
      "%lf %lf %lf %lf %lf %lf %lf\n",
      T, C, C_s, M2, M2_s, E, E_s
    );
  }

//

/** 乱数系 **/

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

//

/** 格子形状系 **/

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

//

/** 初期値系 **/

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

//

/** メトロポリス法系 **/

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
        delta += -2.0 * J * spin[s] * spin[v];
      }
      if (ran() < exp(-beta * delta)) {
        spin[s] = -1 * spin[s];
        *m += (2.0 * spin[s]) / num_sites; 
        *e += delta;
        // printf("# delta: %lf\n", delta);
        // print_spin(spin,L);
      }
    }
    
  }

  // メトロポリス法に従ってモンテカルロ
  void metropolis_montecarlo(
    double T, int L, int* spin, int** lattice,
    int MS, int ignore,
    double** E, double** M
  ) {

    // エネルギー、磁化
    double e, m;

    calc_em(L, spin, lattice, &e, &m);
    
    // エネルギー、エネルギーの自乗、磁化の自乗
    for (int i = 0; i < MS; ++i) {
      // 最初のignoreステップが過ぎたら記録をしていく
      if (i >= ignore) {
        int eq_index = i - ignore;
        vec_elem(*E, eq_index) = e;
        vec_elem(*M, eq_index) = m;
      }
      each_montecarlo_step(
        T, L, spin, lattice,
        &e, &m
      );
      // printf("\r#%d", i);
    }
    // printf("\n");
  }

//

/** 計算結果要約系 **/

  // 物理量について自己相関時間を計算する
  double calc_tau (
    int sample_steps, double* A,
    double A_ave, double A2_ave
  ) {
    double A_ave2 = pow(A_ave, 2);
    double denomi = A2_ave - A_ave2;
    double tau = 0.0;
    for (int t = 0; t < sample_steps; ++t) {
      // t = 0,1,2...sample_steps,
      int i_number = sample_steps - t;
      double C = 0.0;
      for (int i = 0; i < i_number; ++i) {
        C += (1.0 / i_number) * (vec_elem(A, i + t) * vec_elem(A, i) - A_ave2) / denomi;
      }
      if (0 < C && C < 1) {
        tau += (-1.0 / log(C)) * t / sample_steps;
      }
      // printf("\r# t: %d, C: %lf, tau: %lf", t, C, tau);
    }
    // printf("\n");
    return tau;
  }

  // サンプルについて量を自乗した配列を作る。
  double* calc_sq (int sample_steps, double* E) {
    double* E2 = alloc_dvector(sample_steps);
    for (int i = 0; i < sample_steps; ++i) {
      vec_elem(E2, i) = pow(vec_elem(E, i), 2);
    }
    return E2;
  }

  // サンプルの平均
  double calc_ave (int sample_steps, double* E) {
    double E_ave = 0;
    for (int i = 0; i < sample_steps; ++i) {
      E_ave += vec_elem(E, i) / sample_steps;
    }
    return E_ave;
  }

  // サンプルの熱容量
  double calc_C (int N, double T, double E_ave, double E2_ave) {
    return 1.0 / (N * pow(T, 2)) * (E2_ave - pow(E_ave, 2));
  }

  // サンプルの熱容量の誤差
  double calc_C_s(int N, double T, double E_ave, double E_ave_s, double E2_ave_s) {
    return 1.0 / (N * pow(T, 2)) * sqrt(pow(E2_ave_s, 2) + 2 * pow(E_ave * E_ave_s, 2));
  }

  // 計算結果要約
  void calc_summary(
    int sample_steps, double* E_arr, double* M_arr, int N, double T,
    double* C, double* C_s, double* M2, double* M2_s, double* E, double* E_s
  ) {
    // サンプルについて計算
      *E = calc_ave(sample_steps, E_arr);

      double* E2_arr = calc_sq(sample_steps, E_arr);
      double E2 = calc_ave(sample_steps, E2_arr);
      *E_s = sqrt(E2 - pow(*E, 2));

      double* M2_arr = calc_sq(sample_steps, M_arr);
      *M2 = calc_ave(sample_steps, M2_arr);

      double* M4_arr = calc_sq(sample_steps, M2_arr);
      double M4 = calc_ave(sample_steps, M4_arr);
      *M2_s = sqrt(M4 - pow(*M2, 2));

      *C = calc_C(N, T, *E, E2);
      double* E4_arr = calc_sq(sample_steps, E2_arr);
      double E4 = calc_ave(sample_steps, E4_arr);
      double E2_s = sqrt(E4 - pow(E2, 2));
      *C_s = calc_C_s(N, T, *E, *E_s, E2_s);

      /*

      printf("# E:%lf\n",*E);
      printf("# E_s:%lf\n",*E_s);
      printf("# M2:%lf\n",*M2);
      printf("# M2_s:%lf\n",*M2_s);
      printf("# C:%lf\n",*C);
      printf("# C_s:%lf\n",*C_s);

      */
    //

    // 自己相関時間と係数
      double tau = calc_tau(
        sample_steps, E_arr,
        *E, E2
      );
      // printf("# tau:%lf\n",tau);
      double scale = sqrt(1 + 2 * tau);
      // printf("# scale:%lf\n",scale);
    //

    // 再推定
      *E_s = *E_s / scale;
      *M2_s = *M2_s / scale;
      *C_s = *C_s / scale;
    //
  }

//

// 全体の処理
int main(int argc, char** argv) {
  // プログラム初期化
    init_rand();

    int L, MS, ignore, mode;
    double T;
    extract_args(
      argc, argv,
      &L, &MS, &ignore, &T, &mode
    );
    int N = L * L;
  //

  // スピン系初期化
    // スピン初期配置
    int *spin = ising_random_grid(L);
    // イジングモデル用の結合行列
    int **lattice = square_lattice(L);
    // 計算に用いるステップ数
    int sample_steps = MS - ignore;
    // エネルギー、磁場の格納用の配列
    double* E = alloc_dvector(sample_steps);
    double* M = alloc_dvector(sample_steps);
  //

  if (mode == 1) {
    // 温度自動
    double T_start = 2.0;
    double T_end = 3.0;
    int n = 20;
    double T_step = (T_end - T_start) / n;
    // 計算結果の要約
    double C, C_s, M2, M2_s, E_ave, E_ave_s;
    for (int i = 0; i <= n; ++i) {
      int *spin = ising_random_grid(L);
      T = T_start + T_step * i;
      // モンテカルロ法
      metropolis_montecarlo(
        T, L, spin, lattice,
        MS, ignore,
        &E, &M
      );
      calc_summary(
        sample_steps, E, M, N, T,
        &C, &C_s, &M2, &M2_s, &E_ave, &E_ave_s
      );
      // 数々の物理量を表示
      print_values(T, C, C_s, M2, M2_s, E_ave, E_ave_s);
    }
  } else {
    // 温度確定
    // モンテカルロ法
    metropolis_montecarlo(
      T, L, spin, lattice,
      MS, ignore,
      &E, &M
    );
    // モンテカルロ過程を表示
    print_sequence(sample_steps, E, M);
  }
}