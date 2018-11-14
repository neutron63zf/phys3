#include "../../../00_lib/mersenne_twister.h"
#include "../../../00_lib/cmatrix.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define M_PI  3.14159265358979323846 /* pi */

/* uniform random number on (0,1) */
double ran () {
  double x = genrand_real3();
  return x;
}

// 指数分布に従う確率密度関数
// 逆関数法
double p_exp (double lambda) {
  double x = ran();
  double l = - lambda * log(1 - x);
  return l;
}

// sin(\theta)/2に従う確率密度関数
// 0 ~ +piだが、
// 実際は反転もありうるのでそれも考慮して-pi ~ 0も付け加える
// 棄却法
double p_sin () {
  double x,y;
  do {
    x = M_PI * ran();
    y = ran() / 2;
  } while (
    sin(x) / 2 <= y
  );
  double sign = -1.0;
  if (ran() > 1.0 / 2.0) {
    sign = 1.0;
  }
  return x * sign;
}

// 吸収されたかいなか
// 0 or 1
int p_absorbtion (double p_c) {
  double x = ran();
  int ret = 0;
  if (x < p_c) {
    ret = 1;
  }
  return ret;
}

// 次の位置、角度を計算
// 位置、角度は書き換えられる
void calc_next_position (
  // 平均自由行程
  double lambda,
  // <- mutate -> 位置、角度
  double* z, double* theta
) {
  double l = p_exp(lambda);
  *theta += p_sin();
  if (*theta > M_PI) {
    *theta -= 2.0 * M_PI;
  } else if (*theta < -M_PI) {
    *theta += 2.0 * M_PI;
  }
  // printf("# l: %lf, theta: %lf\n", l, *theta);
  // printf("# z: %lf -> ", *z);
  *z += l * cos(*theta);
  // printf("%lf\n", *z);
}

// 一回あたりの試行
void neutron_board (
  // 吸収率、平均自由行程、板の厚さ
  double p_c, double lambda, double D,
  // -> 反射率、吸収率、透過率（0or1だがdoubleで）
  double* reflection, double* absorption, double* through
) {
  // 吸収されたか
  int is_absorbed = 0;
  // 位置、角度
  double z = 0;
  double theta = 0;
  do {
    calc_next_position(
      lambda,
      &z, &theta
    );
    // 吸収されたかを判定
    is_absorbed = p_absorbtion(p_c);
    if (is_absorbed) {
      break;
    }
  } while (
    0 <= z && z <= D
  );

  *reflection = 0.0;
  *absorption = 0.0;
  *through    = 0.0;

  if (is_absorbed) {
    *absorption = 1.0;
  } else if (z < 0) {
    *reflection = 1.0;
  } else if (z > D) {
    *through = 1.0;
  }
}

// 分散を計算
double nb_sigma (double avg, int M) {
  double variance = avg * ( 1.0 - avg ) / M;
  if (variance <= 0) {
    return 0.0;
  }
  return sqrt(
    avg * ( 1.0 - avg ) / M
  );
}

// M回の試行
void montecarlo_neutron_board (
  // 吸収率、平均自由行程、板の厚さ、繰り返し回数
  double p_c, double lambda, double D, int M,
  // -> 反射率、吸収率、透過率
  double* reflection, double* absorption, double* through,
  // -> 反射率の誤差、吸収率の誤差、透過率の誤差
  double* err_ref, double* err_abs, double* err_thr
) {
  *reflection = 0.0;
  *absorption = 0.0;
  *through    = 0.0;
  double r, a, t;
  // #pragma omp parallel for reduction(+:_r,_a,_t)
  for (int i = 0; i < M; ++i) {
    neutron_board (
      p_c, lambda, D,
      &r, &a, &t
    );
    *reflection += r / M;
    *absorption += a / M;
    *through += t / M;
  }
  *err_ref = nb_sigma(*reflection , M);
  *err_abs = nb_sigma(*absorption , M);
  *err_thr = nb_sigma(*through    , M);
}

// 板の厚さを変えながら試行を繰り返す
// gnuplot用出力もする
void change_d_montecarlo_nb (
  // -> 吸収率、平均自由行程、繰り返し回数
  double p_c, double lambda, int M,
  // -> 板の厚さスタート、板の厚さエンド、板の厚さ刻み数
  double D_start, double D_end, int D_step_num
) {
  // 出力ヘッダー
  printf("# row description\n");
  printf("# D ref err_ref abs err_abs thr err_thr\n\n");

  // 反射率、吸収率、透過率
  double reflection, absorption, through;
  // 反射率の誤差、吸収率の誤差、透過率の誤差
  double err_ref, err_abs, err_thr;

  double step = ( D_end - D_start ) / D_step_num;
  double D;
  // #pragma omp parallel for
  for (int i = 0; i <= D_step_num; ++i) {
    // i = 0          -> D = D_start
    // i = D_step_num -> D_end
    D = D_start + step * i;
    montecarlo_neutron_board (
      p_c, lambda, D, M,
      &reflection, &absorption, &through,
      &err_ref, &err_abs, &err_thr
    );
    printf(
      "%lf %lf %lf %lf %lf %lf %lf\n",
      D, reflection, err_ref, absorption, err_abs, through, err_thr
    );
  }
}

// コマンドラインから引数を抜き取る
void extract_args(
  // mainの引数
  int argc, char** argv,
  // -> 吸収率、平均自由行程、繰り返し回数
  double* p_c, double* lambda, int* M,
  // -> 板の厚さスタート、板の厚さエンド、板の厚さ刻み数
  double* D_start, double* D_end, int* D_step_num
) {
  if (argc < 6) {
    fprintf(stderr, "Usage: %s p_c lambda M D_start D_end D_step_num\n", argv[0]);
    exit(1);
  }
  *p_c = atof(argv[1]);
  *lambda = atof(argv[2]);
  *M = atoi(argv[3]);
  *D_start = atof(argv[4]);
  *D_end = atof(argv[5]);
  *D_step_num = atoi(argv[6]);
}

// 実行エントリー関数
int main(int argc, char** argv) {
  // シードで乱数を初期化
  int seed = time(NULL);
  init_genrand(seed);

  // 吸収率、平均自由行程
  double p_c, lambda;
  // 繰り返し回数
  int M;
  // 板の厚さスタート、板の厚さエンド
  double D_start, D_end;
  // 板の厚さ刻み数
  int D_step_num;
  // 引数からデータ抽出
  extract_args(
    argc, argv,
    &p_c, &lambda, &M,
    &D_start, &D_end, &D_step_num
  );

  // 本体
  printf("\n");
  change_d_montecarlo_nb(
    p_c, lambda, M,
    D_start, D_end, D_step_num
  );

  // メタデータ出力
  printf("\n");
  printf("# simulation metadata\n");
  printf("# p_c: %lf\n", p_c);
  printf("# lambda: %lf\n", lambda);
  printf("# M: %d\n", M);
  printf("# D: %lf -> %lf (divide by %d)\n", D_start, D_end, D_step_num);
  printf("\n");

  return 0;
}