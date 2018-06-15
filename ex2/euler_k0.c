#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// 方程式に出現する各種定数
// 型を担保するため「.0」をつける
#define x0 10.0
#define v0 0.0
#define kappa 0.0 // kappaをゼロに
#define k 2.0
#define m 1.0

// 終了時間
#define targetSec 30.0

// 最終的にhを2^(-maxNPower)まで減らす
// #define maxNPower 10

/*
    基本方針：
    位置だけでなく速度に対応する変数も設定し、それぞれについてオイラー法を適用する。
*/

// 加速度を求める関数
double delta_v (double x, double v) {
    return ( -k * x - kappa * v ) / m;
}

// 現在のv,xを元に次のv,xを錬成する関数
// パラメーターhを必要とする
double next_vx (double *x, double *v, double h) {
    double a = delta_v(*x, *v);
    // 現在の速度のままhだけ進んだらどうなるか
    *x = *x + h * *v;
    // 速度を更新(現在の加速度のままhだけ進んだらどうなるか)
    *v = *v + h * a;
    return 0.0;
}

// コマンドライン引数による処理の分岐
int main (int argc, char *argv[]) {
    double t = 0;
    double h = pow(2, -atoi(argv[1]));
    double x = x0;
    double v = v0;
    for (int i = 0; i * h <= targetSec ; i++){
        printf("%lf %lf\n", t, 1.0/2.0*m*v*v+1.0/2.0*k*x*x);
        next_vx(&x, &v, h);
        t += h;
    }
    return 0;
}

