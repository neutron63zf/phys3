#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// 方程式に出現する各種定数
// 型を担保するため「.0」をつける
#define x0 10.0
#define v0 0.0
#define kappa 0.0 // kappaをゼロにする。
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

// 速度を求める関数
// vなのだが、明示的にxの微分を得るために使用する
double delta_x (double x, double v) {
    return v;
}

// 現在のv,xを元に次のv,xを錬成する関数
// パラメーターhを必要とする
// 今回は中点法なので少し複雑です。
double next_vx (double *x, double *v, double h) {
    double k1 = h * delta_v(*x, *v);
    double l1 = h * delta_x(*x, *v);
    double k2 = h * delta_v(
        *x + 1.0 / 2.0 * k1,
        *v + 1.0 / 2.0 * l1
    );
    double l2 = h * delta_x(
        *x + 1.0 / 2.0 * k1,
        *v + 1.0 / 2.0 * l1
    );
    // 位置・速度を更新
    *x = *x + l2;
    *v = *v + k2;
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

