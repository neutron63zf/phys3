#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// 方程式に出現する各種定数
// 型を担保するため「.0」をつける
#define x0 10.0
#define v0 0.0
#define kappa 0.2
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
// vなのだが、RKを明示的にxに適用するために使用する
double delta_x (double x, double v) {
    return v;
}

// 現在のv,xを元に次のv,xを錬成する関数
// パラメーターhを必要とする
// 今回はRK法なのでさらに複雑です。
// x, v双方にRK法を適用しつつ計算をしていく
double next_vx (double *x, double *v, double h) {
    // v計算部分
    double k1 = h * delta_v(*x, *v);
    double l1 = h * delta_x(*x, *v);
    double k2 = h * delta_v(
        *x + l1 / 2.0,
        *v + k1 / 2.0
    );
    double l2 = h * delta_x(
        *x + l1 / 2.0,
        *v + k1 / 2.0
    );
    double k3 = h * delta_v(
        *x + l2 / 2.0,
        *v + k2 / 2.0        
    );
    double l3 = h * delta_x(
        *x + l2 / 2.0,
        *v + k2 / 2.0        
    );
    double k4 = h * delta_v(
        *x + l3,
        *v + k3
    );
    double l4 = h * delta_x(
        *x + l3,
        *v + k3
    );
    double va = *v + k1 / 6.0 + k2 / 3.0 + k3 / 3.0 + k4 / 6.0;
    double xa = *x + l1 / 6.0 + l2 / 3.0 + l3 / 3.0 + l4 / 6.0;

    // 速度・位置を更新
    *x = xa;
    *v = va;
    return 0.0;
}

// 真の値
double true_val (double t) {
    double gamma = kappa / (2.0 * m);
    double w = sqrt(k/m - pow(gamma,2));
    return x0 * exp(-gamma*t) * cos(w*t);
}

// コマンドライン引数による処理の分岐
int main (int argc, char *argv[]) {
    double t = 0;
    double h = pow(2, -atoi(argv[1]));
    double x = x0;
    double v = v0;
    for (int i = 0; i * h <= targetSec ; i++){
        // printf("%lf %lf %lf\n", t, x, v);
        printf("%lf %lf\n", t, fabs(
            x-true_val(t)
        ));
        next_vx(&x, &v, h);
        t += h;
    }
    return 0;
}

