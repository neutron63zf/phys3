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

// 現在のv,xを元に次のv,xを錬成する関数
// パラメーターhを必要とする
// 今回は中点法なので少し複雑です。
double next_vx (double *x, double *v, double h) {
    // 現在時点での加速度
    double a = delta_v(*x, *v);
    // Euler法による、「h/2経過した」ときの状態
    double vc = *v + ( h / 2.0 ) * a;
    double xc = *x + ( h / 2.0 ) * *v;
    // その状態での加速度を求める
    double ac = delta_v(xc, vc);
    // 現在の速度のままhだけ進んだらどうなるか
    *x = *x + h * vc;
    // 速度を更新(現在の加速度のままhだけ進んだらどうなるか)
    *v = *v + h * ac;
    return 0.0;
}

// コマンドライン引数による処理の分岐
int main (int argc, char *argv[]) {
    double t = 0;
    double h = pow(2, -atoi(argv[1]));
    double x = x0;
    double v = v0;
    for (int i = 0; i * h <= targetSec ; i++){
        printf("%lf %lf %lf\n", t, x, v);
        next_vx(&x, &v, h);
        t += h;
    }
    return 0;
}

