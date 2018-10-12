#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// 方程式に出現する各種定数
// 型を担保するため「.0」をつける
#define x0 10.0
#define v0 0.0
#define kappa 0.0 // kappaの値はゼロ
#define k 2.0
#define m 1.0

// 終了時間
#define targetSec 30.0

// 最終的にhを2^(-maxNPower)まで減らす
// #define maxNPower 10

/*
    基本方針：
    リープ・フロッグ法
*/


// リープ・フロッグ法の中に出現するポテンシャルの微分
double round_V (double q) {
    return k*q;
}

// 現在のq,pを元に次のq,pを錬成する関数
// パラメーターhを必要とする
// 今回はRK法なのでさらに複雑です。
// q, p双方にRK法を適用しつつ計算をしていく
double next_pq (double *q, double *p, double h) {
    double p_h_2 = *p - h / 2.0 * round_V(*q);
    *q = *q + h * p_h_2;
    *p = p_h_2 - h / 2.0 * round_V(*q);
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
        next_pq(&x, &v, h);
        t += h;
    }
    return 0;
}

