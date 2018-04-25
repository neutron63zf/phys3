#include <stdio.h>
#include <math.h>

int main () {
    // とりあえず倍精度にしてみる
    // ターゲット
    double x = 0.3 * M_PI;
    // 基準点の値を入れる
    double f_base = sin(x);
    // 微分の刻み幅
    double h = 1;
    // 変化先を見る
    double f_diff = sin(x+h);
    // 微分
    double delta = (f_diff - f_base) / h;
    printf("step => %lf, delta => %lf\n", h, delta);
    return 0;
}