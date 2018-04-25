#include <stdio.h>
#include <math.h>

int main () {
    // とりあえず倍精度にしてみる
    // ターゲット
    double x = 0.3 * M_PI;
    // 精度の限界
    int MAX_COUNT = 25;
    double f_aft,delta,f_bef;
    double h = 1;
    // csvファイルに出力するためにカンマ区切りにする
    printf("step, delta\n");
    for(int i = 0; i < MAX_COUNT; i++){
        f_aft = sin(x+h);
        f_bef = sin(x-h);
        // 微分の定義式
        delta = (f_aft - f_bef) / (2 * h);
        printf("%lf, %lf\n", h, delta);
        h = h / 2;
    }
    return 0;
}