#include <stdio.h>
#include <math.h>

int main () {
    // とりあえず倍精度にしてみる
    // ターゲット
    double x = 0.3 * M_PI;
    // 基準点の値を入れる
    double f_base = sin(x);
    // 精度の限界
    int MAX_COUNT = 25;
    double f_diff,delta;
    double h = 1;
    // csvファイルに出力するためにカンマ区切りにする
    printf("step, delta\n");
    for(int i = 0; i < MAX_COUNT; i++){
        // 変化先を見る
        f_diff = sin(x+h);
        // 微分の定義式
        delta = (f_diff - f_base) / h;
        printf("%lf, %lf\n", h, delta);
        h = h / 2;
    }
    return 0;
}