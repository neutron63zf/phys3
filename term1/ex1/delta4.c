#include <stdio.h>
#include <math.h>

double df(double x, double h) {
    double y1 = sin(x+h);
    double y2 = sin(x-h);
    double y3 = sin(x+2*h);
    double y4 = sin(x-2*h);
    return ( y4 - 8 * y2 + 8 * y1 - y3 ) / ( 12 * h );
}

int main () {
    // とりあえず倍精度にしてみる
    // ターゲット
    double x = 0.3 * M_PI;
    // 精度の限界
    int MAX_COUNT = 21;
    double f_aft,delta,f_bef;
    double h = 1;
    printf("#step delta\n");
    for(int i = 0; i < MAX_COUNT; i++){
        // 複雑になってきたので分離
        delta = df(x, h);
        printf("%lf %lf\n", h, fabs(delta-cos(x)));
        h = h / 2;
    }
    return 0;
}