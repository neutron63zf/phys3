#include <stdio.h>
#include <math.h>

// 繰り返し回数の上限
#define MAX_COUNT 25
// 最初の点
#define START_POINT 1

double input_fn (double x) {
    return pow(x, 3) - 3 * pow(x, 2) + x + 3;
}

// これも数値微分するのはめんどい。というか求められてなさそう。
double diffed_fn (double x) {
    return 3 * pow(x, 2) - 6 * x + 1;
}

double next_xn (double x) {
    return x - input_fn(x) / diffed_fn(x);
}

int main () {
    // 最初のx
    double x = START_POINT;
    for(int i = 0; i < MAX_COUNT; i++){
        x = next_xn(x);
        printf("%lf\n",x);
    }
    return 0;
}