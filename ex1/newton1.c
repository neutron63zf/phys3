#include <stdio.h>
#include <math.h>

// 繰り返し回数の上限
#define MAX_COUNT 25
// 最初の点
#define START_POINT -100

double input_fn (double x) {
    return tanh(x) + 0.2 * x + 0.3;
}

// これも数値微分するのはめんどい。というか求められてなさそう。
double diffed_fn (double x) {
    return 1 / (cosh(x) * cosh(x)) + 0.2;
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