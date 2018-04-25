#include <stdio.h>
#include <math.h>

// 繰り返し回数の上限
#define MAX_COUNT 25
// 最初の点
#define START_POINT 1

double input_fn (double x) {
    return tanh(x) + 0.2 * x + 0.3;
}

// これも数値微分するのはめんどい。というか求められてなさそう。
double diffed_fn (double x) {
    return 1 / (cosh(x) * cosh(x)) + 0.2;
}

int main () {
    // 最初のx
    double x = START_POINT;
    double f_now,df_now;
    f_now = input_fn(x);
    df_now = diffed_fn(x);
    printf("%lf",x / 2);
    return 0;
}