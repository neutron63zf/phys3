#include "../../../00_lib/mersenne_twister.h"
#include "../../../00_lib/cmatrix.h"
#include <stdio.h>
#include <math.h>

/* uniform random number on (0,1) */
double ran() {
    double x = genrand_real3();
    return x;
}

double p (double x) {
    double y;
    if (0 < x && x <= 1.0/2.0) {
        y = 4.0 * x;
    } else if (1.0/2.0 < x && x < 1) {
        y = 4.0 * (1.0 - x);
    } else {
        y = 0;
    }
    return y;
}

int main (int argc, char** argv) {
    // 全体の設定値
    int seed = 12345;
    int samples = pow(2,24);
    int bins = pow(2,6);

    // 乱数をセット
    init_genrand(seed);

    // ヒストグラム用配列と各種数値
    double xmin = 0;
    double xmax = 1;
    double dx = (xmax - xmin) / bins;
    double dh = 1.0 / (samples * dx);
    double *hist = alloc_dvector(bins);
    for (int idx = 0; idx < bins; ++idx) hist[idx] = 0.0;

    // 棄却法に用いる箱の縦横の長さ、左下の座標
    double bh = 2;
    double bw = 1;
    double by = 0;
    double bx = 0;

    // 実行
    int count = 0;
    double x,y;
    double px;
    int idx;
    while (count <= samples) {
        x = bx + bw * ran();
        y = by + bh * ran();
        px = p(x);

        // 棄却
        if (y > px) continue;
        
        ++count;
        idx = (x - xmin) / dx;
        hist[idx] += 1.0;
    }

    double ave,err;
    for (int idx = 0; idx < bins; ++idx) {
        x = (idx + 0.5) * dx;
        ave = hist[idx] * dh;
        err = sqrt(hist[idx]) * dh;
        printf("%d %15.10f %15.10f %15.10f\n", idx, x, ave, err);
    }

    free_dvector(hist);
    return 0;
}