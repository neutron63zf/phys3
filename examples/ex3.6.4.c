#include <stdio.h>

void division(int divident, int divisor, int *quotient, int *residual) {
    *quotient = divident / divisor;
    *residual = divident % divisor;
}

int main() {
    int josuu = 3;
    int hi_josuu = 13;
    int shou, amari;
    division(hi_josuu, josuu, &shou, &amari);
    printf("%d / %d = %d ... %d\n", hi_josuu, josuu, shou, amari);
    return 0;
}