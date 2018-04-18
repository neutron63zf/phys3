#include <stdio.h>
#include <math.h>
/* test sign of cosine values */
int main() {
    int i = 0;
    while (i <= 180) {
        double angle;
        double cosval;
        angle = i*M_PI/180.0;
        cosval = cos(angle);
        printf("cos(%lf) is %lf\n", angle, cosval);
        i = i + 20;
    }
    return 0;
}