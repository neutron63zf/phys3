#include <math.h>
#include <stdio.h>

double circle_area(double r) {
    return r*r*M_PI;
}

int main() {
    double radius = 2.0;
    double area = circle_area(radius);
    printf("Radius: %lf, Area: %lf\n", radius, area);
    return 0;
}