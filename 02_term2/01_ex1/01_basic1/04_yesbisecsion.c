#include <stdio.h>
#include <math.h>

double f(double x) { return pow(x, 3) - 3 * pow(x, 2) + x + 3; }

int main() {
  double a, b, x;
  double delta = -1.0e-1;
  double epsilon = 1.0e-8;
  a = 1;
  b = a + delta;
  printf("# Find initial enclosure\n");
  while (f(a) * f(b) > 0) {
    printf("# f(%le) = %le, f(%le) = %le\n", a, f(a), b, f(b));
    b += delta;
  }
  printf("# Start bisection search\n");
  while (fabs(b - a) > epsilon) {
    printf("# f(%le) = %le, f(%le) = %le\n", a, f(a), b, f(b));
    x = (a + b) / 2;
    if (fabs(f(x)) < epsilon) break;
    if (f(a) * f(x) < 0) {
      b = x;
    } else {
      a = x;
    }
    printf("%le\n", x);
  }
  printf("# Solution = %le\n", x);
  printf("%le\n", x);
}
