#include <stdio.h>
int main() {
    int sum = 0;
    int i;
    for (i = 1; i <= 100; ++i) {
        sum = sum + i;
    }
    printf("sum of integers from 1 to 100 is %d\n", sum);
    return 0;
}