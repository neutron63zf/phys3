#include <stdio.h>
int main() {
    int a = 0; // 一つ前の
    int b = 1; // 最新の
    int c = -1; // スワップ用
    for (int i = 1; i <= 60; i++) {
        // iってのは「最新の」のインデックス
        c = b;
        b = a + b;
        a = c;
        printf("%d => %d\n",i,b);
        // 40くらいでオーバーフローする
    }
    return 0;
}