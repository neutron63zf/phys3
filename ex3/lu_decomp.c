// 
#include "../lib/cmatrix.h"
#include <stdio.h>

/* http://www.netlib.org/lapack/explore-html/d3/d6a/dgetrf_8f.html */
extern void dgetrf_(int *M, int *N, double *A, int *LDA, int*IPIV, int *INFO);

// 行列式を計算する関数
double det (int *n, double **mat) {
    int *ipiv;
    int info;

    // pivot用ベクトルにメモリを割り当てる
    ipiv = alloc_ivector(*n);

    fprint_dmatrix(stdout, *n, *n, mat);

    // LU分解を実行
    dgetrf_(n, n, mat_ptr(mat), n, vec_ptr(ipiv), &info);

    printf("%d", *n);
    // ipiv = alloc_ivector(n);
    // dgetrf_(&n, &n, mat_ptr(a), &n, vec_ptr(ipiv), &info);
    return 0.0;
}

// ファイルから行列を読み込む関数
void read_dmat_from_file (char* filename, double **mat, int *n) {
    FILE *fp;
    int m;

    fp = fopen(filename, "r");
    if (fp == NULL) {
        fprintf(stderr, "Error: file can not open\n");
        exit(1);
    }
    read_dmatrix(fp, &m, n, &mat);
    if (m != *n) {
        fprintf(stderr, "Error: inconsistent number of equations\n");
        exit(1);
    }
    // printf("Matrix (Double):\n");
    // fprint_dmatrix(stdout, *n, *n, mat);
}

int main (int argc, char** argv) {

    char* filename;
    int n;
    double **mat;

    // コマンドライン引数が足りなかったりおかしかったりする場合を弾く
    if (argc < 2) {
        fprintf(stderr, "Usage: %s inputfile\n", argv[0]);
        exit(1);
    }

    // コマンドライン変数からファイル名を取得
    filename = argv[1];
    printf("FILE NAME:\n\t%s \n\n", filename);

    read_dmat_from_file(filename, mat, &n);
    fprint_dmatrix(stdout, n, n, mat);
    // det(&n, &mat);
    return 0;
}
