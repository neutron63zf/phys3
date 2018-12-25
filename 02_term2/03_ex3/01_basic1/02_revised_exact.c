#include "bond_table.h"
#include <stdio.h>
#include <math.h>

// コマンドラインから引数を抜き取る
void extract_args(
  // main の引数
  int argc, char** argv,
  // -> システムサイズ、温度
  int* L, double* T
) {
  if (argc < 2) {
    fprintf(stderr, "Usage: %s L T\n", argv[0]);
    exit(1);
  }
  *L=atoi(argv[1]);
  *T=atof(argv[2]);
}

double sumlog(double loga, double logb) {
  if (loga > logb) {
    return loga + log(1+exp(logb-loga));
  } else {
    return logb + log(1+exp(loga-logb));
  }
}

int main(int argc, char** argv) {
  // double temperature = 2.0;
  // int L = 4;
  double temperature;
  int L;
  extract_args(argc, argv, &L, &temperature);
  int num_sites = L * L;
  int num_bonds = 2 * L * L;
  int num_states = 1 << (L * L);
  int **bond;
  double w[2];
  double weight, sum, logsum;
  int i, j, b, s, si, sj;
  bond = alloc_imatrix(num_bonds, 2);
  init_bond_table(L, bond);
  w[0] = exp(1.0 / temperature);
  w[1] = exp(-1.0 / temperature);
  sum = 0;
  logsum = 0;
  for (s = 0; s < num_states; ++s) {
    weight = 1;
    for (b = 0; b < num_bonds; ++b) {
      i = mat_elem(bond, b, 0);
      j = mat_elem(bond, b, 1);
      si = (s >> i) & 1;
      sj = (s >> j) & 1;
      if (si == sj) {
        weight *= w[0];
      } else {
        weight *= w[1];
      }
    }
    logsum = sumlog(logsum, log(weight));
    sum += weight;
  }
  printf("temperature = %15.10f\n", temperature);
  printf("L = %d\n", L);
  printf("---- ordinary method ----\n");
  printf("Z = %15.10f\n", sum);
  printf("free energy = %15.10f\n", - temperature * log(sum));
  printf("free energy density = %15.10f\n", - temperature * log(sum) / num_sites);
  printf("---- log-sum-exp method ----\n");
  printf("Z = %15.10f\n", exp(logsum));
  printf("free energy = %15.10f\n", - temperature * logsum);
  printf("free energy density = %15.10f\n", - temperature * logsum / num_sites);
}
