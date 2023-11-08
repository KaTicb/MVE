#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#define h 0.2f
#define r 0.3f
#define NXB 15
#define NYB 12
#define lambda 1.7f
#define REP 3000
#define EPS 1.e-4
#define T1 20.0f
#define T2 20.0f

int main(int argc, char **argv, char *env[]) {

  // DEFINITION
  int NX = 3 * NXB + 1;
  int NY = 3 * NYB + 1;
  int N = NY * NX;

  double T[NY][NX];
  double max_rel_diff = 0.0;

  int i1 = NXB, i2 = NXB + NXB, i3 = NXB + NXB + NXB, j1 = NYB, j2 = NYB + NYB,
      j3 = NYB + NYB + NYB;

  int i, j, m = NX + 1, k = 0;
  // DEFINITION

  // EMPTY SPACE
  for (j = 0; j <= j3; ++j)
    for (i = 0; i <= i3; ++i) {
      T[j][i] = 0;
    }
  // EMPTY SPACE

  double alf_1 = -h / r;
  double alf_2 = -r / h;
  double alf_3 = alf_2 * 0.5f;
  double alf_4 = alf_1 * 0.5f;
  double gam_1 = -2.f * (alf_1 + alf_2);
  double gam_2 = -1.5f * (alf_1 + alf_2);
  double gam_3 = -(alf_1 + alf_2);
  double gam_4 = -(alf_3 + alf_4);

  // CONTACTS
  for (j = 0, i = 0; i <= i1; ++i) {
    T[j][i] = T1;
  }

  for (j = 0, i = 0; j <= j1; ++j) {
    T[j][i] = T1;
  }

  for (j = j2, i = i3; j <= j3; ++j) {
    T[j][i] = T2;
  }
  for (j = j3, i = i2; i <= i3; ++i) {
    T[j][i] = T2;
  }
  // CONTACTS
  std::ofstream MDD_write("../dataMDD.bin",
                          std::ios_base::out | std::ios_base::binary);
  double T_1, T_k1, dT, delta;
  // RECOUNT CORRECT TEMPERATURE
  while (k < REP) {
    for (int j = 0; j <= j3; ++j) {
      for (int i = 0; i <= i3; ++i) {

        // ABIK +
        if (i > 0 && i < i3 && j > 0 && j < j1) {
          T_1 = -(alf_1 * T[j - 1][i] + alf_2 * T[j][i - 1] +
                  alf_2 * T[j][i + 1] + alf_1 * T[j + 1][i]) /
                gam_1;
        }
        // BCDE +
        else if (i > 0 && i < i1 && j >= j1 && j < j3) {
          T_1 = -(alf_1 * T[j - 1][i] + alf_2 * T[j][i - 1] +
                  alf_2 * T[j][i + 1] + alf_1 * T[j + 1][i]) /
                gam_1;
        }
        // GFIH +
        else if (i > i2 && i < i3 && j >= j1 && j < j3) {
          T_1 = -(alf_1 * T[j - 1][i] + alf_2 * T[j][i - 1] +
                  alf_2 * T[j][i + 1] + alf_1 * T[j + 1][i]) /
                gam_1;
        }
        // MK +
        else if (i > i1 && i < i3 && j == 0) {
          T_1 = -(alf_3 * T[j][i - 1] + alf_3 * T[j][i + 1] +
                  alf_1 * T[j + 1][i]) /
                gam_3;
        }
        // K +
        else if (j == 0 && i == i3) {
          T_1 = -(alf_3 * T[j][i - 1] + alf_4 * T[j + 1][i]) / gam_4;
        }
        // BC +
        else if (j > j1 && j < j3 && i == 0) {
          T_1 = -(alf_4 * T[j - 1][i] + alf_2 * T[j][i + 1] +
                  alf_4 * T[j + 1][i]) /
                gam_3;
        }
        // C +
        else if (j == j3 && i == 0) {
          T_1 = -(alf_4 * T[j - 1][i] + alf_3 * T[j][i + 1]) / gam_4;
        }
        // CD +
        else if (i > 0 && i < i1 && j == j3) {
          T_1 = -(alf_3 * T[j][i - 1] + alf_3 * T[j][i + 1] +
                  alf_1 * T[j + 1][i]) /
                gam_3;
        }
        // LK +
        else if (i == i3 && j > 0 && j < j2) {
          T_1 = -(alf_4 * T[j - 1][i] + alf_2 * T[j][i - 1] +
                  alf_4 * T[j + 1][i]) /
                gam_3;
        }
        // D +
        else if (i == i1 && j == j3) {
          T_1 = -(alf_4 * T[j - 1][i] + alf_3 * T[j][i - 1]) / gam_4;
        }
        // DE +
        else if (i == i1 && j > j1 && j < j3) {
          T_1 = -(alf_4 * T[j - 1][i] + alf_2 * T[j][i - 1] +
                  alf_4 * T[j + 1][i]) /
                gam_3;
        }
        // E +
        else if (i == i1 && j == j1) {
          T_1 = -(alf_1 * T[j - 1][i] + alf_2 * T[j][i - 1] +
                  alf_3 * T[j][i + 1] + alf_4 * T[j + 1][i]) /
                gam_2;
        }
        // EF +
        else if (j == j1 && i > i1 && i < i2) {
          T_1 = -(alf_3 * T[j][i - 1] + alf_3 * T[j][i + 1] +
                  alf_1 * T[j + 1][i]) /
                gam_3;
        }
        // F +
        else if (i == i2 && j == j1) {
          T_1 = -(alf_1 * T[j - 1][i] + alf_3 * T[j][i - 1] +
                  alf_2 * T[j][i + 1] + alf_4 * T[j + 1][i]) /
                gam_2;
        }
        // FG +
        else if (i == i2 && j > j1 && j < j3) {
          T_1 = -(alf_4 * T[j - 1][i] + alf_2 * T[j][i + 1] +
                  alf_4 * T[j + 1][i]) /
                gam_3;
        } else
          continue;
        // iteration
        dT = lambda * (T_1 - T[j][i]);
        T_k1 = T[j][i] + dT;
        // diff
        delta = fabs(dT / T_k1);
        if (delta > max_rel_diff)
          max_rel_diff = delta;
        // results
        T[j][i] = T_k1;
      }
    }
    if (max_rel_diff < EPS)
      break;

    MDD_write << max_rel_diff << " ";

    max_rel_diff = 0.0f;
    ++k;
  }
  // RECOUNT CORRECT TEMPERATURE

  std::ofstream T_write("../dataT.txt", std::ios_base::out);
  std::ofstream T_write_bin("../dataT.bin",
                            std::ios_base::out | std::ios_base::binary);

  // FILE WRITE
  for (j = NY - 1; j >= 0; --j) {
    for (i = 0; i < NX; ++i) {
      T_write << std::fixed << std::setw(6) << std::setprecision(2) << T[j][i];
      T_write_bin << T[j][i] << " ";
    }
    T_write << std::endl;
    T_write_bin << std::endl;
  }
  // FILE WRITE

  std::cout << "Hi, Wasley:)" << std::endl;

  T_write.close();
  MDD_write.close();
  T_write_bin.close();

  return EXIT_SUCCESS;
}
