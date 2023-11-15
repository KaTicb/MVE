#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#define h 0.2f
#define r 0.3f
#define NXB 15
#define NYB 12
#define REP 30000
#define T1 150.0f
#define T2 200.0f

constexpr std::size_t NX = 3 * NXB + 1;
constexpr std::size_t NY = 3 * NYB + 1;

void setup_contacts(std::array<std::array<double, NX>, NY> &T, size_t i1,
                    size_t i2, size_t i3, size_t j1, size_t j2, size_t j3) {
  std::size_t i;
  std::size_t j;

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
}

void write_to_file(const std::filesystem::path &dir_path,
                   const std::string &filename,
                   const std::array<std::array<double, NX>, NY> &T) {

  if (!std::filesystem::exists(dir_path)) {
    if (!std::filesystem::create_directories(dir_path)) {
      std::cerr << "Failed to create directory " << dir_path << std::endl;
      return;
    }
  }

  const auto file_path = dir_path / filename;
  std::ofstream file(file_path,
                     std::ios::out | std::ios::binary | std::ios::trunc);

  if (!file.is_open()) {
    std::cerr << "Unable to open file " << file_path << std::endl;
    return;
  }

  for (int j = NY - 1; j >= 0; --j) {
    for (int i = 0; i < NX; ++i) {
      file << T[j][i] << " ";
    }
    file << std::endl;
  }

  file.close();
}

int main(int argc, char **argv, char *env[]) {

  std::cout << "Command-line arguments:\n";
  for (int i = 0; i < argc; i++) {
    std::cout << "argv[" << i << "]: " << argv[i] << "\n";
  }

  // DEFINITION
  std::array<std::array<double, NX>, NY> T{};
  double s = h * r;

  std::size_t i1 = NXB;
  std::size_t i2 = NXB + NXB;
  std::size_t i3 = NXB + NXB + NXB;
  std::size_t j1 = NYB;
  std::size_t j2 = NYB + NYB;
  std::size_t j3 = NYB + NYB + NYB;

  std::size_t file_num = 0;

  double alf_1 = -h / r;
  double alf_2 = -r / h;
  double alf_3 = alf_2 * 0.5f;
  double alf_4 = alf_1 * 0.5f;
  double gam_1 = -2.f * (alf_1 + alf_2);
  double gam_2 = -1.5f * (alf_1 + alf_2);
  double gam_3 = -(alf_1 + alf_2);
  double gam_4 = -(alf_3 + alf_4);
  double beta_1 = std::pow(std::min(h, r), 2) / (4 * s); // 4 directs
  double beta_2 = beta_1 * 4 / 3;                        // 3 directs
  double beta_3 = beta_1 * 2;                            // 2 directs
  double beta_4 = beta_1 * 4;                            // 1 direct

  std::cout << "Successfully defined parameters" << std::endl;

  setup_contacts(T, i1, i2, i3, j1, j2, j3);

  std::cout << "Successfully setup contacts" << std::endl;

  // RECOUNT CORRECT TEMPERATURE
  for (std::size_t k = 0; k < REP; k++) {
    for (std::size_t j = 0; j <= j3; ++j) {
      for (std::size_t i = 0; i <= i3; ++i) {

        // ABIK +
        if (i > 0 && i < i3 && j > 0 && j < j1) {
          T[j][i] =
              T[j][i] - beta_1 * (alf_1 * T[j - 1][i] + alf_2 * T[j][i - 1] +
                                  alf_2 * T[j][i + 1] + alf_1 * T[j + 1][i] +
                                  T[j][i] * gam_1);
        }
        // BCDE +
        else if (i > 0 && i < i1 && j >= j1 && j < j3) {
          T[j][i] =
              T[j][i] - beta_1 * (alf_1 * T[j - 1][i] + alf_2 * T[j][i - 1] +
                                  alf_2 * T[j][i + 1] + alf_1 * T[j + 1][i] +
                                  T[j][i] * gam_1);
        }
        // GFIH +
        else if (i > i2 && i < i3 && j >= j1 && j < j3) {
          T[j][i] =
              T[j][i] - beta_1 * (alf_1 * T[j - 1][i] + alf_2 * T[j][i - 1] +
                                  alf_2 * T[j][i + 1] + alf_1 * T[j + 1][i] +
                                  T[j][i] * gam_1);
        }
        // MK +
        else if (i > i1 && i < i3 && j == 0) {
          T[j][i] =
              T[j][i] - beta_3 * (alf_3 * T[j][i - 1] + alf_3 * T[j][i + 1] +
                                  alf_1 * T[j + 1][i] + T[j][i] * gam_3);
        }
        // K +
        else if (j == 0 && i == i3) {
          T[j][i] = T[j][i] - beta_4 * (alf_3 * T[j][i - 1] +
                                        alf_4 * T[j + 1][i] + T[j][i] * gam_4);
        }
        // BC +
        else if (j > j1 && j < j3 && i == 0) {
          T[j][i] =
              T[j][i] - beta_3 * (alf_4 * T[j - 1][i] + alf_2 * T[j][i + 1] +
                                  alf_4 * T[j + 1][i] + T[j][i] * gam_3);
        }
        // C +
        else if (j == j3 && i == 0) {
          T[j][i] = T[j][i] - beta_4 * (alf_4 * T[j - 1][i] +
                                        alf_3 * T[j][i + 1] + T[j][i] * gam_4);
        }
        // CD +
        else if (i > 0 && i < i1 && j == j3) {
          T[j][i] =
              T[j][i] - beta_3 * (alf_3 * T[j][i - 1] + alf_3 * T[j][i + 1] +
                                  alf_1 * T[j + 1][i] + T[j][i] * gam_3);
        }
        // LK +
        else if (i == i3 && j > 0 && j < j2) {
          T[j][i] =
              T[j][i] - beta_3 * (alf_4 * T[j - 1][i] + alf_2 * T[j][i - 1] +
                                  alf_4 * T[j + 1][i] + T[j][i] * gam_3);
        }
        // D +
        else if (i == i1 && j == j3) {
          T[j][i] = T[j][i] - beta_4 * (alf_4 * T[j - 1][i] +
                                        alf_3 * T[j][i - 1] + T[j][i] * gam_4);
        }
        // DE +
        else if (i == i1 && j > j1 && j < j3) {
          T[j][i] =
              T[j][i] - beta_3 * (alf_4 * T[j - 1][i] + alf_2 * T[j][i - 1] +
                                  alf_4 * T[j + 1][i] + T[j][i] * gam_3);
        }
        // E +
        else if (i == i1 && j == j1) {
          T[j][i] =
              T[j][i] - beta_2 * (alf_1 * T[j - 1][i] + alf_2 * T[j][i - 1] +
                                  alf_3 * T[j][i + 1] + alf_4 * T[j + 1][i] +
                                  T[j][i] * gam_2);
        }
        // EF +
        else if (j == j1 && i > i1 && i < i2) {
          T[j][i] =
              T[j][i] - beta_3 * (alf_3 * T[j][i - 1] + alf_3 * T[j][i + 1] +
                                  alf_1 * T[j + 1][i] + T[j][i] * gam_3);
        }
        // F +
        else if (i == i2 && j == j1) {
          T[j][i] =
              T[j][i] - beta_2 * (alf_1 * T[j - 1][i] + alf_3 * T[j][i - 1] +
                                  alf_2 * T[j][i + 1] + alf_4 * T[j + 1][i] +
                                  T[j][i] * gam_2);
        }
        // FG +
        else if (i == i2 && j > j1 && j < j3) {
          T[j][i] =
              T[j][i] - beta_3 * (alf_4 * T[j - 1][i] + alf_2 * T[j][i + 1] +
                                  alf_4 * T[j + 1][i] + T[j][i] * gam_3);
        } else
          continue;
      }
    }

    // FILE WRITE
    if (!(k % 100)) {
      std::string dirPath = "data_1";
      std::string filename = "T" + std::to_string(++file_num) + ".bin";
      write_to_file(dirPath, filename, T);
    }
  }

  return EXIT_SUCCESS;
}
