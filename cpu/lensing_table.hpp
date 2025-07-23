#pragma once

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include "hunt.hpp"
#include "matrix.hpp"

struct LensingTable {
  UniformHunt cos_psi, u, cos_alpha;
  Matrix cos_alpha_of_u_cos_psi, lf_of_u_cos_psi, cdt_over_R_of_u_cos_alpha;

  void
  read_matrix(Matrix& mat, std::string const& fname, int n_rows, int n_cols) {
    std::ifstream in_file(fname);
    if (!in_file) {
      std::cerr << "Could not open lensing table file: " << fname << "\n";
      std::exit(EXIT_FAILURE);
    }
    mat.reset(n_rows, n_cols);
    for (int i = 0; i < n_rows; ++i) {
      for (int j = 0; j < n_cols; ++j) {
        if (!(in_file >> mat(i, j))) {
          std::cerr << "Error reading matrix from file: " << fname << "\n";
          std::exit(EXIT_FAILURE);
        }
      }
    }
    in_file.close();
  }

  LensingTable()
      : cos_psi("lensing_table/cos_psi.txt"),
        u("lensing_table/u.txt"),
        cos_alpha("lensing_table/cos_alpha.txt") {
    read_matrix(cos_alpha_of_u_cos_psi, "lensing_table/cos_alpha_of_u_cos_psi.txt", u.n, cos_psi.n);
    read_matrix(lf_of_u_cos_psi, "lensing_table/lf_of_u_cos_psi.txt", u.n, cos_psi.n);
    read_matrix(
        cdt_over_R_of_u_cos_alpha,
        "lensing_table/cdt_over_R_of_u_cos_alpha.txt",
        u.n,
        cos_alpha.n
    );
  }

  auto
  cal_cos_alpha_lf_of_u_cos_psi(double u_, double cos_psi_) {
    auto [i, a] = u(u_);
    auto [j, b] = cos_psi(cos_psi_);

    double cos_alpha_ = (1. - a) * (1. - b) * cos_alpha_of_u_cos_psi(i, j)
                      + a * (1. - b) * cos_alpha_of_u_cos_psi(i + 1, j)
                      + (1. - a) * b * cos_alpha_of_u_cos_psi(i, j + 1)
                      + a * b * cos_alpha_of_u_cos_psi(i + 1, j + 1);

    double lf_ = (1. - a) * (1. - b) * lf_of_u_cos_psi(i, j)
               + a * (1. - b) * lf_of_u_cos_psi(i + 1, j)
               + (1. - a) * b * lf_of_u_cos_psi(i, j + 1)
               + a * b * lf_of_u_cos_psi(i + 1, j + 1);

    return std::make_tuple(cos_alpha_, lf_);
  }

  double
  cal_cdt_over_R_of_u_cos_alpha(double u_, double cos_alpha_) {
    auto [i, a] = u(u_);
    auto [j, b] = cos_alpha(cos_alpha_);

    return (1. - a) * (1. - b) * cdt_over_R_of_u_cos_alpha(i, j)
         + a * (1. - b) * cdt_over_R_of_u_cos_alpha(i + 1, j)
         + (1. - a) * b * cdt_over_R_of_u_cos_alpha(i, j + 1)
         + a * b * cdt_over_R_of_u_cos_alpha(i + 1, j + 1);
  }
};
