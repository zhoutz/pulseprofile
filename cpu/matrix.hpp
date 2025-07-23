#pragma once

#include <vector>

struct Matrix {
  int n_row = 0;
  int n_col = 0;
  std::vector<double> data_;

  Matrix() = default;
  Matrix(int n_row, int n_col) : n_row(n_row), n_col(n_col), data_(n_row * n_col) {}

  double&
  operator()(int row, int col) {
    return data_[row * n_col + col];
  }

  void
  fill(double value) {
    std::fill(data_.begin(), data_.end(), value);
  }
  void
  reset(int n_row_, int n_col_) {
    n_row = n_row_;
    n_col = n_col_;
    data_.resize(n_row * n_col);
  }
};
