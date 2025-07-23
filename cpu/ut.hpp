#pragma once

#include <vector>

inline std::vector<double>
linspace(double beg, double end, int num_points) {
  std::vector<double> vec(num_points);
  double step = (end - beg) / (num_points - 1);
  for (int i = 0; i < num_points; ++i) {
    vec[i] = beg + i * step;
  }
  return vec;
}

inline std::vector<double>
geomspace(double beg, double end, int num_points) {
  std::vector<double> vec(num_points);
  double log_beg = std::log(beg);
  double log_end = std::log(end);
  double step = (log_end - log_beg) / (num_points - 1);
  for (int i = 0; i < num_points; ++i) {
    vec[i] = std::exp(log_beg + i * step);
  }
  return vec;
}
