#pragma once

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <vector>

struct VecHunt {
  std::vector<double> const& xx;
  int jsav;

  VecHunt(std::vector<double> const& x) : xx(x), jsav(0) {}

  auto
  operator()(double x) {
    if ((x - xx.front()) * (x - xx.back()) > 0) {
      printf(
          "VecHunt::operator(): x = %g is out of bounds [%g, %g]: %d\n",
          x,
          xx.front(),
          xx.back(),
          (int)xx.size()
      );
    }
    int n = xx.size();
    int jl = jsav, jm, ju, inc = 1;
    bool ascnd = (xx.back() >= xx.front());
    if (jl < 0 || jl > n - 1) {
      jl = 0;
      ju = n - 1;
    } else {
      if (x >= xx[jl] == ascnd) {
        for (;;) {
          ju = jl + inc;
          if (ju >= n - 1) {
            ju = n - 1;
            break;
          } else if (x < xx[ju] == ascnd)
            break;
          else {
            jl = ju;
            inc += inc;
          }
        }
      } else {
        ju = jl;
        for (;;) {
          jl = jl - inc;
          if (jl <= 0) {
            jl = 0;
            break;
          } else if (x >= xx[jl] == ascnd)
            break;
          else {
            ju = jl;
            inc += inc;
          }
        }
      }
    }
    while (ju - jl > 1) {
      jm = (ju + jl) >> 1;
      if (x >= xx[jm] == ascnd)
        jl = jm;
      else
        ju = jm;
    }
    jsav = jl;
    int i = jl < 0 ? 0 : (jl > n - 2 ? n - 2 : jl);
    double a = (x - xx[i]) / (xx[i + 1] - xx[i]);
    return std::make_pair(i, a);
  }

  double
  operator[](int i) {
    if (i < 0 || i >= (int)xx.size()) {
      printf("VecHunt::operator[]: Index %d out of bounds [0, %d)\n", i, (int)xx.size());
      std::exit(EXIT_FAILURE);
    }
    return xx[i];
  }
};

struct UniformHunt {
  double x_min, x_max;
  int n;

  UniformHunt(double x_min, double x_max, int n) : x_min(x_min), x_max(x_max), n(n) {}

  UniformHunt(char const* fname) {
    std::ifstream in_file(fname);
    if (!in_file.is_open()) {
      std::printf("UniformHunt::UniformHunt: Could not open file %s\n", fname);
      std::exit(EXIT_FAILURE);
    }
    in_file >> x_min >> x_max >> n;
    if (n < 2) {
      std::printf("UniformHunt::UniformHunt: Error: n must be at least 2, got %d\n", n);
      std::exit(EXIT_FAILURE);
    }
    in_file.close();
  }

  auto
  operator()(double x) {
    int i;
    double a;

    if (x < x_min) {
      printf("UniformHunt::operator(): x = %g is out of bounds [%g, %g]: %d\n", x, x_min, x_max, n);
      i = 0;
      a = 0;
      return std::make_pair(i, a);
    }
    if (x > x_max) {
      printf("UniformHunt::operator(): x = %g is out of bounds [%g, %g]: %d\n", x, x_min, x_max, n);
      i = n - 2;
      a = 1;
      return std::make_pair(i, a);
    }
    double t = (x - x_min) / (x_max - x_min) * (n - 1);
    i = int(t);
    a = t - i;
    return std::make_pair(i, a);
  }

  double
  operator[](int i) {
    if (i < 0 || i >= n) {
      printf("UniformHunt::operator[]: Index %d out of bounds [0, %d)\n", i, n);
      std::exit(EXIT_FAILURE);
    }
    return x_min + (x_max - x_min) * i / (n - 1);
  }
};
