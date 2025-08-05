#pragma once

#include <cmath>
#include <cstdio>
#include <fstream>
#include <vector>

#include "hunt.hpp"

struct Interstellar {
  std::vector<double> E, attenuation;
  VecHunt E_hunt;

  Interstellar(char const* file_path) : E_hunt(E) {
    std::ifstream file(file_path);
    if (!file) {
      std::printf("Error opening interstellar absorption file: %s\n", file_path);
      std::exit(EXIT_FAILURE);
    }

    double e, att, _;
    while (file >> e >> _ >> att) {
      E.push_back(e);
      attenuation.push_back(att);
    }
  }

  double
  operator()(double E_keV) {
    auto [i, a] = E_hunt(E_keV);
    double att = (1 - a) * attenuation[i] + a * attenuation[i + 1];
    return std::pow(att, 5.);
  }
};
