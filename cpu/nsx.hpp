#pragma once

#include <cstdio>
#include <fstream>
#include <vector>

#include "hunt.hpp"

struct nsx_H_v200804 {
  static constexpr auto nsx_fname = "model_data/nsx_H_v200804.out";
  static constexpr int NlogT = 35;
  static constexpr int Nlogg = 14;
  static constexpr int NlogEkT = 166;
  static constexpr int Nmu = 67;
  static constexpr double logT_min = 5.10;
  static constexpr double logT_max = 6.80;
  static constexpr double logg_min = 13.7;
  static constexpr double logg_max = 15.0;
  static constexpr double logEkT_min = -1.30;
  static constexpr double logEkT_max = 2.00;
};

struct nsx_H_v171019 {
  static constexpr auto nsx_fname = "model_data/nsx_H_v171019.out";
  static constexpr int NlogT = 35;
  static constexpr int Nlogg = 11;
  static constexpr int NlogEkT = 166;
  static constexpr int Nmu = 67;
  static constexpr double logT_min = 5.10;
  static constexpr double logT_max = 6.80;
  static constexpr double logg_min = 13.7;
  static constexpr double logg_max = 14.7;
  static constexpr double logEkT_min = -1.30;
  static constexpr double logEkT_max = 2.00;
};

template <typename Config>
struct NSX : Config {
  using Config::logEkT_max;
  using Config::logEkT_min;
  using Config::logg_max;
  using Config::logg_min;
  using Config::logT_max;
  using Config::logT_min;
  using Config::NlogEkT;
  using Config::Nlogg;
  using Config::NlogT;
  using Config::Nmu;
  using Config::nsx_fname;

  std::vector<double> mu_vec;
  VecHunt mu_hunt;
  UniformHunt logEkT_hunt, logT_hunt, logg_hunt;
  std::vector<double> data;

  double&
  operator()(int i_logT, int i_logg, int i_logEkT, int i_mu) {
    if (i_logT < 0
        || i_logT >= NlogT
        || i_logg < 0
        || i_logg >= Nlogg
        || i_logEkT < 0
        || i_logEkT >= NlogEkT
        || i_mu < 0
        || i_mu >= Nmu) {
      printf("NSX::operator() : Index out of bounds\n");
      std::exit(EXIT_FAILURE);
    }
    int index = i_logT * Nlogg * NlogEkT * Nmu  //
              + i_logg * NlogEkT * Nmu          //
              + i_logEkT * Nmu                  //
              + i_mu;
    return data[index];
  }

  NSX()
      : mu_vec(Nmu),
        mu_hunt(mu_vec),
        logEkT_hunt(logEkT_min, logEkT_max, NlogEkT),
        logT_hunt(logT_min, logT_max, NlogT),
        logg_hunt(logg_min, logg_max, Nlogg),
        data(NlogT * Nlogg * NlogEkT * Nmu, 0.0) {
    std::ifstream in_file(nsx_fname);
    // in_file.tie(0)->sync_with_stdio(false);
    if (!in_file) {
      printf("Could not open NSX file: %s\n", nsx_fname);
      std::exit(EXIT_FAILURE);
    }

    printf("Loading NSX data from %s\n", nsx_fname);

    for (int i_logT = 0; i_logT < NlogT; ++i_logT) {
      for (int i_logg = 0; i_logg < Nlogg; ++i_logg) {
        for (int i_logEkT = 0; i_logEkT < NlogEkT; ++i_logEkT) {
          for (int i_mu = 0; i_mu < Nmu; ++i_mu) {
            double logEkT, mu, logInuT3, logT, logg;
            in_file >> logEkT >> mu >> logInuT3 >> logT >> logg;

            if (!in_file) {
              printf(
                  "Error reading NSX data at indices: %d, %d, %d, %d\n",
                  i_logT,
                  i_logg,
                  i_logEkT,
                  i_mu
              );
              std::exit(EXIT_FAILURE);
            }

            if (mu_hunt[i_mu] == 0) {
              // mu_hunt[i_mu] = mu;
              mu_vec[i_mu] = mu;
            } else {
              if (!IsClose(mu_hunt[i_mu], mu)) {
                printf("Inconsistent mu value at index %d: %g vs %g\n", i_mu, mu_hunt[i_mu], mu);
                std::exit(EXIT_FAILURE);
              }
            }

            if (!IsClose(logEkT_hunt[i_logEkT], logEkT)) {
              printf(
                  "Inconsistent logEkT value at index %d: %g vs %g\n",
                  i_logEkT,
                  logEkT_hunt[i_logEkT],
                  logEkT
              );
              std::exit(EXIT_FAILURE);
            }

            if (!IsClose(logT_hunt[i_logT], logT)) {
              printf(
                  "Inconsistent logT value at index %d: %g vs %g\n",
                  i_logT,
                  logT_hunt[i_logT],
                  logT
              );
              std::exit(EXIT_FAILURE);
            }

            if (!IsClose(logg_hunt[i_logg], logg)) {
              printf(
                  "Inconsistent logg value at index %d: %g vs %g\n",
                  i_logg,
                  logg_hunt[i_logg],
                  logg
              );
              std::exit(EXIT_FAILURE);
            }

            // this->operator()(i_logT, i_logg, i_logEkT, i_mu) = logInuT3;
            this->operator()(i_logT, i_logg, i_logEkT, i_mu) = std::pow(10, logInuT3);
          }
        }
      }
    }
    in_file.close();
    printf("NSX data loaded successfully from %s\n", nsx_fname);

    for (int i = 0; i < Nmu - 1; ++i) {
      if (mu_hunt[i] < mu_hunt[i + 1]) {
        printf(
            "mu_hunt is not monotonically decreasing at index %d: %g > %g\n",
            i,
            mu_hunt[i],
            mu_hunt[i + 1]
        );
        std::exit(EXIT_FAILURE);
      }
    }
  }

  bool
  IsClose(double a, double b) {
    constexpr double atol = 10 * std::numeric_limits<double>::epsilon();
    constexpr double rtol = 1e-8;
    double bigger = std::max(std::abs(a), std::abs(b));
    return std::abs(a - b) <= atol + rtol * bigger;
  }

  double
  Interp_logIT3(double logT, double logg, double logEkT, double mu) {
    auto [i_logT, a_logT] = logT_hunt(logT);
    auto [i_logg, a_logg] = logg_hunt(logg);
    auto [i_logEkT, a_logEkT] = logEkT_hunt(logEkT);
    auto [i_mu, a_mu] = mu_hunt(mu);

    double ret = 0;
    for (int j_logT = 0; j_logT <= 1; ++j_logT) {
      for (int j_logg = 0; j_logg <= 1; ++j_logg) {
        for (int j_logEkT = 0; j_logEkT <= 1; ++j_logEkT) {
          for (int j_mu = 0; j_mu <= 1; ++j_mu) {
            double v = this->operator()(
                i_logT + j_logT,
                i_logg + j_logg,
                i_logEkT + j_logEkT,
                i_mu + j_mu
            );
            double weight = (j_logT * a_logT + (1 - j_logT) * (1 - a_logT))
                          * (j_logg * a_logg + (1 - j_logg) * (1 - a_logg))
                          * (j_logEkT * a_logEkT + (1 - j_logEkT) * (1 - a_logEkT))
                          * (j_mu * a_mu + (1 - j_mu) * (1 - a_mu));
            ret += v * weight;
            // if (weight > 1e-10) {
            //   printf("j_logEkT = %d, a_logEkT = %g, ", j_logEkT, a_logEkT);
            //   printf("v = %g, weight = %g\n", v, weight);
            //   printf(
            //       "i_logT = %d, i_logg = %d, i_logEkT = %d, i_mu = %d\n",
            //       i_logT + j_logT,
            //       i_logg + j_logg,
            //       i_logEkT + j_logEkT,
            //       i_mu + j_mu
            //   );
            // }
          }
        }
      }
    }
    // printf("ret = %g\n", ret);
    return ret;
  }
};
