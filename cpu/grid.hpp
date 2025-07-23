#pragma once

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>
#include "unit.hpp"

struct Patch {
  double phi;
  double kT_in_keV;
};
struct Ring {
  double cos_theta;
  double dOmega;
  std::vector<Patch> patches;
};
using Grid = std::vector<Ring>;

inline Grid
PointSource(double spot_theta, double spot_angular_radius, double temp) {
  return Grid{
      Ring{
          .cos_theta = std::cos(spot_theta),
          .dOmega = two_pi * (1. - std::cos(spot_angular_radius)),
          .patches = {
              Patch{
                  .phi = 0,
                  .kT_in_keV = temp,
              },
          },
      },
  };
}

inline Grid
CapSourceUniformTheta(
    double spot_center_theta,
    double spot_angular_radius,
    double kT_in_keV,
    int N_theta,
    int N_phi
) {
  double spot_center_phi = 0.;
  double s1 = std::sin(spot_center_theta);
  double c1 = std::cos(spot_center_theta);
  double cos_angular_radius = std::cos(spot_angular_radius);
  double d_theta = pi / N_theta;
  double d_phi = two_pi / N_phi;
  Grid grid;
  for (int i = 0; i < N_theta; ++i) {
    double theta = (i + .5) * d_theta;
    double s2 = std::sin(theta);
    double c2 = std::cos(theta);
    double dOmega = s2 * d_theta * d_phi;
    Ring ring{
        .cos_theta = c2,
        .dOmega = dOmega,
        .patches = {},
    };
    for (int j = 0; j < N_phi; ++j) {
      double phi = j * d_phi;
      double cos_rho = s1 * s2 * std::cos(phi - spot_center_phi) + c1 * c2;
      if (cos_rho > cos_angular_radius) {
        ring.patches.push_back(
            Patch{
                .phi = phi,
                .kT_in_keV = kT_in_keV,
            }
        );
      }
    }
    if (!ring.patches.empty()) {
      grid.push_back(std::move(ring));
    }
  }
  return grid;
}

inline Grid
CapSourceUniformCosTheta(
    double spot_center_theta,
    double spot_angular_radius,
    double kT_in_keV,
    int N_theta,
    int N_phi
) {
  double spot_center_phi = 0.;
  double s1 = std::sin(spot_center_theta);
  double c1 = std::cos(spot_center_theta);
  double cos_angular_radius = std::cos(spot_angular_radius);
  double d_cos_theta = 2. / N_theta;
  double d_phi = two_pi / N_phi;
  double dOmega = d_cos_theta * d_phi;
  Grid grid;
  for (int i = 0; i < N_theta; ++i) {
    double c2 = -1. + (i + .5) * d_cos_theta;
    double s2 = std::sqrt(1. - c2 * c2);
    Ring ring{
        .cos_theta = c2,
        .dOmega = dOmega,
        .patches = {},
    };
    for (int j = 0; j < N_phi; ++j) {
      double phi = j * d_phi;
      double cos_rho = s1 * s2 * std::cos(phi - spot_center_phi) + c1 * c2;
      if (cos_rho > cos_angular_radius) {
        ring.patches.push_back(
            Patch{
                .phi = phi,
                .kT_in_keV = kT_in_keV,
            }
        );
      }
    }
    if (!ring.patches.empty()) {
      grid.push_back(std::move(ring));
    }
  }
  return grid;
}

inline Grid
CapSourceHP(double spot_center_theta, double spot_angular_radius, double kT_in_keV, int N_side) {
  double spot_center_phi = 0.;
  double s1 = std::sin(spot_center_theta);
  double c1 = std::cos(spot_center_theta);
  double cos_angular_radius = std::cos(spot_angular_radius);
  double dOmega = pi / (3. * N_side * N_side);

  // std::cout << "N_side: " << N_side << "\n";
  // std::cout << "sqrt(12 N_side^2) = " << std::sqrt(12. * N_side * N_side) << "\n";
  // std::cout << "dOmega: " << dOmega << "\n";

  auto cal_z_phi_ = [N_side](int i, int j) {
    if (1 <= i && i < N_side && 1 <= j && j <= 4 * i) {
      double z = 1. - 1. * i * i / (3. * N_side * N_side);
      double s = 1.;
      double phi = pi / (2. * i) * (j - s / 2.);
      return std::make_pair(z, phi);
    } else if (N_side <= i && i <= 2 * N_side && 1 <= j && j <= 4 * N_side) {
      double z = 4. / 3. - 2. * i / (3. * N_side);
      double s = (i - N_side + 1) % 2;
      double phi = pi / (2. * N_side) * (j - s / 2.);
      return std::make_pair(z, phi);
    } else {
      std::printf("In CapSourceHP, i: %d, j: %d is out of range.\n", i, j);
      std::exit(1);
    }
  };

  auto cal_z_phi = [N_side, cal_z_phi_](int i, int j) {
    if (i <= 2 * N_side) {
      return cal_z_phi_(i, j);
    } else {
      i = 4 * N_side - i;
      auto [z, phi] = cal_z_phi_(i, j);
      return std::make_pair(-z, phi);
    }
  };

  Grid grid;
  for (int i = 1; i <= 4 * N_side - 1; ++i) {
    Ring ring{
        .cos_theta = -999,
        .dOmega = dOmega,
        .patches = {},
    };
    int jmax = 4 * std::min({i, 4 * N_side - i, N_side});
    for (int j = 1; j <= jmax; ++j) {
      auto [cos_theta, phi] = cal_z_phi(i, j);
      ring.cos_theta = cos_theta;
      double c2 = cos_theta;
      double s2 = std::sqrt(1. - c2 * c2);
      double cos_rho = s1 * s2 * std::cos(phi - spot_center_phi) + c1 * c2;
      if (cos_rho > cos_angular_radius) {
        ring.patches.push_back(
            Patch{
                .phi = phi,
                .kT_in_keV = kT_in_keV,
            }
        );
      }
    }
    if (!ring.patches.empty()) {
      grid.push_back(std::move(ring));
    }
  }
  return grid;
}

// inline void
// print_grid(Grid const& grid) {
//   for (auto const& ring : grid) {
//     std::cout << "Ring cos_theta: " << ring.cos_theta << ", dOmega: " << ring.dOmega << "\n";
//     for (auto const& patch : ring.patches) {
//       std::cout << "  Patch phi: " << patch.phi << ", temp: " << patch.kT_in_keV << "\n";
//     }
//   }
// }
