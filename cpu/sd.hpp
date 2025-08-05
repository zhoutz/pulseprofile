#pragma once

#include <vector>
#include "grid.hpp"
#include "hunt.hpp"
#include "lensing_table.hpp"
#include "matrix.hpp"
#include "unit.hpp"

inline Matrix
CalculateTotalFluxSD(
    LensingTable& lt,
    Grid const& grid,
    std::vector<double> const& output_phase_grid,
    std::vector<double> const& output_E_grid,
    double M,
    double R,
    double D,
    double frequency_nu,
    double obs_theta
) {
  double u = M / R * schwarzschild_radius_of_sun_in_km;
  double uu = std::sqrt(1. - u);
  double R2 = R * R;
  double D2 = D * D;
  double frequency_Omega = frequency_nu * two_pi;
  double cos_obs_theta = std::cos(obs_theta);
  double sin_obs_theta = std::sin(obs_theta);

  int N_phase = output_phase_grid.size() * 3;

  std::vector<double> fluxes_over_I(N_phase + 1), redshift_factors(N_phase + 1);
  std::vector<double> phase_s(N_phase + 1), phase_o(N_phase + 1);

  VecHunt phase_o_hunt{phase_o};

  Matrix total_flux(output_phase_grid.size(), output_E_grid.size());
  total_flux.fill(0);

  for (auto const& ring : grid) {
    double cos_spot_theta = ring.cos_theta;
    double sin_spot_theta = std::sqrt(1. - cos_spot_theta * cos_spot_theta);

    double dOmega = ring.dOmega;
    double dS = dOmega * R2;

    double cc = cos_obs_theta * cos_spot_theta;
    double ss = sin_obs_theta * sin_spot_theta;

    int invisible_i_phase_min = N_phase + 1;
    int invisible_i_phase_max = -1;
    for (int i_phase = 0; i_phase < N_phase; ++i_phase) {
      bool is_visible = true;
      do {
        double phase = 1. * i_phase / N_phase;
        phase_s[i_phase] = phase;
        double spot_phi = two_pi * phase;
        double cos_spot_phi = std::cos(spot_phi);
        double sin_spot_phi = std::sin(spot_phi);
        double cos_psi = cc + ss * cos_spot_phi;
        if (cos_psi < lt.cos_psi.x_min) {
          is_visible = false;
          break;
        }
        double sin_psi = std::sqrt(1. - cos_psi * cos_psi);
        auto [cos_alpha, lf] = lt.cal_cos_alpha_lf_of_u_cos_psi(u, cos_psi);
        if (cos_alpha <= 0.) {
          is_visible = false;
          break;
        }
        double sin_alpha = std::sqrt(1. - cos_alpha * cos_alpha);
        double sin_alpha_over_sin_psi = cos_psi == 1. ? std::sqrt(lf) : sin_alpha / sin_psi;
        double beta = frequency_Omega * R * sin_spot_theta / uu / c_in_km_s;
        double gamma = 1. / std::sqrt(1. - beta * beta);
        double cos_xi = -sin_alpha_over_sin_psi * sin_obs_theta * sin_spot_phi;
        double delta = 1. / (gamma * (1. - beta * cos_xi));
        double delta3 = delta * delta * delta;

        double cdt_over_R = lt.cal_cdt_over_R_of_u_cos_alpha(u, cos_alpha);
        double dt1 = cdt_over_R * R / c_in_km_s;
        double delta_phase = dt1 * frequency_nu;

        double cos_alpha_prime = cos_alpha * delta;

        // This is the analytic formula used in Bogdanov 2019
        // fluxes_over_I[i_phase] = uu * delta3 * cos_alpha_prime * lf * (dS * gamma) / D2;

        // But we can only reproduce sharon's numerical results using the following formula
        fluxes_over_I[i_phase] = uu * delta3 * cos_alpha * lf * (dS) / D2;

        // we can only reproduce bogdanov's numerical results using the following formula
        // fluxes_over_I[i_phase] = uu * delta3 * cos_alpha * lf * (dS / gamma) / D2;

        redshift_factors[i_phase] = 1. / (delta * uu);
        phase_o[i_phase] = phase + delta_phase;
      } while (0);
      if (is_visible == false) {
        fluxes_over_I[i_phase] = 0;
        redshift_factors[i_phase] = -1;
        phase_o[i_phase] = -1;
        invisible_i_phase_min = std::min(invisible_i_phase_min, i_phase);
        invisible_i_phase_max = std::max(invisible_i_phase_max, i_phase);
      }
    }

    if (invisible_i_phase_min == 0 || invisible_i_phase_max == N_phase - 1) continue;
    if (invisible_i_phase_min <= invisible_i_phase_max) {
      int beg = invisible_i_phase_min - 1;
      int end = invisible_i_phase_max + 1;
      double phase0 = phase_o[beg];
      double phase1 = phase_o[end];
      double rf0 = redshift_factors[beg];
      double rf1 = redshift_factors[end];
      for (int i_phase = invisible_i_phase_min; i_phase <= invisible_i_phase_max; ++i_phase) {
        phase_o[i_phase] = phase0 + (phase1 - phase0) * (double(i_phase - beg) / (end - beg));
        redshift_factors[i_phase] = rf0 + (rf1 - rf0) * (double(i_phase - beg) / (end - beg));
      }
    }

    // https://numpy.org/doc/stable/reference/generated/numpy.gradient.html
    for (int i_phase = 0; i_phase < N_phase; ++i_phase) {
      double hs = i_phase == 0 ? phase_o[0] + 1. - phase_o[N_phase - 1]
                               : phase_o[i_phase] - phase_o[i_phase - 1];
      double hd = i_phase == N_phase - 1 ? phase_o[0] + 1. - phase_o[N_phase - 1]
                                         : phase_o[i_phase + 1] - phase_o[i_phase];
      double df = 1. / N_phase;
      double d_phase_s_d_phase_o = df * (hs * hs + hd * hd) / (hs * hd * (hs + hd));
      fluxes_over_I[i_phase] *= d_phase_s_d_phase_o;
    }

    fluxes_over_I[N_phase] = fluxes_over_I[0];
    redshift_factors[N_phase] = redshift_factors[0];
    phase_s[N_phase] = 1.0;
    phase_o[N_phase] = 1.0 + phase_o[0];

    for (int i_output_phase = 0; i_output_phase < output_phase_grid.size(); ++i_output_phase) {
      for (auto const& patch : ring.patches) {
        double phase_shift = patch.phi / two_pi;
        double patch_phase = std::fmod(output_phase_grid[i_output_phase] + phase_shift + 1., 1.);
        if (patch_phase < phase_o[0]) patch_phase += 1.;
        auto [i, a] = phase_o_hunt(patch_phase);
        double flux_over_I = (1 - a) * fluxes_over_I[i] + a * fluxes_over_I[i + 1];
        double redshift_factor = (1 - a) * redshift_factors[i] + a * redshift_factors[i + 1];
        for (int i_output_E = 0; i_output_E < output_E_grid.size(); ++i_output_E) {
          double E_obs = output_E_grid[i_output_E];
          double E_emit = E_obs * redshift_factor;
          double I = Ibb(E_emit, patch.kT_in_keV);
          total_flux(i_output_phase, i_output_E) += flux_over_I * I;
        }
      }
    }
  }

  return total_flux;
}
