#pragma once

#include <cmath>

// https://pdg.lbl.gov/2024/reviews/rpp2024-rev-astrophysical-constants.pdf
constexpr double schwarzschild_radius_of_sun_in_km = 2.9532501;  // km

// https://en.wikipedia.org/wiki/Solar_mass
constexpr double sqrt_km3_over_s2_GMsun = 2.745011592867327e-6;

constexpr double pi = 3.14159265358979323846;
constexpr double two_pi = 2. * pi;

constexpr double h_in_keV_s = 4.135668e-18;
constexpr double c_in_m_s = 299792458.;
constexpr double c_in_km_s = c_in_m_s / 1e3;
constexpr double c_in_cm_s = c_in_m_s * 1e2;

constexpr double Ibb_constant = 2. / (h_in_keV_s * h_in_keV_s * h_in_keV_s * c_in_cm_s * c_in_cm_s);

constexpr double kpc_in_km = 3.08567758e+16;

constexpr double degree = pi / 180.;

// blackbody spectrum, E in keV, kT in keV
inline double
Ibb(double E, double kT) {
  return Ibb_constant * (E * E * E) / std::expm1(E / kT);
}
