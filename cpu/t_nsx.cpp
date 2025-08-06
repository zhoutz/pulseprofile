#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include "nsx.hpp"
#include "unit.hpp"
#include "ut.hpp"

int
main(int argc, char* argv[]) {
  double logg = 13.8;
  double logT = 6.11;
  double mu = 0.5;

  double T_in_kelvin = std::pow(10., logT);
  double kT_in_keV = T_in_kelvin * 8.617333262E-8;
  printf("kT_in_keV = %g\n", kT_in_keV);
  double T3 = T_in_kelvin * T_in_kelvin * T_in_kelvin;

  NSX<nsx_H_v200804> nsx;
  // NSX<nsx_H_v171019> nsx;
  std::ofstream out_file("nsx_output.txt");

  auto E_obs_in_keV_grid = geomspace(std::pow(10., -2.), std::pow(10., 0.5), 1000);
  for (auto E_obs_in_keV : E_obs_in_keV_grid) {
    double logEkT = std::log10(E_obs_in_keV / kT_in_keV);
    // double logEkT = std::log10(E_obs_in_keV);
    double logIT3 = nsx.Interp_logIT3(logT, logg, logEkT, mu);
    // double I = std::pow(10, logIT3) * T3 / 6.62607015E-27;
    // double I = std::pow(10, logIT3) * std::pow(10, 3*logT);
    double I = logIT3 * std::pow(10, 3 * logT);
    I /= 6.62607015E-27;
    double Ibb_value = Ibb(E_obs_in_keV, kT_in_keV);
    out_file << E_obs_in_keV << " " << I << " " << Ibb_value << "\n";
    // break;
  }
  out_file.close();
}
