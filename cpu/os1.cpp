#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>
#include "grid.hpp"
#include "lensing_table.hpp"
#include "os.hpp"
#include "unit.hpp"
#include "ut.hpp"

double frequency_nu, spot_center_theta, obs_theta, angular_radius;
std::string ss;
char beaming;

int
main(int argc, char* argv[]) {
  if (argc != 6) {
    std::cout
        << "Usage: "
        << argv[0]
        << " <frequency_nu> <spot_center_theta> <obs_theta> <angular_radius> <a-j>\n";
  } else {
    frequency_nu = std::atof(argv[1]);
    spot_center_theta = std::atof(argv[2]) * degree;
    obs_theta = std::atof(argv[3]) * degree;
    angular_radius = std::atof(argv[4]);
    ss = argv[5];

    if ('a' <= ss[0] && ss[0] <= 'f') {
      beaming = 'i';  // isotropic
    } else if (ss[0] == 'g' || ss[0] == 'i') {
      beaming = 'c';  // cos
    } else if (ss[0] == 'h' || ss[0] == 'j') {
      beaming = 's';  // sin
    } else {
      std::cerr << "Beaming must be [a-f] for isotropic, [g,i] for cos, or [h,j] for sin.\n";
      return 1;
    }
  }

  constexpr double M = 1.4;      // M_sun;
  constexpr double Re = 12;      // km
  constexpr double kT = 0.35;    // keV
  constexpr double E_obs = 1.0;  // keV
  constexpr double D = 0.2 * kpc_in_km;

  std::vector<double> output_phase_grid = linspace(0, 0.992187500, 128);
  std::vector<double> output_E_grid{E_obs};

  LensingTable lt;

#define USE_HEALPIX

#ifdef USE_HEALPIX
  int N_side = 256;
#else
  int N_theta = 800;
  int N_phi = N_theta;
#endif

  Grid source = (angular_radius < 0.1)  //
                  ? PointSource(spot_center_theta, angular_radius, kT)
#ifdef USE_HEALPIX
                  : CapSourceHP(spot_center_theta, angular_radius, kT, N_side);
#else
                  : CapSourceUniformTheta(spot_center_theta, angular_radius, kT, N_theta, N_phi);
#endif

  auto time_start = std::chrono::high_resolution_clock::now();
  auto total_flux = CalculateTotalFluxOS(
      lt,
      source,
      output_phase_grid,
      output_E_grid,
      M,
      Re,
      D,
      frequency_nu,
      obs_theta,
      beaming
  );
  auto time_end = std::chrono::high_resolution_clock::now();
  double t_ms = std::chrono::duration<double, std::milli>(time_end - time_start).count();
  std::cout << "Time taken: " << t_ms << " ms\n";

  std::ofstream out_file("os1" + ss + ".txt");
  out_file << std::setprecision(16);

  for (double output_phase : output_phase_grid) {
    out_file << output_phase << " ";
  }
  out_file << "\n";

  for (double output_E : output_E_grid) {
    out_file << output_E << " ";
  }
  out_file << "\n";

  for (int i_phase = 0; i_phase < output_phase_grid.size(); ++i_phase) {
    for (int i_E = 0; i_E < output_E_grid.size(); ++i_E) {
      double E = output_E_grid[i_E];
      out_file << total_flux(i_phase, i_E) / E << " ";
    }
    out_file << "\n";
  }
}
