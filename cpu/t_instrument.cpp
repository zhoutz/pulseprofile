#include "instrument.hpp"

int
main() {
  Instrument instrument(
      "model_data/nicer-consim135p-teamonly-array50_arf.txt",
      "model_data/nicer-rmf6s-teamonly-array50_full_matrix.txt",
      0,     // Eo_beg
      1800,  // Eo_end
      30,    // ch_beg
      300    // ch_end
  );
  double E_beg = instrument.Eo.front();
  double E_end = instrument.Eo.back();
  int n_E = instrument.Eo.size();
  printf("E: [%g,%g] keV, n_E=%d\n", E_beg, E_end, n_E);
}
