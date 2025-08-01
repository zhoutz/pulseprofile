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
}
