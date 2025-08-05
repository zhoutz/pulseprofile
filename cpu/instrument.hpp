#pragma once

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <vector>
#include "matrix.hpp"

struct Instrument {
  std::vector<double> Eo;
  Matrix RSP;

  Instrument(
      char const* arf_file_path,
      char const* rmf_file_path,
      int Eo_beg,
      int Eo_end,
      int ch_beg,
      int ch_end,
      int arf_skiprows = 3,
      int rmf_skiprows = 3
  ) {
    constexpr auto streamsize = std::numeric_limits<std::streamsize>::max();

    int n_ch = ch_end - ch_beg;
    int n_Eo = Eo_end - Eo_beg;

    RSP.reset(n_ch, n_Eo);

    // Read ARF file
    int arf_lineno;
    {
      std::ifstream arf_file(arf_file_path);
      if (!arf_file) {
        printf("Error opening ARF file: %s\n", arf_file_path);
        std::exit(EXIT_FAILURE);
      }
      Eo.resize(n_Eo);
      while (arf_skiprows--) arf_file.ignore(streamsize, '\n');
      for (arf_lineno = 0;; ++arf_lineno) {
        double Eo1, Eo2, effective_area;
        arf_file >> Eo1 >> Eo2 >> effective_area;
        if (arf_file.eof()) break;
        if (Eo_beg <= arf_lineno && arf_lineno < Eo_end) {
          Eo[arf_lineno - Eo_beg] = 0.5 * (Eo1 + Eo2);
          for (int ch = 0; ch < n_ch; ++ch) {
            RSP(ch, arf_lineno - Eo_beg) = effective_area;
          }
        }
      }
      arf_file.close();
      printf(
          "ARF file \"%s\":\n  contains %d lines, using line [%d,%d).\n",
          arf_file_path,
          arf_lineno,
          Eo_beg,
          Eo_end
      );
    }

    // Read RMF file
    {
      std::ifstream rmf_file(rmf_file_path);
      if (!rmf_file) {
        printf("Error opening RMF file: %s\n", rmf_file_path);
        std::exit(EXIT_FAILURE);
      }
      while (rmf_skiprows--) rmf_file.ignore(streamsize, '\n');

      for (int i_Eo = 0; i_Eo < arf_lineno; ++i_Eo) {
        double a, b;
        int c, d, n_total_channels;
        rmf_file >> a >> b >> c >> d >> n_total_channels;

        if (n_total_channels != 1501) {
          printf(
              "Warning: Expected 1501 channels, but got %d channels in RMF file.\n",
              n_total_channels
          );
        }

        for (int ch = 0; ch < n_total_channels; ++ch) {
          double value;
          rmf_file >> value;

          if (ch_beg <= ch && ch < ch_end && Eo_beg <= i_Eo && i_Eo < Eo_end) {
            RSP(ch - ch_beg, i_Eo - Eo_beg) *= value;
          }
        }
      }
      rmf_file.close();
    }

    // print the RSP matrix for debugging
    // printf("RSP matrix:\n");
    // for (int i = 0; i < n_ch; ++i) {
    //   for (int j = 0; j < n_Eo; ++j) {
    //     printf("RSP[%d,%d]=%f\n", i, j, RSP(i, j));
    //   }
    // }
  }
};
