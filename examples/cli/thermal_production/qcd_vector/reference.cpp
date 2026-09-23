// Reuse the author's full vector-channel assembly from Zenodo 3478144.
// Rename its example entry point; the master integrals and Rho_V are unchanged.
#include <cstdlib>

#define main spectral_example_main
#include "examples/dileptons.cpp"
#undef main

int main(int argc, char **argv) {
  if (argc != 4) {
    std::cerr << "usage: reference k0_over_T k_over_T absolute_tolerance\n";
    return 1;
  }
  k0 = std::atof(argv[1]);
  k = std::atof(argv[2]);
  tol = std::atof(argv[3]);
  Rho_V rho;
  rho();
  std::cout << std::setprecision(17) << k0 << ' ' << k << ' ' << tol << ' '
            << rho.lo << ' ' << rho.nlo << '\n';
}
