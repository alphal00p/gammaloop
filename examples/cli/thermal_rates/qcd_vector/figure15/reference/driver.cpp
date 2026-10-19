// Reuse the exact vector and density assemblies from Zenodo3478144 v0.9-spf.
#include <chrono>
#include <cstdlib>
#include <string>
#define main spectral_example_main
#include "examples/dileptons.cpp"
#undef main

int main(int argc, char **argv) {
  if (argc != 5) {
    std::cerr << "usage: reference E_over_T k_over_T tolerance vector|density\n";
    return 1;
  }
  k0 = std::atof(argv[1]);
  k = std::atof(argv[2]);
  tol = std::atof(argv[3]);
  const std::string component = argv[4];
  if (k <= 0 || k == k0 || tol <= 0 ||
      (component != "vector" && component != "density")) {
    std::cerr << "require k>0, E!=k, tolerance>0 and a valid component\n";
    return 2;
  }
  const auto start = std::chrono::steady_clock::now();
  double lo, nlo;
  if (component == "vector") {
    Rho_V rho;
    rho();
    lo = rho.lo;
    nlo = rho.nlo;
  } else {
    Rho_00 rho;
    rho();
    lo = rho.lo;
    nlo = rho.nlo;
  }
  const double seconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
          .count();
  std::cout << std::setprecision(17) << k0 << ',' << k << ',' << tol << ','
            << lo << ',' << nlo << ',' << seconds << '\n';
}
