#include "seed.hpp"
#include "spqr.hpp"
#include <iostream>

int main() {
  std::string line;
  while (std::getline(std::cin, line)) {
    try {
      const auto input = ec::Json::parse(line);
      const auto result =
          input.contains("diagram")
              ? ec::initialize(input.at("diagram"), input.value("scale", 2.4),
                               input.value("external_sides", true))
              : ec::decompose(input);
      std::cout << result.dump() << '\n';
    } catch (const std::exception &error) {
      std::cout << ec::Json({{"error", error.what()}}).dump() << '\n';
    }
    std::cout.flush();
  }
  return std::cin.bad() ? 1 : 0;
}
