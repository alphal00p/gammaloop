#include <iostream>
#include <math.h>
#include <quadmath.h>
int main() {
  char buffer[128];
  __float128 two=2;
  __float128 dp=sqrt(2.);
  __float128 qp=sqrtq(two);
  quadmath_snprintf(buffer, sizeof(buffer), "%.32Qe", dp*dp);
  std::cout << buffer << std::endl;
  quadmath_snprintf(buffer, sizeof(buffer), "%.32Qe", qp*qp);
  std::cout << buffer << std::endl;
  return 0;
}
