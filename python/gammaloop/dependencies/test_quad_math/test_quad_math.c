#include <stdio.h>
#include <math.h>
#include <quadmath.h>

int main() {
    char buffer[128];
    __float128 two = 2;
    __float128 dp = sqrt(2.); // Using double precision sqrt
    __float128 qp = sqrtq(two); // Using quad precision sqrt
    
    // Print the square of dp using quadmath_snprintf and printf
    quadmath_snprintf(buffer, sizeof(buffer), "%.32Qe", dp * dp);
    printf("%s\n", buffer);
    
    // Print the square of qp using quadmath_snprintf and printf
    quadmath_snprintf(buffer, sizeof(buffer), "%.32Qe", qp * qp);
    printf("%s\n", buffer);
    
    return 0;
}
