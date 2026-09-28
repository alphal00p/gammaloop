= Zero polynomial in multivariate partial fractions

Run `python reproduce.py` with Symbolica revision
`939c4de8b8da5546101b937a425eda443f409e8a`. No Spenso, Idenso or HEP
objects are involved. The expected result of `0.apart(x, y)` is zero.

The observed panic originates in `MultivariatePolynomial::from_coefficient_list`:
it divides the exponent-array length by the coefficient count before handling
an empty coefficient list. Multivariate partial fractioning reaches it through
`to_polynomial_in`. Retaining an empty polynomial with its declared variables
requires handling zero terms before that division.

Upstream main commit `6a96c9d77b217ea394770524c119e0613d4d9ef5`
fixes zero-polynomial construction by taking the variable count from the
variable list. It also accepts zero denominator factors in the univariate
partial-fraction shortcut. These are the only two source changes since the
failing revision above.

Keep this script and the tensor metadata suite's original zero-partial-fraction
assertions as regressions. Final qualification uses the fixed dependency;
verification of the updated executable is recorded separately from the original
panic receipt.
