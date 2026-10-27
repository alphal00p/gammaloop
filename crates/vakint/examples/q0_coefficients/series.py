"""Exact, factorwise Laurent jets of q0 scalar master-product coefficients.

No host is imported and no work runs at import time.  Public Atom iterators and
polynomial conversion expose individual polynomial factors; no complete scalar
series, cancellation, common denominator, or rational-polynomial constructor is
used.  Python Fraction recurrences are checked against independent native Series
arithmetic on pole-stripped polynomial units.  Numerical evaluation is absent.
"""

import hashlib
import json
import math
import re
from fractions import Fraction


class FactorwiseSeriesError(ValueError):
    """Unsupported scalar structure or disagreement of exact arithmetic paths."""


class _NotPolynomial(FactorwiseSeriesError):
    pass


def _require(condition, message):
    if not condition:
        raise FactorwiseSeriesError(message)


def _fraction(expression, context):
    text = expression.to_canonical_string()
    if not re.fullmatch(r"[+-]?\d+(?:/[+-]?\d+)?", text):
        raise _NotPolynomial(
            f"{context}: expected a rational literal, found {text[:160]!r}"
        )
    try:
        return Fraction(text)
    except (ValueError, ZeroDivisionError) as error:
        raise _NotPolynomial(f"{context}: invalid rational literal") from error


def _integer(expression, context):
    number = _fraction(expression, context)
    _require(
        number.denominator == 1, f"{context}: only integer factor powers are supported"
    )
    return number.numerator


def _trim(polynomial):
    return {
        degree: coefficient for degree, coefficient in polynomial.items() if coefficient
    }


def _multiply(left, right, depth):
    result = {}
    for i, a in left.items():
        for j, b in right.items():
            if i + j <= depth:
                result[i + j] = result.get(i + j, Fraction(0)) + a * b
    return _trim(result)


def _power(polynomial, exponent, depth):
    _require(
        type(exponent) is int and exponent >= 0, "nonnegative polynomial power required"
    )
    result = {0: Fraction(1)}
    base = {degree: value for degree, value in polynomial.items() if degree <= depth}
    while exponent:
        if exponent & 1:
            result = _multiply(result, base, depth)
        exponent >>= 1
        if exponent:
            base = _multiply(base, base, depth)
    return result


def _divide_jet(p, q, depth):
    """Solve q*r=p modulo ep**(depth+1), then check the residual separately."""
    q0 = q.get(0, Fraction(0))
    _require(q0 != 0, "pole-stripped denominator has zero constant coefficient")
    result = {}
    for degree in range(depth + 1):
        previous = sum(
            (
                coefficient * result.get(degree - j, Fraction(0))
                for j, coefficient in q.items()
                if 0 < j <= degree
            ),
            Fraction(0),
        )
        result[degree] = (p.get(degree, Fraction(0)) - previous) / q0
    residual = _multiply(q, result, depth)
    _require(
        all(
            residual.get(degree, Fraction(0)) == p.get(degree, Fraction(0))
            for degree in range(depth + 1)
        ),
        "Fraction P-Q*jet truncated-ring certificate failed",
    )
    return _trim(result)


def _atom(core, polynomial, variable):
    result = core.E("0")
    for degree, coefficient in sorted(polynomial.items()):
        result += core.E(str(coefficient)) * variable**degree
    return result


def _gamma_atom(core, polynomial, gamma):
    return _atom(core, polynomial, gamma)


def _add_gamma(target, power, gamma_power, coefficient):
    if coefficient:
        polynomial = target.setdefault(power, {})
        polynomial[gamma_power] = polynomial.get(gamma_power, Fraction(0)) + coefficient


def _clean_gamma(values):
    return {
        power: cleaned
        for power, polynomial in values.items()
        if (cleaned := _trim(polynomial))
    }


class _Context:
    def __init__(self, core, ep, gamma, cutoff):
        self.core, self.ep, self.gamma, self.cutoff = core, ep, gamma, cutoff
        self.zero = core.E("0")
        self.one = core.E("1")
        self.normalization = core.E("𝑒") ** (-5 * gamma * ep)
        self.polynomials = {}
        self.nonpolynomials = set()
        self.native_powers = {}
        self.normalization_series = {}
        self.polynomial_catalog = {}

    def polynomial(self, expression):
        if expression in self.polynomials:
            return self.polynomials[expression]
        if expression in self.nonpolynomials:
            raise _NotPolynomial(
                "factor was already found to be nonpolynomial over Q(ep)"
            )
        try:
            # The public polynomial API performs no rational-function GCD.
            # It can introduce additional variables, so inspect EVERY resulting
            # coefficient: only exact Fraction literals may survive.
            converted = expression.to_polynomial(vars=[self.ep])
            values = {}
            for exponents, coefficient in converted.coefficient_list(self.ep):
                _require(
                    len(exponents) == 1
                    and type(exponents[0]) is int
                    and exponents[0] >= 0,
                    "invalid univariate polynomial exponent",
                )
                degree = exponents[0]
                value = _fraction(coefficient.to_expression(), "polynomial coefficient")
                values[degree] = values.get(degree, Fraction(0)) + value
            values = _trim(values)
        except _NotPolynomial:
            self.nonpolynomials.add(expression)
            raise
        except FactorwiseSeriesError:
            raise
        except Exception as error:
            self.nonpolynomials.add(expression)
            raise _NotPolynomial(
                "factor cannot be converted to a verified polynomial over Q(ep)"
            ) from error
        encoded = json.dumps(
            [[degree, str(value)] for degree, value in sorted(values.items())],
            separators=(",", ":"),
        ).encode()
        identity = hashlib.sha256(encoded).hexdigest()
        order = min(values) if values else None
        record = {
            "values": values,
            "order": order,
            "unit": {degree - order: value for degree, value in values.items()}
            if values
            else {},
            "identity": identity,
        }
        self.polynomials[expression] = record
        self.polynomial_catalog.setdefault(
            identity,
            {
                "sha256": identity,
                "degree": max(values) if values else None,
                "leading_order": order,
                "nonzero_coefficients": len(values),
                "coefficient_domain": "Q; all coefficients verified as exact Fraction literals",
            },
        )
        return record

    def branches(self, factors):
        """Distribute only immediate nonpolynomial Add factors, never a full AST."""
        for index, factor in enumerate(factors):
            if factor.get_type() != self.core.AtomType.Add:
                continue
            try:
                self.polynomial(factor)
            except _NotPolynomial:
                for child in factor:
                    children = (
                        list(child)
                        if child.get_type() == self.core.AtomType.Mul
                        else [child]
                    )
                    yield from self.branches(
                        factors[:index] + children + factors[index + 1 :]
                    )
                return
        yield factors

    def factors(self, term):
        factors = list(term) if term.get_type() == self.core.AtomType.Mul else [term]
        yield from self.branches(factors)

    def unit_power(self, polynomial, exponent, depth):
        key = (tuple(sorted(polynomial["unit"].items())), exponent, depth)
        if key not in self.native_powers:
            unit = {
                degree: value
                for degree, value in polynomial["unit"].items()
                if degree <= depth
            }
            _require(unit.get(0, 0) != 0, "native unit Series would start with zero")
            base = _atom(self.core, unit, self.ep).series(
                self.ep, 0, depth, depth_is_absolute=True
            )
            _require(
                base.get_ramification() == 1, "unexpected unit-Series ramification"
            )
            # Zero powers are skipped by the caller: native rpow(0) is undefined.
            self.native_powers[key] = base**exponent
        return self.native_powers[key]

    def analytic(self, depth):
        if depth not in self.normalization_series:
            series = self.normalization.series(
                self.ep, 0, depth, depth_is_absolute=True
            )
            for degree in range(depth + 1):
                wanted = (
                    self.core.E(str(Fraction((-5) ** degree, math.factorial(degree))))
                    * self.gamma**degree
                )
                _require(
                    (series[degree] - wanted).expand() == self.zero,
                    f"native normalization Series differs from closed formula at order {degree}",
                )
            self.normalization_series[depth] = series
        return self.normalization_series[depth]


def compute_scalar_jet(core, scalar, ep, gamma, cutoff=5):
    """Return ``(exact_jet_Atom, JSON_proof)`` through an absolute epsilon order.

    Every nonzero scalar term must contain exactly the saved exp(-5*gamma*ep)
    factor.  Other factors are integer powers of polynomials over Q(ep).
    Top-level sums and immediate sums within products are supported.  A rational
    sum underneath a power is deliberately rejected instead of invoking a broad
    CAS fallback.  No numerical retry, approximation, or master substitution is
    performed.  Remainders are coefficientwise in the formal master-product basis.
    """
    _require(type(cutoff) is int, "integer absolute Laurent cutoff required")
    context = _Context(core, ep, gamma, cutoff)
    native_total = context.zero.series(ep, 0, cutoff, depth_is_absolute=True)
    expected = {}
    rational_total = {}
    certificates = []
    minimum_term_valuation = None
    eligible_terms = 0
    exact_zero_terms = 0
    source_terms = list(scalar) if scalar.get_type() == core.AtomType.Add else [scalar]
    for source_index, source_term in enumerate(source_terms):
        if source_term == context.zero:
            exact_zero_terms += 1
            continue
        for factors in context.factors(source_term):
            normalization_count = sum(
                factor == context.normalization for factor in factors
            )
            _require(
                normalization_count == 1,
                f"source term {source_index}: expected one exact exp(-5*gamma*ep) factor, found {normalization_count}",
            )
            data = []
            is_zero = False
            valuation = 0
            for factor in factors:
                if factor == context.normalization:
                    continue
                if factor.get_type() == core.AtomType.Pow:
                    base, exponent_atom = list(factor)
                    exponent = _integer(
                        exponent_atom, f"source term {source_index} factor exponent"
                    )
                else:
                    base, exponent = factor, 1
                if exponent == 0:
                    continue
                try:
                    polynomial = context.polynomial(base)
                except _NotPolynomial as error:
                    raise FactorwiseSeriesError(
                        f"source term {source_index}: an integer-powered factor is not a polynomial over Q(ep); "
                        "unexpected constants, normalization, or rational sum underneath a power"
                    ) from error
                if polynomial["order"] is None:
                    _require(
                        exponent > 0,
                        f"source term {source_index}: zero denominator factor",
                    )
                    is_zero = True
                else:
                    valuation += exponent * polynomial["order"]
                data.append((polynomial, exponent))
            if is_zero:
                exact_zero_terms += 1
                continue
            minimum_term_valuation = (
                valuation
                if minimum_term_valuation is None
                else min(minimum_term_valuation, valuation)
            )
            bindings = [
                {"polynomial_sha256": p["identity"], "integer_power": n}
                for p, n in data
            ]
            certificate = {
                "source_term": source_index,
                "factor_bindings": bindings,
                "exact_term_valuation": valuation,
                "remainder_order_at_least": cutoff + 1,
            }
            if valuation > cutoff:
                certificate.update(
                    status="ABOVE_CUTOFF",
                    independent_native_series="not needed; nonzero unit denominators establish the exact valuation",
                )
                certificates.append(certificate)
                continue
            depth = cutoff - valuation
            p, q = {0: Fraction(1)}, {0: Fraction(1)}
            native = context.one.series(ep, 0, depth, depth_is_absolute=True)
            for polynomial, exponent in data:
                unit = {
                    degree: value
                    for degree, value in polynomial["unit"].items()
                    if degree <= depth
                }
                powered = _power(unit, abs(exponent), depth)
                if exponent > 0:
                    p = _multiply(p, powered, depth)
                else:
                    q = _multiply(q, powered, depth)
                try:
                    native = native * context.unit_power(polynomial, exponent, depth)
                except Exception as error:
                    raise FactorwiseSeriesError(
                        f"source term {source_index}: independent native unit-factor Series failed"
                    ) from error
            rational = _divide_jet(p, q, depth)
            for degree in range(depth + 1):
                _require(
                    _fraction(native[degree], "native rational unit coefficient")
                    == rational.get(degree, 0),
                    f"source term {source_index}: native unit Series and Fraction recurrence disagree at order {degree}",
                )
            normalized = native * context.analytic(depth)
            _require(
                normalized.get_ramification() == 1,
                "unexpected normalized Series ramification",
            )
            shifted = normalized.shift(valuation)
            order_n, order_d = shifted.get_absolute_order()
            _require(
                Fraction(order_n, order_d) >= cutoff + 1,
                "native shifted term loses requested absolute order",
            )
            native_total = native_total + shifted
            term_expected = {}
            for degree, value in rational.items():
                epsilon_power = valuation + degree
                rational_total[epsilon_power] = (
                    rational_total.get(epsilon_power, Fraction(0)) + value
                )
                for j in range(cutoff - epsilon_power + 1):
                    coefficient = value * Fraction((-5) ** j, math.factorial(j))
                    _add_gamma(term_expected, epsilon_power + j, j, coefficient)
                    _add_gamma(expected, epsilon_power + j, j, coefficient)
            term_expected = _clean_gamma(term_expected)
            for power in range(valuation, cutoff + 1):
                wanted = _gamma_atom(core, term_expected.get(power, {}), gamma)
                _require(
                    (shifted[power] - wanted).expand() == context.zero,
                    f"source term {source_index}: normalized native/Fraction jet disagreement at epsilon^{power}",
                )
            certificate.update(
                status="PASS",
                required_unit_depth=depth,
                unit_denominator_constant=str(q[0]),
                exact_unit_residual="Ptilde-Qtilde*jet == 0 modulo epsilon^(required_unit_depth+1)",
                independent_native_series="all pole-stripped unit factors and closed normalization compared exactly; then shifted without re-expansion",
            )
            certificates.append(certificate)
            eligible_terms += 1
    expected = _clean_gamma(expected)
    rational_total = _trim(rational_total)
    first = min(expected) if expected else None
    observed = {}
    for exponent, coefficient in native_total:
        power = _integer(exponent, "native summed Series exponent")
        _require(power <= cutoff, "native summed jet includes an unrequested order")
        observed[power] = coefficient
    for power in sorted(set(expected) | set(observed)):
        wanted = _gamma_atom(core, expected.get(power, {}), gamma)
        _require(
            (observed.get(power, context.zero) - wanted).expand() == context.zero,
            f"summed native/Fraction jet disagreement at epsilon^{power}",
        )
    order_n, order_d = native_total.get_absolute_order()
    _require(
        Fraction(order_n, order_d) >= cutoff + 1,
        "summed native Series has insufficient absolute order",
    )
    jet = context.zero
    for power, polynomial in sorted(expected.items()):
        jet += _gamma_atom(core, polynomial, gamma) * ep**power
    encoded = json.dumps(
        [
            [
                power,
                [[degree, str(value)] for degree, value in sorted(polynomial.items())],
            ]
            for power, polynomial in sorted(expected.items())
        ],
        separators=(",", ":"),
    ).encode()
    structural_zero = not certificates
    leading_cancellation = (
        minimum_term_valuation is not None
        and minimum_term_valuation <= cutoff
        and (first is None or first > minimum_term_valuation)
    )
    proof = {
        "schema": "vakint.q0-factorwise-exact-scalar-jet.v1",
        "status": "PASS",
        "cutoff_inclusive": cutoff,
        "remainder_order_at_least": cutoff + 1,
        "remainder_scope": "scalar coefficient of one formal master-product; masters are not expanded",
        "normalization": "one exp(-5*gamma*ep) per nonzero rational term; never applied again",
        "normalization_depth_required": max(0, cutoff - minimum_term_valuation)
        if minimum_term_valuation is not None
        else 0,
        "minimum_term_valuation": minimum_term_valuation,
        "observed_valuation_through_cutoff": first,
        "valuation_lower_bound": first if first is not None else cutoff + 1,
        "identically_zero_proven": structural_zero,
        "leading_cancellation_detected": leading_cancellation,
        "source_Add_terms": len(source_terms),
        "expanded_rational_terms": len(certificates),
        "native_compared_terms": eligible_terms,
        "structurally_zero_terms": exact_zero_terms,
        "term_certificates": certificates,
        "polynomial_catalog": [
            context.polynomial_catalog[key]
            for key in sorted(context.polynomial_catalog)
        ],
        "symbolica_series_absolute_order": {
            "numerator": order_n,
            "denominator": order_d,
        },
        "rational_jet": [
            [power, str(value)] for power, value in sorted(rational_total.items())
        ],
        "exact_normalized_coefficient_sha256": hashlib.sha256(encoded).hexdigest(),
        "exact_normalized_nonzero_coefficients": len(expected),
        "whole_scalar_series_calls": 0,
        "global_GCD_or_cancel_calls": 0,
        "rational_polynomial_constructions": 0,
        "numerical_evaluation_performed": False,
        "proof": "individually verified polynomial factors; exact unit P-Q*jet truncated-ring residual; independent native unit-factor Series and normalization; finite sum of coefficientwise remainders",
    }
    return jet, proof


def probe(core):
    """Exact synthetic probes; called only by the separately guarded toy job."""
    ep = core.S("q0_factorwise_probe::ep")
    gamma = core.E("symbolica::γ")
    one, zero = core.E("1"), core.E("0")
    n = core.E("𝑒") ** (-5 * gamma * ep)
    cases = {
        "geometric": n / (one - ep),
        "sixth_order_pole": n * ep**-6,
        "high_degree_rational": n
        * ep**-6
        * (one + 2 * ep + 3 * ep**2 + ep**120)
        / (one - ep + 2 * ep**2 + ep**100) ** 3,
        "leading_cancellation": n * (one / (ep * (one - ep)) - one / ep),
        "cancelled_through_cutoff": n
        * (one / (ep * (one - ep)) - one / ep - one / (one - ep)),
        "above_cutoff": n * ep**6,
        "structural_zero": zero,
    }
    results, jets = {}, {}
    for name, expression in cases.items():
        jet, proof = compute_scalar_jet(core, expression, ep, gamma)
        # Whole-expression Series is used ONLY on these tiny synthetic controls,
        # never on an actual saved scalar or in compute_scalar_jet.
        oracle = expression.series(ep, 0, 5, depth_is_absolute=True).to_expression()
        _require(
            (jet - oracle).expand() == zero,
            f"tiny control differs from ordinary Series: {name}",
        )
        results[name], jets[name] = proof, jet
    _require(
        results["sixth_order_pole"]["normalization_depth_required"] == 11,
        "pole depth is not eleven",
    )
    _require(
        results["leading_cancellation"]["leading_cancellation_detected"],
        "leading cancellation was not recorded",
    )
    _require(
        (jets["leading_cancellation"] - jets["geometric"]).expand() == zero,
        "leading-cancellation identity failed",
    )
    _require(
        jets["cancelled_through_cutoff"] == zero
        and results["cancelled_through_cutoff"]["valuation_lower_bound"] >= 6
        and not results["cancelled_through_cutoff"]["identically_zero_proven"],
        "finite cancellation was incorrectly promoted into an all-order zero proof",
    )
    _require(
        results["structural_zero"]["identically_zero_proven"],
        "structural zero was not retained",
    )
    _require(
        not results["above_cutoff"]["leading_cancellation_detected"]
        and not results["above_cutoff"]["identically_zero_proven"]
        and results["above_cutoff"]["valuation_lower_bound"] >= 6,
        "an above-cutoff term was incorrectly reported as cancelled or identically zero",
    )
    rejected = []
    for label, bad in (
        ("missing_normalization", one / (one - ep)),
        ("unexpected_constant", gamma * n / (one - ep)),
        ("fractional_power", n * (one + ep) ** core.E("1/2")),
    ):
        try:
            compute_scalar_jet(core, bad, ep, gamma)
        except FactorwiseSeriesError:
            rejected.append(label)
        else:
            raise FactorwiseSeriesError(
                f"unsupported tiny structure was accepted: {label}"
            )
    return {
        "status": "PASS",
        "scope": "exact synthetic factorwise controls only",
        "cases": results,
        "rejected_structures": rejected,
        "numerical_evaluations": 0,
        "cache_reads": 0,
    }
