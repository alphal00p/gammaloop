"""Lossless, offline reconstruction of saved formal-master coefficients.

Only public Symbolica Float/ComplexFloat/Expression interfaces are used.  This
module never evaluates master integrals or numerical formulas.  The exact-ratio
JSON and native Atom are authoritative; a decimal ``.expr`` is a human display.
"""

import copy
import hashlib
import json
import re
from fractions import Fraction
from pathlib import Path
from typing import ClassVar


class CoefficientArtifactError(ValueError):
    """A saved coefficient artifact is missing, inconsistent, or corrupt."""


def _require(condition, message):
    if not condition:
        raise CoefficientArtifactError(message)


def _integer(value, label):
    _require(type(value) is int, f"{label} must be an integer")
    return value


def _ratio(record, label):
    try:
        numerator_text = record["numerator"]
        denominator_text = record["denominator"]
        _require(
            type(numerator_text) is str
            and re.fullmatch(r"-?[0-9]+", numerator_text) is not None,
            f"invalid integer numerator for {label}",
        )
        _require(
            type(denominator_text) is str
            and re.fullmatch(r"[0-9]+", denominator_text) is not None,
            f"invalid integer denominator for {label}",
        )
        numerator = int(numerator_text)
        denominator = int(denominator_text)
    except (KeyError, TypeError, ValueError) as error:
        raise CoefficientArtifactError(f"invalid exact ratio for {label}") from error
    _require(denominator > 0, f"nonpositive denominator for {label}")
    return Fraction(numerator, denominator)


def _component(api, record, requested, label):
    center = _ratio(record["center_exact_ratio"], label + " center")
    radius = _ratio(
        record["absolute_error_exact_ratio_including_export_rounding"],
        label + " error radius",
    )
    _require(radius >= 0, f"negative error radius for {label}")
    bits = _integer(record["float_binary_precision"], label + " binary precision")
    _require(bits > 0, f"nonpositive binary precision for {label}")
    exact_zero = record["exact_zero"]
    _require(type(exact_zero) is bool, f"invalid exact-zero flag for {label}")
    if exact_zero:
        _require(center == 0 and radius == 0, f"invalid exact zero for {label}")
    else:
        achieved = _integer(
            record["exported_reliable_decimal_digits_at_least"],
            label + " achieved precision",
        )
        tracked = _integer(
            record["tracked_significant_digits"], label + " tracked precision"
        )
        _require(
            achieved >= requested and tracked >= requested,
            f"insufficient recorded precision for {label}",
        )
        _require(
            abs(center) > radius and radius * 10**requested <= abs(center) - radius,
            f"error radius does not establish requested digits for {label}",
        )
    # No decimal formatting, native float conversion, or precision inference.
    value = api.Float.from_ratio(center.numerator, center.denominator, precision=bits)
    _require(value.precision == bits, f"binary precision changed for {label}")
    _require(
        Fraction(*value.as_integer_ratio()) == center,
        f"center is not exactly representable at its recorded bits: {label}",
    )
    return value


def reconstruct_coefficient_table(table, *, api):
    """Return an Atom without evaluating any formula or formal master.

    ``api`` is the already imported canonical Symbolica module.  Passing the
    existing native host in development exercises exactly the same public API;
    this function imports no second Symbolica kernel.  The caller retains
    ``table``: its component errors, error quality, warnings and provenance are
    preserved as metadata and are not silently encoded as extra accurate digits.
    The exported centers are reconstructed exactly; discarded working guard
    digits are not recovered.
    """
    try:
        _require(
            table["schema"] == "vakint.q0-formal-master-product-coefficients.v1",
            "unsupported coefficient-table schema",
        )
        _require(
            table["unknown_master_values_evaluated"] is False,
            "expected formal, unevaluated masters",
        )
        _require(
            table["normalization_applied_again"] is False,
            "unexpected additional normalization",
        )
        cutoff = table["cutoff"]
        _require(
            cutoff["included_maximum_power"] == 5
            and cutoff["remainder_starts_at"] == 6,
            "unsupported coefficientwise Laurent cutoff",
        )
        expression = api.E("0")
        epsilon = api.E("vakint::ε")
        coupling = api.E(table["common_coupling_prefactor"])
        _require(
            coupling == api.E("UFO::GC_11") ** 10,
            "unexpected common coupling prefactor",
        )
        products_seen = set()
        for row in table["products"]:
            product_key = row["master_product"]
            _require(
                product_key not in products_seen, "duplicate formal master product"
            )
            products_seen.add(product_key)
            _require(
                row["remainder"]["coefficientwise"] is True
                and row["remainder"]["epsilon_power"] == 6,
                f"incorrect remainder for {product_key}",
            )
            requested = _integer(
                row["requested_decimal_precision"], "requested decimal precision"
            )
            _require(requested > 0, "nonpositive requested precision")
            coefficient = api.E("0")
            powers_seen = set()
            for number in row["coefficients"]:
                power = _integer(number["epsilon_power"], "epsilon power")
                _require(
                    power <= 5 and power not in powers_seen,
                    f"invalid or duplicate epsilon power for {product_key}",
                )
                powers_seen.add(power)
                label = f"{product_key}, epsilon^{power}"
                real = _component(api, number["real"], requested, label + " real")
                imaginary = _component(
                    api, number["imaginary"], requested, label + " imaginary"
                )
                value = api.ComplexFloat(real, imaginary)
                _require(
                    value.real.precision == real.precision
                    and value.imag.precision == imaginary.precision,
                    f"complex construction changed component precision: {label}",
                )
                coefficient += api.Expression.num(value) * epsilon**power
            expression += api.E(product_key) * coefficient
        return coupling * expression
    except (KeyError, TypeError, OverflowError) as error:
        raise CoefficientArtifactError("malformed coefficient table") from error


def require_raw_equal(actual, expected, label):
    """Compare component storage precision as well as values in the current host."""
    _require(bytes(actual.__getstate__()) == bytes(expected.__getstate__()), label)


class SavedQ0Result:
    """Authenticated saved expressions and their coefficient-precision sidecar.

    Supply the application's existing Symbolica API. No host is imported here,
    and loading never reduces or numerically evaluates an integral. Native and
    JSON formats are selected explicitly; an error never selects another format.
    """

    _NAMES: ClassVar[dict[str, str]] = {
        "retained": "q0-retained",
        "symbolic": "q0-coefficients-exact-through-epsilon5",
        "numerical": "q0-coefficients-300digits-through-epsilon5",
    }
    _BUNDLED_PROVENANCE_SHA256 = (
        "941d5a2e888557e8cabd8e8ab118935e95631be0a16110bd28a41e49436db36a"
    )

    def __init__(self, directory, api, provenance, table):
        self._directory = directory
        self._api = api
        self._provenance = provenance
        self._table = table
        self._json_expression = reconstruct_coefficient_table(table, api=api)

    @classmethod
    def load(cls, directory=None, *, api, expected_provenance_sha256=None):
        """Validate local payloads and construct one owned, offline inspection."""
        bundled = directory is None
        directory = (
            Path(__file__).resolve().parents[1] / "data/five_loop_q0"
            if bundled
            else Path(directory)
        )
        directory = directory.resolve()
        try:
            payload = (directory / "provenance.json").read_bytes()
            expected = (
                cls._BUNDLED_PROVENANCE_SHA256
                if bundled and expected_provenance_sha256 is None
                else expected_provenance_sha256
            )
            if expected is not None:
                _require(
                    hashlib.sha256(payload).hexdigest() == expected,
                    "provenance.json checksum mismatch",
                )
            provenance = json.loads(payload)
            _require(
                provenance["schema"] == "vakint.q0-coefficient-artifact-provenance.v1",
                "unsupported artifact provenance schema",
            )
            _require(
                provenance["tracked_requested_decimal_digits"] == 300
                and provenance["fixed_guard_decimal_digits"] == 32
                and provenance["automatic_retries"] == 0,
                "unsupported precision contract",
            )
            required = {
                base + extension
                for base in cls._NAMES.values()
                for extension in (".atom", ".expr")
            } | {"coefficient-table.json"}
            manifest = provenance["artifacts"]
            _require(
                required.issubset(manifest), "required artifact missing from manifest"
            )
            for name, reference in manifest.items():
                cls._checked_payload(directory, name, reference)
            table = json.loads(
                cls._checked_payload(
                    directory,
                    "coefficient-table.json",
                    manifest["coefficient-table.json"],
                )
            )
            _require(
                all(
                    row["requested_decimal_precision"] == 300
                    and row["fixed_guard_decimal_digits"] == 32
                    and row["automatic_retries"] == 0
                    for row in table["products"]
                ),
                "coefficient precision contract differs from provenance",
            )
            result = cls(directory, api, provenance, table)
            summary = result.summary()
            _require(
                summary["master_products"] == provenance["master_product_count"]
                and summary["master_symbols"] == provenance["individual_master_count"],
                "master inventory differs from provenance",
            )
            return result
        except (OSError, json.JSONDecodeError, KeyError, TypeError) as error:
            raise CoefficientArtifactError("cannot load saved q0 artifacts") from error

    @staticmethod
    def _checked_payload(directory, name, reference):
        _require(
            type(name) is str
            and Path(name).name == name
            and name not in ("", ".", ".."),
            "invalid artifact filename",
        )
        path = directory / name
        _require(
            path.is_file() and not path.is_symlink(),
            f"regular artifact required: {name}",
        )
        payload = path.read_bytes()
        _require(
            type(reference["bytes"]) is int and len(payload) == reference["bytes"],
            f"{name} length mismatch",
        )
        _require(
            hashlib.sha256(payload).hexdigest() == reference["sha256"],
            f"{name} checksum mismatch",
        )
        return payload

    def _payload(self, name):
        return self._checked_payload(
            self._directory, name, self._provenance["artifacts"][name]
        )

    @property
    def directory(self):
        return self._directory

    @property
    def provenance(self):
        return copy.deepcopy(self._provenance)

    @property
    def coefficient_table(self):
        """Independent metadata copy: errors and warnings are never remeasured."""
        return copy.deepcopy(self._table)

    def expression(self, variant="retained", *, format="native"):
        """Load native/exact text, or reconstruct numerical values from ratios."""
        _require(variant in self._NAMES, f"unknown expression variant: {variant}")
        if format == "json":
            _require(
                variant == "numerical",
                "JSON contains numerical scalar coefficients only",
            )
            return self._json_expression
        _require(format in ("native", "text"), f"unknown expression format: {format}")
        _require(
            not (format == "text" and variant == "numerical"),
            "numerical .expr is a human display; choose native or JSON for an Atom",
        )
        name = self._NAMES[variant] + (".atom" if format == "native" else ".expr")
        payload = self._payload(name)
        if format == "text":
            return self._api.E(payload.decode("utf-8"))
        expression = self._api.Expression.load(str(self._directory / name))
        self._payload(name)
        return expression

    def display_text(self, variant="numerical"):
        """Return recorded text; numerical decimal text is not a Float transport."""
        _require(variant in self._NAMES, f"unknown expression variant: {variant}")
        return self._payload(self._NAMES[variant] + ".expr").decode("utf-8")

    def summary(self):
        products = self._table["products"]
        coefficients = [number for row in products for number in row["coefficients"]]
        components = [
            number[part] for number in coefficients for part in ("real", "imaginary")
        ]
        nonzero = [component for component in components if not component["exact_zero"]]
        identities = set(
            re.findall(
                r"vakint::\{\}::RustRedMaster\([^()]*\)",
                " ".join(row["master_product"] for row in products),
            )
        )
        return {
            "master_symbols": len(identities),
            "master_products": len(products),
            "scalar_coefficients": len(coefficients),
            "retained_native_file_bytes": self._provenance["artifacts"][
                "q0-retained.atom"
            ]["bytes"],
            "retained_text_bytes": self._provenance["artifacts"]["q0-retained.expr"][
                "bytes"
            ],
            "epsilon_minimum": min(number["epsilon_power"] for number in coefficients),
            "epsilon_maximum": max(number["epsilon_power"] for number in coefficients),
            "remainder_starts_at": self._table["cutoff"]["remainder_starts_at"],
            "remainder_scope": "each scalar coefficient of a formal master product",
            "minimum_exported_nonzero_decimal_digits": min(
                component["exported_reliable_decimal_digits_at_least"]
                for component in nonzero
            ),
            "minimum_tracked_nonzero_decimal_digits": min(
                component["tracked_significant_digits"] for component in nonzero
            ),
            "component_storage_bits": sorted(
                {component["float_binary_precision"] for component in components}
            ),
            "precision_warnings": sum(len(row["warnings"]) for row in products),
            "master_values_evaluated": False,
        }

    def verify(self):
        """Check files, exact texts and numerical precision without reevaluation."""
        for name in self._provenance["artifacts"]:
            self._payload(name)
        require_raw_equal(
            self.expression("retained"),
            self.expression("retained", format="text"),
            "exact retained native/text mismatch",
        )
        require_raw_equal(
            self.expression("symbolic"),
            self.expression("symbolic", format="text"),
            "exact symbolic coefficient native/text mismatch",
        )
        native = self.expression("numerical")
        require_raw_equal(
            native,
            self.expression("numerical", format="json"),
            "lossless numerical JSON/native mismatch",
        )
        retained_name = self._NAMES["retained"] + ".expr"
        _require(
            hashlib.sha256(self._payload(retained_name)).hexdigest()
            == self._provenance["input"]["sha256"],
            "retained exact source differs from authenticated input",
        )
        coupling = self._api.E("UFO::GC_11")
        terms = native.coefficient_list(coupling)
        _require(
            len(terms) == 1 and terms[0][0] == coupling**10,
            "common coupling is not GC_11^10",
        )
        return {"status": "PASS", **self.summary(), "numerical_evaluations": 0}

    def __repr__(self):
        summary = self.summary()
        return (
            f"SavedQ0Result({summary['master_symbols']} masters, "
            f"{summary['master_products']} products, "
            f"{summary['scalar_coefficients']} scalar coefficients)"
        )
