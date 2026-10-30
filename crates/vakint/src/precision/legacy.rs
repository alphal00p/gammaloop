//! Finite FORM table entries and input floats remain identifiable until complete coefficients
//! reach the shared tracked evaluator. The leaves contain exact metadata and
//! need neither a process-global value map nor a mutable engine registry.

use std::collections::BTreeMap;

use sha2::{Digest, Sha256};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    coefficient::CoefficientView,
    domains::{
        float::{Float, SingleFloat},
        rational::Rational,
    },
    id::Replacement,
    symbol,
};

use crate::{
    VakintError, get_integer_from_atom,
    inspection::FiniteSource,
    utils::vakint_macros::{vk_parse, vk_symbol},
};

fn failure(detail: impl Into<String>) -> VakintError {
    VakintError::NumericalPrecision(detail.into())
}

pub(crate) fn is_source(part: AtomView) -> bool {
    matches!(part, AtomView::Fun(function) if function.get_symbol() == vk_symbol!("NumericalSource"))
}

pub(crate) fn coefficient_power(monomial: AtomView, epsilon: Symbol) -> Result<i64, VakintError> {
    match monomial {
        AtomView::Var(variable) if variable.get_symbol() == epsilon => Ok(1),
        AtomView::Pow(power) if power.get_base() == Atom::var(epsilon).as_view() => {
            get_integer_from_atom(power.get_exp())
                .ok_or_else(|| failure("noninteger FORM table epsilon order"))
        }
        _ if monomial == Atom::one().as_view() => Ok(0),
        _ => Err(failure(format!(
            "FORM numerical table has a nonmonomial epsilon factor: {monomial}"
        ))),
    }
}

fn source_hash(
    backend: &str,
    identity: &str,
    order: i64,
    component: i64,
    ordinal: usize,
    center: &Atom,
    bits: u32,
) -> String {
    let mut digest = Sha256::new();
    for field in [
        backend.to_owned(),
        identity.to_owned(),
        order.to_string(),
        component.to_string(),
        ordinal.to_string(),
        center.to_canonical_string(),
        bits.to_string(),
    ] {
        digest.update(field.as_bytes());
        digest.update([0]);
    }
    format!("{:x}", digest.finalize())
}

fn staged_float(
    value: &Float,
    source: AtomView,
    backend: &str,
    order: i64,
    component: i64,
    ordinal: usize,
) -> Atom {
    let limited;
    let value = if backend == "PySecDec" {
        limited = {
            let mut value = value.clone();
            value.set_prec(value.prec().min(53));
            value
        };
        &limited
    } else {
        value
    };
    if value.is_zero() {
        return Atom::Zero;
    }
    let center = Atom::num(value.to_rational());
    let identity = source.to_canonical_string();
    let encoded = identity
        .as_bytes()
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect::<String>();
    let digest = source_hash(
        backend,
        &identity,
        order,
        component,
        ordinal,
        &center,
        value.prec(),
    );
    let key = Atom::var(symbol!(format!(
        "vakint_sources::{backend}__{encoded}__{ordinal}__{digest}"
    )));
    vk_symbol!("NumericalSource").call((
        key,
        Atom::num(order),
        Atom::num(component),
        center,
        Atom::num(value.prec()),
    ))
}

/// Stage every finite literal in an existing Laurent table without changing
/// its exact coefficients, symbolic tails, epsilon orders, or source accuracy.
pub(crate) fn stage_table(
    source: AtomView,
    target: AtomView,
    epsilon: Symbol,
    backend: &str,
) -> Result<Atom, VakintError> {
    let mut has_float = false;
    let mut non_finite = false;
    target.visitor(&mut |part| {
        if let AtomView::Num(number) = part
            && let CoefficientView::Float(real, imaginary) = number.get_coeff_view()
        {
            has_float = true;
            non_finite |= !real.to_float().is_finite() || !imaginary.to_float().is_finite();
        }
        !non_finite && !is_source(part)
    });
    if non_finite {
        return Err(failure(format!(
            "{backend} source {source} contains a non-finite floating coefficient"
        )));
    }
    if !has_float {
        return Ok(target.to_owned());
    }
    if backend == "Input" {
        let mut ordinal = 0;
        return Ok(target.replace_map(|part, _context, out| {
            if is_source(part) {
                **out = part.to_owned();
            } else if let AtomView::Num(number) = part
                && let CoefficientView::Float(real, imaginary) = number.get_coeff_view()
            {
                let real = staged_float(&real.to_float(), source, backend, 0, 0, ordinal);
                let imaginary = staged_float(&imaginary.to_float(), source, backend, 0, 1, ordinal);
                **out = real + Atom::i() * imaginary;
                ordinal += 1;
            }
        }));
    }
    let mut result = Atom::Zero;
    for (monomial, coefficient) in target.coefficient_list::<i32>(&[Atom::var(epsilon)]) {
        let order = coefficient_power(monomial.as_view(), epsilon)?;
        if backend == "PySecDec" {
            let mut ordinal = 0;
            let staged = coefficient.replace_map(|part, _context, out| {
                if is_source(part) {
                    **out = part.to_owned();
                } else if let AtomView::Num(number) = part
                    && let CoefficientView::Float(real, imaginary) = number.get_coeff_view()
                {
                    let real = staged_float(&real.to_float(), source, backend, order, 0, ordinal);
                    let imaginary =
                        staged_float(&imaginary.to_float(), source, backend, order, 1, ordinal);
                    **out = real + Atom::i() * imaginary;
                    ordinal += 1;
                }
            });
            result += staged * monomial;
            continue;
        }
        let mut replacements = BTreeMap::new();
        let mut ordinal = 0;
        coefficient.visitor(&mut |part| {
            if let AtomView::Num(number) = part
                && let CoefficientView::Float(real, imaginary) = number.get_coeff_view()
            {
                let real = staged_float(&real.to_float(), source, backend, order, 0, ordinal);
                let imaginary =
                    staged_float(&imaginary.to_float(), source, backend, order, 1, ordinal);
                replacements
                    .entry(part.to_owned())
                    .or_insert_with(|| real + Atom::i() * imaginary);
                ordinal += 1;
            }
            !is_source(part)
        });
        let replacements = replacements
            .into_iter()
            .map(|(from, to)| Replacement::new(from, to))
            .collect::<Vec<_>>();
        result += coefficient.replace_multiple(&replacements) * monomial;
    }
    Ok(result)
}

struct DecodedSource {
    original: Atom,
    backend: String,
    ordinal: usize,
    order: i64,
    component: i64,
    center: Rational,
    bits: u32,
    digest: String,
}

fn decode(part: AtomView) -> Result<Option<DecodedSource>, VakintError> {
    let AtomView::Fun(function) = part else {
        return Ok(None);
    };
    if function.get_symbol() != vk_symbol!("NumericalSource") {
        return Ok(None);
    }
    let arguments = function.iter().collect::<Vec<_>>();
    let [identity, order, component, center, bits] = arguments.as_slice() else {
        return Err(failure(
            "invalid numerical source: expected five metadata arguments",
        ));
    };
    let AtomView::Var(identity) = identity else {
        return Err(failure("invalid numerical source identity"));
    };
    if identity.get_symbol().get_namespace() != "vakint_sources" {
        return Err(failure("invalid numerical source namespace"));
    }
    let identity_symbol = identity.get_symbol();
    let name = identity_symbol.get_name();
    let name = name.rsplit("::").next().unwrap();
    let fields = name.split("__").collect::<Vec<_>>();
    let [backend, encoded, ordinal, digest] = fields.as_slice() else {
        return Err(failure("invalid numerical source identity fields"));
    };
    if !["FMFT", "MATAD", "AlphaLoop", "Input", "PySecDec"].contains(backend)
        || encoded.len() % 2 != 0
    {
        return Err(failure("unsupported numerical table source"));
    }
    // Byte-indexed pairs are safe only after validating the complete encoding.
    // An even UTF-8 byte length alone does not imply character boundaries.
    if !encoded.bytes().all(|byte| byte.is_ascii_hexdigit()) {
        return Err(failure("invalid encoded numerical source"));
    }
    let bytes = (0..encoded.len())
        .step_by(2)
        .map(|start| u8::from_str_radix(&encoded[start..start + 2], 16))
        .collect::<Result<Vec<_>, _>>()
        .map_err(|_| failure("invalid encoded numerical source"))?;
    let original_text =
        String::from_utf8(bytes).map_err(|_| failure("invalid numerical source UTF-8"))?;
    let original = vk_parse!(&original_text).map_err(failure)?;
    let ordinal = ordinal
        .parse::<usize>()
        .map_err(|_| failure("invalid numerical source literal index"))?;
    let order = get_integer_from_atom(*order)
        .ok_or_else(|| failure("invalid numerical source epsilon order"))?;
    let component = get_integer_from_atom(*component)
        .filter(|component| [0, 1].contains(component))
        .ok_or_else(|| failure("invalid numerical source component"))?;
    let bits = get_integer_from_atom(*bits)
        .and_then(|bits| u32::try_from(bits).ok())
        .filter(|bits| *bits >= 1)
        .ok_or_else(|| {
            failure(format!(
                "invalid numerical source precision: {backend} {original}, epsilon^{order}, component {component}, bits {bits}, center {center}"
            ))
        })?;
    let expected = source_hash(
        backend,
        &original_text,
        order,
        component,
        ordinal,
        &center.to_owned(),
        bits,
    );
    if *digest != expected {
        return Err(failure(
            "numerical source identity does not match its payload",
        ));
    }
    let center = Rational::try_from(*center)
        .map_err(|_| failure("numerical source center is not an exact real rational"))?;
    Ok(Some(DecodedSource {
        original,
        backend: (*backend).into(),
        ordinal,
        order,
        component,
        center,
        bits,
        digest: (*digest).into(),
    }))
}

/// Recover a named primitive so an authenticated active registry may supersede
/// the historical fallback. Identity metadata never contains the original Atom.
pub(crate) fn original_source(part: AtomView) -> Result<Option<Atom>, VakintError> {
    Ok(decode(part)?.map(|source| source.original))
}

/// The authenticated table position of a leaf: which literal of which table
/// entry, at which epsilon order and complex component, it denotes.
pub(crate) struct SourceIdentity {
    pub backend: String,
    pub original: Atom,
    pub ordinal: usize,
    pub order: i64,
    pub component: i64,
}

pub(crate) fn source_identity(part: AtomView) -> Result<Option<SourceIdentity>, VakintError> {
    Ok(decode(part)?.map(|source| SourceIdentity {
        backend: source.backend,
        original: source.original,
        ordinal: source.ordinal,
        order: source.order,
        component: source.component,
    }))
}

/// Runtime literals retain their binary accuracy, including less than one
/// complete decimal digit after cancellation in the caller's Float algebra.
pub(super) fn input_value(part: AtomView) -> Result<Option<Float>, VakintError> {
    let Some(source) = decode(part)? else {
        return Ok(None);
    };
    if source.backend != "Input" {
        return Ok(None);
    }
    let value = source.center.to_multi_prec_float(source.bits);
    // Staging omits literal zeros, and every retained center must be exactly
    // representable at its declared precision. Never round forged metadata or
    // reinterpret an uncertain zero as a declared exact-zero input.
    if value.is_zero() || value.to_rational() != source.center {
        return Err(failure(format!(
            "runtime numerical source {} does not encode a nonzero Float at {} bits",
            source.original, source.bits
        )));
    }
    Ok(Some(value))
}

pub(crate) fn collect_sources(
    expression: AtomView,
) -> Result<BTreeMap<Atom, FiniteSource>, VakintError> {
    let mut sources = BTreeMap::new();
    let mut error = None;
    expression.visitor(&mut |part| {
        if error.is_some() {
            return false;
        }
        match decode(part) {
            Ok(Some(source)) => {
                let mut significant_digits =
                    ((u128::from(source.bits.saturating_sub(1)) * 301_029_995_663_981_195)
                        / 1_000_000_000_000_000_000) as u32;
                // The archived calculations were independently checked through
                // 20,000 decimal digits. Their internal guard digits are not
                // additional validated source precision.
                if ["MATAD", "AlphaLoop"].contains(&source.backend.as_str()) {
                    significant_digits = significant_digits.min(20_000);
                }
                let value = source.center.to_multi_prec_float(source.bits).to_string();
                let master = format!("{} {}", source.backend, source.original.to_plain_string());
                sources.insert(
                    part.to_owned(),
                    FiniteSource {
                        value,
                        significant_digits,
                        identity: format!(
                            "{master} epsilon^{} literal {} ({})",
                            source.order,
                            source.ordinal,
                            if source.component == 0 {
                                "real"
                            } else {
                                "imaginary"
                            }
                        ),
                        master_id: (!["Input", "PySecDec"].contains(&source.backend.as_str()))
                            .then_some(master),
                        epsilon_power: (source.backend != "Input").then_some(source.order),
                        data_version: match source.backend.as_str() {
                            "Input" => "runtime inputs".into(),
                            "PySecDec" => "PySecDec estimates".into(),
                            _ => format!("built-in {} tables", source.backend),
                        },
                        source_sha256: source.digest,
                        absolute_error: None,
                        error_is_bounded: false,
                    },
                );
                false
            }
            Ok(None) => true,
            Err(failure) => {
                error = Some(failure);
                false
            }
        }
    });
    if let Some(error) = error {
        return Err(error);
    }
    Ok(sources)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        EvaluationMethod, EvaluationOrder, FMFTOptions, LoopNormalizationFactor, MATADOptions,
        Vakint, VakintSettings,
    };
    use std::collections::HashMap;
    use symbolica::domains::float::Real;

    fn staged(name: &str, value: &str) -> Atom {
        Vakint::initialize_vakint_symbols();
        stage_table(
            vk_parse!(name).unwrap().as_view(),
            vk_parse!(value).unwrap().as_view(),
            vk_symbol!("ep"),
            "FMFT",
        )
        .unwrap()
    }

    #[test]
    fn table_leaves_preserve_orders_source_caps_and_exact_identity() {
        let source = staged("test_master", "1.2345`30/ep^2+ep*(2.3456`40+3.4567`25*𝑖)");
        source.visitor(&mut |part| {
            if let AtomView::Num(number) = part {
                assert!(!matches!(
                    number.get_coeff_view(),
                    CoefficientView::Float(_, _)
                ));
            }
            true
        });
        let sources = collect_sources(source.as_view()).unwrap();
        assert_eq!(sources.len(), 3);
        assert!(
            sources
                .values()
                .any(|source| source.epsilon_power == Some(-2))
        );
        assert!(
            sources
                .values()
                .any(|source| source.epsilon_power == Some(1))
        );
        assert!(
            sources
                .values()
                .all(|source| (24..=40).contains(&source.significant_digits))
        );
        let replaced = source
            .replace(vk_parse!("test_master").unwrap().to_pattern())
            .with(Atom::num(17));
        assert_eq!(
            source, replaced,
            "table substitutions cannot rewrite source identity metadata"
        );
        let cancelled = (Atom::num(3) * source.clone() - source * Atom::num(3)).expand();
        assert!(cancelled.is_zero());
        assert!(collect_sources(cancelled.as_view()).unwrap().is_empty());
    }

    #[test]
    fn equal_midpoints_from_independent_master_rows_remain_uncertain() {
        let left = staged("first_master", "1.23456789`30");
        let right = staged("second_master", "1.23456789`30");
        let difference = left - right;
        assert!(!difference.is_zero());
        assert_eq!(collect_sources(difference.as_view()).unwrap().len(), 2);
        let settings = VakintSettings {
            run_time_decimal_precision: 50,
            ..Default::default()
        };
        let error =
            super::super::materialize(difference.as_view(), &settings, &HashMap::new(), None)
                .unwrap_err();
        assert!(error.to_string().contains("no reliable"), "{error}");
    }

    #[test]
    fn corrupted_source_metadata_does_not_silently_claim_more_digits() {
        let source = staged("test_master", "1.2345`30");
        assert_eq!(
            original_source(source.as_view()).unwrap(),
            Some(vk_parse!("test_master").unwrap())
        );
        let AtomView::Fun(function) = source.as_view() else {
            panic!("expected source leaf")
        };
        let mut args = function
            .iter()
            .map(|argument| argument.to_owned())
            .collect::<Vec<_>>();
        args[4] = Atom::num(10_000);
        let corrupted = function.get_symbol().call(&args);
        assert!(
            collect_sources(corrupted.as_view())
                .unwrap_err()
                .to_string()
                .contains("does not match")
        );
    }

    #[test]
    fn high_precision_table_centers_round_trip_without_f64() {
        let source = staged("high_precision_master", "1.234567890123456789`20000");
        let sources = collect_sources(source.as_view()).unwrap();
        let data = sources.values().next().unwrap();
        assert!(data.significant_digits >= 19_999);
        assert!(data.value.starts_with("1.234567890123456789"));
        assert_eq!(
            original_source(source.as_view()).unwrap(),
            Some(vk_parse!("high_precision_master").unwrap())
        );
    }

    #[test]
    fn input_source_metadata_survives_form_scalar_algebra() {
        let vakint = Vakint::new().unwrap();
        let settings = VakintSettings::default();
        let source = stage_table(
            vk_parse!("input_coefficient").unwrap().as_view(),
            vk_parse!("1.23456789`30").unwrap().as_view(),
            vk_symbol!("ep"),
            "Input",
        )
        .unwrap();
        let (header, body, indices) = vakint
            .prepare_expression_for_form(&settings, source.clone(), false, &[])
            .unwrap();
        assert!(header.contains("NumericalSource"));
        let script = format!(
            "Off Statistics;\n{header}\nLocal F = 3*({body})-2*({body});\n.sort\n#write <out.txt> \"%E\", F\n.end\n"
        );
        let output = vakint
            .run_form(
                &settings,
                &[],
                ("source.frm".into(), script),
                Vec::new(),
                true,
                settings.temporary_directory.clone(),
            )
            .unwrap();
        let restored = vakint
            .process_form_output(&settings, output, indices, BTreeMap::new())
            .unwrap();
        assert_eq!(source, restored);
        assert_eq!(collect_sources(source.as_view()).unwrap().len(), 1);
    }

    #[test]
    fn input_float_sources_preserve_nonpolynomial_epsilon_factors() {
        Vakint::initialize_vakint_symbols();
        let input = vk_parse!("1.2345`30/(1-2*ep)").unwrap();
        let staged = stage_table(
            vk_parse!("input_coefficient").unwrap().as_view(),
            input.as_view(),
            vk_symbol!("ep"),
            "Input",
        )
        .unwrap();
        assert_eq!(collect_sources(staged.as_view()).unwrap().len(), 1);
        assert!(staged.contains_symbol(vk_symbol!("ep")));
        let exact = vk_parse!("1/(1-2*ep)").unwrap();
        assert_eq!(
            stage_table(
                vk_parse!("input_coefficient").unwrap().as_view(),
                exact.as_view(),
                vk_symbol!("ep"),
                "Input"
            )
            .unwrap(),
            exact
        );
    }

    #[test]
    fn input_literal_occurrences_are_independent_after_symbol_substitution() {
        Vakint::initialize_vakint_symbols();
        for input in ["1.2345`30*x+1.2345`30*y", "1.2345`30*x-1.2345`30*y"] {
            let original = vk_parse!(input).unwrap();
            let staged = stage_table(
                vk_parse!("input_coefficient").unwrap().as_view(),
                original.as_view(),
                vk_symbol!("ep"),
                "Input",
            )
            .unwrap();
            assert_eq!(collect_sources(staged.as_view()).unwrap().len(), 2);
            let y_value = if input.contains("+1.2345") {
                -Atom::one()
            } else {
                Atom::one()
            };
            let cancelled_centers = staged.replace_multiple(&[
                Replacement::new(vk_parse!("x").unwrap(), Atom::one()),
                Replacement::new(vk_parse!("y").unwrap(), y_value),
            ]);
            assert!(!cancelled_centers.is_zero());
            let settings = VakintSettings {
                run_time_decimal_precision: 40,
                ..Default::default()
            };
            let error = super::super::materialize(
                cancelled_centers.as_view(),
                &settings,
                &HashMap::new(),
                None,
            )
            .unwrap_err();
            assert!(error.to_string().contains("no reliable"), "{error}");
            assert_eq!(
                stage_table(
                    vk_parse!("another_input").unwrap().as_view(),
                    staged.as_view(),
                    vk_symbol!("ep"),
                    "Input"
                )
                .unwrap(),
                staged,
                "existing identities must survive repeated entry boundaries"
            );
        }
    }

    #[test]
    fn archived_guard_digits_are_not_extra_validated_precision() {
        Vakint::initialize_vakint_symbols();
        let value = vk_parse!("1.25`21000").unwrap();
        for backend in ["MATAD", "AlphaLoop", "Input"] {
            let staged = stage_table(
                vk_parse!("precision_source").unwrap().as_view(),
                value.as_view(),
                vk_symbol!("ep"),
                backend,
            )
            .unwrap();
            let sources = collect_sources(staged.as_view()).unwrap();
            let digits = sources.values().next().unwrap().significant_digits;
            if backend == "Input" {
                assert!(digits > 20_000);
            } else {
                assert_eq!(digits, 20_000);
            }
        }
    }

    #[test]
    fn runtime_inputs_preserve_one_to_four_binary_bits() {
        Vakint::initialize_vakint_symbols();
        // The actual ghost_2loop/d4 residue was imaginary 2^-331 at two bits.
        let center = Rational::from((1, 2)).pow(331);
        let settings = VakintSettings {
            run_time_decimal_precision: 50,
            ..Default::default()
        };
        for bits in 1..=4 {
            let value = center.to_multi_prec_float(bits);
            let stage = |name: &str| {
                stage_table(
                    vk_parse!(name).unwrap().as_view(),
                    (Atom::i() * Atom::num(value.clone())).as_view(),
                    vk_symbol!("ep"),
                    "Input",
                )
                .unwrap()
            };
            let left = stage("low_bit_left");
            let right = stage("low_bit_right");
            let sources = collect_sources(left.as_view()).unwrap();
            let (leaf, source) = sources.first_key_value().unwrap();
            assert_eq!(source.significant_digits, 0);
            let reconstructed = input_value(leaf.as_view()).unwrap().unwrap();
            assert_eq!(reconstructed.prec(), bits);
            assert_eq!(reconstructed.to_rational(), center);

            let evaluate = |expression: &Atom| {
                super::super::evaluate_coefficient(
                    expression.as_view(),
                    &BTreeMap::new(),
                    &HashMap::new(),
                    &settings,
                )
                .unwrap()
            };
            let (small, warnings) = evaluate(&left);
            assert!(small.precision().real().exact_zero());
            assert_eq!(small.value().im.to_rational(), center);
            let binary_radius = &center * &Rational::from((1, 2)).pow(u64::from(bits - 1));
            assert!(small.precision().imaginary().absolute_error().to_rational() >= binary_radius);
            assert!(small.to_atom().is_err());
            assert!(
                warnings
                    .iter()
                    .any(|warning| warning.sources.iter().any(|source| {
                        source.available_decimal_precision == Some(0)
                            && source.identity.contains("low_bit_left")
                    }))
            );

            let (combined, _) = evaluate(&(Atom::i() + &left));
            assert!(combined.to_atom().is_ok());
            let discrepancy =
                (combined.value().im.to_rational() - (&center + &Rational::from(1))).abs();
            assert!(
                discrepancy
                    <= combined
                        .precision()
                        .imaginary()
                        .absolute_error()
                        .to_rational()
            );

            let (cancelled, _) = evaluate(&(&left - right));
            assert!(cancelled.value().im.is_zero());
            assert!(!cancelled.precision().imaginary().exact_zero());
            assert!(cancelled.to_atom().is_err());
            let identical_source = left.clone();
            let (exact_zero, warnings) = evaluate(&(left - identical_source));
            assert_eq!(exact_zero.to_atom().unwrap(), Atom::Zero);
            assert!(warnings.is_empty());
        }
    }

    #[test]
    fn runtime_input_metadata_requires_a_representable_nonzero_float() {
        Vakint::initialize_vakint_symbols();
        let original = vk_parse!("invalid_runtime_input")
            .unwrap()
            .to_canonical_string();
        let encoded = original
            .as_bytes()
            .iter()
            .map(|byte| format!("{byte:02x}"))
            .collect::<String>();
        for (center, bits) in [
            (Rational::from(1), 0),
            (Rational::from((1, 3)), 2),
            (Rational::from(0), 2),
        ] {
            let center = Atom::num(center);
            let digest = source_hash("Input", &original, 0, 0, 0, &center, bits);
            let key = Atom::var(symbol!(format!(
                "vakint_sources::Input__{encoded}__0__{digest}"
            )));
            let source = vk_symbol!("NumericalSource").call((
                key,
                Atom::Zero,
                Atom::Zero,
                center,
                Atom::num(bits),
            ));
            assert!(input_value(source.as_view()).is_err());
        }
    }

    #[test]
    fn malformed_unicode_source_encoding_returns_an_error() {
        Vakint::initialize_vakint_symbols();
        for encoded in ["a€", "€a", "λ", "ag"] {
            let identity = Atom::var(symbol!(format!(
                "vakint_sources::Input__{encoded}__0__invalid"
            )));
            let source = vk_symbol!("NumericalSource").call((
                identity,
                Atom::Zero,
                Atom::Zero,
                Atom::one(),
                Atom::num(53),
            ));
            let error = collect_sources(source.as_view()).unwrap_err();
            assert!(
                error
                    .to_string()
                    .contains("invalid encoded numerical source")
            );
        }
    }

    #[test]
    fn runtime_binary_accuracy_matches_parameters_without_relaxing_table_caps() {
        Vakint::initialize_vakint_symbols();
        let value = Float::parse("1.23456789", Some(53)).unwrap();
        let variable = vk_parse!("binary_parameter").unwrap();
        let source = stage_table(
            variable.as_view(),
            Atom::num(value.clone()).as_view(),
            vk_symbol!("ep"),
            "Input",
        )
        .unwrap();
        let settings = VakintSettings {
            run_time_decimal_precision: 80,
            ..Default::default()
        };
        let (input, warnings) = super::super::evaluate_coefficient(
            source.as_view(),
            &BTreeMap::new(),
            &HashMap::new(),
            &settings,
        )
        .unwrap();
        let parameters = HashMap::from([(
            variable,
            symbolica::domains::float::Complex::new(value, Float::new(53)),
        )]);
        let (parameter, _) = super::super::evaluate_coefficient(
            parameters.keys().next().unwrap().as_view(),
            &BTreeMap::new(),
            &parameters,
            &settings,
        )
        .unwrap();
        assert_eq!(input.value(), parameter.value());
        assert_eq!(
            input.precision().real().absolute_error(),
            parameter.precision().real().absolute_error()
        );
        assert!(warnings.iter().any(|warning| {
            warning
                .sources
                .iter()
                .any(|source| source.available_decimal_precision == Some(15))
        }));

        let table = staged("fixed_decimal_table", "1.23456789`30");
        assert!(input_value(table.as_view()).unwrap().is_none());
        // The native 30-digit literal stores 100 bits. The unchanged table
        // policy advertises floor((100 - 1) * log10(2)) = 29 complete digits.
        assert_eq!(
            collect_sources(table.as_view())
                .unwrap()
                .values()
                .next()
                .unwrap()
                .significant_digits,
            29
        );
        let (coefficient, warnings) = super::super::evaluate_coefficient(
            table.as_view(),
            &BTreeMap::new(),
            &HashMap::new(),
            &settings,
        )
        .unwrap();
        assert!(coefficient.precision().real().significant_digits().unwrap() <= 29);
        assert!(warnings.iter().any(|warning| {
            warning
                .sources
                .iter()
                .any(|source| source.available_decimal_precision == Some(29))
        }));
    }

    #[test]
    fn pysecdec_outputs_have_f64_caps_and_independent_occurrences() {
        Vakint::initialize_vakint_symbols();
        let input = vk_parse!("1.23456789`100*x+1.23456789`100*y").unwrap();
        let staged = stage_table(
            vk_parse!("pysecdec_result(7)").unwrap().as_view(),
            input.as_view(),
            vk_symbol!("ep"),
            "PySecDec",
        )
        .unwrap();
        let sources = collect_sources(staged.as_view()).unwrap();
        assert_eq!(sources.len(), 2);
        assert!(
            sources
                .values()
                .all(|source| source.significant_digits <= 15
                    && source.data_version == "PySecDec estimates"
                    && source.master_id.is_none())
        );
        let cancelled_centers = staged.replace_multiple(&[
            Replacement::new(vk_parse!("x").unwrap(), Atom::one()),
            Replacement::new(vk_parse!("y").unwrap(), -Atom::one()),
        ]);
        assert!(!cancelled_centers.is_zero());
        let error = super::super::materialize(
            cancelled_centers.as_view(),
            &VakintSettings::default(),
            &HashMap::new(),
            None,
        )
        .unwrap_err();
        assert!(error.to_string().contains("no reliable"), "{error}");

        let laurent = vk_parse!("1.25`100/ep+1.25`100*ep").unwrap();
        let staged = stage_table(
            vk_parse!("pysecdec_result(8)").unwrap().as_view(),
            laurent.as_view(),
            vk_symbol!("ep"),
            "PySecDec",
        )
        .unwrap();
        let sources = collect_sources(staged.as_view()).unwrap();
        assert_eq!(sources.len(), 2);
        assert_eq!(
            sources
                .values()
                .filter_map(|source| source.epsilon_power)
                .collect::<std::collections::BTreeSet<_>>(),
            std::collections::BTreeSet::from([-1, 1])
        );
    }

    #[test]
    fn pysecdec_numerators_are_not_given_private_numerical_sources() {
        let vakint = Vakint::new().unwrap();
        let settings = form_settings(
            EvaluationMethod::PySecDec(crate::PySecDecOptions::default()),
            LoopNormalizationFactor::MSbar,
            1,
        );
        let mut term = crate::VakintTerm {
            integral: vk_parse!("topo(I1L(muvsq,1))").unwrap(),
            numerator: vk_parse!("1.23456789`30").unwrap(),
            vectors: Vec::new(),
        };
        term.stage_numerical_inputs(&vakint, &settings).unwrap();
        assert!(
            !term
                .numerator
                .contains_symbol(vk_symbol!("NumericalSource"))
        );
        let (_, serialized, _) = vakint
            .prepare_expression_for_form(&settings, term.numerator, false, &[])
            .unwrap();
        assert!(!serialized.contains("NumericalSource"));
    }

    #[test]
    fn direct_form_backend_wrappers_preserve_uncertain_input_cancellations() {
        let vakint = Vakint::new().unwrap();
        let numerator = vk_parse!("1.25`30*dot(k(1),k(1))-5/4*muvsq").unwrap();
        for (method, topology, loops) in [
            (
                EvaluationMethod::MATAD(MATADOptions::default()),
                "topo(I1L(muvsq,1))",
                1,
            ),
            (
                EvaluationMethod::FMFT(FMFTOptions::default()),
                "topo(I4L_H(muvsq,1,1,1,1,0,0,0,0,0))",
                4,
            ),
        ] {
            let settings =
                form_settings(method.clone(), LoopNormalizationFactor::FMFTandMATAD, loops);
            let mut rules = vakint
                .match_topologies_for_settings(&settings, vk_parse!(topology).unwrap().as_view())
                .unwrap()
                .unwrap();
            rules.apply_replacement_rules().unwrap();
            let direct = match &method {
                EvaluationMethod::MATAD(options) => {
                    vakint.matad_evaluate(&settings, numerator.as_view(), &rules, options)
                }
                EvaluationMethod::FMFT(options) => {
                    vakint.fmft_evaluate(&settings, numerator.as_view(), &rules, options)
                }
                _ => unreachable!(),
            };
            let error = direct.unwrap_err();
            assert!(
                error.to_string().contains("no reliable"),
                "{method}: {error}"
            );
            let error = method
                .evaluate_integral(&vakint, &settings, numerator.as_view(), &rules)
                .unwrap_err();
            assert!(
                error.to_string().contains("no reliable"),
                "{method}: {error}"
            );
        }
    }

    fn form_settings(
        method: EvaluationMethod,
        normalization: LoopNormalizationFactor,
        loops: i64,
    ) -> VakintSettings {
        VakintSettings {
            evaluation_order: EvaluationOrder(vec![method]),
            integral_normalization_factor: normalization,
            number_of_terms_in_epsilon_expansion: loops + 1,
            run_time_decimal_precision: 32,
            use_dot_product_notation: true,
            ..Default::default()
        }
    }

    #[test]
    fn cold_fmft_nonunit_mass_input_precision_and_normalization() {
        let vakint = Vakint::new().unwrap();
        for normalization in [
            LoopNormalizationFactor::FMFTandMATAD,
            LoopNormalizationFactor::MSbar,
        ] {
            let settings = form_settings(
                EvaluationMethod::FMFT(FMFTOptions::default()),
                normalization.clone(),
                4,
            );
            let input = vk_parse!("1.23456789`30*topo(I4L_H(muvsq,1,1,1,1,1,1,1,1,1))").unwrap();
            let evaluated = vakint.evaluate(&settings, input.as_view()).unwrap();
            assert!(!evaluated.contains_symbol(vk_symbol!("NumericalSource")));
            let parameters =
                ahash::HashMap::from_iter([("muvsq".into(), settings.real_to_prec("1.3"))]);
            let (result, _) = vakint
                .numerical_evaluation(
                    &settings,
                    evaluated.as_view(),
                    &parameters,
                    &ahash::HashMap::default(),
                    None,
                )
                .unwrap();
            let coefficient = result.coefficients().get(&0).unwrap();
            let expected = &crate::fmft_numerics::MASTERS_NUMERIC_SUBSTITUTIONS
                [&vk_parse!("PR12ep0").unwrap()]
                .0
                * vk_parse!("10/13*123456789/100000000").unwrap();
            let expected = match normalization {
                LoopNormalizationFactor::MSbar => expected * vk_parse!("(16*𝜋^2)^-4").unwrap(),
                _ => expected,
            };
            let expected = expected
                .evaluate_with_prec(
                    &ahash::HashMap::<Atom, symbolica::domains::float::Complex<Float>>::default(),
                    240,
                )
                .unwrap();
            let delta = (coefficient.value().re.clone() - &expected.re).norm();
            assert!(
                delta < expected.re.norm() * Float::parse("1e-24", Some(240)).unwrap(),
                "{normalization}: {delta}"
            );
            assert!(coefficient.precision().imaginary().exact_zero());
            assert!(coefficient.precision().real().significant_digits().unwrap() >= 20);
            assert!(coefficient.precision().real().significant_digits().unwrap() <= 30);
        }
    }

    #[test]
    fn cold_fmft_mass_related_terms_cancel_before_finite_master_substitution() {
        let vakint = Vakint::new().unwrap();
        for normalization in [
            LoopNormalizationFactor::FMFTandMATAD,
            LoopNormalizationFactor::MSbar,
        ] {
            let settings = form_settings(
                EvaluationMethod::FMFT(FMFTOptions::default()),
                normalization,
                4,
            );
            // H starts at epsilon^0 and scales as one inverse mass squared at
            // that order. Different mass symbols force separate FORM calls;
            // keeping the mass ratio exact also checks global algebra.
            let input = vk_parse!(
                "topo(I4L_H(mass_a_sq,1,1,1,1,1,1,1,1,1))-mass_b_sq/mass_a_sq*topo(I4L_H(mass_b_sq,1,1,1,1,1,1,1,1,1))"
            )
            .unwrap();
            let evaluated = vakint.evaluate(&settings, input.as_view()).unwrap();
            assert!(evaluated.is_zero(), "{evaluated}");
        }
    }

    #[test]
    fn cold_form_custom_normalization_tracks_independent_epsilon_coefficients() {
        let vakint = Vakint::new().unwrap();
        for (method, loops, topology, polynomial) in [
            (
                EvaluationMethod::MATAD(MATADOptions::default()),
                1,
                "topo(I1L(muvsq,1))",
                "1.25`30-1.25`30*eps",
            ),
            (
                EvaluationMethod::FMFT(FMFTOptions::default()),
                4,
                "topo(I4L_H(muvsq,1,1,1,1,0,0,0,0,0))",
                "1.25`30-5.0`30*eps",
            ),
        ] {
            // K1 has equal leading and next Laurent coefficients; K1^4 has
            // ratio four. Equal floating centers in the custom polynomial
            // are independent inputs, so the next coefficient is uncertain
            // zero rather than an algebraically certified zero.
            let normalization = LoopNormalizationFactor::Custom(format!(
                "({polynomial})*({})",
                LoopNormalizationFactor::FMFTandMATAD.to_expression()
            ));
            let settings = form_settings(method.clone(), normalization, loops);
            let error = vakint
                .evaluate(&settings, vk_parse!(topology).unwrap().as_view())
                .unwrap_err();
            assert!(
                error.to_string().contains("no reliable"),
                "{method}: {error}"
            );
        }
    }

    #[test]
    fn cold_matad_symbolic_and_direct_tables_agree_at_nonunit_mass() {
        let vakint = Vakint::new().unwrap();
        let input = vk_parse!("topo(I3L(muvsq,1,1,1,1,1,1))").unwrap();
        for normalization in [
            LoopNormalizationFactor::FMFTandMATAD,
            LoopNormalizationFactor::MSbar,
        ] {
            let mut outputs = Vec::new();
            for direct_numerical_substition in [false, true] {
                let settings = form_settings(
                    EvaluationMethod::MATAD(MATADOptions {
                        direct_numerical_substition,
                        ..Default::default()
                    }),
                    normalization.clone(),
                    3,
                );
                let evaluated = vakint.evaluate(&settings, input.as_view()).unwrap();
                assert!(!evaluated.contains_symbol(vk_symbol!("NumericalSource")));
                let parameters = ahash::HashMap::from_iter([
                    (settings.mu_r_sq_symbol.clone(), settings.real_to_prec("1")),
                    ("muvsq".into(), settings.real_to_prec("1.3")),
                ]);
                let (result, _) = vakint
                    .numerical_evaluation(
                        &settings,
                        evaluated.as_view(),
                        &parameters,
                        &ahash::HashMap::default(),
                        None,
                    )
                    .unwrap();
                outputs.push(result);
            }
            for (order, reference) in outputs[0].coefficients() {
                let actual = outputs[1].coefficients().get(order).unwrap();
                for (left, right) in [
                    (&reference.value().re, &actual.value().re),
                    (&reference.value().im, &actual.value().im),
                ] {
                    if left.is_zero() && right.is_zero() {
                        continue;
                    }
                    let delta = (left.clone() - right).norm();
                    assert!(
                        delta < left.norm() * Float::parse("1e-24", Some(240)).unwrap(),
                        "{normalization}, epsilon^{order}: {delta}"
                    );
                }
            }
        }
    }
}
