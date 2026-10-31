//! Five-loop retained post-processing proxy built from the saved q0 result.
//!
//! The saved historical q0 coefficients are exactly the complete Q(d) master
//! product coefficients of the five-loop zero-momentum reduction. Rebuilding
//! them as native coefficients in `d` and assembling the caller expression the
//! way the saved application does reproduces the shape of the retained
//! catalogue combination without running the multi-hour native reduction.

use std::time::Instant;

use symbolica::prelude::{IntegerRing, RationalPolynomial};

use super::*;
use crate::utils::vakint_macros::vk_parse;

/// The previous whole-term canonicalization, kept as an exactness reference.
pub(crate) fn reference_combine_exact(expression: Atom) -> Atom {
    let expression = expression.expand();
    if MasterCatalog::has_floating_coefficients(&expression) {
        expression
    } else {
        expression.together()
    }
}

/// The previous retained catalogue combination, kept as an exactness reference.
pub(crate) fn reference_combine_master_coefficients(raw: Atom) -> Atom {
    let raw = raw.expand();
    if MasterCatalog::has_floating_coefficients(&raw) {
        return raw;
    }
    let symbols: Vec<_> = MasterCatalog::referenced_masters(&raw)
        .into_iter()
        .collect();
    let terms = raw.coefficient_list::<i32>(&symbols);
    if terms
        .iter()
        .any(|(_, coefficient)| !MasterCatalog::referenced_masters(coefficient).is_empty())
    {
        return reference_combine_exact(raw);
    }
    Atom::add_many(
        terms
            .into_iter()
            .map(|(monomial, coefficient)| monomial * reference_combine_exact(coefficient)),
    )
}

/// Workers of the timed combination: the default `n_cores` (`VakintSettings::default()`).
fn proxy_threads() -> usize {
    super::super::postprocessing_threads(&VakintSettings::default())
}

fn proxy_products() -> usize {
    std::env::var("VAKINT_POSTPROC_PROXY_PRODUCTS")
        .ok()
        .map(|value| value.parse().expect("product count"))
        .unwrap_or(usize::MAX)
}

/// The caller expression handed to the retained catalogue at five loops:
/// native coefficients in `d` mapped to `4-2*ε`, master products, the common
/// spectator and mass power, and the five-loop MSbar measure normalization.
pub(crate) fn five_loop_proxy_raw(products: usize) -> Atom {
    crate::Vakint::initialize_vakint_symbols();
    let bytes = include_bytes!("../../tests/resources/rustred/q0-retained.raw.atom");
    let mut stream = bytes.as_slice();
    let original = Atom::import(&mut stream, None).unwrap();
    let symbols: Vec<_> = MasterCatalog::referenced_masters(&original)
        .into_iter()
        .collect();
    let epsilon = vk_parse!("ε").unwrap();
    // The saved coefficients carry the historical measure factor exp(-5γε);
    // native coefficients are rational functions of d alone.
    let euler_gamma = symbolica::parse!("symbolica::γ");
    let coupling = vk_parse!("UFO::GC_11").unwrap();
    let mass = vk_parse!("mUV").unwrap();
    let dimension = vk_parse!("rustred::d").unwrap();
    let physical = Atom::num(4) - Atom::num(2) * &epsilon;
    let mut raw = Atom::Zero;
    for (monomial, coefficient) in original
        .coefficient_list::<i32>(&symbols)
        .into_iter()
        .take(products)
    {
        let scalar = coefficient
            .replace(euler_gamma.to_pattern())
            .with(Atom::Zero)
            .replace(coupling.to_pattern())
            .with(Atom::num(1))
            .replace(mass.to_pattern())
            .with(Atom::num(1))
            .replace(epsilon.to_pattern())
            .with(((Atom::num(4) - &dimension) / Atom::num(2)).to_pattern());
        let native: RationalPolynomial<IntegerRing, u16> =
            scalar.to_rational_polynomial(&Q, &Z, None);
        assert_eq!(
            native.get_variables().len(),
            1,
            "native coefficients depend on d only"
        );
        raw += native
            .to_expression()
            .replace(dimension.to_pattern())
            .with(physical.to_pattern())
            * monomial;
    }
    let mass_power = Atom::num((5, 2)) * &physical + Atom::num(-14);
    let settings = VakintSettings {
        integral_normalization_factor: crate::LoopNormalizationFactor::MSbar,
        ..VakintSettings::default()
    };
    raw * coupling.pow(Atom::num(10))
        * mass.pow(Atom::num(2)).pow(mass_power)
        * super::super::five_loop::SavedApplication::normalization(&settings, 5).unwrap()
}

/// Times the retained combination on the proxy. With
/// `VAKINT_POSTPROC_PROXY_REFERENCE=1` the previous implementation also runs
/// and both results must be identical expressions.
#[test]
#[ignore = "five-loop retained post-processing timing proxy; run explicitly"]
fn five_loop_retained_post_processing_proxy() {
    let started = Instant::now();
    let raw = five_loop_proxy_raw(proxy_products());
    println!(
        "proxy construction {:.3} s, {} bytes",
        started.elapsed().as_secs_f64(),
        raw.as_view().get_byte_size()
    );
    let started = Instant::now();
    let simplified = crate::utils::simplify_real(raw.as_view());
    println!("simplify_real {:.3} s", started.elapsed().as_secs_f64());
    let threads = proxy_threads();
    let started = Instant::now();
    let combined = MasterCatalog::combine_master_coefficients(simplified.clone(), threads);
    println!(
        "combine_master_coefficients ({threads} threads) {:.3} s, {} bytes, canonical sha256 {:x}",
        started.elapsed().as_secs_f64(),
        combined.as_view().get_byte_size(),
        Sha256::digest(combined.to_canonical_string().as_bytes())
    );
    if threads > 1 {
        let started = Instant::now();
        let sequential = MasterCatalog::combine_master_coefficients(simplified.clone(), 1);
        println!(
            "combine_master_coefficients (1 thread) {:.3} s, identical {}",
            started.elapsed().as_secs_f64(),
            sequential == combined
        );
        assert_eq!(sequential, combined);
    }
    if std::env::var("VAKINT_POSTPROC_PROXY_REFERENCE").as_deref() == Ok("1") {
        let started = Instant::now();
        let reference = reference_combine_master_coefficients(simplified);
        println!(
            "reference combine_master_coefficients {:.3} s, identical {}",
            started.elapsed().as_secs_f64(),
            reference == combined
        );
        assert_eq!(reference, combined);
    }
}

/// Per-coefficient timings of the retained combination on the proxy.
#[test]
#[ignore = "five-loop retained post-processing diagnostics; run explicitly"]
fn five_loop_retained_post_processing_breakdown() {
    let raw = crate::utils::simplify_real(five_loop_proxy_raw(proxy_products()).as_view());
    let started = Instant::now();
    let expanded = raw.expand();
    println!("expand {:.3} s", started.elapsed().as_secs_f64());
    let symbols: Vec<_> = MasterCatalog::referenced_masters(&expanded)
        .into_iter()
        .collect();
    let started = Instant::now();
    let terms = expanded.coefficient_list::<i32>(&symbols);
    println!(
        "coefficient_list {:.3} s, {} monomials",
        started.elapsed().as_secs_f64(),
        terms.len()
    );
    for (index, (_, coefficient)) in terms.iter().enumerate() {
        let started = Instant::now();
        let inner = coefficient.expand();
        let expand_seconds = started.elapsed().as_secs_f64();
        let shared = together_over_shared_remainder(&inner);
        let terms = match inner.as_view() {
            AtomView::Add(add) => add.get_nargs(),
            _ => 1,
        };
        println!(
            "coefficient {index}: {terms} terms, expand {expand_seconds:.3} s, shared {} {:.3} s",
            shared.is_some(),
            started.elapsed().as_secs_f64() - expand_seconds
        );
        if shared.is_none() {
            let distinct: BTreeSet<_> = match inner.as_view() {
                AtomView::Add(add) => add
                    .iter()
                    .map(|term| {
                        SplitTerm::new(term)
                            .remainder
                            .iter()
                            .map(|f| f.to_owned())
                            .collect::<Vec<_>>()
                    })
                    .collect(),
                _ => BTreeSet::new(),
            };
            println!("  distinct remainders {}", distinct.len());
            for remainder in distinct.iter().take(3) {
                println!(
                    "  remainder: {}",
                    Atom::mul_many(remainder.iter())
                        .to_string()
                        .chars()
                        .take(400)
                        .collect::<String>()
                );
            }
        }
    }
}

/// Prints the structure of the first expanded term at each retained stage.
#[test]
#[ignore = "five-loop retained post-processing structure diagnostics; run explicitly"]
fn five_loop_retained_post_processing_structure() {
    fn first_term(label: &str, expression: &Atom, symbols: &[Atom]) {
        let terms = expression.coefficient_list::<i32>(symbols);
        let (monomial, coefficient) = terms
            .iter()
            .find(|(monomial, _)| {
                monomial == &symbols[0] || monomial.to_string().contains(&symbols[0].to_string())
            })
            .unwrap_or(&terms[0]);
        let first = match coefficient.as_view() {
            AtomView::Add(add) => add.iter().next().unwrap().to_owned(),
            view => view.to_owned(),
        };
        let factors: Vec<_> = match first.as_view() {
            AtomView::Mul(mul) => mul
                .iter()
                .map(|f| f.to_string().chars().take(120).collect::<String>())
                .collect(),
            view => vec![view.to_string()],
        };
        println!(
            "{label}: monomial {}\n  {}",
            monomial.to_string().chars().take(200).collect::<String>(),
            factors.join("\n  ")
        );
    }
    let raw = five_loop_proxy_raw(proxy_products());
    let symbols: Vec<_> = MasterCatalog::referenced_masters(&raw)
        .into_iter()
        .collect();
    first_term("raw", &raw, &symbols);
    let simplified = crate::utils::simplify_real(raw.as_view());
    first_term("simplified", &simplified, &symbols);
    first_term("raw expanded", &raw.expand(), &symbols);
    first_term("simplified expanded", &simplified.expand(), &symbols);
}
