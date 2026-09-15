//! Exact terminal data in the canonical unit-mass K6 parent slots.
//!
//! These records neither construct a catalog nor grant artifact authority.
//! The existing loader pairs them with genuine K6 ClosedArtifact bytes;
//! TerminalCatalog::compile must still bind the schema,
//! algorithm, family fingerprint and complete typed key set.
//!
//! The 38 keys were printed and individually checked against their sectors in
//! RustRed's canonical checked-program regression (623 rules, 2026-09-15).
//! Denominators are k1²-1, k2²-1, k3²-1, (k3-k1)²-1, (k1-k2)²-1,
//! (k2-k3)²-1. These are already parent-slot keys: do not apply the earlier
//! source-reference slot exchange a second time.
//!
//! RustRed's authenticated S4 and unit-Jacobian product proofs give six
//! sector orbits but only five distinct values. The three irreducible anchors
//! are the existing k6_matad_oracle_tests.rs corner records, independently
//! rerun with the preserved MATAD/FORM executable. The five-line record keeps
//! its entire outer minus sign; it is not simply miD5. Constants use Vakint's
//! internal ep and pre-normalization Minkowski MATAD convention at m²=1.
//! No numerical truncation or extra independent master evaluation is needed.

use super::TerminalSource;
#[cfg(test)]
use super::{TerminalValueSource, compile_value};
#[cfg(test)]
use crate::utils::vakint_macros::{vk_parse, vk_symbol};
#[cfg(test)]
use rustred::family::IntegralKey;
#[cfg(test)]
use std::collections::{BTreeMap, BTreeSet};
#[cfg(test)]
use symbolica::atom::{Atom, AtomCore};

const TADPOLE_CUBED: &str = "(-Gam(1,1)/(ep*(ep-1)))^3";
const SUNSET_TIMES_TADPOLE: &str = "(-Gam(1,1)/(ep*(ep-1)))*(-miT111)";
const FOUR_LINE: &str = "miBN";
const SIX_LINE: &str = "miD6";
// Only the common mass was set to one in the complete raw MATAD record.
// Preserve the original expression and signs; Symbolica handles arithmetic.
const FIVE_LINE: &str = "(((-2*ep+4)*-2208+(-2*ep+4)^2*1736+(-2*ep+4)^3*-718+(-2*ep+4)^4*165+(-2*ep+4)^5*-20+(-2*ep+4)^6+1152)^(-1)*Gam(1,1)^3*-16+((-2*ep+4)*-3+12)*((-2*ep+4)*2+-6)^(-1)*miD5+((-2*ep+4)*-7+(-2*ep+4)^2+12)^(-1)*-4*Gam(1,1)*miT111+((-2*ep+4)*3+-8)*((-2*ep+4)*8+-24)^(-1)*miBN)*-1";

pub(in crate::rustred_evaluation) const SOURCES: [TerminalSource<'static>; 38] = [
    TerminalSource::exact_matad_basis(&[0, 0, 1, 0, 1, 1], TADPOLE_CUBED),
    TerminalSource::exact_matad_basis(&[0, 0, 1, 1, 0, 1], TADPOLE_CUBED),
    TerminalSource::exact_matad_basis(&[0, 0, 1, 1, 1, 0], TADPOLE_CUBED),
    TerminalSource::exact_matad_basis(&[0, 0, 1, 1, 1, 1], SUNSET_TIMES_TADPOLE),
    TerminalSource::exact_matad_basis(&[0, 1, 0, 0, 1, 1], TADPOLE_CUBED),
    TerminalSource::exact_matad_basis(&[0, 1, 0, 1, 0, 1], TADPOLE_CUBED),
    TerminalSource::exact_matad_basis(&[0, 1, 0, 1, 1, 0], TADPOLE_CUBED),
    TerminalSource::exact_matad_basis(&[0, 1, 0, 1, 1, 1], SUNSET_TIMES_TADPOLE),
    TerminalSource::exact_matad_basis(&[0, 1, 1, 0, 1, 0], TADPOLE_CUBED),
    TerminalSource::exact_matad_basis(&[0, 1, 1, 0, 1, 1], SUNSET_TIMES_TADPOLE),
    TerminalSource::exact_matad_basis(&[0, 1, 1, 1, 0, 0], TADPOLE_CUBED),
    TerminalSource::exact_matad_basis(&[0, 1, 1, 1, 0, 1], SUNSET_TIMES_TADPOLE),
    TerminalSource::exact_matad_basis(&[0, 1, 1, 1, 1, 0], FOUR_LINE),
    TerminalSource::exact_matad_basis(&[0, 1, 1, 1, 1, 1], FIVE_LINE),
    TerminalSource::exact_matad_basis(&[1, 0, 0, 0, 1, 1], TADPOLE_CUBED),
    TerminalSource::exact_matad_basis(&[1, 0, 0, 1, 0, 1], TADPOLE_CUBED),
    TerminalSource::exact_matad_basis(&[1, 0, 0, 1, 1, 0], TADPOLE_CUBED),
    TerminalSource::exact_matad_basis(&[1, 0, 0, 1, 1, 1], SUNSET_TIMES_TADPOLE),
    TerminalSource::exact_matad_basis(&[1, 0, 1, 0, 0, 1], TADPOLE_CUBED),
    TerminalSource::exact_matad_basis(&[1, 0, 1, 0, 1, 0], TADPOLE_CUBED),
    TerminalSource::exact_matad_basis(&[1, 0, 1, 0, 1, 1], FOUR_LINE),
    TerminalSource::exact_matad_basis(&[1, 0, 1, 1, 0, 1], SUNSET_TIMES_TADPOLE),
    TerminalSource::exact_matad_basis(&[1, 0, 1, 1, 1, 0], SUNSET_TIMES_TADPOLE),
    TerminalSource::exact_matad_basis(&[1, 0, 1, 1, 1, 1], FIVE_LINE),
    TerminalSource::exact_matad_basis(&[1, 1, 0, 0, 0, 1], TADPOLE_CUBED),
    TerminalSource::exact_matad_basis(&[1, 1, 0, 0, 1, 1], SUNSET_TIMES_TADPOLE),
    TerminalSource::exact_matad_basis(&[1, 1, 0, 1, 0, 0], TADPOLE_CUBED),
    TerminalSource::exact_matad_basis(&[1, 1, 0, 1, 0, 1], FOUR_LINE),
    TerminalSource::exact_matad_basis(&[1, 1, 0, 1, 1, 0], SUNSET_TIMES_TADPOLE),
    TerminalSource::exact_matad_basis(&[1, 1, 0, 1, 1, 1], FIVE_LINE),
    TerminalSource::exact_matad_basis(&[1, 1, 1, 0, 0, 0], TADPOLE_CUBED),
    TerminalSource::exact_matad_basis(&[1, 1, 1, 0, 0, 1], SUNSET_TIMES_TADPOLE),
    TerminalSource::exact_matad_basis(&[1, 1, 1, 0, 1, 0], SUNSET_TIMES_TADPOLE),
    TerminalSource::exact_matad_basis(&[1, 1, 1, 0, 1, 1], FIVE_LINE),
    TerminalSource::exact_matad_basis(&[1, 1, 1, 1, 0, 0], SUNSET_TIMES_TADPOLE),
    TerminalSource::exact_matad_basis(&[1, 1, 1, 1, 0, 1], FIVE_LINE),
    TerminalSource::exact_matad_basis(&[1, 1, 1, 1, 1, 0], FIVE_LINE),
    TerminalSource::exact_matad_basis(&[1, 1, 1, 1, 1, 1], SIX_LINE),
];

#[test]
fn candidate_keys_match_the_checked_program_and_exact_symmetry_orbits() {
    // Collision-free 64-bit characteristic set of the 38 printed six-bit
    // corner keys, with the first physical denominator as the highest bit.
    const CHECKED_KEY_SET: u64 = 0xfffaeee8fce8e800;
    // Images of denominator slots under the two independently authenticated
    // momentum maps in RustRed's three_loop/symmetry.rs. This is a structural
    // data check, not another graph matcher or a source of symmetry authority.
    const GENERATORS: [[usize; 6]; 2] = [[0, 4, 3, 2, 1, 5], [1, 2, 0, 4, 5, 3]];
    let mut values = BTreeMap::new();
    let mut key_set = 0_u64;
    for source in &SOURCES {
        let key = IntegralKey::try_new(source.powers.iter().copied()).unwrap();
        let powers: [i64; 6] = key.powers().try_into().unwrap();
        assert!(powers.iter().all(|&power| power == 0 || power == 1));
        let bits = powers.iter().fold(0_u32, |bits, &power| {
            (bits << 1) | u32::try_from(power).unwrap()
        });
        key_set |= 1_u64 << bits;
        let TerminalValueSource::ExactMatadBasis(value) = source.value else {
            panic!("the candidate must not introduce numerical terminal tables");
        };
        assert!(
            values.insert(powers, value).is_none(),
            "duplicate {powers:?}"
        );
    }
    assert_eq!(values.len(), 38);
    assert_eq!(key_set, CHECKED_KEY_SET);

    let mut orbit_counts = BTreeMap::new();
    for (&powers, &value) in &values {
        let mut orbit = BTreeSet::from([powers]);
        let mut pending = vec![powers];
        while let Some(key) = pending.pop() {
            for permutation in GENERATORS {
                let mut image = [0_i64; 6];
                for axis in 0..6 {
                    image[permutation[axis]] = key[axis];
                }
                assert_eq!(values.get(&image), Some(&value), "orbit of {powers:?}");
                if orbit.insert(image) {
                    assert!(orbit.len() <= 24);
                    pending.push(image);
                }
            }
        }
        *orbit_counts.entry(*orbit.first().unwrap()).or_insert(0) += 1;
    }
    assert_eq!(
        orbit_counts,
        BTreeMap::from([
            ([0, 0, 1, 0, 1, 1], 12),
            ([0, 0, 1, 1, 0, 1], 4),
            ([0, 0, 1, 1, 1, 1], 12),
            ([0, 1, 1, 1, 1, 0], 3),
            ([0, 1, 1, 1, 1, 1], 6),
            ([1, 1, 1, 1, 1, 1], 1),
        ])
    );
}

#[test]
fn candidate_anchors_equal_the_recorded_matad_values_at_unit_mass() {
    // Existing raw oracle records, with only namespace spelling and the
    // external epsilon name mapped to this module's internal ep convention.
    // Keeping their mass dependence checks that unit-scale extraction does
    // not accidentally ship the later per-loop normalization a second time.
    let records = [
        ([0, 1, 1, 1, 1, 0], "(muvsq^(-1))^(3*ep)*miBN*muvsq^2"),
        (
            [0, 1, 1, 1, 1, 1],
            "(((-2*ep+4)*-2208+(-2*ep+4)^2*1736+(-2*ep+4)^3*-718+(-2*ep+4)^4*165+(-2*ep+4)^5*-20+(-2*ep+4)^6+1152)^(-1)*Gam(1,1)^3*-16*muvsq+((-2*ep+4)*-3+12)*((-2*ep+4)*2+-6)^(-1)*miD5*muvsq+((-2*ep+4)*-7+(-2*ep+4)^2+12)^(-1)*-4*Gam(1,1)*miT111*muvsq+((-2*ep+4)*3+-8)*((-2*ep+4)*8+-24)^(-1)*miBN*muvsq)*(muvsq^(-1))^(3*ep)*-1",
        ),
        ([1, 1, 1, 1, 1, 1], "(muvsq^(-1))^(3*ep)*miD6"),
    ];
    let mass = vk_parse!("muvsq").unwrap();
    let one = Atom::num(1);
    for (powers, raw) in records {
        let source = SOURCES
            .iter()
            .find(|source| source.powers == powers)
            .unwrap();
        let key = IntegralKey::try_new(powers).unwrap();
        let value = compile_value(&key, &source.value)
            .unwrap()
            .expression(false);
        let expected = vk_parse!(raw)
            .unwrap()
            .replace(mass.to_pattern())
            .with(one.to_pattern())
            .together();
        assert!((value - expected).together().is_zero(), "{powers:?}");
    }
}

#[test]
fn candidate_products_preserve_the_existing_minkowski_tadpole_signs() {
    // These are the same K1 and K3 unit-mass anchors already loaded in
    // artifact.rs. Product proofs have determinant ±1 and normalization one.
    let tadpole = vk_parse!("-Gam(1,1)/(ep*(ep-1))").unwrap();
    let sunset = vk_parse!("-miT111").unwrap();
    let expected = [
        ([0, 0, 1, 0, 1, 1], tadpole.clone().pow(Atom::num(3))),
        ([0, 0, 1, 1, 0, 1], tadpole.clone().pow(Atom::num(3))),
        ([0, 0, 1, 1, 1, 1], tadpole * sunset),
    ];
    for (powers, product) in expected {
        let source = SOURCES
            .iter()
            .find(|source| source.powers == powers)
            .unwrap();
        let key = IntegralKey::try_new(powers).unwrap();
        let value = compile_value(&key, &source.value)
            .unwrap()
            .expression(false);
        assert!((value - product).together().is_zero(), "{powers:?}");
    }
}

#[test]
fn candidate_values_are_exact_internal_basis_expressions_without_mass_or_tails() {
    for source in &SOURCES {
        let key = IntegralKey::try_new(source.powers.iter().copied()).unwrap();
        let compiled = compile_value(&key, &source.value).unwrap();
        let raw = compiled.expression(false);
        assert_eq!(raw, compiled.expression(true));
        for forbidden in ["muvsq", "mursq", "ε", "Oep", "RustRedMaster", "d"] {
            assert!(
                !raw.contains_symbol(vk_symbol!(forbidden)),
                "{:?}: {forbidden}",
                source.powers
            );
        }
        assert!(!raw.is_zero());
    }
}

#[test]
fn candidate_master_substitution_reaches_the_finite_part_at_nonunit_mass_without_form() {
    use crate::{MATADOptions, VakintSettings, matad::MATAD};

    let settings = VakintSettings {
        form_exe_path: "/definitely/not/a/form/executable".into(),
        number_of_terms_in_epsilon_expansion: 4,
        run_time_decimal_precision: 80,
        use_dot_product_notation: true,
        ..VakintSettings::default()
    };
    let mass = vk_parse!("13/10").unwrap();
    let finalizer = MATAD::with_settings(settings);
    let mut checked_values = BTreeSet::new();
    for source in &SOURCES {
        let TerminalValueSource::ExactMatadBasis(value) = source.value else {
            panic!("candidate K6 data must remain exact");
        };
        // Equal orbit values also have equal dimensional mass powers. Test
        // each genuinely different master expression once, not 38 copies of
        // the same native series expansion.
        let mass_power = 6 - source.powers.iter().sum::<i64>();
        if !checked_values.insert((value, mass_power)) {
            continue;
        }
        let key = IntegralKey::try_new(source.powers.iter().copied()).unwrap();
        let raw = compile_value(&key, &source.value)
            .unwrap()
            .expression(false)
            * mass.clone().pow(Atom::num(mass_power));
        let result = finalizer
            .finalize_reduced_masters(raw, 3, None, &mass, &MATADOptions::default())
            .unwrap_or_else(|error| panic!("{:?}: {error}", source.powers));
        assert!(!result.is_zero(), "{:?}", source.powers);
        for forbidden in ["Gam", "miT111", "miBN", "miD5", "miD6", "Oep"] {
            assert!(
                !result.contains_symbol(vk_symbol!(forbidden)),
                "{:?}: unmaterialized {forbidden}",
                source.powers
            );
        }
    }
    assert_eq!(checked_values.len(), 5);
}

#[test]
#[ignore = "offline FORM/MATAD check of all 38 candidate values and symbolic mass restoration; not a RustRed artifact test"]
fn candidate_all_38_values_match_matad_with_symbolic_mass_restoration() {
    use crate::{
        EvaluationOrder, LoopNormalizationFactor, MATADOptions, Vakint, VakintSettings,
        matad::MATAD,
    };

    let options = MATADOptions {
        expand_masters: false,
        ..MATADOptions::default()
    };
    let oracle_settings = VakintSettings {
        form_exe_path: std::env::var("VAKINT_K6_ORACLE_FORM_PATH")
            .or_else(|_| std::env::var("FORM_PATH"))
            .expect("offline K6 oracle requires VAKINT_K6_ORACLE_FORM_PATH or FORM_PATH"),
        evaluation_order: EvaluationOrder::matad_only(Some(options.clone())),
        integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD,
        run_time_decimal_precision: 80,
        number_of_terms_in_epsilon_expansion: 4,
        use_dot_product_notation: true,
        allow_unknown_integrals: false,
        ..VakintSettings::default()
    };
    let candidate_settings = VakintSettings {
        form_exe_path: "/definitely/not/a/form/executable".into(),
        ..oracle_settings.clone()
    };
    let finalizer = MATAD::with_settings(candidate_settings);
    let vakint = Vakint::new().unwrap();
    let mass = vk_parse!("muvsq").unwrap();
    for source in &SOURCES {
        let [n1, n2, n3, n4, n5, n6] = source.powers else {
            panic!("a canonical K6 key must have six coordinates");
        };
        let input = vk_parse!(&format!("topo(I3L(muvsq,{n1},{n2},{n3},{n4},{n5},{n6}))")).unwrap();
        let oracle = vakint
            .evaluate_integral(&oracle_settings, input.as_view())
            .unwrap_or_else(|error| panic!("{:?}: {error}", source.powers));
        let key = IntegralKey::try_new(source.powers.iter().copied()).unwrap();
        // The stored constants are unit-mass, pre-normalization Minkowski
        // values. Restore only the classical power here; the existing native
        // finalizer supplies the common -3*epsilon mass dependence exactly
        // once. The real reducer will obtain its additional power shifts from
        // the existing homogeneity witness, not from this candidate-data test.
        let raw = compile_value(&key, &source.value)
            .unwrap()
            .expression(false)
            * mass
                .clone()
                .pow(Atom::num(6 - source.powers.iter().sum::<i64>()));
        let candidate = finalizer
            .finalize_reduced_masters(raw, 3, None, &mass, &options)
            .unwrap();
        assert!(
            (candidate - oracle).together().is_zero(),
            "candidate value or symbolic mass restoration differs at {:?}",
            source.powers
        );
    }
}
