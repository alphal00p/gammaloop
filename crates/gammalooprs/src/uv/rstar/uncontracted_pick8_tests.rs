//! Pick 8 (`crown_legs_at_node`) at a POSITIVE Taylor degree: helper VH's
//! caution on the `100000 + e` leg of a paired crown edge.
//!
//! The record is the K2 child record of L5 006961 id 5 that the gate binary
//! of 0484ddd9f6 wrote (`testdata/u7a/l5-006961-id5-paired-leg.record.txt`):
//! crown edge 13 has both ends in the child, its legs are 13 (out of vertex
//! 4) and 100013 (into the `ct:` vertex 3), and the insertion factor at
//! vertex 3 names `Q(100013, ..)` (the rename of pick 8). Its Taylor degree
//! is 0, as for every record with a paired crown edge (a child or IR hard
//! part with one has at least four legs, `local_dod <= 0`); the variants
//! below raise it to 1 and 2 to test what a positive degree would do:
//!
//! * the leg is scaled exactly as its partner: the record with the two legs'
//!   roles swapped (13 into vertex 3, 100013 out of vertex 4, every `Q(13)`
//!   <-> `Q(100013)`, `external 3 edge=13` kept) has the same value, and so
//!   has the record whose `scaled_edges` also lists 100013 or lists 100013
//!   instead of 13 (vakint scales the momentum `p(j)` that a leg carries,
//!   `scaled_edges` only has to cover every `p(j)`);
//! * an undeclared `Q(100000 + e)` (the leg's line dropped, or the factor
//!   naming an id without a line) is REFUSED, at degree 0 and above: every
//!   record the emitter writes has a `taylor` line, and vakint's closure
//!   validation refuses a constant external vector with a Taylor step;
//! * the pre-pick-8 form `Q(13)` in the insertion at vertex 3 is refused
//!   (`not a boundary momentum of the insertion`).

use eyre::{Result, ensure, eyre};
use symbolica::atom::{Atom, AtomCore};
use vakint::{
    Vakint,
    p_integral::factorized::{FactorizedKind, FactorizedNumerator, evaluate_factorized_closure},
};

use crate::initialisation::test_initialise;

use super::{recursive_tests::settings, uncontracted_oracle_tests::atoms_equal};

const RECORD: &str = include_str!("testdata/u7a/l5-006961-id5-paired-leg.record.txt");

fn lines_map(text: &str, f: impl Fn(&str) -> Option<String>) -> String {
    let mut out = String::new();
    for line in text.lines() {
        if let Some(line) = f(line) {
            out.push_str(&line);
            out.push('\n');
        }
    }
    out
}

fn with_degree(text: &str, degree: u32) -> String {
    let out = text.replace("taylor\tdegree=0\t", &format!("taylor\tdegree={degree}\t"));
    assert!(out.contains(&format!("taylor\tdegree={degree}\t")));
    out
}

fn with_scaled(text: &str, scaled: &str) -> String {
    let out = text.replace("\tscaled_edges=0,1,13\t", &format!("\tscaled_edges={scaled}\t"));
    assert!(out.contains(&format!("\tscaled_edges={scaled}\t")));
    out
}

/// The legs' roles swapped: 13 enters vertex 3, 100013 leaves vertex 4, every
/// momentum renamed with them; `external 3 edge=13` now names the incoming leg.
fn swapped(text: &str) -> String {
    let out = text
        .replace("::Q(13,", "::Q(@@,")
        .replace("::Q(100013,", "::Q(13,")
        .replace("::Q(@@,", "::Q(100013,")
        .replace(
            "edge\t13\tsource=4\ttarget=ext\thedges=24,-\t",
            "edge\t13\tsource=ext\ttarget=3\thedges=-,25\t",
        )
        .replace(
            "edge\t100013\tsource=ext\ttarget=3\thedges=-,25\t",
            "edge\t100013\tsource=4\ttarget=ext\thedges=24,-\t",
        );
    assert!(out.contains("edge\t13\tsource=ext\ttarget=3\thedges=-,25\t"));
    assert!(out.contains("edge\t100013\tsource=4\ttarget=ext\thedges=24,-\t"));
    out
}

/// The insertion factor's `Q(100013, ..)` renamed to `Q(id, ..)`.
fn insertion_names(text: &str, id: usize) -> String {
    lines_map(text, |line| {
        Some(if line.starts_with("insertion_factor\t") {
            line.replace("::Q(100013,", &format!("::Q({id},"))
        } else {
            line.to_owned()
        })
    })
}

/// The leg's line dropped (the record keeps `Q(100013)` in its insertion).
fn leg_dropped(text: &str) -> String {
    lines_map(text, |line| {
        if line.starts_with("edge\t100013\t") || line.starts_with("edge_factor\t100013\t") {
            None
        } else if line == "edges\t10" {
            Some("edges\t9".to_owned())
        } else if line == "lines\t61" {
            // the record's own line count, so that the refusal is not the count's
            Some("lines\t59".to_owned())
        } else {
            Some(line.to_owned())
        }
    })
}

fn read(text: &str) -> std::result::Result<FactorizedNumerator, String> {
    let value = FactorizedNumerator::from_text(text).map_err(|e| e.to_string())?;
    value.validate().map_err(|e| e.to_string())?;
    Ok(value)
}

/// vakint's raw closure value of a record text (`Err`: the reader's, the
/// validation's or the evaluation's error text).
fn evaluate(vakint: &Vakint, text: &str) -> std::result::Result<Atom, String> {
    let value = read(text)?;
    let FactorizedKind::Closure(spec) = &value.kind else {
        return Err("not a closure record".into());
    };
    let mut vakint_settings = settings().vakint.true_settings();
    vakint_settings.number_of_terms_in_epsilon_expansion =
        i64::try_from(spec.depth.expect("depth")).expect("depth");
    let prefactor = value.prefactor.total();
    let evaluated = evaluate_factorized_closure(
        vakint,
        &vakint_settings,
        std::slice::from_ref(&value),
        std::slice::from_ref(&prefactor),
    )
    .map_err(|e| e.to_string())?;
    evaluated
        .records
        .into_iter()
        .next()
        .map(|record| record.value)
        .ok_or_else(|| "no record value".into())
}

fn refused(what: &str, result: std::result::Result<Atom, String>, needle: &str) -> Result<()> {
    match result {
        Ok(value) => Err(eyre!("{what}: evaluated, expected a refusal: {value:.300}")),
        Err(error) => {
            eprintln!("PICK8 {what}: refused: {error:.400}");
            ensure!(
                error.contains(needle),
                "{what}: refused, but not by `{needle}`: {error}"
            );
            Ok(())
        }
    }
}

#[test]
fn pick8_paired_leg_scales_like_its_partner_at_a_positive_degree() -> Result<()> {
    test_initialise()?;
    let vakint = Vakint::new()?;
    ensure!(RECORD.contains("taylor\tdegree=0\tcompensate_mass=1\tscaled_edges=0,1,13\t"));
    ensure!(RECORD.contains("edge\t100013\tsource=ext\ttarget=3\thedges=-,25\t"));
    ensure!(RECORD.contains("insertion\t3\thedges=5,18,25\t"));
    // degree 0 (the record as written): the legs swapped agree too
    let base = evaluate(&vakint, RECORD).map_err(|e| eyre!("degree 0: {e}"))?;
    let base_swapped = evaluate(&vakint, &swapped(RECORD)).map_err(|e| eyre!("{e}"))?;
    eprintln!("PICK8 degree 0: zero={} len={}", base.is_zero(), base.to_string().len());
    ensure!(atoms_equal(&base, &base_swapped), "degree 0: swapped legs differ");
    let mut nonzero = false;
    for degree in [1, 2] {
        let text = with_degree(RECORD, degree);
        let value = evaluate(&vakint, &text).map_err(|e| eyre!("degree {degree}: {e}"))?;
        let zero = value.is_zero();
        nonzero |= !zero;
        let differs_from_0 = !atoms_equal(&value, &base);
        eprintln!(
            "PICK8 degree {degree}: zero={zero} len={} differs_from_degree_0={differs_from_0}",
            value.to_string().len()
        );
        for (what, variant) in [
            ("legs swapped", swapped(&text)),
            ("scaled_edges=0,1,13,100013", with_scaled(&text, "0,1,13,100013")),
            ("scaled_edges=0,1,100013", with_scaled(&text, "0,1,100013")),
        ] {
            let other = evaluate(&vakint, &variant).map_err(|e| eyre!("degree {degree} {what}: {e}"))?;
            let equal = atoms_equal(&value, &other);
            eprintln!("PICK8 degree {degree} {what}: equal={equal}");
            ensure!(equal, "degree {degree}: {what} differs");
        }
    }
    ensure!(nonzero, "every positive-degree value is zero: the test shows nothing");
    // an undeclared Q(100000 + e): refused at every degree
    for degree in [0, 1, 2] {
        let text = with_degree(RECORD, degree);
        refused(
            &format!("degree {degree} factor names Q(200013) (no line)"),
            evaluate(&vakint, &insertion_names(&text, 200_013)),
            "constant external vectors",
        )?;
        refused(
            &format!("degree {degree} leg line 100013 dropped"),
            evaluate(&vakint, &leg_dropped(&text)),
            "",
        )?;
        refused(
            &format!("degree {degree} pre-pick-8 Q(13) at vertex 3"),
            evaluate(&vakint, &insertion_names(&text, 13)),
            "not a boundary momentum of the insertion",
        )?;
    }
    Ok(())
}
