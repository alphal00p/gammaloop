//! D6 (helper VH): the Forcer cross-check of the open-index p-integrals. The
//! open-index value of RustRed (FeynKit tensor reduction with the external
//! momentum declared, then the p lowering and RustRed) closed afterwards with
//! every metric/p closing tensor equals FORCER's value of the input contracted
//! with that tensor (the scalar p-integral, FORM 5), exactly (`Symbolic`
//! masters, G-scheme), at one, two and three loops. Test-only: FORM 5 comes
//! from the environment (`VAKINT_TEST_FORM5`, `FORM_PATH` or `form` on the
//! PATH, e.g. the effort's devenv); without it the test reports a skip, and
//! fails only under `VAKINT_REQUIRE_FORM5=1` (the gate's setting).
#![cfg(feature = "feynkit-ingress")]

#[path = "support/open_indices.rs"]
mod open_indices;
#[path = "support/p_backends.rs"]
mod p_backends;
#[path = "support/p_integrals.rs"]
mod p_integrals;
#[path = "support/p_rustred.rs"]
mod p_rustred;

use open_indices::*;
use p_backends::*;
use p_integrals::parse;
use p_rustred::rustred_p_settings;
use vakint::{LoopNormalizationFactor, RustRedMasterMode, vakint_symbol};

const CASES: [(&str, &str, usize, &[&str]); 3] = [
    (
        "k(1,mu)*k(1,nu)",
        "topo(prop(1,edge(1,2),k(1),0,1)*prop(2,edge(2,1),k(1)-p(1),0,2))",
        1,
        &["mu", "nu"],
    ),
    (
        "k(2,mu)*k(2,nu)*k(1,rho)",
        "topo(P2L(p(1),1,1,1,1,1))",
        2,
        &["mu", "nu", "rho"],
    ),
    (
        "k(1,mu)*k(3,nu)",
        "topo(P3L_LA(p(1),1,1,1,1,1,1,1,1))",
        3,
        &["mu", "nu"],
    ),
];

#[test]
fn closed_open_index_values_equal_forcer_on_the_contracted_scalars() {
    let Some(form) = live_form("the D6 open-index values against Forcer") else {
        return;
    };
    let engine = p_battery_engine();
    for (numerator, topology, loops, indices) in CASES {
        let terms = battery_terms(loops);
        let mode = RustRedMasterMode::Symbolic;
        let epsilon = vakint_symbol!(
            rustred_p_settings(mode, LoopNormalizationFactor::GScheme, terms)
                .epsilon_symbol
                .as_str()
        );
        let open = rustred_value(&parse(&format!("({numerator})*{topology}")), mode, terms)
            .unwrap_or_else(|error| panic!("{numerator}: {error}"));
        let mut closings = 0;
        for closing in closing_tensors(indices) {
            let contracted = parse(&format!("({numerator})*({closing})*{topology}"));
            let forcer = forcer_side_value(engine, &form, &contracted, mode, terms)
                .unwrap_or_else(|error| panic!("{numerator} closed with {closing}: {error}"));
            assert_closed_case(
                &open,
                &parse(&closing),
                &forcer,
                epsilon,
                &format!("{numerator} {topology} vs Forcer"),
            );
            closings += 1;
        }
        eprintln!("D6 Forcer cross-check {loops}L {numerator}: {closings} closings equal");
    }
}
