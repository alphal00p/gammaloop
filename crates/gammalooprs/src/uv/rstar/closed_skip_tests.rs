//! `GAMMALOOP_CLOSED_SKIP_SUMMAND_SIMPLIFY` spot checks on small closed
//! remainders: on the closing-root path every ordinary-UV remainder of the
//! colour-factored two- and three-loop background diagrams must be certified
//! a closed scalar by idenso when it is contracted, and its outer Taylor
//! polynomial must be a byte-for-byte fixed point of the summand
//! simplification the skip omits.

use eyre::Result;
use symbolica::atom::AtomCore;
use vakint::Vakint;

use crate::{
    feyngen::feynkit::FeynmanDiagramGammaLoopExt, initialisation::test_initialise,
    uv::approx::streamed::with_production_piece_limit,
};

use super::{
    RStarPlan,
    expression::{RStarIntegrand, closed_dot_scalar},
    recursive_tests::{colour_factored, load, settings},
};

#[test]
fn closed_root_taylor_polynomials_are_certified_fixed_points_of_the_summand_simplification()
-> Result<()> {
    // Production piece limit: polynomials below it are simplified whole, as
    // in production (the test limit of one byte would stream every summand).
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let vakint = Vakint::new()?;
        let settings = settings();
        let mut checked = 0;
        let mut uncertified = 0;
        let mut with_children = 0;
        let mut certified_loops = std::collections::BTreeSet::new();
        for (name, dot) in [
            ("L2/0", include_str!("testdata/l2-d0.dot")),
            ("L2/6", include_str!("testdata/l2-d6.dot")),
            ("L2/7", include_str!("testdata/l2-d7.dot")),
            ("L2/8", include_str!("testdata/l2-d8.dot")),
            ("L2/12", include_str!("testdata/l2-d12.dot")),
            ("L3/166", include_str!("testdata/l3-d166.dot")),
            ("L3/185", include_str!("testdata/l3-d185.dot")),
            ("L3/241", include_str!("testdata/l3-d241.dot")),
        ] {
            let diagram = colour_factored(&load(dot)?)?;
            let graph = diagram.to_gamma_loop_graph(None, true)?;
            let model = diagram.model();
            let plan = RStarPlan::with_automatic_edge(&graph, model, &settings)?;
            let polynomials =
                plan.closing_root_outer_taylor_polynomials(&graph, model, &vakint, &settings)?;
            assert!(!polynomials.is_empty(), "{name}: no closed remainder");
            with_children += plan
                .ordinary_uv
                .iter()
                .filter(|term| !term.counterterms.is_empty())
                .count();
            for (index, (loops, contracted, certified, polynomial)) in
                polynomials.into_iter().enumerate()
            {
                if !certified {
                    // No producer certificate: production simplifies this
                    // polynomial as before. Report what it is.
                    eprintln!(
                        "{name} remainder {index}: uncertified, loops={loops} contracted_bytes={} closed_dot_scalar={} contracted={}",
                        contracted.as_view().get_byte_size(),
                        closed_dot_scalar(contracted.as_view()),
                        contracted
                            .to_canonical_string()
                            .chars()
                            .take(300)
                            .collect::<String>()
                    );
                    uncertified += 1;
                    continue;
                }
                assert!(
                    closed_dot_scalar(polynomial.as_view()),
                    "{name} remainder {index}: a certified outer Taylor polynomial fails the syntactic closed-scalar test"
                );
                let simplified = RStarIntegrand::simplify_summands(polynomial.as_view())?;
                // Atom equality compares the normalized byte representation.
                assert!(
                    simplified == polynomial,
                    "{name} remainder {index}: the summand simplification changed a certified polynomial\nbefore: {polynomial}\nafter: {simplified}"
                );
                checked += 1;
                certified_loops.insert(loops);
            }
        }
        eprintln!(
            "closed-skip spot checks: {checked} certified polynomials (remainder loops {certified_loops:?}), {uncertified} uncertified, {with_children} terms with children"
        );
        assert!(
            checked >= 8 && with_children > 0,
            "{checked} {with_children}"
        );
        Ok(())
    })
}
