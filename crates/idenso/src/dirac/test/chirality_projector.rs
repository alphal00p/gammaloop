//! Regression coverage for the chirality-projector expansion
//! (`ℙ± = ½(𝟙 ± γ5)`) used to collapse closed fermion loops.
//!
//! The expansion is intrinsically four-dimensional: every downstream `gamma5`
//! rule is `FourDimensional`, so away from four dimensions an expanded
//! projector degrades into an irreducible `gamma5` polynomial. These tests pin
//! both halves of that contract — it fires in 4D, and it must not fire in D.
use super::*;

use spenso::shadowing::IntoAtom;
use symbolica::atom::{FunctionBuilder, Symbol};

fn projector(symbol: Symbol, a: Atom, b: Atom) -> Atom {
    FunctionBuilder::new(symbol).add_arg(a).add_arg(b).finish()
}

fn projp(a: Atom, b: Atom) -> Atom {
    projector(AGS.projp, a, b)
}

fn projm(a: Atom, b: Atom) -> Atom {
    projector(AGS.projm, a, b)
}

// ------------------------------------------------------------------------
// Dimension guard: a D-dimensional projector must be left alone.
// ------------------------------------------------------------------------

#[test]
fn d_dimensional_projector_is_not_expanded() {
    let r = test_initialize();
    let expr = projp(
        slot!(r.bis_d, a).into_atom(),
        slot!(r.bis_d, b).into_atom(),
    );

    assert_snapshot!(expr.simplify_gamma().to_bare_ordered_string(), @"chain(bis(d,a),bis(d,b),projp(in,out))");
}

#[test]
fn d_dimensional_projector_pair_stays_projector_valued() {
    let r = test_initialize();
    // ℙ₊ℙ₊ = ℙ₊ is a 4D identity, so in D dimensions the product simply stays
    // inert. It must *not* become a `chain(gamma5,gamma5)` that no rule can
    // collapse.
    let expr = projp(slot!(r.bis_d, a).into_atom(), slot!(r.bis_d, b).into_atom())
        * projp(slot!(r.bis_d, b).into_atom(), slot!(r.bis_d, c).into_atom());

    assert_snapshot!(expr.simplify_gamma().to_bare_ordered_string(), @"chain(bis(d,a),bis(d,c),projp(in,out),projp(in,out))");
}

#[test]
fn d_dimensional_open_projector_chain_stays_inert() {
    let r = test_initialize();
    // γ^μ ℙ₊ γ^ν on an open D-dimensional spin line. An injected γ5 could never
    // be anticommuted out, so it would permanently block the two gammas from
    // meeting.
    let expr = gamma!(slot!(r.mink_d, mu), slot!(r.bis_d, a), slot!(r.bis_d, b))
        * projp(slot!(r.bis_d, b).into_atom(), slot!(r.bis_d, c).into_atom())
        * gamma!(slot!(r.mink_d, nu), slot!(r.bis_d, c), slot!(r.bis_d, e));

    assert_snapshot!(expr.simplify_gamma().to_bare_ordered_string(), @"chain(bis(d,a),bis(d,e),gamma(in,out,mink(d,mu)),projp(in,out),gamma(in,out,mink(d,nu)))");
}

#[test]
fn d_dimensional_closed_loop_with_projector_stays_a_single_trace() {
    let r = test_initialize();
    // Tr(γ^μ ℙ₊ γ^ν) in D dimensions: `TRACE_GAMMA5_RECURSION` is 4D-only, so
    // expanding here would only split one inert trace into two.
    let expr = gamma!(slot!(r.mink_d, mu), slot!(r.bis_d, a), slot!(r.bis_d, b))
        * projp(slot!(r.bis_d, b).into_atom(), slot!(r.bis_d, c).into_atom())
        * gamma!(slot!(r.mink_d, nu), slot!(r.bis_d, c), slot!(r.bis_d, a));

    assert_snapshot!(expr.simplify_gamma().to_bare_ordered_string(), @"trace(bis(d),cyclic(gamma(in,out,mink(d,nu)),gamma(in,out,mink(d,mu)),projp(in,out)))");
}

// ------------------------------------------------------------------------
// Four-dimensional behaviour the expansion exists for.
// ------------------------------------------------------------------------

#[test]
fn four_dimensional_open_projector_chain_expands() {
    let r = test_initialize();
    let expr = gamma!(slot!(r.mink4, mu), slot!(r.bis4, a), slot!(r.bis4, b))
        * projp(slot!(r.bis4, b).into_atom(), slot!(r.bis4, c).into_atom())
        * gamma!(slot!(r.mink4, nu), slot!(r.bis4, c), slot!(r.bis4, e));

    assert_snapshot!(expr.simplify_gamma().to_bare_ordered_string(), @"-1/2*chain(bis(4,a),bis(4,e),gamma(in,out,mink(4,mu)),gamma(in,out,mink(4,nu)),gamma5(in,out))+1/2*chain(bis(4,a),bis(4,e),gamma(in,out,mink(4,mu)),gamma(in,out,mink(4,nu)))");
}

#[test]
fn four_dimensional_projector_pair_contracts() {
    let r = test_initialize();
    // ℙ₊ℙ₋ = 0 must still be reached before any expansion turns the pair into
    // a γ5 polynomial.
    let expr = projp(slot!(r.bis4, a).into_atom(), slot!(r.bis4, b).into_atom())
        * projm(slot!(r.bis4, b).into_atom(), slot!(r.bis4, c).into_atom());

    assert_snapshot!(expr.simplify_gamma().to_bare_ordered_string(), @"0");
}

#[test]
fn pre_chained_projector_closed_loop_reduces_through_trace() {
    let r = test_initialize();
    // Same loop as `projector_closed_loop_reduces_through_trace`, but handed to
    // the simplifier already chained — the shape `reduce_bridge` used to build.
    // It must reach the identical Lorentz structure.
    let expr = gamma!(slot!(r.mink4, mu), slot!(r.bis4, 1), slot!(r.bis4, 2))
        * projp(slot!(r.bis4, 2).into_atom(), slot!(r.bis4, 3).into_atom())
        * gamma!(slot!(r.mink4, nu), slot!(r.bis4, 3), slot!(r.bis4, 4))
        * gamma!(slot!(r.mink4, rho), slot!(r.bis4, 4), slot!(r.bis4, 5))
        * gamma!(slot!(r.mink4, sigma), slot!(r.bis4, 5), slot!(r.bis4, 1));

    let direct = expr.simplify_gamma().expand().to_bare_ordered_string();
    let pre_chained = expr
        .collect_gamma_chains()
        .simplify_gamma()
        .expand()
        .to_bare_ordered_string();

    assert_eq!(direct, pre_chained);
    assert_snapshot!(pre_chained, @"-2*epsilon(mink(4,mu),mink(4,nu),mink(4,rho),mink(4,sigma))+-2*g(mink(4,mu),mink(4,rho))*g(mink(4,nu),mink(4,sigma))+2*g(mink(4,mu),mink(4,nu))*g(mink(4,rho),mink(4,sigma))+2*g(mink(4,mu),mink(4,sigma))*g(mink(4,nu),mink(4,rho))");
}

#[test]
fn projector_simplification_is_idempotent() {
    let r = test_initialize();
    // A second `simplify_gamma()` must be a no-op: `reduce_bridge` runs after
    // `serialization` has already gamma-simplified, so a non-idempotent pass
    // would silently keep rewriting the numerator.
    let expr = gamma!(slot!(r.mink4, mu), slot!(r.bis4, a), slot!(r.bis4, b))
        * projp(slot!(r.bis4, b).into_atom(), slot!(r.bis4, c).into_atom())
        * gamma!(slot!(r.mink4, nu), slot!(r.bis4, c), slot!(r.bis4, e));

    let once = expr.simplify_gamma();
    let twice = once.simplify_gamma();

    assert_eq!(
        once.to_bare_ordered_string(),
        twice.to_bare_ordered_string()
    );
}

// ------------------------------------------------------------------------
// `without_trace_evaluation()` must keep closed loops inert.
// ------------------------------------------------------------------------

#[test]
fn projector_closed_loop_is_inert_without_trace_evaluation() {
    let r = test_initialize();
    // Tr(γ^μ ℙ₊ γ^ν) in 4D. A caller who asked for no trace evaluation must not
    // get a `trace(...)` evaluated behind its back, nor a projector expansion
    // that only exists to feed the trace evaluator.
    let expr = gamma!(slot!(r.mink4, mu), slot!(r.bis4, a), slot!(r.bis4, b))
        * projp(slot!(r.bis4, b).into_atom(), slot!(r.bis4, c).into_atom())
        * gamma!(slot!(r.mink4, nu), slot!(r.bis4, c), slot!(r.bis4, a));

    assert_snapshot!(
        expr.simplify_gamma_with(
            GammaSimplifySettings::repeated_pairs().without_trace_evaluation()
        )
        .to_bare_ordered_string(),
        @"trace(bis(4),cyclic(gamma(in,out,mink(4,nu)),gamma(in,out,mink(4,mu)),projp(in,out)))"
    );
}

#[test]
fn projector_position_flips_the_epsilon_sign() {
    let r = test_initialize();
    // γ^μ anticommutes past γ5, so γ^μ ℙ₊ = ℙ₋ γ^μ: moving the projector across
    // a single gamma flips the sign of the ε term. This pins the convention the
    // `projector_closed_loop_reduces_through_trace` label depends on.
    let projector_first = projp(slot!(r.bis4, 1).into_atom(), slot!(r.bis4, 2).into_atom())
        * gamma!(slot!(r.mink4, mu), slot!(r.bis4, 2), slot!(r.bis4, 3))
        * gamma!(slot!(r.mink4, nu), slot!(r.bis4, 3), slot!(r.bis4, 4))
        * gamma!(slot!(r.mink4, rho), slot!(r.bis4, 4), slot!(r.bis4, 5))
        * gamma!(slot!(r.mink4, sigma), slot!(r.bis4, 5), slot!(r.bis4, 1));

    // Tr(ℙ₊ γ^μ γ^ν γ^ρ γ^σ) = 2(g^{μν}g^{ρσ} − g^{μρ}g^{νσ} + g^{μσ}g^{νρ}) + 2 ε^{μνρσ}.
    assert_snapshot!(
        projector_first.simplify_gamma().expand().to_bare_ordered_string(),
        @"-2*g(mink(4,mu),mink(4,rho))*g(mink(4,nu),mink(4,sigma))+2*epsilon(mink(4,mu),mink(4,nu),mink(4,rho),mink(4,sigma))+2*g(mink(4,mu),mink(4,nu))*g(mink(4,rho),mink(4,sigma))+2*g(mink(4,mu),mink(4,sigma))*g(mink(4,nu),mink(4,rho))"
    );
}
