use spenso::{
    chain,
    network::library::symbolic::ETS,
    s, slot,
    structure::{
        abstract_index::AbstractIndex,
        representation::{Minkowski, RepName, Representation},
    },
    symbolica_atom::IntoAtom,
    trace,
};
use symbolica::{atom::Atom, function, symbol};

use crate::{
    color::{CS, ColorSimplifier},
    dirac::{GammaSimplifier, GammaSimplifySettings},
    representations::{Bispinor, ColorAdjoint, ColorFundamental, initialize},
};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ReferenceDomain {
    Dirac4,
    ColorSu3,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ReferenceSimplification {
    GammaDefault,
    GammaCanonical,
    ColorDefault,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ReferenceValidation {
    Enabled,
    Blocked(&'static str),
}

impl ReferenceValidation {
    pub fn is_enabled(self) -> bool {
        matches!(self, Self::Enabled)
    }
}

pub struct ReferenceCase {
    pub name: &'static str,
    pub domain: ReferenceDomain,
    pub expression: fn() -> Atom,
    pub simplification: ReferenceSimplification,
    pub validation: ReferenceValidation,
}

impl ReferenceCase {
    pub fn expression(&self) -> Atom {
        (self.expression)()
    }

    pub fn simplify(&self, expression: &Atom) -> Atom {
        match self.simplification {
            ReferenceSimplification::GammaDefault => expression.simplify_gamma(),
            ReferenceSimplification::GammaCanonical => {
                expression.simplify_gamma_with(GammaSimplifySettings::canonical())
            }
            ReferenceSimplification::ColorDefault => expression.simplify_color(),
        }
    }

    pub fn simplified(&self) -> Atom {
        let expression = self.expression();
        self.simplify(&expression)
    }
}

struct ConcreteReps {
    mink4: Representation<Minkowski>,
    bis4: Representation<Bispinor>,
    cof3: Representation<ColorFundamental>,
    coad8: Representation<ColorAdjoint>,
}

impl ConcreteReps {
    fn new() -> Self {
        initialize();
        Self {
            mink4: Minkowski {}.new_rep(4),
            bis4: Bispinor {}.new_rep(4),
            cof3: ColorFundamental {}.new_rep(3),
            coad8: ColorAdjoint {}.new_rep(8),
        }
    }
}

fn metric(left: impl IntoAtom, right: impl IntoAtom) -> Atom {
    function!(ETS.metric, left.into_atom(), right.into_atom())
}

fn mom(name: &str, rep: &Representation<Minkowski>) -> Atom {
    function!(symbol!(format!("spenso::{name}")), rep.to_symbolic([]))
}

fn p4(r: &ConcreteReps) -> Atom {
    mom("p", &r.mink4)
}

fn q4(r: &ConcreteReps) -> Atom {
    mom("q", &r.mink4)
}

fn color_t(r: &ConcreteReps, index: symbolica::atom::Symbol) -> Atom {
    crate::color_t!(r.coad8.slot::<AbstractIndex, _>(index))
}

fn color_f(
    r: &ConcreteReps,
    a: symbolica::atom::Symbol,
    b: symbolica::atom::Symbol,
    c: symbolica::atom::Symbol,
) -> Atom {
    crate::color_f!(
        r.coad8.slot::<AbstractIndex, _>(a),
        r.coad8.slot::<AbstractIndex, _>(b),
        r.coad8.slot::<AbstractIndex, _>(c)
    )
}

fn color_chain(
    r: &ConcreteReps,
    start: symbolica::atom::Symbol,
    end: symbolica::atom::Symbol,
    indices: &[symbolica::atom::Symbol],
) -> Atom {
    let factors = indices
        .iter()
        .copied()
        .map(|index| color_t(r, index))
        .collect::<Vec<_>>();
    chain!(
        r.cof3.slot::<AbstractIndex, _>(start),
        r.cof3.dual().slot::<AbstractIndex, _>(end);
        factors
    )
}

fn color_trace(r: &ConcreteReps, indices: &[symbolica::atom::Symbol]) -> Atom {
    let factors = indices
        .iter()
        .copied()
        .map(|index| color_t(r, index))
        .collect::<Vec<_>>();
    trace!(&r.cof3; factors)
}

fn dirac_form_two_gamma_trace() -> Atom {
    initialize();
    crate::gamma!(mu, a, b) * crate::gamma!(nu, b, a)
}

fn dirac_form_odd_gamma_trace() -> Atom {
    initialize();
    crate::gamma!(mu, a, b) * crate::gamma!(nu, b, c) * crate::gamma!(rho, c, a)
}

fn dirac_form_four_gamma_trace() -> Atom {
    initialize();
    crate::gamma!(mu, a, b)
        * crate::gamma!(nu, b, c)
        * crate::gamma!(rho, c, d)
        * crate::gamma!(sigma, d, a)
}

fn dirac_feyncalc_open_chain_chisholm_id2() -> Atom {
    let r = ConcreteReps::new();
    chain!(
        slot!(r.bis4, a),
        slot!(r.bis4, b),
        crate::gamma!(slot!(r.mink4, mu)),
        crate::gamma!(slot!(r.mink4, nu)),
        crate::gamma!(slot!(r.mink4, rho)),
        crate::gamma!(slot!(r.mink4, sigma)),
        crate::gamma!(slot!(r.mink4, mu)),
    )
}

fn dirac_feyncalc_open_chain_chisholm_id3() -> Atom {
    let r = ConcreteReps::new();
    chain!(
        slot!(r.bis4, i),
        slot!(r.bis4, j),
        crate::gamma!(slot!(r.mink4, mu)),
        crate::gamma!(slot!(r.mink4, alpha)),
        crate::gamma!(slot!(r.mink4, beta)),
        crate::gamma!(slot!(r.mink4, rho)),
        crate::gamma!(slot!(r.mink4, sigma)),
        crate::gamma!(slot!(r.mink4, mu)),
    ) / 2
}

fn dirac_feyncalc_slash_sandwich_id4() -> Atom {
    let r = ConcreteReps::new();
    let p = p4(&r);
    let q = q4(&r);
    Atom::var(s!(m))
        * crate::gamma!(p.clone(), slot!(r.bis4, i), slot!(r.bis4, a))
        * crate::gamma!(p.clone(), slot!(r.bis4, a), slot!(r.bis4, j))
        - crate::gamma!(p.clone(), slot!(r.bis4, i), slot!(r.bis4, a))
            * crate::gamma!(q.clone(), slot!(r.bis4, a), slot!(r.bis4, b))
            * crate::gamma!(p.clone(), slot!(r.bis4, b), slot!(r.bis4, j))
}

fn dirac_feyncalc_gamma5_anticommutes_id5() -> Atom {
    let r = ConcreteReps::new();
    chain!(
        slot!(r.bis4, i),
        slot!(r.bis4, j),
        crate::gamma5!(),
        crate::gamma!(slot!(r.mink4, mu)),
    )
}

fn dirac_feyncalc_order_anticommutator_id30() -> Atom {
    let r = ConcreteReps::new();
    let p = p4(&r);
    chain!(
        slot!(r.bis4, i),
        slot!(r.bis4, j),
        crate::gamma!(slot!(r.mink4, nu)),
        crate::gamma!(p.clone()),
    ) + chain!(
        slot!(r.bis4, i),
        slot!(r.bis4, j),
        crate::gamma!(p.clone()),
        crate::gamma!(slot!(r.mink4, nu)),
    ) - Atom::num(2) * metric(slot!(r.mink4, nu), p) * chain!(slot!(r.bis4, i), slot!(r.bis4, j))
}

fn dirac_gamma0_square() -> Atom {
    let r = ConcreteReps::new();
    chain!(
        slot!(r.bis4, a),
        slot!(r.bis4, b),
        crate::gamma0!(),
        crate::gamma0!(),
    )
}

fn dirac_gamma0_conjugates_gamma5() -> Atom {
    let r = ConcreteReps::new();
    chain!(
        slot!(r.bis4, a),
        slot!(r.bis4, b),
        crate::gamma0!(),
        crate::gamma5!(),
        crate::gamma0!(),
    )
}

fn dirac_gamma5_four_gamma_trace_is_epsilon() -> Atom {
    let r = ConcreteReps::new();
    trace!(
        &r.bis4,
        crate::gamma5!(),
        crate::gamma!(slot!(r.mink4, mu)),
        crate::gamma!(slot!(r.mink4, nu)),
        crate::gamma!(slot!(r.mink4, rho)),
        crate::gamma!(slot!(r.mink4, sigma)),
    )
}

fn color_form_one_generator_trace() -> Atom {
    let r = ConcreteReps::new();
    color_trace(&r, &[s!(a)])
}

fn color_form_two_generator_trace() -> Atom {
    let r = ConcreteReps::new();
    color_trace(&r, &[s!(a), s!(b)])
}

fn color_form_fierz_generator_contraction() -> Atom {
    initialize();
    CS.t_pattern(3, 8, s!(a), s!(i), s!(j)) * CS.t_pattern(3, 8, s!(a), s!(k), s!(l))
}

fn color_feyncalc_open_chain_separated_casimir_id2() -> Atom {
    let r = ConcreteReps::new();
    color_chain(&r, s!(i), s!(j), &[s!(a), s!(b), s!(a)])
}

fn color_feyncalc_structure_loop_to_adjoint_delta_id3() -> Atom {
    let r = ConcreteReps::new();
    color_f(&r, s!(a), s!(c), s!(d)) * color_f(&r, s!(b), s!(c), s!(d))
}

fn color_feyncalc_structure_times_open_chain_id4() -> Atom {
    let r = ConcreteReps::new();
    color_f(&r, s!(a), s!(b), s!(c)) * color_chain(&r, s!(i), s!(j), &[s!(b), s!(c)])
}

fn color_feyncalc_delta_closes_doubled_two_generator_chain_id18() -> Atom {
    let r = ConcreteReps::new();
    let line = color_chain(&r, s!(a), s!(d), &[s!(i)]) * color_chain(&r, s!(d), s!(b), &[s!(j)]);
    metric(slot!(r.cof3, b), slot!(r.cof3.dual(), a)) * (line.clone() + line)
}

fn color_feyncalc_three_generator_trace_terminal_id30() -> Atom {
    let r = ConcreteReps::new();
    color_trace(&r, &[s!(i), s!(j), s!(k)])
}

fn color_feyncalc_repeated_trace_pair_id52() -> Atom {
    let r = ConcreteReps::new();
    color_trace(&r, &[s!(i1), s!(i2), s!(i1), s!(i2)])
}

fn color_form_four_generator_trace_terminal() -> Atom {
    let r = ConcreteReps::new();
    color_trace(&r, &[s!(a), s!(b), s!(c), s!(d)])
}

pub static REFERENCE_CASES: &[ReferenceCase] = &[
    ReferenceCase {
        name: "dirac_form_two_gamma_trace",
        domain: ReferenceDomain::Dirac4,
        expression: dirac_form_two_gamma_trace,
        simplification: ReferenceSimplification::GammaDefault,
        validation: ReferenceValidation::Enabled,
    },
    ReferenceCase {
        name: "dirac_form_odd_gamma_trace",
        domain: ReferenceDomain::Dirac4,
        expression: dirac_form_odd_gamma_trace,
        simplification: ReferenceSimplification::GammaDefault,
        validation: ReferenceValidation::Enabled,
    },
    ReferenceCase {
        name: "dirac_form_four_gamma_trace",
        domain: ReferenceDomain::Dirac4,
        expression: dirac_form_four_gamma_trace,
        simplification: ReferenceSimplification::GammaDefault,
        validation: ReferenceValidation::Enabled,
    },
    ReferenceCase {
        name: "dirac_feyncalc_open_chain_chisholm_id2",
        domain: ReferenceDomain::Dirac4,
        expression: dirac_feyncalc_open_chain_chisholm_id2,
        simplification: ReferenceSimplification::GammaDefault,
        validation: ReferenceValidation::Enabled,
    },
    ReferenceCase {
        name: "dirac_feyncalc_open_chain_chisholm_id3",
        domain: ReferenceDomain::Dirac4,
        expression: dirac_feyncalc_open_chain_chisholm_id3,
        simplification: ReferenceSimplification::GammaDefault,
        validation: ReferenceValidation::Enabled,
    },
    ReferenceCase {
        name: "dirac_feyncalc_slash_sandwich_id4",
        domain: ReferenceDomain::Dirac4,
        expression: dirac_feyncalc_slash_sandwich_id4,
        simplification: ReferenceSimplification::GammaDefault,
        validation: ReferenceValidation::Enabled,
    },
    ReferenceCase {
        name: "dirac_feyncalc_gamma5_anticommutes_id5",
        domain: ReferenceDomain::Dirac4,
        expression: dirac_feyncalc_gamma5_anticommutes_id5,
        simplification: ReferenceSimplification::GammaDefault,
        validation: ReferenceValidation::Enabled,
    },
    ReferenceCase {
        name: "dirac_feyncalc_order_anticommutator_id30",
        domain: ReferenceDomain::Dirac4,
        expression: dirac_feyncalc_order_anticommutator_id30,
        simplification: ReferenceSimplification::GammaCanonical,
        validation: ReferenceValidation::Enabled,
    },
    ReferenceCase {
        name: "dirac_gamma0_square",
        domain: ReferenceDomain::Dirac4,
        expression: dirac_gamma0_square,
        simplification: ReferenceSimplification::GammaDefault,
        validation: ReferenceValidation::Enabled,
    },
    ReferenceCase {
        name: "dirac_gamma0_conjugates_gamma5",
        domain: ReferenceDomain::Dirac4,
        expression: dirac_gamma0_conjugates_gamma5,
        simplification: ReferenceSimplification::GammaDefault,
        validation: ReferenceValidation::Enabled,
    },
    ReferenceCase {
        name: "dirac_gamma5_four_gamma_trace_is_epsilon",
        domain: ReferenceDomain::Dirac4,
        expression: dirac_gamma5_four_gamma_trace_is_epsilon,
        simplification: ReferenceSimplification::GammaDefault,
        validation: ReferenceValidation::Enabled,
    },
    ReferenceCase {
        name: "color_form_one_generator_trace",
        domain: ReferenceDomain::ColorSu3,
        expression: color_form_one_generator_trace,
        simplification: ReferenceSimplification::ColorDefault,
        validation: ReferenceValidation::Enabled,
    },
    ReferenceCase {
        name: "color_form_two_generator_trace",
        domain: ReferenceDomain::ColorSu3,
        expression: color_form_two_generator_trace,
        simplification: ReferenceSimplification::ColorDefault,
        validation: ReferenceValidation::Enabled,
    },
    ReferenceCase {
        name: "color_form_fierz_generator_contraction",
        domain: ReferenceDomain::ColorSu3,
        expression: color_form_fierz_generator_contraction,
        simplification: ReferenceSimplification::ColorDefault,
        validation: ReferenceValidation::Enabled,
    },
    ReferenceCase {
        name: "color_feyncalc_open_chain_separated_casimir_id2",
        domain: ReferenceDomain::ColorSu3,
        expression: color_feyncalc_open_chain_separated_casimir_id2,
        simplification: ReferenceSimplification::ColorDefault,
        validation: ReferenceValidation::Blocked(
            "hep-lib parser cannot yet materialize repeated open color-chain adjoint contractions",
        ),
    },
    ReferenceCase {
        name: "color_feyncalc_structure_loop_to_adjoint_delta_id3",
        domain: ReferenceDomain::ColorSu3,
        expression: color_feyncalc_structure_loop_to_adjoint_delta_id3,
        simplification: ReferenceSimplification::ColorDefault,
        validation: ReferenceValidation::Enabled,
    },
    ReferenceCase {
        name: "color_feyncalc_structure_times_open_chain_id4",
        domain: ReferenceDomain::ColorSu3,
        expression: color_feyncalc_structure_times_open_chain_id4,
        simplification: ReferenceSimplification::ColorDefault,
        validation: ReferenceValidation::Blocked(
            "hep-lib parser cannot yet materialize repeated open color-chain adjoint contractions",
        ),
    },
    ReferenceCase {
        name: "color_feyncalc_delta_closes_doubled_two_generator_chain_id18",
        domain: ReferenceDomain::ColorSu3,
        expression: color_feyncalc_delta_closes_doubled_two_generator_chain_id18,
        simplification: ReferenceSimplification::ColorDefault,
        validation: ReferenceValidation::Blocked(
            "hep-lib parser cannot yet materialize repeated open color-chain endpoint contractions",
        ),
    },
    ReferenceCase {
        name: "color_feyncalc_three_generator_trace_terminal_id30",
        domain: ReferenceDomain::ColorSu3,
        expression: color_feyncalc_three_generator_trace_terminal_id30,
        simplification: ReferenceSimplification::ColorDefault,
        validation: ReferenceValidation::Blocked(
            "hep-lib parser cannot yet materialize length-three color traces reliably",
        ),
    },
    ReferenceCase {
        name: "color_feyncalc_repeated_trace_pair_id52",
        domain: ReferenceDomain::ColorSu3,
        expression: color_feyncalc_repeated_trace_pair_id52,
        simplification: ReferenceSimplification::ColorDefault,
        validation: ReferenceValidation::Blocked(
            "hep-lib parser cannot yet materialize repeated color-trace adjoint contractions",
        ),
    },
    ReferenceCase {
        name: "color_form_four_generator_trace_terminal",
        domain: ReferenceDomain::ColorSu3,
        expression: color_form_four_generator_trace_terminal,
        simplification: ReferenceSimplification::ColorDefault,
        validation: ReferenceValidation::Blocked(
            "hep-lib parser cannot yet materialize higher color traces reliably",
        ),
    },
];

pub fn reference_cases() -> &'static [ReferenceCase] {
    REFERENCE_CASES
}

#[cfg(test)]
mod tests {
    use insta::assert_snapshot;
    use spenso::shadowing::symbolica_utils::AtomCoreExt;

    use super::*;

    #[test]
    fn reference_case_simplification_snapshots() {
        let snapshot = reference_cases()
            .iter()
            .map(|case| {
                let validation = match case.validation {
                    ReferenceValidation::Enabled => "enabled".to_string(),
                    ReferenceValidation::Blocked(reason) => format!("blocked: {reason}"),
                };
                format!(
                    "{} [{:?}, {:?}, {}]\n{}",
                    case.name,
                    case.domain,
                    case.simplification,
                    validation,
                    case.simplified().to_bare_ordered_string(),
                )
            })
            .collect::<Vec<_>>()
            .join("\n\n");

        assert_snapshot!("reference_case_simplifications", snapshot);
    }
}
