use std::ops::Neg;

use color_eyre::Result;
use eyre::eyre;
use gammaloop_tracing_filter::debug_instrument;
use idenso::{
    color::ColorSimplifier,
    dirac::GammaSimplifier,
    representations::Bispinor,
    shorthands::{
        chain::Chain,
        metric::MetricSimplifier,
        schoonschip::{Schoonschip, SchoonschipSettings},
    },
};

use linnet::half_edge::subgraph::{Inclusion, SuBitGraph, SubGraphLike, SubSetLike};
use spenso::shadowing::TensorCollectExt;
use symbolica::{
    atom::{AtomCore, AtomView},
    prelude::*,
};

use crate::{
    debug_tags,
    graph::{FourDDenominator, Graph, LMBext, LoopMomentumBasis},
    numerator::aind::Aind,
    utils::{GS, W_},
    uv::{
        ApproximationType, UltravioletGraph,
        approx::{ForestNodeLike, Rooted, UVCtx, integrated::IntegratedCts},
        marker::{UvMarker, UvOperation},
    },
};
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) struct Local4dCts(Atom);
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) struct Full4dCts(Atom);

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct FourDTerm {
    pub(crate) numerator: Atom,
    pub(crate) denominators: Vec<FourDDenominator>,
}

impl Full4dCts {
    pub(crate) fn atom(&self) -> &Atom {
        &self.0
    }

    /// Pole-part subtraction inserts only completed poles; MUV keeps the local
    /// counterterm together with its integrated finite contribution.
    pub(crate) fn recursion_input(
        local: &Local4dCts,
        integrated: &IntegratedCts,
        scheme: ApproximationType,
        is_root: bool,
    ) -> Result<Self> {
        match scheme {
            ApproximationType::MUV => Ok(Self(&local.0 + integrated.finite_counterterm_atom())),
            ApproximationType::PolePart if !is_root => Ok(Self(integrated.pole_atom())),
            ApproximationType::PolePart => Ok(Self(&local.0 + integrated.pole_atom())),
            scheme => Err(eyre!("No recursive counterterm projection for {scheme}")),
        }
    }

    /// A disconnected local value already contains every component's recursive
    /// projection, so it must not be projected again when used by an outer limit.
    pub(crate) fn from_factorized_local(local: &Local4dCts) -> Self {
        Self(local.0.clone())
    }

    /// Complete one typed local term with the propagators outside its owning
    /// spinney. Numerator factors outside the spinney are deliberately not
    /// attached here; final assembly grows those under each exact source-local
    /// energy map after the complete 4D forest has been built.
    pub(crate) fn with_cograph<S: SubGraphLike>(
        local: &Local4dCts,
        graph: &Graph,
        cograph: &S,
    ) -> Self {
        Self(&local.0 * graph.denominator(cograph, |_| -1))
    }

    pub(crate) fn from_coefficient<S: SubGraphLike>(
        coefficient: &Atom,
        graph: &Graph,
        cograph: &S,
    ) -> Self {
        Self(coefficient * graph.denominator(cograph, |_| -1))
    }

    pub(crate) fn terms(&self) -> Result<Vec<FourDTerm>> {
        FourDTerm::from_view(self.0.as_view())
    }
}

impl Neg for Local4dCts {
    type Output = Local4dCts;
    fn neg(self) -> Self::Output {
        Local4dCts(-&self.0)
    }
}

impl Local4dCts {
    pub(crate) fn atom(&self) -> &Atom {
        &self.0
    }

    pub(crate) fn from_full_product(factors: impl IntoIterator<Item = Full4dCts>) -> Self {
        Self(
            factors
                .into_iter()
                .fold(Atom::one(), |product, factor| product * factor.0),
        )
    }

    /// Enumerate only the additive terms needed by the 3D projection while
    /// retaining every non-denominator factor as a Symbolica atom. In
    /// particular, numerator products and powers are not materialized into a
    /// parallel expression tree or expanded polynomial.
    #[cfg(test)]
    pub(crate) fn terms(&self) -> Result<Vec<FourDTerm>> {
        FourDTerm::from_view(self.0.as_view())
    }
}

impl FourDTerm {
    fn numerator(numerator: Atom) -> Self {
        Self {
            numerator,
            denominators: Vec::new(),
        }
    }

    fn product(mut left: Self, right: Self) -> Self {
        left.numerator *= right.numerator;
        left.denominators.extend(right.denominators);
        left
    }

    fn from_view(view: AtomView<'_>) -> Result<Vec<Self>> {
        match view {
            AtomView::Add(add) => {
                let terms = add
                    .iter()
                    .map(Self::from_view)
                    .collect::<Result<Vec<_>>>()
                    .map(|terms| terms.into_iter().flatten().collect::<Vec<_>>())?;
                if terms.iter().all(|term| term.denominators.is_empty()) {
                    Ok(vec![Self::numerator(view.to_owned())])
                } else {
                    Ok(terms)
                }
            }
            AtomView::Mul(mul) => {
                let mut terms = vec![Self::numerator(Atom::one())];
                for factor in mul.iter() {
                    let factor_terms = Self::from_view(factor)?;
                    terms = terms
                        .into_iter()
                        .flat_map(|left| {
                            factor_terms
                                .iter()
                                .cloned()
                                .map(move |right| Self::product(left.clone(), right))
                        })
                        .collect();
                }
                Ok(terms)
            }
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                let Some(denominator) = FourDDenominator::from_view(base)? else {
                    return Ok(vec![Self::numerator(view.to_owned())]);
                };
                let Ok(exponent) = i64::try_from(exponent) else {
                    return Err(eyre!(
                        "4D denominator has non-integer power `{}`",
                        exponent.to_owned()
                    ));
                };
                if exponent >= 0 {
                    return Ok(vec![Self::numerator(view.to_owned())]);
                }
                let multiplicity = usize::try_from(exponent.unsigned_abs())
                    .map_err(|_| eyre!("4D denominator multiplicity does not fit in memory"))?;
                Ok(vec![Self {
                    numerator: Atom::one(),
                    denominators: std::iter::repeat_n(denominator, multiplicity).collect(),
                }])
            }
            _ => Ok(vec![Self::numerator(view.to_owned())]),
        }
    }
}

impl FourDDenominator {
    fn from_view(view: AtomView<'_>) -> Result<Option<Self>> {
        let AtomView::Fun(function) = view else {
            return Ok(None);
        };
        if function.get_symbol() != GS.den {
            return Ok(None);
        }
        if function.get_nargs() != 4 {
            return Err(eyre!(
                "expected a 4D denominator wrapper with four arguments, found {}",
                function.get_nargs()
            ));
        }
        let source_edge = linnet::half_edge::involution::EdgeIndex(
            usize::try_from(function.get(0)).map_err(|_| {
                eyre!(
                    "4D denominator wrapper has non-integer edge id `{}`",
                    function.get(0).to_owned()
                )
            })?,
        );
        Ok(Some(Self {
            source_edge,
            momentum: function.get(1).to_owned(),
            mass_squared: function.get(2).to_owned(),
            full_expr: function.get(3).to_owned(),
        }))
    }
}

impl Rooted for Local4dCts {
    fn root() -> Self {
        Self(Atom::one())
    }
}

pub struct Local4d;

impl Graph {
    pub(crate) fn uv_rescaled(
        &self,
        replacement_subgraph: &SuBitGraph,
        n_loops: usize,
        lmb: &LoopMomentumBasis,
        atom: &Atom,
    ) -> Atom {
        // only apply replacements for edges in the reduced graph
        let mom_reps = self.uv_wrapped_replacement(replacement_subgraph, lmb, &[W_.x___]);
        let mut atomarg = atom.replace_multiple(&mom_reps);

        // rescale the loop momenta in the whole subgraph, including previously expanded cycles
        for edge in &lmb.loop_edges {
            atomarg = atomarg
                .replace(GS.emr_mom(*edge, W_.x___))
                .with(GS.emr_mom(*edge, W_.x___) / GS.rescale);
        }

        // Free `mUVexp` occurrences are left untouched here: with the
        // inverse loop-momentum expansion, soft dependence is generated by
        // the shifted hard momenta and denominator expansion rather than by
        // an explicit soft rescaling pass.

        // Free `mUV` occurrences come from previously integrated CT coefficients.
        // They are hard vacuum scales in the parent UV expansion and therefore
        // carry the same weight as the parent loop momentum. The denominator
        // rewrite below introduces the unscaled vacuum mass for the parent
        // propagator basis.
        atomarg = atomarg
            .replace(GS.m_uv_vacuum)
            .with(Atom::var(GS.m_uv_vacuum) / GS.rescale);

        // The stored scale restores consumed loop measures independently of mUV. It is set to
        // one when the next integration consumes it.
        atomarg = atomarg
            .replace(GS.integrated_loop_scale)
            .with(Atom::var(GS.integrated_loop_scale) * GS.rescale);

        let tsquare = Atom::var(GS.rescale).pow(2);
        let m_uv_expansion_sq = Atom::var(GS.m_uv_expansion).pow(2);
        let m_uv_vacuum_sq = Atom::var(GS.m_uv_vacuum).pow(2);

        debug_tags!(#uv, #integrated, #inspect;
            log.res = atomarg,
            "Rescaled momenta expanded"
        );
        atomarg = atomarg
            .replace(GS.den(W_.a_, W_.mom_, W_.mass_, W_.prop_))
            .with_map(move |m| {
                let edge = m.get(W_.a_).unwrap().to_atom();
                let momentum = m.get(W_.mom_).unwrap().to_atom();
                let mass = m.get(W_.mass_).unwrap().to_atom();
                let propagator = m.get(W_.prop_).unwrap().to_atom();

                let (mass, propagator) = if mass == m_uv_expansion_sq {
                    // The MUV deformation was already introduced by an inner UV limit.
                    // Rescale it without adding the vacuum mass a second time.
                    (mass, propagator * &tsquare)
                } else {
                    (
                        mass * &tsquare + &m_uv_expansion_sq,
                        propagator * &tsquare + &m_uv_expansion_sq * &tsquare - &m_uv_vacuum_sq,
                    )
                };

                GS.den(edge, momentum, mass, propagator) / &tsquare
            })
            .replace(function!(GS.den, W_.a_, W_.mom_, W_.a___))
            .with_map(move |m| {
                let mut f = symbolica::atom::FunctionBuilder::new(GS.den);
                f = f.add_arg(m.get(W_.a_).unwrap().to_atom());
                f = f.add_arg(
                    (m.get(W_.mom_).unwrap().to_atom() * GS.rescale)
                        .expand()
                        .replace(GS.rescale)
                        .with(Atom::Zero),
                );
                f = f.add_arg(m.get(W_.a___).unwrap().to_atom());
                f.finish()
            });

        atomarg *= Atom::var(GS.rescale).pow(-4 * n_loops as i64);
        atomarg
    }
}

/// Add the numerator of the reduced subgraph, (without given), to the integrand.
/// Then, 4d -> d-dim on minkowski indices
#[debug_instrument(
        current = %current.log_display(),
        given = %given.log_display(),
    )]
fn grow<S: super::ForestNodeLike>(
    integrand: &Full4dCts,
    ctx: &UVCtx<'_>,
    current: &S,
    given: &S,
) -> Result<Atom> {
    let reduced = current.reduced_subgraph(given);
    let graph = ctx.graph;

    let mut t_arg = ctx
        .graph
        .numerator(&reduced, given.subgraph())
        .to_d_dim(GS.dim)
        .get_single_atom()
        .unwrap();

    t_arg /= graph.denominator(&reduced, |_| 1);
    let integrand = (integrand.atom() * t_arg).simplify_metrics();

    debug_tags!(#uv, #integrated, #algebra, #start; log.integrand = integrand, reduced = %reduced.string_label());

    Ok(integrand)
}

#[debug_instrument(
        current = %current.log_display(),
        given = %given.log_display(),
    )]
fn t<S: super::ForestNodeLike>(
    integrand: &Atom,
    ctx: &UVCtx<'_>,
    current: &S,
    given: &S,
) -> Result<Local4dCts> {
    let graph = ctx.graph;
    let reduced = current.reduced_subgraph(given);
    let n_loops = graph.n_loops(current.subgraph());

    // Keep the child loops as a subset of the outer basis so integrating the child
    // commutes with the outer Taylor expansion.
    let is_compatible = |candidate: &LoopMomentumBasis| {
        given
            .lmb()
            .loop_edges
            .iter()
            .all(|edge| candidate.loop_edges.contains(edge))
            && candidate
                .loop_edges
                .iter()
                .filter(|edge| given.subgraph().includes(&graph[*edge].1))
                .count()
                == given.lmb().loop_edges.len()
    };
    let generated_lmb;
    let lmb = if given.subgraph().is_empty() || is_compatible(current.lmb()) {
        current.lmb()
    } else {
        generated_lmb = graph
            .generate_loop_momentum_bases_of(current.subgraph())
            .into_iter()
            .find(is_compatible)
            .ok_or_else(|| eyre!("no loop momentum basis compatible with nested UV subgraph"))?;
        &generated_lmb
    };

    let rescaled = graph.uv_rescaled(&reduced, n_loops, lmb, integrand);
    debug_tags!(#uv,#integrated,#rescaled;log.res = rescaled, n_loops=%n_loops,"Rescaled expanded");

    let series = rescaled
        .series(GS.rescale, Atom::Zero, 0)
        .unwrap()
        .to_atom();
    debug_tags!(#uv,#integrated, #series;log.res = series, "Series expanded");

    let evalutated = series.replace(GS.rescale).with(Atom::num(1));
    debug_tags!(#uv,#integrated,#series;log.res = evalutated, "Evaluated at t = 1");

    let collected = evalutated
        .simplify_metrics()
        .collect_rep((Bispinor {}).into())
        .collect_gamma_chains();
    debug_tags!(#uv,#integrated,#collect;log.expr = collected, "After gamma chain collection");

    let schoonschip = collected
        .schoonschip_with_settings(&SchoonschipSettings {
            simplify_chain_like_functions: true,
            schoonschip_rank1_tensors: true,
            ..Default::default()
        })
        .normalize_chains();
    debug_tags!(#uv, #integrated, #profile, #trace, #start, #collect;
        log.expr = schoonschip,
        "After gamma schoonschip"
    );
    let collected = schoonschip
        .collect_chains_and_traces()
        .simplify_metrics()
        .collect_gamma_chains()
        .collect_color()
        .collect_factors();
    debug_tags!(#uv, #integrated, #profile, #trace, #start, #collect;
        log.expr = collected,
        "After gamma collection"
    );

    let simplified = collected.simplify_gamma();
    debug_tags!(#uv, #integrated, #vakint, #profile, #trace, #start, #gamma;
        log.expr = simplified,
        "After gamma simplification"
    );
    let schoonschipped = simplified.schoonschip_net::<Aind>();
    debug_tags!(#uv, #integrated, #vakint, #profile, #trace,#schoonschip, #start;
        log.expr = schoonschipped,
        "After Schoonschip net"
    );
    let dotted = schoonschipped.to_dots().normalize_dots();
    debug_tags!(#uv, #integrated, #vakint, #profile, #trace, #dots;
        log.expr = dotted,
        "After dots"
    );

    Ok(Local4dCts(dotted))
}

pub(crate) fn uv_limit<S: ForestNodeLike, M: ForestNodeLike>(
    integrand: &Full4dCts,
    ctx: &UVCtx<'_>,
    current: &S,
    given: &S,
    marker_current: &M,
    marker_given: &M,
) -> Result<Local4dCts> {
    match current.renormalization_scheme() {
        ApproximationType::MUV | ApproximationType::PolePart => {
            let grown = grow(integrand, ctx, current, given)?;
            let result = -t(&grown, ctx, current, given)?;
            Ok(Local4dCts(UvMarker::new(ctx.settings).apply(
                UvOperation::Approx,
                marker_current.subgraph(),
                marker_given.subgraph(),
                result.atom(),
            )))
        }
        atype => Err(eyre!("Not yet implemented {:?}", atype)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::initialisation::test_initialise;
    use symbolica::{function, symbol};

    #[test]
    fn term_projection_keeps_factorized_numerator_atoms() -> Result<()> {
        test_initialise()?;
        let affine = symbol!("local_4d_test::Affine");
        let k0 = symbol!("local_4d_test::K0");
        let k1 = symbol!("local_4d_test::K1");
        let sum = function!(affine, Atom::var(k0)) + function!(affine, Atom::var(k1));
        let numerator = sum.pow(2) * (Atom::var(k0) + Atom::var(k1));

        let terms = Local4dCts(numerator.clone()).terms()?;
        assert_eq!(terms.len(), 1);
        assert_eq!(terms[0].numerator, numerator);
        assert!(terms[0].denominators.is_empty());
        Ok(())
    }

    #[test]
    fn term_projection_preserves_dots_and_distinct_same_edge_expressions() -> Result<()> {
        test_initialise()?;
        let first_full = Atom::var(symbol!("local_4d_test::first_full"));
        let second_full = Atom::var(symbol!("local_4d_test::second_full"));
        let first = GS.den(
            0,
            FunctionBuilder::new(GS.emr_mom).add_arg(0).finish(),
            0,
            &first_full,
        );
        let second = GS.den(
            0,
            FunctionBuilder::new(GS.emr_mom).add_arg(0).finish(),
            0,
            &second_full,
        );
        let source = Full4dCts(first.pow(-2) * second.pow(-1));

        let terms = source.terms()?;
        assert_eq!(terms.len(), 1);
        assert_eq!(terms[0].denominators.len(), 3);
        assert_eq!(
            terms[0]
                .denominators
                .iter()
                .filter(|denominator| denominator.full_expr == first_full)
                .count(),
            2
        );
        assert_eq!(
            terms[0]
                .denominators
                .iter()
                .filter(|denominator| denominator.full_expr == second_full)
                .count(),
            1
        );
        Ok(())
    }
}
