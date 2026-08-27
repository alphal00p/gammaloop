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
pub(crate) struct FourDSector {
    atom: Atom,
    frozen_lmbs: Vec<LoopMomentumBasis>,
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
struct FourDSectors {
    active: Vec<FourDSector>,
    recursive_completion: Vec<FourDSector>,
    atom: Atom,
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) struct Local4dCts(FourDSectors);
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) struct Full4dCts(FourDSectors);

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct FourDTerm {
    pub(crate) numerator: Atom,
    pub(crate) denominators: Vec<FourDDenominator>,
}

impl FourDSector {
    fn new(atom: Atom, frozen_lmbs: Vec<LoopMomentumBasis>) -> Self {
        Self { atom, frozen_lmbs }
    }

    pub(crate) fn atom(&self) -> &Atom {
        &self.atom
    }

    pub(crate) fn frozen_lmbs(&self) -> &[LoopMomentumBasis] {
        &self.frozen_lmbs
    }

    pub(crate) fn physical_terms(&self) -> Result<Vec<FourDTerm>> {
        let physical_atom = self
            .atom
            .replace(GS.integrated_loop_scale)
            .with(Atom::one());
        FourDTerm::from_view(physical_atom.as_view())
    }
}

impl FourDSectors {
    fn new(active: Vec<FourDSector>, recursive_completion: Vec<FourDSector>) -> Self {
        let atom = active
            .iter()
            .chain(&recursive_completion)
            .fold(Atom::Zero, |sum, sector| sum + &sector.atom);
        Self::with_atom(active, recursive_completion, atom)
    }

    fn with_atom(
        mut active: Vec<FourDSector>,
        mut recursive_completion: Vec<FourDSector>,
        atom: Atom,
    ) -> Self {
        active.retain(|sector| !sector.atom.is_zero());
        recursive_completion.retain(|sector| !sector.atom.is_zero());
        Self {
            active,
            recursive_completion,
            atom,
        }
    }

    fn active_atom(atom: Atom) -> Self {
        Self::new(vec![FourDSector::new(atom, Vec::new())], Vec::new())
    }

    fn all(&self) -> impl Iterator<Item = &FourDSector> {
        self.active.iter().chain(&self.recursive_completion)
    }
}

impl Full4dCts {
    /// Pole-part subtraction inserts only completed poles; MUV keeps the local
    /// counterterm together with its integrated finite contribution.
    pub(crate) fn recursion_input(
        local: &Local4dCts,
        integrated: &IntegratedCts,
        scheme: ApproximationType,
        is_root: bool,
        lmb: &LoopMomentumBasis,
    ) -> Result<Self> {
        let completed = |atom| FourDSector::new(atom, vec![lmb.clone()]);
        match scheme {
            ApproximationType::MUV => {
                let integrated = integrated.finite_counterterm_atom();
                Ok(Self(FourDSectors::with_atom(
                    local.0.active.clone(),
                    local
                        .0
                        .recursive_completion
                        .iter()
                        .cloned()
                        .chain([completed(integrated.clone())])
                        .collect(),
                    &local.0.atom + integrated,
                )))
            }
            ApproximationType::PolePart if !is_root => {
                let integrated = integrated.pole_atom();
                Ok(Self(FourDSectors::with_atom(
                    Vec::new(),
                    vec![completed(integrated.clone())],
                    integrated,
                )))
            }
            ApproximationType::PolePart => {
                let integrated = integrated.pole_atom();
                Ok(Self(FourDSectors::with_atom(
                    local.0.active.clone(),
                    local
                        .0
                        .recursive_completion
                        .iter()
                        .cloned()
                        .chain([completed(integrated.clone())])
                        .collect(),
                    &local.0.atom + integrated,
                )))
            }
            scheme => Err(eyre!("No recursive counterterm projection for {scheme}")),
        }
    }

    /// A disconnected local value already contains every component's recursive
    /// projection, so it must not be projected again when used by an outer limit.
    pub(crate) fn from_factorized_local(local: &Local4dCts) -> Self {
        Self(local.0.clone())
    }

    /// Complete every still-local typed sector with the propagators outside its
    /// owning spinney. Fully integrated recursive completions are deliberately
    /// excluded: final assembly adds that coefficient exactly once on its
    /// separately localized branch. Numerator factors outside the spinney are
    /// grown later under each exact source-local energy map.
    pub(crate) fn with_cograph<S: SubGraphLike>(
        local: &Local4dCts,
        graph: &Graph,
        cograph: &S,
    ) -> Self {
        let cograph = graph.denominator(cograph, |_| -1);
        let active = local
            .0
            .active
            .iter()
            .map(|sector| FourDSector::new(&sector.atom * &cograph, sector.frozen_lmbs.clone()))
            .collect();
        let recursive_completion = local
            .0
            .recursive_completion
            .iter()
            .fold(Atom::Zero, |sum, sector| sum + &sector.atom);
        Self(FourDSectors::with_atom(
            active,
            Vec::new(),
            (&local.0.atom - recursive_completion) * cograph,
        ))
    }

    #[cfg(test)]
    pub(crate) fn from_coefficient<S: SubGraphLike>(
        coefficient: &Atom,
        graph: &Graph,
        cograph: &S,
    ) -> Self {
        Self(FourDSectors::active_atom(
            coefficient * graph.denominator(cograph, |_| -1),
        ))
    }

    pub(crate) fn from_frozen_coefficient<S: SubGraphLike>(
        coefficient: &Atom,
        graph: &Graph,
        cograph: &S,
        lmb: &LoopMomentumBasis,
    ) -> Self {
        Self(FourDSectors::new(
            vec![FourDSector::new(
                coefficient * graph.denominator(cograph, |_| -1),
                vec![lmb.clone()],
            )],
            Vec::new(),
        ))
    }

    #[cfg(test)]
    pub(crate) fn terms(&self) -> Result<Vec<FourDTerm>> {
        FourDTerm::from_view(self.0.atom.as_view())
    }

    pub(crate) fn sectors(&self) -> impl Iterator<Item = &FourDSector> {
        self.0.all()
    }
}

impl Neg for Local4dCts {
    type Output = Local4dCts;
    fn neg(self) -> Self::Output {
        let atom = -self.0.atom;
        Local4dCts(FourDSectors::with_atom(
            self.0
                .active
                .into_iter()
                .map(|sector| FourDSector::new(-sector.atom, sector.frozen_lmbs))
                .collect(),
            self.0
                .recursive_completion
                .into_iter()
                .map(|sector| FourDSector::new(-sector.atom, sector.frozen_lmbs))
                .collect(),
            atom,
        ))
    }
}

impl Local4dCts {
    pub(crate) fn atom(&self) -> &Atom {
        &self.0.atom
    }

    pub(crate) fn from_full_product(factors: impl IntoIterator<Item = Full4dCts>) -> Self {
        let mut products = vec![(FourDSector::new(Atom::one(), Vec::new()), false)];
        let mut atom = Atom::one();
        for factor in factors {
            atom *= &factor.0.atom;
            let factor_sectors = factor
                .0
                .active
                .into_iter()
                .map(|sector| (sector, true))
                .chain(
                    factor
                        .0
                        .recursive_completion
                        .into_iter()
                        .map(|sector| (sector, false)),
                )
                .collect::<Vec<_>>();
            products = products
                .into_iter()
                .flat_map(|(left, left_active)| {
                    factor_sectors
                        .iter()
                        .cloned()
                        .map(move |(right, right_active)| {
                            let mut frozen_lmbs = left.frozen_lmbs.clone();
                            frozen_lmbs.extend(right.frozen_lmbs);
                            (
                                FourDSector::new(left.atom.clone() * right.atom, frozen_lmbs),
                                left_active || right_active,
                            )
                        })
                })
                .collect();
        }
        let (active, recursive_completion) = products
            .into_iter()
            .partition::<Vec<_>, _>(|(_, is_active)| *is_active);
        Self(FourDSectors::with_atom(
            active.into_iter().map(|(sector, _)| sector).collect(),
            recursive_completion
                .into_iter()
                .map(|(sector, _)| sector)
                .collect(),
            atom,
        ))
    }

    /// Enumerate only the additive terms needed by the 3D projection while
    /// retaining every non-denominator factor as a Symbolica atom. In
    /// particular, numerator products and powers are not materialized into a
    /// parallel expression tree or expanded polynomial.
    #[cfg(test)]
    pub(crate) fn terms(&self) -> Result<Vec<FourDTerm>> {
        FourDTerm::from_view(self.0.atom.as_view())
    }

    #[cfg(test)]
    pub(crate) fn active_sectors(&self) -> &[FourDSector] {
        &self.0.active
    }

    #[cfg(test)]
    pub(crate) fn recursive_completion(&self) -> &[FourDSector] {
        &self.0.recursive_completion
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
        Self::from_view_with_provenance(view)?
            .into_iter()
            .map(Self::cancel_positive_denominator_factors)
            .collect()
    }

    fn from_view_with_provenance(view: AtomView<'_>) -> Result<Vec<Self>> {
        match view {
            AtomView::Add(add) => {
                // A collected Taylor coefficient can contain provenance
                // wrappers in its additive numerator which cancel different
                // powers of the surrounding common denominator. Distribute
                // only such sums; ordinary numerator sums remain factorized.
                if !view.to_owned().contains_symbol(GS.den) {
                    return Ok(vec![Self::numerator(view.to_owned())]);
                }
                add.iter()
                    .map(Self::from_view_with_provenance)
                    .collect::<Result<Vec<_>>>()
                    .map(|terms| terms.into_iter().flatten().collect())
            }
            AtomView::Mul(mul) => {
                let mut terms = vec![Self::numerator(Atom::one())];
                for factor in mul.iter() {
                    let factor_terms = Self::from_view_with_provenance(factor)?;
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

    fn cancel_positive_denominator_factors(self) -> Result<Self> {
        let Self {
            numerator,
            mut denominators,
        } = self;
        let factors = match numerator.as_view() {
            AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
            factor => vec![factor],
        };
        let mut retained_factors = Vec::new();
        for factor in factors {
            let positive_denominator = match factor {
                AtomView::Fun(_) => {
                    FourDDenominator::from_view(factor)?.map(|denominator| (denominator, 1))
                }
                AtomView::Pow(power) => {
                    let (base, exponent) = power.get_base_exp();
                    let Some(denominator) = FourDDenominator::from_view(base)? else {
                        retained_factors.push(factor.to_owned());
                        continue;
                    };
                    let Ok(exponent) = i64::try_from(exponent) else {
                        return Err(eyre!(
                            "4D denominator has non-integer numerator power `{}`",
                            exponent.to_owned()
                        ));
                    };
                    if exponent > 0 {
                        Some((
                            denominator,
                            usize::try_from(exponent).map_err(|_| {
                                eyre!("4D denominator numerator power does not fit in memory")
                            })?,
                        ))
                    } else {
                        None
                    }
                }
                _ => None,
            };
            let Some((positive_denominator, multiplicity)) = positive_denominator else {
                retained_factors.push(factor.to_owned());
                continue;
            };
            for _ in 0..multiplicity {
                if let Some(position) = denominators
                    .iter()
                    .position(|denominator| denominator == &positive_denominator)
                {
                    denominators.remove(position);
                } else {
                    retained_factors.push(positive_denominator.full_expr.clone());
                }
            }
        }
        let numerator = retained_factors
            .into_iter()
            .fold(Atom::one(), |product, factor| product * factor)
            .replace(GS.den(W_.a_, W_.b_, W_.c_, W_.d_))
            .with(W_.d_);
        Ok(Self {
            numerator,
            denominators,
        })
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
        Self(FourDSectors::active_atom(Atom::one()))
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
    integrand: &Atom,
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
    let integrand = (integrand * t_arg).simplify_metrics();

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
) -> Result<Atom> {
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

    Ok(dotted)
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
            let marker = UvMarker::new(ctx.settings);
            let project = |atom: &Atom| -> Result<Atom> {
                if atom.is_zero() {
                    return Ok(Atom::Zero);
                }
                let grown = grow(atom, ctx, current, given)?;
                let result = -t(&grown, ctx, current, given)?;
                Ok(marker.apply(
                    UvOperation::Approx,
                    marker_current.subgraph(),
                    marker_given.subgraph(),
                    &result,
                ))
            };
            let sectors = integrand
                .sectors()
                .map(|sector| {
                    Ok(FourDSector::new(
                        project(sector.atom())?,
                        sector.frozen_lmbs().to_vec(),
                    ))
                })
                .collect::<Result<Vec<_>>>()?;
            let atom = if sectors.len() == 1 {
                sectors[0].atom.clone()
            } else {
                project(&integrand.0.atom)?
            };
            Ok(Local4dCts(FourDSectors::with_atom(
                sectors,
                Vec::new(),
                atom,
            )))
        }
        atype => Err(eyre!("Not yet implemented {:?}", atype)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        cff::generation::ExactCffGenerationCache,
        dot,
        graph::{Graph, GraphThreeDSource, parse::IntoGraph},
        initialisation::test_initialise,
        uv::{Spinney, UVgenerationSettings, hedge_poset::OwnedForestNode},
    };
    use linnet::half_edge::involution::EdgeIndex;
    use linnet::half_edge::subgraph::{InternalSubGraph, SubSetOps};
    use symbolica::{function, symbol};

    #[test]
    fn term_projection_keeps_factorized_numerator_atoms() -> Result<()> {
        test_initialise()?;
        let affine = symbol!("local_4d_test::Affine");
        let k0 = symbol!("local_4d_test::K0");
        let k1 = symbol!("local_4d_test::K1");
        let sum = function!(affine, Atom::var(k0)) + function!(affine, Atom::var(k1));
        let numerator = sum.pow(2) * (Atom::var(k0) + Atom::var(k1));

        let terms = Local4dCts(FourDSectors::active_atom(numerator.clone())).terms()?;
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
        let source = Full4dCts(FourDSectors::active_atom(first.pow(-2) * second.pow(-1)));

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

    #[test]
    fn term_projection_cancels_provenance_without_expanding_ordinary_sums() -> Result<()> {
        test_initialise()?;
        let first_full = Atom::var(symbol!("local_4d_test::first_cancellable_full"));
        let second_full = Atom::var(symbol!("local_4d_test::second_cancellable_full"));
        let first = GS.den(
            0,
            FunctionBuilder::new(GS.emr_mom).add_arg(0).finish(),
            0,
            &first_full,
        );
        let second = GS.den(
            1,
            FunctionBuilder::new(GS.emr_mom).add_arg(1).finish(),
            0,
            &second_full,
        );
        let left = Atom::var(symbol!("local_4d_test::left"));
        let right = Atom::var(symbol!("local_4d_test::right"));
        let ordinary_sum = left + right;
        let collected = (&ordinary_sum * &first + &second) * first.pow(-1) * second.pow(-1);

        let mut terms = FourDTerm::from_view(collected.as_view())?;
        terms.sort_by_key(|term| usize::from(term.denominators[0].source_edge));
        assert_eq!(terms.len(), 2);
        assert_eq!(
            terms[0].denominators,
            vec![FourDDenominator::from_view(first.as_view())?.unwrap()]
        );
        assert_eq!(terms[0].numerator, Atom::one());
        assert_eq!(
            terms[1].denominators,
            vec![FourDDenominator::from_view(second.as_view())?.unwrap()]
        );
        assert_eq!(terms[1].numerator, ordinary_sum);
        assert!(
            terms
                .iter()
                .all(|term| !term.numerator.contains_symbol(GS.den))
        );
        Ok(())
    }

    #[test]
    fn dod_one_triangle_taylor_terms_keep_source_owners_and_reuse_two_cffs() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph exact_uv_triangle_taylor {
            edge [num=1 mass=1]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]

            incoming -> v5 [id=0]
            v0 -> v2 [id=1 lmb_id=0]
            v3 -> v0 [id=2]
            v0 -> v5 [id=3]
            v2 -> v1 [id=4]
            v1 -> v3 [id=5 lmb_id=1]
            v1 -> v4 [id=6]
            v2 -> v3 [id=7]
            v4 -> outgoing [id=8]
        })?;
        let owners = [EdgeIndex(4), EdgeIndex(5), EdgeIndex(7)];
        let uv_filter = owners
            .into_iter()
            .map(|edge| graph.get_edge_subgraph(edge))
            .reduce(|left, right| left.union(&right))
            .expect("the UV triangle has three source edges");
        let uv_subgraph =
            InternalSubGraph::cleaned_filter_optimist(uv_filter.clone(), graph.as_ref());
        let current = OwnedForestNode {
            spinney: Spinney::with_scheme(
                uv_subgraph,
                &graph,
                &graph.loop_momentum_basis,
                ApproximationType::MUV,
                1,
            )
            .expect("the source triangle has a compatible loop-momentum basis"),
            topo_order: 0,
        };
        let given = OwnedForestNode {
            spinney: Spinney::empty(&graph),
            topo_order: 0,
        };
        let numerator = owners.into_iter().fold(Atom::one(), |product, edge| {
            product * GS.emr_mom(edge, GS.cind(0))
        });
        let integrand = numerator / graph.denominator(&uv_filter, |_| 1);
        let settings = UVgenerationSettings::default();
        let expanded = t(&integrand, &UVCtx::new(&graph, &settings), &current, &given)?;
        let cograph = graph
            .get_edge_subgraph(EdgeIndex(1))
            .union(&graph.get_edge_subgraph(EdgeIndex(2)));
        let terms = Full4dCts::from_coefficient(&expanded, &graph, &cograph).terms()?;

        let mut owner_multiplicities = terms
            .iter()
            .map(|term| {
                owners.map(|owner| {
                    term.denominators
                        .iter()
                        .filter(|denominator| denominator.source_edge == owner)
                        .count()
                })
            })
            .collect::<Vec<_>>();
        owner_multiplicities.sort_unstable();
        assert_eq!(
            owner_multiplicities,
            vec![[1, 1, 1], [1, 1, 1], [1, 1, 2], [2, 1, 1]],
            "the real Taylor coefficient must retain one base topology and one dotted topology for each shifted source owner",
        );

        let options = graph.denominator_only_cff_3d_expression_options();
        let mut cache = ExactCffGenerationCache::default();
        let mut exact_bounds = Vec::new();
        let active_terms = terms
            .iter()
            .map(|term| {
                term.denominators
                    .iter()
                    .map(|denominator| {
                        let is_uv = uv_filter.includes(&graph[&denominator.source_edge].1);
                        Ok(denominator
                            .depends_on_loop(&graph, is_uv)?
                            .then(|| denominator.clone()))
                    })
                    .collect::<Result<Vec<_>>>()
                    .map(|denominators| {
                        denominators
                            .into_iter()
                            .flatten()
                            .collect::<Vec<FourDDenominator>>()
                    })
            })
            .collect::<Result<Vec<_>>>()?;
        for (term, active_denominators) in terms.iter().zip(&active_terms) {
            let source = GraphThreeDSource::from_exact_denominators_in_uv_edges(
                &graph,
                active_denominators,
                owners,
            )?;
            graph.register_3d_expression_for_4d_term(
                &source,
                &options,
                &term.numerator,
                &mut cache,
            )?;
        }
        for (term, active_denominators) in terms.iter().zip(&active_terms) {
            let source = GraphThreeDSource::from_exact_denominators_in_uv_edges(
                &graph,
                active_denominators,
                owners,
            )?;
            let (_, _, plan, _) = graph.generate_3d_expression_for_4d_term(
                &source,
                &options,
                &term.numerator,
                Some(&mut cache),
            )?;
            exact_bounds.push(plan.energy_degree_bounds().to_vec());
        }
        assert_eq!(
            cache.len(),
            2,
            "the actual numerator-hit and denominator-hit Taylor terms must generate only the base and dotted canonical CFF topologies; exact bounds: {exact_bounds:?}",
        );
        Ok(())
    }

    #[test]
    fn factorized_product_separates_active_and_completed_sectors() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(
            digraph G {
                edge [particle="scalar_1"];
                node [num=1];
                a -> a [id=0];
                b -> b [id=1];
            },
            "scalars"
        )?;
        let lmb_a = graph
            .generate_loop_momentum_bases_of(&graph.get_edge_subgraph(EdgeIndex(0)))
            .into_iter()
            .next()
            .expect("the first self-loop has a loop-momentum basis");
        let lmb_b = graph
            .generate_loop_momentum_bases_of(&graph.get_edge_subgraph(EdgeIndex(1)))
            .into_iter()
            .next()
            .expect("the second self-loop has a loop-momentum basis");
        let local_a = Atom::var(symbol!("local_4d_test::local_a"));
        let finite_a = Atom::var(symbol!("local_4d_test::finite_a"));
        let local_b = Atom::var(symbol!("local_4d_test::local_b"));
        let finite_b = Atom::var(symbol!("local_4d_test::finite_b"));
        let component = |local: Atom, finite: Atom, lmb: LoopMomentumBasis| {
            Full4dCts(FourDSectors::new(
                vec![FourDSector::new(local, Vec::new())],
                vec![FourDSector::new(finite, vec![lmb])],
            ))
        };

        let product = Local4dCts::from_full_product([
            component(local_a.clone(), finite_a.clone(), lmb_a.clone()),
            component(local_b.clone(), finite_b.clone(), lmb_b.clone()),
        ]);
        assert_eq!(product.active_sectors().len(), 3);
        assert_eq!(product.recursive_completion().len(), 1);

        let sector = |atom: Atom| {
            product
                .active_sectors()
                .iter()
                .find(|sector| sector.atom == atom)
                .expect("the factorized active sector is present")
        };
        assert!(sector(&local_a * &local_b).frozen_lmbs.is_empty());
        assert_eq!(
            sector(&finite_a * &local_b).frozen_lmbs,
            vec![lmb_a.clone()]
        );
        assert_eq!(
            sector(&local_a * &finite_b).frozen_lmbs,
            vec![lmb_b.clone()]
        );
        assert_eq!(
            product.recursive_completion()[0],
            FourDSector::new(&finite_a * &finite_b, vec![lmb_a, lmb_b])
        );

        let cograph = graph.empty_subgraph::<SuBitGraph>();
        let projected_source = Full4dCts::with_cograph(&product, &graph, &cograph);
        assert_eq!(projected_source.sectors().count(), 3);
        let recursive_source = Full4dCts::from_factorized_local(&product);
        assert_eq!(recursive_source.sectors().count(), 4);
        assert_eq!(
            product.atom(),
            &((&local_a + &finite_a) * (&local_b + &finite_b))
        );
        Ok(())
    }
}
