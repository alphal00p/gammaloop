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

use linnet::half_edge::{
    involution::EdgeIndex,
    subgraph::{Inclusion, InternalSubGraph, SuBitGraph, SubGraphLike, SubSetLike},
};
use spenso::shadowing::TensorCollectExt;
use symbolica::{
    atom::{AtomCore, AtomView},
    prelude::*,
};

#[cfg(test)]
use crate::utils::symbols::UvMomentumProvenanceRole;
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
    // Each active component retains its denominator-owner set, the full source
    // scope whose omitted prefix is contracted by exact CFF reconstruction, and
    // the quotient LMB in which its energy residues must be generated. The
    // complete Taylor LMB records the compatible coordinates used by the latest
    // enclosing T.
    // Both are distinct from frozen LMBs whose integrations have already been
    // completed and which only own localization factors.
    pub(crate) active_components: Vec<(SuBitGraph, SuBitGraph, LoopMomentumBasis)>,
    pub(crate) taylor_lmb: Option<LoopMomentumBasis>,
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
    pub(crate) fn new(
        atom: Atom,
        active_components: Vec<(SuBitGraph, SuBitGraph, LoopMomentumBasis)>,
        taylor_lmb: Option<LoopMomentumBasis>,
        frozen_lmbs: Vec<LoopMomentumBasis>,
    ) -> Self {
        Self {
            atom,
            active_components,
            taylor_lmb,
            frozen_lmbs,
        }
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
        Self::new(
            vec![FourDSector::new(atom, Vec::new(), None, Vec::new())],
            Vec::new(),
        )
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
        let completed = |atom| FourDSector::new(atom, Vec::new(), None, vec![lmb.clone()]);
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
            .map(|sector| {
                FourDSector::new(
                    &sector.atom * &cograph,
                    sector.active_components.clone(),
                    sector.taylor_lmb.clone(),
                    sector.frozen_lmbs.clone(),
                )
            })
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

    #[cfg(test)]
    pub(crate) fn from_frozen_coefficient<S: SubGraphLike>(
        coefficient: &Atom,
        graph: &Graph,
        cograph: &S,
        lmb: &LoopMomentumBasis,
    ) -> Self {
        Self(FourDSectors::new(
            vec![FourDSector::new(
                coefficient * graph.denominator(cograph, |_| -1),
                Vec::new(),
                None,
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
                .map(|sector| {
                    FourDSector::new(
                        -sector.atom,
                        sector.active_components,
                        sector.taylor_lmb,
                        sector.frozen_lmbs,
                    )
                })
                .collect(),
            self.0
                .recursive_completion
                .into_iter()
                .map(|sector| {
                    FourDSector::new(
                        -sector.atom,
                        sector.active_components,
                        sector.taylor_lmb,
                        sector.frozen_lmbs,
                    )
                })
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
        let mut products = vec![(
            FourDSector::new(Atom::one(), Vec::new(), None, Vec::new()),
            false,
        )];
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
                            let left_has_active_components = !left.active_components.is_empty();
                            let right_has_active_components = !right.active_components.is_empty();
                            let mut active_components = left.active_components.clone();
                            active_components.extend(right.active_components);
                            let taylor_lmb =
                                match (left_has_active_components, right_has_active_components) {
                                    (true, false) => left.taylor_lmb.clone(),
                                    (false, true) => right.taylor_lmb.clone(),
                                    // Two independent active factors have no single
                                    // enclosing Taylor frame. Once ambiguous, later
                                    // products must not spuriously restore one.
                                    _ => None,
                                };
                            let mut frozen_lmbs = left.frozen_lmbs.clone();
                            frozen_lmbs.extend(right.frozen_lmbs);
                            (
                                FourDSector::new(
                                    left.atom.clone() * right.atom,
                                    active_components,
                                    taylor_lmb,
                                    frozen_lmbs,
                                ),
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
        match view {
            AtomView::Add(add) => {
                // The outer additive shell of a completed Taylor coefficient can
                // contain terms with different denominator topologies. Split only
                // those distinct topologies. A common topology retains its additive
                // numerator as one factorized atom, even when it contains positive
                // typed denominator factors for CFF to pinch internally.
                let terms = add
                    .iter()
                    .map(Self::from_view)
                    .collect::<Result<Vec<_>>>()?
                    .into_iter()
                    .flatten()
                    .collect::<Vec<_>>();
                let common_denominators = terms[0].denominators.clone();
                if terms
                    .iter()
                    .all(|term| term.denominators == common_denominators)
                {
                    Ok(vec![Self {
                        numerator: terms
                            .into_iter()
                            .fold(Atom::Zero, |sum, term| sum + term.numerator),
                        denominators: common_denominators,
                    }])
                } else {
                    Ok(terms)
                }
            }
            AtomView::Mul(mul) => {
                let mut products = vec![Self::numerator(Atom::one())];
                for factor in mul.iter() {
                    let factor_terms = Self::from_view(factor)?;
                    products = products
                        .into_iter()
                        .flat_map(|left| {
                            factor_terms
                                .iter()
                                .cloned()
                                .map(move |right| Self::product(left.clone(), right))
                        })
                        .collect();
                }
                Ok(products)
            }
            _ => Ok(vec![Self::from_factorized_term(view)?]),
        }
    }

    fn from_factorized_term(view: AtomView<'_>) -> Result<Self> {
        match view {
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                let Some(denominator) = FourDDenominator::from_view(base)? else {
                    return Ok(Self::numerator(view.to_owned()));
                };
                let Ok(exponent) = i64::try_from(exponent) else {
                    return Err(eyre!(
                        "4D denominator has non-integer power `{}`",
                        exponent.to_owned()
                    ));
                };
                if exponent >= 0 {
                    return Ok(Self::numerator(view.to_owned()));
                }
                let multiplicity = usize::try_from(exponent.unsigned_abs())
                    .map_err(|_| eyre!("4D denominator multiplicity does not fit in memory"))?;
                Ok(Self {
                    numerator: Atom::one(),
                    denominators: std::iter::repeat_n(denominator, multiplicity).collect(),
                })
            }
            _ => Ok(Self::numerator(view.to_owned())),
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
            momentum: GS.erase_uv_momentum_provenance(&function.get(1).to_owned()),
            mass_squared: function.get(2).to_owned(),
            full_expr: GS.erase_uv_momentum_provenance(&function.get(3).to_owned()),
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
        expansion_subgraph: &SuBitGraph,
        n_loops: usize,
        lmb: &LoopMomentumBasis,
        atom: &Atom,
    ) -> Atom {
        let scaled_edges = self
            .underlying
            .iter_edges_of(expansion_subgraph)
            .filter_map(|(pair, edge, _)| {
                (matches!(
                    pair,
                    linnet::half_edge::involution::HedgePair::Paired { .. }
                ) && lmb.edge_signatures[edge]
                    .internal
                    .iter()
                    .any(|coefficient| *coefficient != crate::momentum::SignOrZero::Zero))
                .then_some(edge)
            })
            .collect::<Vec<_>>();
        // Keep the denominator owner in the momentum tag before the Taylor
        // operator acts. The derivative of `den` traverses only its fourth
        // argument, so a Q carrying this tag can later leave the wrapper while
        // still identifying the source line which produced it.
        let tagged_edges = scaled_edges.clone();
        let tagged = atom
            .replace(GS.den(W_.a_, W_.mom_, W_.mass_, W_.prop_))
            .with_map(move |matched| {
                let edge = matched.get(W_.a_).unwrap().to_atom();
                let is_scaled = usize::try_from(edge.as_view()).is_ok_and(|edge| {
                    tagged_edges.contains(&linnet::half_edge::involution::EdgeIndex(edge))
                });
                if !is_scaled {
                    return GS.den(
                        edge,
                        matched.get(W_.mom_).unwrap().to_atom(),
                        matched.get(W_.mass_).unwrap().to_atom(),
                        matched.get(W_.prop_).unwrap().to_atom(),
                    );
                }
                let source_momentum = function!(GS.emr_mom, edge.as_view());
                // Role two is transient and distinguishes a denominator which
                // belongs to this Taylor operation from a persistent role-zero
                // or role-one tag left by a nested child operation.
                let tag = GS.uv_momentum_provenance.call_args([
                    edge.as_view(),
                    Atom::num(2).as_view(),
                    source_momentum.as_view(),
                ]);
                let tag_momenta = |value: Atom| {
                    value.replace_map(|view, _, output| {
                        let AtomView::Fun(momentum) = view else {
                            return;
                        };
                        if momentum.get_symbol() != GS.emr_mom || momentum.get_nargs() == 0 {
                            return;
                        }
                        if matches!(
                            momentum.get(0),
                            AtomView::Fun(provenance)
                                if provenance.get_symbol() == GS.uv_momentum_provenance
                        ) {
                            // Existing child provenance is immutable metadata.
                            // Re-emitting the complete momentum also makes this
                            // node opaque to the top-down traversal, so the outer
                            // operation cannot tag a literal Q inside its stored
                            // hard-momentum payload.
                            **output = Atom::from(momentum.to_owned());
                            return;
                        }
                        if momentum.get(0) != edge.as_view() {
                            return;
                        }
                        let mut tagged = FunctionBuilder::new(GS.emr_mom).add_arg(tag.clone());
                        for index in momentum.iter().skip(1) {
                            tagged = tagged.add_arg(index.to_owned());
                        }
                        **output = tagged.finish();
                    })
                };
                GS.den(
                    edge.clone(),
                    tag_momenta(matched.get(W_.mom_).unwrap().to_atom()),
                    matched.get(W_.mass_).unwrap().to_atom(),
                    tag_momenta(matched.get(W_.prop_).unwrap().to_atom()),
                )
            });

        // Scale the hard part of every momentum in the complete current UV
        // subgraph while retaining both its immutable source owner and its
        // complete child-LMB hard projection in the Q tag. Thus
        // Q_e(t) = H_e/t + S_e. Keeping H_e as one tagged rank-one momentum is
        // essential: a derivative of D_e must not lose its owner when the LMB
        // happens to spell H_e with another physical edge ID.
        let inverse_rescale = Atom::one() / GS.rescale;
        let mut derivative_replacements = Vec::new();
        let mut fixed_replacements = Vec::new();
        for edge in scaled_edges {
            let edge_atom = Atom::num(usize::from(edge) as i64);
            let soft = lmb.ext_atom(edge, GS.emr_mom, &[W_.x___], true);
            let hard = lmb.loop_atom::<Atom>(edge, GS.emr_mom, &[], true);
            let ordinary = function!(GS.emr_mom, edge_atom.as_view(), W_.x___);
            let fixed_tag =
                GS.uv_momentum_provenance_tag(edge_atom.as_view(), false, hard.as_view());
            let fixed_hard = function!(GS.emr_mom, fixed_tag.as_view(), W_.x___);
            fixed_replacements.push(Replacement::new(
                ordinary.to_pattern(),
                (&fixed_hard * &inverse_rescale + &soft).to_pattern(),
            ));

            let preliminary_tag = GS.uv_momentum_provenance.call_args([
                edge_atom.as_view(),
                Atom::num(2).as_view(),
                function!(GS.emr_mom, edge_atom.as_view()).as_view(),
            ]);
            let derivative = function!(GS.emr_mom, preliminary_tag.as_view(), W_.x___);
            let derivative_tag =
                GS.uv_momentum_provenance_tag(edge_atom.as_view(), true, hard.as_view());
            let derivative_hard = function!(GS.emr_mom, derivative_tag.as_view(), W_.x___);
            derivative_replacements.push(Replacement::new(
                derivative.to_pattern(),
                (&derivative_hard * &inverse_rescale + &soft).to_pattern(),
            ));
        }
        // A nested child tag already stores its complete frozen hard
        // projection. The compatible outer LMB scales every carrier of that
        // projection, so the whole tagged momentum is homogeneous and acquires
        // one scalar 1/t without any soft shift or metadata rewrite.
        let persistent = function!(
            GS.emr_mom,
            function!(GS.uv_momentum_provenance, W_.a_, W_.b_, W_.mom_),
            W_.x___
        );
        derivative_replacements.push(Replacement::new(
            persistent.to_pattern(),
            (&persistent * &inverse_rescale).to_pattern(),
        ));
        derivative_replacements.extend(fixed_replacements);
        let mut atomarg = tagged.replace_multiple(&derivative_replacements);

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
    retained_loop_edges: &[EdgeIndex],
    required_prefix_loop_count: usize,
) -> Result<(Atom, LoopMomentumBasis)> {
    let graph = ctx.graph;
    let n_loops = graph.n_loops(current.subgraph());

    // Keep the child loops as a subset of the outer basis so integrating the child
    // commutes with the outer Taylor expansion.
    let is_compatible = |candidate: &LoopMomentumBasis| {
        (if retained_loop_edges.is_empty() {
            given
                .lmb()
                .loop_edges
                .iter()
                .all(|edge| candidate.loop_edges.contains(edge))
        } else {
            retained_loop_edges
                .iter()
                .all(|edge| candidate.loop_edges.contains(edge))
        }) && candidate
            .loop_edges
            .iter()
            .filter(|edge| given.subgraph().includes(&graph[*edge].1))
            .count()
            == required_prefix_loop_count
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

    let rescaled = graph.uv_rescaled(current.subgraph(), n_loops, lmb, integrand);
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

    Ok((dotted, lmb.clone()))
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
            let reduced_subgraph = current.reduced_subgraph(given);
            let crown = ctx
                .graph
                .dummy_stripped_external_flows_of(current.subgraph());
            let expected_component_loops =
                ctx.graph.n_loops(current.subgraph()) - ctx.graph.n_loops(given.subgraph());
            let project = |sector: &FourDSector| -> Result<FourDSector> {
                if sector.atom.is_zero() {
                    return Ok(FourDSector::new(
                        Atom::Zero,
                        Vec::new(),
                        None,
                        sector.frozen_lmbs.clone(),
                    ));
                }
                let grown = grow(&sector.atom, ctx, current, given)?;
                let mut retained_loop_edges = sector
                    .active_components
                    .iter()
                    .flat_map(|(_, _, lmb)| lmb.loop_edges.iter().copied())
                    .chain(
                        sector
                            .frozen_lmbs
                            .iter()
                            .flat_map(|lmb| lmb.loop_edges.iter().copied()),
                    )
                    .collect::<Vec<_>>();
                retained_loop_edges.sort();
                retained_loop_edges.dedup();
                let expected_prefix_loops = ctx.graph.n_loops(given.subgraph());
                // Fix the quotient carrier before applying T.  Otherwise an
                // enclosing basis may choose an equivalent physical edge whose
                // full-graph momentum differs by a crown shift; after Taylor
                // projection that edge denotes the hard dummy coordinate, not
                // its shifted production momentum.  This canonical component
                // LMB keeps the completed local-4D Taylor coefficient in one
                // deterministic child frame through exact-source reconstruction
                // and later CFF projection.  It is projected-local4D routing
                // data; direct local3D instead applies T to the complete CFF
                // expression after loop-energy integration.
                let canonical_component_lmb = if given.subgraph().is_empty() {
                    None
                } else {
                    let given_internal =
                        InternalSubGraph::try_new(given.subgraph().clone(), ctx.graph.as_ref())
                            .ok_or_else(|| {
                                eyre!("nested UV prefix is not a valid internal subgraph")
                            })?;
                    Some(ctx.graph.shrunken_sub_lmb(
                        current.subgraph(),
                        &given_internal,
                        crown.clone(),
                        None,
                    )?)
                };
                if let Some(component_lmb) = &canonical_component_lmb {
                    retained_loop_edges.extend(component_lmb.loop_edges.iter().copied());
                    retained_loop_edges.sort();
                    retained_loop_edges.dedup();
                }
                let expected_retained_loops = expected_prefix_loops
                    + canonical_component_lmb
                        .as_ref()
                        .map_or(0, |lmb| lmb.loop_edges.len());
                if retained_loop_edges.len() != expected_retained_loops {
                    return Err(eyre!(
                        "local 4D Taylor sector retains {} prefix-plus-quotient carriers, expected {expected_retained_loops}",
                        retained_loop_edges.len()
                    ));
                }
                let (result, coordinate_lmb) = t(
                    &grown,
                    ctx,
                    current,
                    given,
                    &retained_loop_edges,
                    expected_prefix_loops,
                )?;
                // Exact residues must use the very coordinates in which T
                // produced their hard momenta. Rebuilding a canonical quotient
                // LMB here can choose an equivalent but differently spelled
                // carrier and thereby reintroduce a physical crown shift. Use
                // graphic-matroid contraction with the enclosing Taylor LMB as
                // its compatibility guide, so prefix directions disappear
                // instead of becoming artificial external coordinates.
                let component_lmb = if given.subgraph().is_empty() {
                    coordinate_lmb.clone()
                } else {
                    let given_internal =
                        InternalSubGraph::try_new(given.subgraph().clone(), ctx.graph.as_ref())
                            .ok_or_else(|| {
                                eyre!("nested UV prefix is not a valid internal subgraph")
                            })?;
                    ctx.graph.shrunken_sub_lmb(
                        current.subgraph(),
                        &given_internal,
                        crown.clone(),
                        Some(&coordinate_lmb),
                    )?
                };
                if component_lmb.loop_edges.len() != expected_component_loops {
                    return Err(eyre!(
                        "local 4D Taylor component has {} loop generators after demoting its nested prefix, expected L(current)-L(given) = {expected_component_loops}",
                        component_lmb.loop_edges.len()
                    ));
                }
                if let Some(canonical_component_lmb) = &canonical_component_lmb
                    && component_lmb.loop_edges != canonical_component_lmb.loop_edges
                {
                    return Err(eyre!(
                        "contracting the local 4D Taylor basis changed the preselected quotient carriers from {:?} to {:?}",
                        canonical_component_lmb.loop_edges,
                        component_lmb.loop_edges,
                    ));
                }
                let mut active_components = sector.active_components.clone();
                active_components.push((
                    reduced_subgraph.clone(),
                    current.subgraph().clone(),
                    component_lmb.clone(),
                ));
                debug_tags!(#generation, #uv, #local, #four_d, #trace;
                    current = %current.log_display(),
                    given = %given.log_display(),
                    retained_loop_edges = ?retained_loop_edges,
                    component_lmb = ?component_lmb,
                    coordinate_lmb = ?coordinate_lmb,
                    active_components = ?active_components,
                    "Framed local-4D Taylor sector"
                );
                Ok(FourDSector::new(
                    marker.apply(
                        UvOperation::Approx,
                        marker_current.subgraph(),
                        marker_given.subgraph(),
                        &-result,
                    ),
                    active_components,
                    Some(coordinate_lmb),
                    sector.frozen_lmbs.clone(),
                ))
            };
            let sectors = integrand
                .sectors()
                .map(project)
                .collect::<Result<Vec<_>>>()?;
            // T is linear. Summing the individually framed sectors preserves the
            // aggregate compatibility atom without inventing one coordinate LMB
            // for a disconnected product.
            Ok(Local4dCts(FourDSectors::new(sectors, Vec::new())))
        }
        atype => Err(eyre!("Not yet implemented {:?}", atype)),
    }
}

#[cfg(test)]
mod tests {
    use std::collections::{BTreeMap, BTreeSet};

    use super::*;
    use crate::{
        cff::generation::ExactCffGenerationCache,
        dot,
        graph::{Graph, GraphThreeDSource, parse::IntoGraph},
        initialisation::test_initialise,
        numerator::energy_degree::EnergyPowerAnalyzer,
        uv::{Spinney, UVgenerationSettings, hedge_poset::OwnedForestNode},
    };
    use linnet::half_edge::involution::EdgeIndex;
    use linnet::half_edge::subgraph::{InternalSubGraph, SubSetOps};
    use spenso::structure::representation::{LibraryRep, Minkowski, RepName};
    use symbolica::{domains::rational::Rational, function, symbol};

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
    fn term_projection_exposes_denominators_below_an_outer_product() -> Result<()> {
        test_initialise()?;
        let coefficient = Atom::var(symbol!("local_4d_test::coefficient"));
        let first_numerator = Atom::var(symbol!("local_4d_test::first_numerator"));
        let second_numerator = Atom::var(symbol!("local_4d_test::second_numerator"));
        let first_full = Atom::var(symbol!("local_4d_test::first_full"));
        let second_full = Atom::var(symbol!("local_4d_test::second_full"));
        let first = GS.den(
            1,
            FunctionBuilder::new(GS.emr_mom).add_arg(1).finish(),
            GS.m_uv_expansion.pow(2),
            &first_full,
        );
        let second = GS.den(
            2,
            -FunctionBuilder::new(GS.emr_mom).add_arg(1).finish(),
            GS.m_uv_expansion.pow(2),
            &second_full,
        );
        let expression = &coefficient
            * (&first_numerator / (&first * &second)
                + &second_numerator / (&first * second.pow(2)));

        let mut terms = FourDTerm::from_view(expression.as_view())?;
        terms.sort_by_key(|term| term.denominators.len());
        assert_eq!(terms.len(), 2);
        assert_eq!(
            terms
                .iter()
                .map(|term| term.denominators.len())
                .collect::<Vec<_>>(),
            vec![2, 3],
        );
        assert_eq!(terms[0].numerator, &coefficient * &first_numerator);
        assert_eq!(terms[1].numerator, &coefficient * &second_numerator);
        assert!(
            terms
                .iter()
                .all(|term| !matches!(term.numerator.as_view(), AtomView::Add(_))),
            "enumerating denominator topologies must not expand either factorized numerator",
        );

        let reconstructed = terms.into_iter().fold(Atom::Zero, |sum, term| {
            let denominators =
                term.denominators
                    .into_iter()
                    .fold(Atom::one(), |product, denominator| {
                        product
                            * GS.den(
                                usize::from(denominator.source_edge),
                                denominator.momentum,
                                denominator.mass_squared,
                                denominator.full_expr,
                            )
                    });
            sum + term.numerator / denominators
        });
        assert!((expression - reconstructed).together().is_zero());
        Ok(())
    }

    #[test]
    fn nested_uv_rescaling_keeps_child_momentum_provenance_immutable() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph nested_provenance_rescaling {
            edge [num=1 mass=1]
            node [num=1]

            a -> b [id=0 lmb_id=0]
            b -> c [id=1]
            c -> a [id=2]
        })?;
        let filter = (0..3)
            .map(|edge| graph.get_edge_subgraph(EdgeIndex(edge)))
            .reduce(|left, right| left.union(&right))
            .expect("the vacuum triangle has three edges");
        let owner = EdgeIndex(1);
        let carrier = graph
            .loop_momentum_basis
            .loop_edges
            .iter()
            .next()
            .copied()
            .unwrap();
        let hard = -FunctionBuilder::new(GS.emr_mom)
            .add_arg(usize::from(carrier))
            .finish();
        let tagged = |denominator_derived| {
            FunctionBuilder::new(GS.emr_mom)
                .add_arg(GS.uv_momentum_provenance_tag(
                    Atom::num(usize::from(owner) as i64).as_view(),
                    denominator_derived,
                    hard.as_view(),
                ))
                .add_arg(GS.cind(0))
                .finish()
        };
        let fixed = tagged(false);
        let derived = tagged(true);
        let factorized = (&fixed + Atom::one()) * (&derived + Atom::num(2));

        let rescaled = graph.uv_rescaled(&filter, 0, &graph.loop_momentum_basis, &factorized);
        let expected = (&fixed / GS.rescale + Atom::one()) * (&derived / GS.rescale + Atom::num(2));
        assert_eq!(rescaled, expected);
        assert_eq!(
            GS.erase_uv_momentum_provenance(&rescaled),
            GS.erase_uv_momentum_provenance(&graph.uv_rescaled(
                &filter,
                0,
                &graph.loop_momentum_basis,
                &GS.erase_uv_momentum_provenance(&factorized),
            )),
            "erasing provenance must commute with a compatible outer UV rescaling",
        );
        assert!(matches!(rescaled.as_view(), AtomView::Mul(_)));
        Ok(())
    }

    #[test]
    fn uv_taylor_provenance_erasure_matches_plain_child_lmb_expansion() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(
            digraph uv_taylor_provenance_oracle {
                num = 1
                edge [particle="scalar_1" num=1]
                node [num=1]

                ext0 [style=invis is_cut=0]
                v2 -> ext0 [id=0]
                ext0 -> v3
                v0 -> v1 [id=1 lmb_id=0 num="Q(1,spenso::mink(4,1))"]
                v0 -> v1 [id=2]
                v3 -> v0 [id=3 lmb_id=1 num="Q(3,spenso::mink(4,1))"]
                v1 -> v2 [id=4]
                v2 -> v3 [id=5]
            },
            "scalars"
        )?;
        let first = EdgeIndex(1);
        let second = EdgeIndex(2);
        let uv_filter = graph
            .get_edge_subgraph(first)
            .union(&graph.get_edge_subgraph(second));
        let uv_subgraph = InternalSubGraph::cleaned_filter_optimist(uv_filter, graph.as_ref());
        let current = OwnedForestNode {
            spinney: Spinney::with_scheme(
                uv_subgraph,
                &graph,
                &graph.loop_momentum_basis,
                ApproximationType::MUV,
                1,
            )
            .expect("the self-energy bubble has a compatible child LMB"),
            topo_order: 0,
        };
        let given = OwnedForestNode {
            spinney: Spinney::empty(&graph),
            topo_order: 0,
        };
        let reduced = current.reduced_subgraph(&given);
        let mut input = graph
            .numerator(&reduced, given.subgraph())
            .to_d_dim(GS.dim)
            .get_single_atom()
            .unwrap();
        input /= graph.denominator(&reduced, |_| 1);

        let first_hard = current
            .lmb()
            .loop_atom::<Atom>(first, GS.emr_mom, &[], true);
        let second_hard = current
            .lmb()
            .loop_atom::<Atom>(second, GS.emr_mom, &[], true);
        let first_soft = current.lmb().ext_atom::<Atom>(first, GS.emr_mom, &[], true);
        let second_soft = current
            .lmb()
            .ext_atom::<Atom>(second, GS.emr_mom, &[], true);
        let carrier = FunctionBuilder::new(GS.emr_mom).add_arg(1).finish();
        assert_eq!(first_hard, carrier);
        assert_eq!(second_hard, -first_hard.clone());
        assert_eq!(first_soft, Atom::Zero);
        assert_eq!(
            second_soft,
            FunctionBuilder::new(GS.emr_mom).add_arg(4).finish()
        );

        let taylor = |rescaled: Atom| -> Result<Atom> {
            Ok(rescaled
                .series(GS.rescale, Atom::Zero, 0)?
                .to_atom()
                .replace(GS.rescale)
                .with(Atom::one())
                .simplify_metrics()
                .to_dots()
                .normalize_dots())
        };
        let tagged = taylor(graph.uv_rescaled(
            current.subgraph(),
            graph.n_loops(current.subgraph()),
            current.lmb(),
            &input,
        ))?;

        // Independent child-sub-LMB oracle. This duplicates only the scalar
        // rescaling law inside the test and never calls `uv_rescaled`, `t`, or
        // any exact-source reconstruction API.
        let mut plain_rescaled = input.replace_multiple(graph.uv_wrapped_replacement(
            &reduced,
            current.lmb(),
            &[W_.x___],
        ));
        let rescale = Atom::var(GS.rescale);
        for loop_edge in current.lmb().loop_edges.iter() {
            let loop_momentum = function!(GS.emr_mom, usize::from(*loop_edge), W_.x___);
            plain_rescaled = plain_rescaled
                .replace(loop_momentum.clone())
                .with(loop_momentum / &rescale);
        }
        let rescale_squared = rescale.pow(2);
        let expansion_mass_squared = Atom::var(GS.m_uv_expansion).pow(2);
        let vacuum_mass_squared = Atom::var(GS.m_uv_vacuum).pow(2);
        plain_rescaled = plain_rescaled
            .replace(GS.den(W_.a_, W_.mom_, W_.mass_, W_.prop_))
            .with_map({
                let rescale = rescale.clone();
                let rescale_squared = rescale_squared.clone();
                move |matched| {
                    let edge = matched.get(W_.a_).unwrap().to_atom();
                    let momentum = matched.get(W_.mom_).unwrap().to_atom();
                    let mass = matched.get(W_.mass_).unwrap().to_atom();
                    let propagator = matched.get(W_.prop_).unwrap().to_atom();
                    let hard = (momentum * &rescale)
                        .expand()
                        .replace(GS.rescale)
                        .with(Atom::Zero);
                    GS.den(
                        edge,
                        hard,
                        mass * &rescale_squared + &expansion_mass_squared,
                        propagator * &rescale_squared + &expansion_mass_squared * &rescale_squared
                            - &vacuum_mass_squared,
                    ) / &rescale_squared
                }
            });
        plain_rescaled *= rescale.pow(-4);
        let plain = taylor(plain_rescaled)?;
        let erased = GS.erase_uv_momentum_provenance(&tagged);
        assert!(!erased.contains_symbol(GS.uv_momentum_provenance));
        assert!(
            (erased - plain).expand().together().is_zero(),
            "tag erasure must recover the independent child-sub-LMB Taylor coefficient",
        );

        let mut tags = BTreeMap::new();
        let mut malformed = false;
        let _ = tagged.replace_map(|view, _, _| {
            let AtomView::Fun(momentum) = view else {
                return;
            };
            if momentum.get_symbol() != GS.emr_mom || momentum.get_nargs() < 1 {
                return;
            }
            if let Some((owner, role, hard)) = GS.uv_momentum_provenance_data(momentum.get(0)) {
                let denominator_derived = role == UvMomentumProvenanceRole::DenominatorDerived;
                if let Some(previous) = tags.insert((owner, denominator_derived), hard.clone()) {
                    assert_eq!(previous, hard);
                }
            } else if matches!(
                momentum.get(0),
                AtomView::Fun(provenance)
                    if provenance.get_symbol() == GS.uv_momentum_provenance
            ) {
                malformed = true;
            }
        });
        assert!(
            !malformed,
            "no transient provenance role may survive Taylor expansion"
        );
        assert_eq!(
            tags,
            BTreeMap::from([
                ((first, false), first_hard.clone()),
                ((first, true), first_hard.clone()),
                ((second, true), second_hard.clone()),
            ])
        );

        // Expand and cancel only inside this ownership oracle. Production
        // retains the common denominator and the factorized Taylor numerator.
        let expanded = tagged.expand();
        let expanded_terms = match expanded.as_view() {
            AtomView::Add(add) => add.iter().map(|term| term.to_owned()).collect(),
            _ => vec![expanded],
        };
        let mut leaves = Vec::new();
        for expanded_term in expanded_terms {
            leaves.extend(FourDTerm::from_view(expanded_term.cancel().as_view())?);
        }
        leaves.sort_by_key(|leaf| {
            leaf.denominators
                .iter()
                .filter(|denominator| denominator.source_edge == second)
                .count()
        });
        assert_eq!(leaves.len(), 2);
        assert_eq!(
            leaves
                .iter()
                .map(|leaf| {
                    [first, second].map(|owner| {
                        leaf.denominators
                            .iter()
                            .filter(|denominator| denominator.source_edge == owner)
                            .count()
                    })
                })
                .collect::<Vec<_>>(),
            vec![[1, 1], [1, 2]],
        );
        let numerator_tags = leaves
            .iter()
            .map(|leaf| {
                let mut numerator_tags = BTreeMap::new();
                let _ = leaf.numerator.replace_map(|view, _, _| {
                    let AtomView::Fun(momentum) = view else {
                        return;
                    };
                    if momentum.get_symbol() == GS.emr_mom
                        && momentum.get_nargs() >= 1
                        && let Some((owner, role, hard)) =
                            GS.uv_momentum_provenance_data(momentum.get(0))
                    {
                        let denominator_derived =
                            role == UvMomentumProvenanceRole::DenominatorDerived;
                        numerator_tags.insert((owner, denominator_derived), hard);
                    }
                });
                numerator_tags
            })
            .collect::<Vec<_>>();
        assert_eq!(
            numerator_tags[0],
            BTreeMap::from([((first, false), first_hard.clone())]),
            "the original numerator must remain fixed on its source edge",
        );
        assert_eq!(
            numerator_tags[1],
            BTreeMap::from([((first, false), first_hard), ((second, true), second_hard),]),
            "only the new denominator-derivative momentum may dispatch on the raised line",
        );
        Ok(())
    }

    #[test]
    fn gl24_dod_two_q1_quartic_taylor_keeps_owner_local_energy_families() -> Result<()> {
        test_initialise()?;
        // This is the GL24 skeleton and generation LMB used by the scalar LU
        // acceptance test. The quartic factor is local to e1. Algebraic leaf
        // extraction below is test-only; production retains one factorized
        // common-denominator Taylor numerator.
        let graph: Graph = dot!(
            digraph gl24_dod_two_taylor {
                edge [particle="scalar_0" num=1]
                node [num=1]
                incoming [style=invis]
                outgoing [style=invis]

                incoming -> v1 [id=0 particle="scalar_1"]
                v0 -> v3 [id=1 lmb_id=0 particle="scalar_1" num="Q(1,spenso::mink(4,1))^2*Q(1,spenso::mink(4,2))^2"]
                v0 -> v5 [id=2 particle="scalar_2"]
                v1 -> v2 [id=3 lmb_id=1]
                v1 -> v3 [id=4]
                v2 -> v4 [id=5 lmb_id=2]
                v2 -> v4 [id=6]
                v3 -> v5 [id=7 particle="scalar_1"]
                v4 -> v5 [id=8 particle="scalar_1"]
                v0 -> outgoing [id=9 particle="scalar_1"]
            },
            "scalars"
        )?;
        let owners = [EdgeIndex(1), EdgeIndex(2), EdgeIndex(7)];
        let uv_filter = owners
            .into_iter()
            .map(|edge| graph.get_edge_subgraph(edge))
            .reduce(|left, right| left.union(&right))
            .expect("the GL24 UV triangle has three edges");
        let uv_subgraph =
            InternalSubGraph::cleaned_filter_optimist(uv_filter.clone(), graph.as_ref());
        let current = OwnedForestNode {
            spinney: Spinney::new(uv_subgraph, &graph, &graph.loop_momentum_basis)
                .expect("the GL24 UV triangle has a compatible child LMB"),
            topo_order: 0,
        };
        assert_eq!(current.spinney.dod, 2);

        let minkowski = LibraryRep::from(Minkowski {});
        let q1_first = GS.emr_mom(owners[0], minkowski.to_symbolic([Atom::num(1)]));
        let q1_second = GS.emr_mom(owners[0], minkowski.to_symbolic([Atom::num(2)]));
        let numerator = q1_first.pow(2) * q1_second.pow(2);
        let integrand = &numerator / graph.denominator(&uv_filter, |_| 1);
        // DOD two starts at t^-2, so the t^0 Laurent coefficient is exactly
        // the second-order Taylor layer, without its leading and linear peers.
        let expanded = graph
            .uv_rescaled(
                current.subgraph(),
                graph.n_loops(current.subgraph()),
                current.lmb(),
                &integrand,
            )
            .series(GS.rescale, Atom::Zero, 0)?
            .coefficient(Rational::from(0))
            .simplify_metrics()
            .to_dots()
            .normalize_dots()
            .together();

        let common_terms = FourDTerm::from_view(expanded.as_view())?;
        assert_eq!(common_terms.len(), 1);
        let common = &common_terms[0];
        let owner_multiplicities = owners.map(|owner| {
            common
                .denominators
                .iter()
                .filter(|denominator| denominator.source_edge == owner)
                .count()
        });
        assert_eq!(
            owner_multiplicities,
            [2, 3, 3],
            "the DOD-two Taylor coefficient must retain the e1^2 e2^3 e7^3 common denominator",
        );

        let excluded_edges = graph
            .underlying
            .iter_edges()
            .filter_map(|(pair, edge, edge_data)| {
                (pair.is_paired() && !edge_data.data.is_dummy && !owners.contains(&edge))
                    .then_some(edge)
            })
            .collect::<Vec<_>>();
        let physical_bounds = graph
            .automatic_numerator_energy_degree_bounds_in_atom_excluding_with_min_degree(
                &common.numerator,
                excluded_edges.iter().copied(),
                1,
            )?;
        assert_eq!(physical_bounds, vec![(1, 6), (2, 4), (7, 4)]);

        let source = GraphThreeDSource::from_exact_denominators_in_uv_edges(
            &graph,
            &common.denominators,
            owners,
        )?;
        let occurrence_owners = source
            .physical_energy_edge_index_map()
            .expect("the exact GL24 source exposes occurrence owners")
            .internal;
        let candidates = source
            .exact_source_energy_mapper()
            .expect("the exact GL24 source exposes an energy mapper")
            .equivalent_energy_candidates(owners)?;
        let plan = graph.plan_numerator_energy_assignment_in_atom_excluding(
            &common.numerator,
            excluded_edges,
            &candidates,
        )?;
        let degree_by_occurrence = plan
            .energy_degree_bounds()
            .iter()
            .copied()
            .collect::<BTreeMap<_, _>>();
        let mut loads_by_owner = BTreeMap::<usize, Vec<(usize, usize)>>::new();
        for (occurrence, owner) in occurrence_owners {
            loads_by_owner.entry(owner).or_default().push((
                occurrence,
                degree_by_occurrence.get(&occurrence).copied().unwrap_or(0),
            ));
        }
        for loads in loads_by_owner.values_mut() {
            loads.sort_unstable();
        }
        assert_eq!(
            loads_by_owner[&1]
                .iter()
                .map(|(_, degree)| *degree)
                .collect::<Vec<_>>(),
            vec![4, 2],
            "the fixed q1^4 load stays on the canonical e1 occurrence",
        );
        for owner in [2, 7] {
            let mut loads = loads_by_owner[&owner]
                .iter()
                .map(|(_, degree)| *degree)
                .collect::<Vec<_>>();
            loads.sort_unstable();
            assert_eq!(
                loads,
                vec![0, 2, 2],
                "the degree-four derived load has minimax degree two on owner e{owner}",
            );
        }

        let provenance = |atom: &Atom| {
            let mut tags = BTreeSet::new();
            let _ = atom.replace_map(|view, _, _| {
                let AtomView::Fun(momentum) = view else {
                    return;
                };
                if momentum.get_symbol() == GS.emr_mom
                    && momentum.get_nargs() >= 1
                    && let Some((owner, role, _)) = GS.uv_momentum_provenance_data(momentum.get(0))
                {
                    let denominator_derived = role == UvMomentumProvenanceRole::DenominatorDerived;
                    tags.insert((usize::from(owner), denominator_derived));
                }
            });
            tags
        };
        assert_eq!(
            provenance(&common.numerator),
            BTreeSet::from([(1, false), (1, true), (2, true), (7, true)]),
        );

        // Reduce only the test copy to its natural Taylor leaves. Regrouping by
        // denominator multiplicity restores each factorized C_i after the
        // test-only expansion of its soft-shift and mass pieces.
        let expanded_for_oracle = expanded.expand();
        let additive_terms = match expanded_for_oracle.as_view() {
            AtomView::Add(add) => add.iter().map(|term| term.to_owned()).collect(),
            _ => vec![expanded_for_oracle],
        };
        let analyzer = EnergyPowerAnalyzer::for_physical_emr_edges(owners);
        let mut reduced_numerators = BTreeMap::<[usize; 3], Atom>::new();
        for additive_term in additive_terms {
            for leaf in FourDTerm::from_view(additive_term.cancel().as_view())? {
                if leaf.numerator.is_zero() {
                    continue;
                }
                let multiplicities = owners.map(|owner| {
                    leaf.denominators
                        .iter()
                        .filter(|denominator| denominator.source_edge == owner)
                        .count()
                });
                let reduced_numerator = reduced_numerators
                    .entry(multiplicities)
                    .or_insert(Atom::Zero);
                *reduced_numerator = &*reduced_numerator + &leaf.numerator;
            }
        }
        let reduced_leaves = reduced_numerators
            .iter()
            .map(|(multiplicities, numerator)| {
                Ok((
                    (
                        *multiplicities,
                        analyzer.analyze_atom(numerator)?.into_generation_bounds(),
                    ),
                    provenance(numerator),
                ))
            })
            .collect::<Result<BTreeMap<_, _>>>()?;
        let fixed_e1 = BTreeSet::from([(1, false)]);
        assert_eq!(
            reduced_leaves.len(),
            6,
            "the second-order coefficient has exactly three derivative-hit and three constant leaves",
        );
        for (key, expected_provenance) in [
            (
                ([1, 3, 1], vec![(1, 4), (2, 2)]),
                BTreeSet::from([(1, false), (2, true)]),
            ),
            (
                ([1, 2, 2], vec![(1, 4), (2, 1), (7, 1)]),
                BTreeSet::from([(1, false), (2, true), (7, true)]),
            ),
            (
                ([1, 1, 3], vec![(1, 4), (7, 2)]),
                BTreeSet::from([(1, false), (7, true)]),
            ),
            (([2, 1, 1], vec![(1, 4)]), fixed_e1.clone()),
            (([1, 2, 1], vec![(1, 4)]), fixed_e1.clone()),
            (([1, 1, 2], vec![(1, 4)]), fixed_e1.clone()),
        ] {
            assert_eq!(
                reduced_leaves.get(&key),
                Some(&expected_provenance),
                "missing or malformed second-order GL24 Taylor leaf {key:?}",
            );
        }

        let expansion_mass_squared = Atom::var(GS.m_uv_expansion).pow(2);
        for (position, owner) in owners.into_iter().enumerate() {
            let mut constant_leaf = [1, 1, 1];
            constant_leaf[position] += 1;
            let coefficient = &reduced_numerators[&constant_leaf];
            let physical_mass = graph.underlying[owner].particle.mass_atom();
            let AtomView::Var(physical_mass_variable) = physical_mass.as_view() else {
                panic!("the mass-sensitive GL24 fixture must use symbolic owner masses");
            };
            let without_mass_constants = coefficient
                .replace(GS.m_uv_expansion)
                .with(Atom::Zero)
                .replace(physical_mass_variable.get_symbol())
                .with(Atom::Zero);
            let actual_mass_part =
                GS.erase_uv_momentum_provenance(&(coefficient - without_mass_constants).expand());
            let expected_mass_part =
                (&numerator * (physical_mass.pow(2) - &expansion_mass_squared)).expand();
            assert_eq!(
                actual_mass_part,
                expected_mass_part,
                "the C_{} leaf must contain m_{}^2 - mUVexp^2 with zero EMR degree",
                position + 1,
                usize::from(owner),
            );
        }
        let vacuum_mass_squared = Atom::var(GS.m_uv_vacuum).pow(2);
        for coefficient in [expansion_mass_squared, vacuum_mass_squared]
            .into_iter()
            .chain(owners.map(|edge| graph.underlying[edge].particle.mass_atom().pow(2)))
        {
            assert!(
                analyzer.analyze_atom(&coefficient)?.is_empty(),
                "a mass-only coefficient must not request an EMR candidate: {coefficient}",
            );
        }
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
    fn term_projection_keeps_typed_numerator_factors_uncancelled() -> Result<()> {
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
        let ordinary_sum = left.clone() + right.clone();
        let factorized_numerator = &ordinary_sum * &first + &second;
        let numerator_terms = FourDTerm::from_view(factorized_numerator.as_view())?;
        assert_eq!(numerator_terms.len(), 1);
        assert!(numerator_terms[0].denominators.is_empty());
        assert_eq!(numerator_terms[0].numerator, factorized_numerator);

        let shared_topology = left.clone() * first.pow(-1) + right.clone() * first.pow(-1);
        let shared_terms = FourDTerm::from_view(shared_topology.as_view())?;
        assert_eq!(shared_terms.len(), 1);
        assert_eq!(shared_terms[0].denominators.len(), 1);
        assert_eq!(shared_terms[0].numerator, ordinary_sum);

        let distinct_topologies = left * first.pow(-1) + right * second.pow(-1);
        let distinct_terms = FourDTerm::from_view(distinct_topologies.as_view())?;
        assert_eq!(distinct_terms.len(), 2);
        assert!(
            distinct_terms
                .iter()
                .all(|term| term.denominators.len() == 1)
        );

        let collected = &factorized_numerator * first.pow(-1) * second.pow(-1);

        let terms = FourDTerm::from_view(collected.as_view())?;
        assert_eq!(terms.len(), 1);
        assert_eq!(terms[0].denominators.len(), 2);
        assert!(
            terms[0]
                .denominators
                .contains(&FourDDenominator::from_view(first.as_view())?.unwrap())
        );
        assert!(
            terms[0]
                .denominators
                .contains(&FourDDenominator::from_view(second.as_view())?.unwrap())
        );
        assert_eq!(terms[0].numerator, &ordinary_sum * &first + &second);
        assert!(terms[0].numerator.contains_symbol(GS.den));
        Ok(())
    }

    #[test]
    fn dod_one_triangle_keeps_one_uncancelled_exact_source() -> Result<()> {
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
        let (expanded, _) = t(
            &integrand,
            &UVCtx::new(&graph, &settings),
            &current,
            &given,
            &[],
            0,
        )?;
        let cograph = graph
            .get_edge_subgraph(EdgeIndex(1))
            .union(&graph.get_edge_subgraph(EdgeIndex(2)));
        let terms = Full4dCts::from_coefficient(&expanded, &graph, &cograph).terms()?;

        assert_eq!(terms.len(), 1);
        let owner_multiplicities = owners.map(|owner| {
            terms[0]
                .denominators
                .iter()
                .filter(|denominator| denominator.source_edge == owner)
                .count()
        });
        assert_eq!(
            owner_multiplicities,
            [2, 1, 2],
            "the common-denominator Taylor coefficient must retain every raised occurrence",
        );
        let has_additive_factor = match terms[0].numerator.as_view() {
            AtomView::Add(_) => true,
            AtomView::Mul(product) => product
                .iter()
                .any(|factor| matches!(factor, AtomView::Add(_))),
            _ => false,
        };
        assert!(
            has_additive_factor && terms[0].numerator.contains_symbol(GS.den),
            "the additive Taylor numerator and its positive typed denominator factors must remain factorized and uncancelled: {}",
            terms[0].numerator,
        );

        let options = graph.denominator_only_cff_3d_expression_options();
        let mut cache = ExactCffGenerationCache::default();
        let active_denominators = terms[0]
            .denominators
            .iter()
            .map(|denominator| {
                let is_uv = uv_filter.includes(&graph[&denominator.source_edge].1);
                Ok(denominator
                    .depends_on_loop(&graph, is_uv)?
                    .then(|| denominator.clone()))
            })
            .collect::<Result<Vec<_>>>()?
            .into_iter()
            .flatten()
            .collect::<Vec<FourDDenominator>>();
        let source = GraphThreeDSource::from_exact_denominators_in_uv_edges(
            &graph,
            &active_denominators,
            owners,
        )?;
        graph.register_3d_expression_for_4d_term(
            &source,
            &options,
            &terms[0].numerator,
            &mut cache,
        )?;
        let (_, _, plan, _) = graph.generate_3d_expression_for_4d_term(
            &source,
            &options,
            &terms[0].numerator,
            Some(&mut cache),
        )?;
        assert_eq!(
            cache.len(),
            1,
            "the uncancelled Taylor coefficient must generate one common-denominator canonical CFF topology; exact bounds: {:?}",
            plan.energy_degree_bounds(),
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
            let active_subgraph = graph.get_edge_subgraph(
                *lmb.loop_edges
                    .first()
                    .expect("the test component has one loop carrier"),
            );
            Full4dCts(FourDSectors::new(
                vec![FourDSector::new(
                    local,
                    vec![(active_subgraph.clone(), active_subgraph, lmb.clone())],
                    Some(lmb.clone()),
                    Vec::new(),
                )],
                vec![FourDSector::new(finite, Vec::new(), None, vec![lmb])],
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
            FourDSector::new(&finite_a * &finite_b, Vec::new(), None, vec![lmb_a, lmb_b])
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
