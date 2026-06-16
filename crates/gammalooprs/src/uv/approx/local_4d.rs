use std::ops::Neg;

use color_eyre::Result;
use eyre::eyre;
use gammaloop_tracing_filter::debug_instrument;
use idenso::{
    color::ColorSimplifier,
    dirac::{AGS, GammaSimplifier},
    representations::Bispinor,
    shorthands::{
        chain::Chain,
        metric::MetricSimplifier,
        schoonschip::{Schoonschip, SchoonschipSettings},
    },
};

#[cfg(test)]
use linnet::half_edge::subgraph::SubGraphLike;
use linnet::half_edge::{
    involution::EdgeIndex,
    subgraph::{Inclusion, InternalSubGraph, SuBitGraph, SubSetLike},
};
use spenso::shadowing::TensorCollectExt;
use symbolica::{
    atom::{AtomCore, AtomView},
    prelude::*,
};

use crate::utils::symbols::UvMomentumProvenanceRole;
use crate::{
    debug_tags,
    graph::{FourDDenominator, Graph, LMBext, LoopMomentumBasis},
    numerator::{aind::Aind, ufo::UFO},
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
    // the quotient LMB in which its energy residues must be generated. These
    // independent component frames are distinct from frozen LMBs whose
    // integrations have already completed and only own localization factors.
    pub(crate) active_components: Vec<(SuBitGraph, SuBitGraph, LoopMomentumBasis)>,
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
        frozen_lmbs: Vec<LoopMomentumBasis>,
    ) -> Self {
        Self {
            atom,
            active_components,
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
            vec![FourDSector::new(atom, Vec::new(), Vec::new())],
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
        let completed = |atom| FourDSector::new(atom, Vec::new(), vec![lmb.clone()]);
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
    #[cfg(test)]
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

    /// Enumerate only the additive terms needed by the 3D projection while
    /// retaining every non-denominator factor as a Symbolica atom. In
    /// particular, numerator products and powers are not materialized into a
    /// parallel expression tree or expanded polynomial.
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
                    FourDSector::new(-sector.atom, sector.active_components, sector.frozen_lmbs)
                })
                .collect(),
            self.0
                .recursive_completion
                .into_iter()
                .map(|sector| {
                    FourDSector::new(-sector.atom, sector.active_components, sector.frozen_lmbs)
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
        let mut products = vec![(FourDSector::new(Atom::one(), Vec::new(), Vec::new()), false)];
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
                            let mut active_components = left.active_components.clone();
                            active_components.extend(right.active_components);
                            let mut frozen_lmbs = left.frozen_lmbs.clone();
                            frozen_lmbs.extend(right.frozen_lmbs);
                            (
                                FourDSector::new(
                                    left.atom.clone() * right.atom,
                                    active_components,
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

    pub(super) fn from_view(view: AtomView<'_>) -> Result<Vec<Self>> {
        // A denominator-free subtree is an opaque numerator. Inspect its symbols
        // once instead of recursively rebuilding its sums/products just to prove
        // that every child has the same empty denominator topology.
        if matches!(view, AtomView::Add(_) | AtomView::Mul(_)) && !view.contains_symbol(GS.den) {
            return Ok(vec![Self::numerator(view.to_owned())]);
        }
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

impl Graph {
    pub(crate) fn uv_rescaled(
        &self,
        expansion_subgraph: &SuBitGraph,
        n_loops: usize,
        lmb: &LoopMomentumBasis,
        reference_lmb: &LoopMomentumBasis,
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
                // belongs to this Taylor operation from persistent role-zero or
                // role-one tags retained by a nested child, role-three tags
                // fixing physical-source momenta, or role-four tags retaining
                // the literal carrier of a new soft crown momentum.
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
                            // Existing child provenance retains its owner and
                            // role. Re-emitting the complete momentum makes this
                            // node opaque to the top-down traversal, so the outer
                            // operation cannot tag a literal Q inside its stored
                            // hard-momentum payload before affine transport.
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
        // subgraph while retaining its immutable source owner. The compatible
        // LMB supplies the carrier coordinates, while the current UV node fixes
        // the soft routing: Q_e(t) = (Q_e - S_current(e))/t + S_current(e).
        // Keeping the hard part as one tagged rank-one momentum is essential:
        // a derivative of D_e must not lose its owner when the LMB happens to
        // spell its hard momentum with another physical edge ID.
        let inverse_rescale = Atom::one() / GS.rescale;
        let mut atomarg = tagged.replace_map(|view, _, output| {
            let AtomView::Fun(momentum) = view else {
                return;
            };
            if momentum.get_symbol() != GS.emr_mom || momentum.get_nargs() == 0 {
                return;
            }
            let (owner, role, routed) = if let Some((owner, role, mut routed)) =
                GS.uv_momentum_provenance_data(momentum.get(0))
            {
                // A hard child tag retains its complete frozen hard projection;
                // a soft child tag retains its literal crown carrier. Express
                // that literal carrier in the compatible enclosing LMB, just
                // like an ordinary Q, before the hard/soft split. Otherwise
                // a newly hard numerator can retain an independent Q absent
                // from the vacuum denominator's momentum system.
                if role == UvMomentumProvenanceRole::DenominatorDerivedSoft {
                    if !scaled_edges.contains(&owner) {
                        **output = Atom::from(momentum.to_owned());
                        return;
                    }
                    routed = lmb.loop_atom::<Atom>(owner, GS.emr_mom, &[], true)
                        + lmb.ext_atom::<Atom>(owner, GS.emr_mom, &[], true);
                }
                // Retained hard payloads can have a soft shift in this LMB.
                // Existing hard ownership stays fixed until assignment to CFF.
                // A previously soft carrier receives a new fixed hard tag below,
                // independently of the child's denominator.
                (owner, role, routed)
            } else {
                let (owner, derived) = match momentum.get(0) {
                    AtomView::Fun(provenance)
                        if provenance.get_symbol() == GS.uv_momentum_provenance
                            && provenance.get_nargs() == 3
                            && provenance.get(1) == Atom::num(2).as_view() =>
                    {
                        (usize::try_from(provenance.get(0)), true)
                    }
                    edge => (usize::try_from(edge), false),
                };
                let Ok(owner) = owner.map(EdgeIndex) else {
                    return;
                };
                if !scaled_edges.contains(&owner) {
                    return;
                }
                (
                    owner,
                    derived.into(),
                    lmb.loop_atom::<Atom>(owner, GS.emr_mom, &[], true)
                        + lmb.ext_atom::<Atom>(owner, GS.emr_mom, &[], true),
                )
            };
            let soft = routed.replace_map(|view, _, output| {
                if let AtomView::Fun(momentum) = view
                    && momentum.get_symbol() == GS.emr_mom
                    && momentum.get_nargs() == 1
                    && let Ok(edge) = usize::try_from(momentum.get(0))
                    && !reference_lmb.ext_edges.contains(&EdgeIndex(edge))
                {
                    // External-coordinate carriers are fixed literally. A
                    // paired crown edge can have a zero row outside this UV
                    // subgraph even though it names one of these coordinates.
                    **output =
                        reference_lmb.ext_atom::<Atom>(EdgeIndex(edge), GS.emr_mom, &[], true);
                }
            });
            let hard = (&routed - &soft).expand();
            // A child coefficient's soft carrier can become hard in an
            // enclosing Taylor limit. It is then fixed to that carrier, not a
            // derivative-created copy of the child's original denominator.
            let hard_role = match role {
                UvMomentumProvenanceRole::DenominatorDerivedSoft => {
                    UvMomentumProvenanceRole::TaylorFixed
                }
                role => role,
            };
            let tag = GS.uv_momentum_provenance_tag(
                Atom::num(usize::from(owner) as i64).as_view(),
                hard_role,
                hard.as_view(),
            );
            let mut tagged_hard = FunctionBuilder::new(GS.emr_mom).add_arg(tag);
            for index in momentum.iter().skip(1) {
                tagged_hard = tagged_hard.add_arg(index.to_owned());
            }
            // In particular, a child soft carrier can remain entirely
            // external to an enclosing limit. Do not manufacture a rank-one
            // tagged momentum for its exactly vanishing hard projection.
            let hard_component = if hard.is_zero() {
                Atom::Zero
            } else {
                tagged_hard.finish()
            };
            let soft = soft.replace_map(|view, _, output| {
                if let AtomView::Fun(soft_momentum) = view
                    && soft_momentum.get_symbol() == GS.emr_mom
                    && soft_momentum.get_nargs() == 1
                {
                    let source = if matches!(
                        role,
                        UvMomentumProvenanceRole::DenominatorDerived
                            | UvMomentumProvenanceRole::DenominatorDerivedSoft
                    ) {
                        GS.uv_momentum_provenance_tag(
                            soft_momentum.get(0),
                            UvMomentumProvenanceRole::DenominatorDerivedSoft,
                            view,
                        )
                    } else {
                        soft_momentum.get(0).to_owned()
                    };
                    let mut component = FunctionBuilder::new(GS.emr_mom).add_arg(source);
                    for index in momentum.iter().skip(1) {
                        component = component.add_arg(index.to_owned());
                    }
                    **output = component.finish();
                }
            });
            **output = hard_component * &inverse_rescale + soft;
        });

        // Free `mUVexp` occurrences are left untouched here: with the
        // inverse loop-momentum expansion, soft dependence is generated by
        // the shifted hard momenta and denominator expansion rather than by
        // an explicit soft rescaling pass.

        // Vacuum masses in active Taylor terms and previously integrated
        // coefficients carry the same hard weight as the parent loop momenta.
        // This preserves the forest's prescribed Taylor order even when the
        // remaining quotient would converge with its coefficient held fixed.
        // The denominator rewrite below introduces the unscaled vacuum mass
        // for the parent propagator basis.
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

    let rescaled = graph.uv_rescaled(current.subgraph(), n_loops, lmb, current.lmb(), integrand);
    debug_tags!(#uv,#integrated,#rescaled;log.res = rescaled, n_loops=%n_loops,"Rescaled expanded");

    let series = rescaled
        .series(GS.rescale, Atom::Zero, 0)
        .map_err(|error| {
            eyre!(
                "local 4D Taylor expansion failed for {}: {error}",
                current.subgraph().string_label()
            )
        })?
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
    // Keep each Taylor term's propagator powers. Collecting factors across the
    // sum clears denominators and manufactures higher-rank numerator factors.
    let collected = schoonschip
        .collect_chains_and_traces()
        .simplify_metrics()
        .collect_gamma_chains()
        .collect_color();
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
            if ctx.settings.generate_integrated {
                // Inspect the original UV-subgraph factors before multiplication,
                // zero-sector pruning, or Dirac simplification can hide gamma5
                // pairs. Chiral projectors contain gamma5 implicitly. Cograph
                // factors and external-state projectors are outside this integral.
                let forbidden = [
                    AGS.gamma5, AGS.projm, AGS.projp, UFO.gamma5, UFO.projm, UFO.projp,
                ];
                let edges = ctx
                    .graph
                    .underlying
                    .iter_edges_of(current.subgraph())
                    .filter(|(pair, _, _)| pair.is_paired())
                    .map(|(_, edge, data)| ("edge", edge.0, &data.data.num.value));
                let vertices = ctx
                    .graph
                    .underlying
                    .iter_nodes_of(current.subgraph())
                    .map(|(vertex, _, data)| ("vertex", vertex.0, &data.num.value));
                for (kind, index, numerator) in edges.chain(vertices) {
                    if let Some(symbol) = forbidden
                        .iter()
                        .find(|symbol| numerator.contains_symbol(**symbol))
                    {
                        return Err(eyre!(
                            "Cannot analytically integrate UV subgraph {} of graph '{}': {kind} {index} contains {symbol}. Gamma5 and chiral projectors are unsupported in d-dimensional analytic UV numerator algebra; no gamma5 prescription is implemented.",
                            current.subgraph().string_label(),
                            ctx.graph.name,
                        ));
                    }
                }
            }
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
                // Store the forest factor -T here. Later CFF projection acts
                // on this signed coefficient without another subtraction minus.
                Ok(FourDSector::new(
                    marker.apply(
                        UvOperation::Approx,
                        marker_current.subgraph(),
                        marker_given.subgraph(),
                        &-result,
                    ),
                    active_components,
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
    fn analytic_uv_rejects_gamma5_before_simplification_with_subgraph_scope() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph gamma5_uv_scope {
            edge [num=1 mass=1]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]
            incoming -> a [id=0]
            a -> b [id=1 lmb_id=0]
            b -> a [id=2]
            b -> c [id=3]
            c -> d [id=4 lmb_id=1]
            d -> c [id=5]
            d -> outgoing [id=6]
        })?;
        let filter = graph
            .get_edge_subgraph(EdgeIndex(1))
            .union(&graph.get_edge_subgraph(EdgeIndex(2)));
        let sibling = graph
            .get_edge_subgraph(EdgeIndex(4))
            .union(&graph.get_edge_subgraph(EdgeIndex(5)));
        let vertex = graph.underlying.iter_nodes_of(&filter).next().unwrap().0;
        let sibling_vertex = graph.underlying.iter_nodes_of(&sibling).next().unwrap().0;
        let mut current = OwnedForestNode {
            spinney: Spinney::with_scheme(
                InternalSubGraph::cleaned_filter_optimist(filter, graph.as_ref()),
                &graph,
                &graph.loop_momentum_basis,
                ApproximationType::MUV,
                0,
            )
            .expect("the scalar bubble has a compatible UV loop-momentum basis"),
            topo_order: 0,
        };
        let given = OwnedForestNode {
            spinney: Spinney::empty(&graph),
            topo_order: 0,
        };
        let settings = UVgenerationSettings::default();
        let input = Full4dCts(FourDSectors::active_atom(Atom::one()));
        let run = |graph: &Graph, current: &OwnedForestNode, settings: &UVgenerationSettings| {
            uv_limit(
                &input,
                &UVCtx::new(graph, settings),
                current,
                &given,
                current,
                &given,
            )
        };
        let start = idenso::bis!(4, Atom::from(Aind::Normal(17)));
        let end = idenso::bis!(4, Atom::from(Aind::Normal(18)));
        let open = idenso::gamma5!(&start, &end);
        let closed = spenso::trace!(Bispinor {}.new_rep(4).to_symbolic([]), idenso::gamma5!());
        let pair = spenso::chain!(&start, &end, idenso::gamma5!(), idenso::gamma5!());
        assert!(!pair.simplify_gamma().contains_symbol(AGS.gamma5));
        let cases = [
            ("open", open.clone()),
            ("closed", closed),
            ("pair", pair.clone()),
            ("left projector", function!(AGS.projm, &start, &end)),
            ("right projector", function!(AGS.projp, &start, &end)),
            ("UFO gamma5", function!(UFO.gamma5, 1, 2)),
            ("UFO left projector", function!(UFO.projm, 1, 2)),
            ("UFO right projector", function!(UFO.projp, 1, 2)),
        ];
        for scheme in [ApproximationType::MUV, ApproximationType::PolePart] {
            current.spinney.renormalization_scheme = scheme;
            for (label, numerator) in &cases {
                for on_vertex in [false, true] {
                    if on_vertex {
                        graph.underlying[vertex].num.value = numerator.clone();
                    } else {
                        graph.underlying[EdgeIndex(1)].num.value = numerator.clone();
                    }
                    let error = run(&graph, &current, &settings).unwrap_err().to_string();
                    assert!(
                        error.contains("d-dimensional analytic UV numerator algebra"),
                        "{label}: {error}"
                    );
                    assert!(error.contains(&graph.name), "{label}: {error}");
                    assert!(
                        error.contains(&current.subgraph().string_label()),
                        "{label}: {error}"
                    );
                    let source = if on_vertex {
                        format!("vertex {}", vertex.0)
                    } else {
                        "edge 1".to_owned()
                    };
                    assert!(error.contains(&source), "{label}: {error}");
                    graph.underlying[vertex].num.value = Atom::one();
                    graph.underlying[EdgeIndex(1)].num.value = Atom::one();
                }
            }
        }

        let expected = run(&graph, &current, &settings)?;
        assert!(!expected.atom().is_zero());
        // A sibling loop, its vertices, and global external-state projectors
        // belong to the cograph of this integration and must remain untouched.
        graph.underlying[EdgeIndex(4)].num.value = open.clone();
        graph.underlying[sibling_vertex].num.value = open.clone();
        graph.global_prefactor.projector = open;
        assert_eq!(run(&graph, &current, &settings)?.atom(), expected.atom());
        assert!(
            graph.underlying[EdgeIndex(4)]
                .num
                .value
                .contains_symbol(AGS.gamma5)
        );
        assert!(
            graph.underlying[sibling_vertex]
                .num
                .value
                .contains_symbol(AGS.gamma5)
        );
        assert!(graph.global_prefactor.projector.contains_symbol(AGS.gamma5));

        graph.underlying[EdgeIndex(1)].num.value = pair;
        let local_only = UVgenerationSettings {
            generate_integrated: false,
            ..settings
        };
        assert!(run(&graph, &current, &local_only).is_ok());
        // Even a vanishing incoming sector must not bypass the source check.
        let zero = Full4dCts(FourDSectors::active_atom(Atom::Zero));
        assert!(
            uv_limit(
                &zero,
                &UVCtx::new(&graph, &UVgenerationSettings::default()),
                &current,
                &given,
                &current,
                &given,
            )
            .unwrap_err()
            .to_string()
            .contains("Gamma5")
        );
        Ok(())
    }

    #[test]
    fn term_projection_preserves_complete_factorized_values() -> Result<()> {
        test_initialise()?;
        let first_full = Atom::var(symbol!("local_4d_test::first_full"));
        let second_full = Atom::var(symbol!("local_4d_test::second_full"));
        let first = GS.den(0, function!(GS.emr_mom, 0), 0, &first_full);
        let second = GS.den(0, function!(GS.emr_mom, 0), 0, &second_full);
        let left = Atom::var(symbol!("local_4d_test::left"));
        let right = Atom::var(symbol!("local_4d_test::right"));
        let affine = symbol!("local_4d_test::Affine");
        let numerator =
            (function!(affine, &left) + function!(affine, &right)).pow(2) * (&left + &right);
        let typed_numerator = (&left + &right) * &first + &second;
        for expression in [
            numerator.clone(),
            &numerator * (&left / (&first * &second) + &right / (&first * second.pow(2))),
            first.pow(-2) * second.pow(-1),
            typed_numerator.clone(),
            &left / &first + &right / &first,
            &left / &first + &right / &second,
            typed_numerator / (&first * &second),
        ] {
            let terms = FourDTerm::from_view(expression.as_view())?;
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
            assert!(
                (&expression - reconstructed).expand().is_zero(),
                "projection must preserve the complete value and distinct denominators on one owner: {expression}",
            );
        }
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

        let rescaled = graph.uv_rescaled(
            &filter,
            0,
            &graph.loop_momentum_basis,
            &graph.loop_momentum_basis,
            &factorized,
        );
        let expected = (&fixed / GS.rescale + Atom::one()) * (&derived / GS.rescale + Atom::num(2));
        assert_eq!(rescaled, expected);
        assert_eq!(
            GS.erase_uv_momentum_provenance(&rescaled),
            GS.erase_uv_momentum_provenance(&graph.uv_rescaled(
                &filter,
                0,
                &graph.loop_momentum_basis,
                &graph.loop_momentum_basis,
                &GS.erase_uv_momentum_provenance(&factorized),
            )),
            "erasing provenance must commute with a compatible outer UV rescaling",
        );

        // The child-soft carrier is a literal edge, not a previously frozen
        // hard expression. In this directed vacuum triangle Q(1) = Q(0);
        // the enclosing vacuum denominator system uses only carrier Q(0).
        // Check the raw payload, before any later graph identity can hide an
        // uneliminated numerator-only Q(1) from analytic integration.
        let carrier_momentum = FunctionBuilder::new(GS.emr_mom)
            .add_arg(usize::from(carrier))
            .finish();
        assert_ne!(owner, carrier);
        assert_eq!(
            graph
                .loop_momentum_basis
                .loop_atom::<Atom>(owner, GS.emr_mom, &[], true),
            carrier_momentum,
        );
        let child_soft = FunctionBuilder::new(GS.emr_mom)
            .add_arg(
                GS.uv_momentum_provenance_tag(
                    Atom::num(usize::from(owner) as i64).as_view(),
                    UvMomentumProvenanceRole::DenominatorDerivedSoft,
                    FunctionBuilder::new(GS.emr_mom)
                        .add_arg(usize::from(owner))
                        .finish(),
                ),
            )
            .add_arg(GS.cind(0))
            .finish();
        let expected = FunctionBuilder::new(GS.emr_mom)
            .add_arg(GS.uv_momentum_provenance_tag(
                Atom::num(usize::from(owner) as i64).as_view(),
                UvMomentumProvenanceRole::TaylorFixed,
                carrier_momentum,
            ))
            .add_arg(GS.cind(0))
            .finish()
            / GS.rescale;
        assert_eq!(
            graph.uv_rescaled(
                &filter,
                0,
                &graph.loop_momentum_basis,
                &graph.loop_momentum_basis,
                &child_soft
            ),
            expected,
        );
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
            (erased.collect_factors() - plain.collect_factors()).is_zero(),
            "tag erasure must recover the independent child-sub-LMB Taylor coefficient",
        );

        let mut tags = BTreeMap::new();
        let mut soft_carriers = BTreeSet::new();
        let mut malformed = false;
        let _ = tagged.replace_map(|view, _, _| {
            let AtomView::Fun(momentum) = view else {
                return;
            };
            if momentum.get_symbol() != GS.emr_mom || momentum.get_nargs() < 1 {
                return;
            }
            if let Some((owner, role, hard)) = GS.uv_momentum_provenance_data(momentum.get(0)) {
                if role == UvMomentumProvenanceRole::DenominatorDerivedSoft {
                    soft_carriers.insert(owner);
                    return;
                }
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
        assert_eq!(soft_carriers, BTreeSet::from([EdgeIndex(4)]));
        assert_eq!(
            tags,
            BTreeMap::from([
                ((first, false), first_hard.clone()),
                ((first, true), first_hard.clone()),
                ((second, true), second_hard.clone()),
            ])
        );

        // Read the natural Taylor leaves directly in this ownership oracle.
        // Production retains the same topologies and factorized numerators;
        // reconstruct their complete value without fixing how many leaves or
        // which order the projection chooses. Numerator provenance stays intact,
        // while only denominator kinematics erase their transient hard tags.
        let expected = tagged
            .replace(GS.den(W_.a_, W_.mom_, W_.mass_, W_.prop_))
            .with_map(|matched| {
                GS.den(
                    matched.get(W_.a_).unwrap().to_atom(),
                    GS.erase_uv_momentum_provenance(&matched.get(W_.mom_).unwrap().to_atom()),
                    matched.get(W_.mass_).unwrap().to_atom(),
                    GS.erase_uv_momentum_provenance(&matched.get(W_.prop_).unwrap().to_atom()),
                )
            });
        let reconstructed =
            FourDTerm::from_view(tagged.as_view())?
                .into_iter()
                .fold(Atom::Zero, |sum, term| {
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
        assert!(
            (reconstructed.collect_factors() - expected.collect_factors()).is_zero(),
            "projection must retain the full original and denominator-derived numerator on each physical owner",
        );
        Ok(())
    }

    #[test]
    fn gl24_dod_two_q1_quartic_taylor_keeps_owner_local_energy_families() -> Result<()> {
        test_initialise()?;
        // This is the GL24 skeleton and generation LMB used by the scalar LU
        // acceptance test. The quartic factor is local to e1. The deliberately
        // collected common denominator below stresses ownership while each
        // numerator stays factorized; production keeps Taylor topologies separate.
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
                current.lmb(),
                &integrand,
            )
            .series(GS.rescale, Atom::Zero, 0)?
            .coefficient(Rational::from(0))
            .simplify_metrics()
            .to_dots()
            .normalize_dots();

        let natural_terms = FourDTerm::from_view(expanded.as_view())?;
        // Typed denominator kinematics erase hard-momentum tags. A positive
        // denominator factor in the numerator must instead retain the original
        // wrapper, so e2/e7 keep their own derivative families even when their
        // hard momenta share the e1 carrier in this child LMB.
        let mut tagged_denominators = Vec::<(FourDDenominator, Atom)>::new();
        for matched in expanded.pattern_match(
            &GS.den(W_.a_, W_.mom_, W_.mass_, W_.prop_).to_pattern(),
            None,
            None,
        ) {
            let tagged = GS.den(
                &matched[&W_.a_],
                &matched[&W_.mom_],
                &matched[&W_.mass_],
                &matched[&W_.prop_],
            );
            let denominator = FourDDenominator::from_view(tagged.as_view())?.unwrap();
            if let Some((_, previous)) = tagged_denominators
                .iter()
                .find(|(candidate, _)| candidate == &denominator)
            {
                assert_eq!(previous, &tagged);
            } else {
                tagged_denominators.push((denominator, tagged));
            }
        }
        let mut common_denominators = Vec::new();
        for term in &natural_terms {
            let mut available = common_denominators.clone();
            for denominator in &term.denominators {
                if let Some(index) = available.iter().position(|item| item == denominator) {
                    available.remove(index);
                } else {
                    common_denominators.push(denominator.clone());
                }
            }
        }
        let denominator_product = |denominators: &[FourDDenominator]| {
            denominators
                .iter()
                .fold(Atom::one(), |product, denominator| {
                    product
                        * &tagged_denominators
                            .iter()
                            .find(|(candidate, _)| candidate == denominator)
                            .expect("every typed denominator has its original tagged wrapper")
                            .1
                })
        };
        let common_denominator = denominator_product(&common_denominators);
        let common_numerator = natural_terms.iter().fold(Atom::Zero, |sum, term| {
            let mut missing = common_denominators.clone();
            for denominator in &term.denominators {
                let index = missing.iter().position(|item| item == denominator).unwrap();
                missing.remove(index);
            }
            // Only multiply the existing numerator by missing denominator
            // factors. Certify each contribution by exact factor cancellation,
            // without forming or distributing a common polynomial numerator.
            let contribution = &term.numerator * denominator_product(&missing);
            assert_eq!(
                (&contribution / &common_denominator).collect_factors(),
                (&term.numerator / denominator_product(&term.denominators)).collect_factors(),
            );
            sum + contribution
        });
        let common = FourDTerm {
            numerator: common_numerator,
            denominators: common_denominators,
        };
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
            .automatic_numerator_energy_degree_bounds_in_atoms_excluding_with_min_degree(
                [&common.numerator],
                excluded_edges.iter().copied(),
                1,
            )?;
        assert_eq!(physical_bounds, vec![(1, 6), (2, 4), (7, 4)]);

        let provenance = |atom: &Atom| {
            let mut tags = BTreeSet::new();
            let _ = atom.replace_map(|view, _, _| {
                let AtomView::Fun(momentum) = view else {
                    return;
                };
                if momentum.get_symbol() == GS.emr_mom
                    && momentum.get_nargs() >= 1
                    && let Some((owner, role, _)) = GS.uv_momentum_provenance_data(momentum.get(0))
                    && role != UvMomentumProvenanceRole::DenominatorDerivedSoft
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

        // Read the same natural Taylor leaves before collecting their common
        // denominator. Regrouping by denominator multiplicity preserves each
        // factorized C_i and its soft-shift and mass pieces.
        let analyzer = EnergyPowerAnalyzer::for_physical_emr_edges(owners);
        let mut reduced_numerators = BTreeMap::<[usize; 3], Atom>::new();
        for leaf in natural_terms {
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

        // Exercise the production tensor/Taylor pipeline as well as the
        // diagnostic common form. Only the rational Taylor shell may split;
        // the quartic numerator stays factorized and fixed on e1.
        let given = OwnedForestNode {
            spinney: Spinney::empty(&graph),
            topo_order: 0,
        };
        let settings = UVgenerationSettings::default();
        let (production, _) = t(
            &integrand,
            &UVCtx::new(&graph, &settings),
            &current,
            &given,
            &[],
            0,
        )?;
        let scalar_series_oracle = graph
            .uv_rescaled(
                current.subgraph(),
                graph.n_loops(current.subgraph()),
                current.lmb(),
                current.lmb(),
                &integrand,
            )
            .series(GS.rescale, Atom::Zero, 0)?
            .to_atom()
            .replace(GS.rescale)
            .with(Atom::one())
            .simplify_metrics()
            .to_dots()
            .normalize_dots();
        assert!(
            (production.collect_factors() - scalar_series_oracle.collect_factors()).is_zero(),
            "production T must equal the complete scalar Taylor series, including leading and linear layers",
        );
        let mut production_numerators = BTreeMap::<[usize; 3], Atom>::new();
        for term in FourDTerm::from_view(production.as_view())? {
            assert!(
                !term.numerator.contains_symbol(GS.den),
                "natural Taylor numerators must not contain denominator-clearing factors",
            );
            let multiplicities = owners.map(|owner| {
                term.denominators
                    .iter()
                    .filter(|denominator| denominator.source_edge == owner)
                    .count()
            });
            assert_eq!(
                term.denominators.len(),
                multiplicities.iter().sum::<usize>()
            );
            *production_numerators
                .entry(multiplicities)
                .or_insert(Atom::Zero) += term.numerator;
        }
        let production_bounds = production_numerators
            .iter()
            .map(|(powers, numerator)| {
                Ok((
                    *powers,
                    analyzer.analyze_atom(numerator)?.into_generation_bounds(),
                ))
            })
            .collect::<Result<BTreeMap<_, _>>>()?;
        assert_eq!(
            production_bounds,
            BTreeMap::from([
                ([1, 1, 1], vec![(1, 4)]),
                ([1, 1, 2], vec![(1, 4), (7, 1)]),
                ([1, 1, 3], vec![(1, 4), (7, 2)]),
                ([1, 2, 1], vec![(1, 4), (2, 1)]),
                ([1, 2, 2], vec![(1, 4), (2, 1), (7, 1)]),
                ([1, 3, 1], vec![(1, 4), (2, 2)]),
                ([2, 1, 1], vec![(1, 4)]),
            ]),
            "production T must retain seven natural topologies with fixed e1 rank four and only derivative-local rank increases",
        );

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
            let actual_mass_part = GS
                .erase_uv_momentum_provenance(&(coefficient - without_mass_constants))
                .collect_factors()
                // Normalize numerical coefficients within the mass sum only;
                // momentum factors and products of sums stay factorized.
                .expand_num();
            let expected_mass_part = (&numerator
                * (physical_mass.pow(2) - &expansion_mass_squared))
                .collect_factors()
                .expand_num();
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
    fn dod_one_triangle_keeps_separate_denominator_topologies() -> Result<()> {
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
        // Numerator Taylor layers can share a denominator topology without
        // being regrouped into one additive term.
        owner_multiplicities.dedup();
        assert_eq!(owner_multiplicities, vec![[1, 1, 1], [1, 1, 2], [2, 1, 1]]);
        assert!(
            terms
                .iter()
                .all(|term| !term.numerator.contains_symbol(GS.den)),
            "Taylor terms must not acquire denominator-clearing numerator factors",
        );
        let independent = graph
            .uv_rescaled(
                current.subgraph(),
                graph.n_loops(current.subgraph()),
                current.lmb(),
                current.lmb(),
                &integrand,
            )
            .series(GS.rescale, Atom::Zero, 0)?
            .to_atom()
            .replace(GS.rescale)
            .with(Atom::one())
            .simplify_metrics()
            .to_dots()
            .normalize_dots();
        assert!((expanded.collect_factors() - independent.collect_factors()).is_zero());

        let options = graph.denominator_only_cff_3d_expression_options();
        let mut cache = ExactCffGenerationCache::default();
        let active_denominators = terms
            .iter()
            .map(|term| {
                term.denominators
                    .iter()
                    .filter_map(|denominator| {
                        let is_uv = uv_filter.includes(&graph[&denominator.source_edge].1);
                        denominator
                            .depends_on_loop(&graph, is_uv)
                            .map(|active| active.then(|| denominator.clone()))
                            .transpose()
                    })
                    .collect::<Result<Vec<_>, _>>()
            })
            .collect::<Result<Vec<_>, _>>()?;
        let sources = active_denominators
            .iter()
            .map(|denominators| {
                GraphThreeDSource::from_exact_denominators_in_uv_edges(&graph, denominators, owners)
            })
            .collect::<Result<Vec<_>, _>>()?;
        for (source, term) in sources.iter().zip(&terms) {
            graph.generate_3d_expression_for_4d_term(
                source,
                &options,
                &term.numerator,
                Some(&mut cache),
            )?;
        }
        // Undotted and dotted sources have different occurrence counts;
        // compatible dotted owner relabellings may share a canonical CFF.
        // Regardless of reuse, each request must preserve its complete residue.
        let cutset = crate::graph::cuts::CutSet::empty(graph.n_hedges());
        let mut has_nonzero_residue = false;
        for (denominators, term) in active_denominators.iter().zip(&terms) {
            let mut values = Vec::new();
            for generation_cache in [Some(&mut cache), None] {
                let (cff, _) = graph.clone().cff_from_4d_denominators_in_uv_edges(
                    denominators,
                    owners,
                    &cutset,
                    &options,
                    &term.numerator,
                    generation_cache,
                )?;
                let mut value = Atom::Zero;
                for term in cff.terms.values() {
                    for orientation in &term.orientations {
                        value += &orientation.expression
                            * term.map_exact_source_numerator(&orientation.orientation)?;
                    }
                }
                values.push(value * Atom::num(cff.production_prefactor_factor()));
            }
            has_nonzero_residue |= !values[0].collect_factors().is_zero();
            assert!(
                (values[0].collect_factors() - values[1].collect_factors())
                    .collect_factors()
                    .is_zero(),
                "batching must preserve the complete value of every Taylor term"
            );
        }

        assert!(
            has_nonzero_residue,
            "the Taylor terms must exercise nonzero residues"
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
                    Vec::new(),
                )],
                vec![FourDSector::new(finite, Vec::new(), vec![lmb])],
            ))
        };

        let product = Local4dCts::from_full_product([
            component(local_a.clone(), finite_a.clone(), lmb_a.clone()),
            component(local_b.clone(), finite_b.clone(), lmb_b.clone()),
        ]);
        for sector in product
            .active_sectors()
            .iter()
            .chain(product.recursive_completion())
        {
            let expected_frozen = [(EdgeIndex(0), &finite_a), (EdgeIndex(1), &finite_b)]
                .into_iter()
                .filter_map(|(owner, factor)| {
                    (sector.atom.replace(factor.clone()).with(0) != sector.atom).then_some(owner)
                })
                .collect::<BTreeSet<_>>();
            let frozen_owners = sector
                .frozen_lmbs()
                .iter()
                .flat_map(|lmb| lmb.loop_edges.iter().copied())
                .collect::<BTreeSet<_>>();
            assert_eq!(
                frozen_owners, expected_frozen,
                "only completed physical components remain inert under later Taylor operations"
            );
        }
        let cograph = graph.empty_subgraph::<SuBitGraph>();
        let projected_source = Full4dCts::with_cograph(&product, &graph, &cograph);
        let expected_active = &local_a * &local_b + &finite_a * &local_b + &local_a * &finite_b;
        assert!(
            (projected_source
                .sectors()
                .fold(Atom::Zero, |sum, sector| sum + &sector.atom)
                - expected_active)
                .expand()
                .is_zero(),
            "projecting the local product excludes the fully completed contribution"
        );
        let recursive_source = Full4dCts::from_factorized_local(&product);
        assert!(
            (recursive_source
                .sectors()
                .fold(Atom::Zero, |sum, sector| sum + &sector.atom)
                - product.atom())
            .expand()
            .is_zero(),
            "recursive replay includes every completed contribution exactly once"
        );
        assert_eq!(
            product.atom(),
            &((&local_a + &finite_a) * (&local_b + &finite_b))
        );
        Ok(())
    }

    #[test]
    fn affine_uv_rescaling_preserves_enclosing_chart_and_owner() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph affine_enclosing_taylor_chart {
            edge [num=1 mass=1]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]
            a -> b [id=0]
            a -> b [id=1 lmb_id=0]
            a -> b [id=2 lmb_id=1]
            incoming -> a [id=3]
            b -> outgoing [id=4]
        })?;
        let filter = (0..3)
            .map(|edge| graph.get_edge_subgraph(EdgeIndex(edge)))
            .reduce(|left, right| left.union(&right))
            .unwrap();
        let reference = &graph.loop_momentum_basis;
        let chosen = graph
            .generate_loop_momentum_bases_of(&filter)
            .into_iter()
            .find(|lmb| {
                lmb.loop_edges.contains(&EdgeIndex(0)) && lmb.loop_edges.contains(&EdgeIndex(1))
            })
            .expect("the retained child carrier has a compatible enclosing basis");
        assert!(
            !reference
                .ext_atom::<Atom>(EdgeIndex(0), GS.emr_mom, &[], true)
                .is_zero()
        );
        assert_ne!(chosen.loop_edges, reference.loop_edges);
        for index in [GS.cind(0), GS.cind(1)] {
            let indices = std::slice::from_ref(&index);
            for owner in [EdgeIndex(0), EdgeIndex(1), EdgeIndex(2)] {
                let original = GS.emr_mom(owner, &index);
                let expected = reference.loop_atom(owner, GS.emr_mom, indices, true) / GS.rescale
                    + reference.ext_atom(owner, GS.emr_mom, indices, true);
                for role in [None, Some(false), Some(true)] {
                    let input = role.map_or_else(
                        || original.clone(),
                        |derived| {
                            FunctionBuilder::new(GS.emr_mom)
                                .add_arg(
                                    GS.uv_momentum_provenance_tag(
                                        Atom::num(usize::from(owner) as i64).as_view(),
                                        derived,
                                        FunctionBuilder::new(GS.emr_mom)
                                            .add_arg(usize::from(owner))
                                            .finish()
                                            .as_view(),
                                    ),
                                )
                                .add_arg(index.as_view())
                                .finish()
                        },
                    );
                    let actual = graph.uv_rescaled(&filter, 0, &chosen, reference, &input);
                    let erased = GS.erase_uv_momentum_provenance(&actual).replace_multiple(
                        graph.uv_wrapped_replacement(&filter, reference, indices),
                    );
                    assert!(
                        (erased - &expected).expand().is_zero(),
                        "ordinary and retained hard momenta must implement the same enclosing Taylor operator"
                    );
                    let _ = actual.replace_map(|view, _, _| {
                        if let AtomView::Fun(momentum) = view
                            && momentum.get_symbol() == GS.emr_mom
                            && momentum.get_nargs() > 0
                            && let Some((actual_owner, actual_role, _)) =
                                GS.uv_momentum_provenance_data(momentum.get(0))
                            && actual_role != UvMomentumProvenanceRole::DenominatorDerivedSoft
                        {
                            assert_eq!(actual_owner, owner);
                            assert_eq!(actual_role, role.unwrap_or(false).into());
                        }
                    });
                }
                let soft_input = FunctionBuilder::new(GS.emr_mom)
                    .add_arg(
                        GS.uv_momentum_provenance_tag(
                            Atom::num(usize::from(owner) as i64).as_view(),
                            UvMomentumProvenanceRole::DenominatorDerivedSoft,
                            FunctionBuilder::new(GS.emr_mom)
                                .add_arg(usize::from(owner))
                                .finish()
                                .as_view(),
                        ),
                    )
                    .add_arg(index.as_view())
                    .finish();
                let transported = graph.uv_rescaled(&filter, 0, &chosen, reference, &soft_input);
                let ordinary = graph.uv_rescaled(&filter, 0, &chosen, reference, &original);
                assert_eq!(
                    GS.erase_uv_momentum_provenance(&transported),
                    GS.erase_uv_momentum_provenance(&ordinary),
                    "a child-soft carrier must enter the same raw enclosing chart as ordinary Q before any later graph substitution",
                );
                assert!(
                    (GS.erase_uv_momentum_provenance(&transported)
                        .replace_multiple(
                            graph.uv_wrapped_replacement(&filter, reference, indices)
                        )
                        - &expected)
                        .expand()
                        .is_zero()
                );
                let _ = transported.replace_map(|view, _, _| {
                    if let AtomView::Fun(momentum) = view
                        && momentum.get_symbol() == GS.emr_mom
                        && momentum.get_nargs() > 0
                        && let Some((actual_owner, actual_role, _)) =
                            GS.uv_momentum_provenance_data(momentum.get(0))
                        && actual_role != UvMomentumProvenanceRole::DenominatorDerivedSoft
                    {
                        assert_eq!(actual_owner, owner);
                        assert_eq!(actual_role, UvMomentumProvenanceRole::TaylorFixed,
                            "a child soft carrier's new hard part is not an outer denominator derivative");
                    }
                });
            }
        }
        // A proper child also has paired crown edges. Their literal external
        // coordinates remain soft even when their out-of-domain LMB rows vanish.
        let child_filter = graph
            .get_edge_subgraph(EdgeIndex(1))
            .union(&graph.get_edge_subgraph(EdgeIndex(2)));
        let child_subgraph =
            InternalSubGraph::cleaned_filter_optimist(child_filter.clone(), graph.as_ref());
        let child_lmb = graph.try_compatible_sub_lmb(
            &child_subgraph,
            graph.dummy_stripped_external_flows_of(&child_subgraph),
            reference,
        )?;
        let crown = child_lmb
            .ext_edges
            .iter()
            .copied()
            .find(|edge| {
                graph[edge].1.is_paired()
                    && child_lmb
                        .ext_atom::<Atom>(*edge, GS.emr_mom, &[], true)
                        .is_zero()
            })
            .expect("the bubble has an external paired carrier with a zero local row");
        for index in [GS.cind(0), GS.cind(1)] {
            let payload = FunctionBuilder::new(GS.emr_mom)
                .add_arg(usize::from(crown))
                .finish();
            let soft = FunctionBuilder::new(GS.emr_mom)
                .add_arg(GS.uv_momentum_provenance_tag(
                    Atom::num(usize::from(crown) as i64).as_view(),
                    UvMomentumProvenanceRole::DenominatorDerivedSoft,
                    payload.as_view(),
                ))
                .add_arg(index)
                .finish();
            assert_eq!(
                graph.uv_rescaled(&child_filter, 0, &child_lmb, &child_lmb, &soft),
                soft,
                "an entirely external soft carrier must not acquire a synthetic zero hard factor",
            );
        }
        for owner in [EdgeIndex(1), EdgeIndex(2)] {
            let expected = child_lmb.loop_atom(owner, GS.emr_mom, &[GS.cind(0)], true) / GS.rescale
                + child_lmb.ext_atom(owner, GS.emr_mom, &[GS.cind(0)], true);
            let ordinary = GS.emr_mom(owner, GS.cind(0));
            let actual = graph.uv_rescaled(&child_filter, 0, &child_lmb, &child_lmb, &ordinary);
            assert!(
                (GS.erase_uv_momentum_provenance(&actual) - &expected)
                    .expand()
                    .is_zero()
            );
            for derived in [false, true] {
                let hard = child_lmb.loop_atom::<Atom>(owner, GS.emr_mom, &[], true)
                    + child_lmb.ext_atom::<Atom>(owner, GS.emr_mom, &[], true)
                    - FunctionBuilder::new(GS.emr_mom)
                        .add_arg(usize::from(crown))
                        .finish();
                let tagged = FunctionBuilder::new(GS.emr_mom)
                    .add_arg(GS.uv_momentum_provenance_tag(
                        Atom::num(usize::from(owner) as i64).as_view(),
                        derived,
                        hard.as_view(),
                    ))
                    .add_arg(GS.cind(0))
                    .finish();
                let actual = graph.uv_rescaled(&child_filter, 0, &child_lmb, &child_lmb, &tagged);
                assert!(
                    (GS.erase_uv_momentum_provenance(&actual) - &expected
                        + GS.emr_mom(crown, GS.cind(0)))
                    .expand()
                    .is_zero()
                );
            }
        }
        Ok(())
    }
}
