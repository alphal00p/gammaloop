use std::{
    collections::{BTreeMap, BTreeSet},
    time::{Duration, Instant},
};

use color_eyre::Result;
use eyre::{WrapErr, eyre};
use gammaloop_tracing_filter::debug_instrument;
use idenso::{
    dirac::AGS,
    shorthands::{metric::MetricSimplifier, schoonschip::Schoonschip},
};

#[cfg(test)]
use linnet::half_edge::subgraph::SubGraphLike;
use linnet::half_edge::{
    involution::EdgeIndex,
    subgraph::{Inclusion, InternalSubGraph, SuBitGraph, SubSetLike, SubSetOps},
};
use serde::{Deserialize, Serialize};
use symbolica::{
    atom::{AtomCore, AtomView},
    prelude::*,
};

use crate::utils::symbols::{UvDenominatorClassId, UvMomentumProvenanceRole};
use crate::{
    debug_tags,
    graph::{FourDDenominator, Graph, LMBext, LoopMomentumBasis},
    momentum::sample::LoopIndex,
    numerator::{symbolica_ext::NumeratorAtomExt, ufo::UFO},
    utils::{GS, W_},
    uv::{
        ApproximationType, UltravioletGraph,
        approx::{ForestNodeLike, Rooted, UVCtx, integrated::IntegratedCts},
        marker::{UvMarker, UvOperation},
    },
};
use three_dimensional_reps::MomentumSignature;

#[derive(Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub(crate) struct Local4dRouteSignature {
    pub edge: usize,
    pub loop_signature: String,
    pub external_signature: String,
}

#[derive(Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub(crate) struct Local4dComponentProvenance {
    pub component: String,
    pub scheme: ApproximationType,
    pub dod: i32,
    pub route_loop_edges: Vec<usize>,
    pub route_external_edges: Vec<usize>,
    pub route_signatures: Vec<Local4dRouteSignature>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub(crate) enum Local4dBranch {
    #[serde(rename = "U")]
    U,
    #[serde(rename = "S")]
    S,
    #[serde(rename = "US")]
    US,
    Combined,
}

#[derive(Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub(crate) struct Local4dBranchProvenance {
    pub branch: Local4dBranch,
    /// Size at the raw operator-composition boundary, before the single forest
    /// sign, marker, and metric simplification are applied. Tensor/gamma
    /// normalization remains deferred to integration or numerical evaluation.
    pub byte_size: usize,
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) struct FourDSector {
    pub(crate) atom: Atom,
    /// A nonzero S branch occurred in this contribution's ancestry. Retained
    /// soft coordinates require canonical routing and restrict supported parent
    /// schemes, even after later cancellations. Parent operators act once on `atom`;
    /// ordinary-scheme comparisons are constructed independently in tests.
    has_soft_ancestry: bool,
    /// Loop generators fixed by completed Taylor operations. Containing
    /// operations extend this basis without redefining its coordinates.
    /// Boundary momenta remain routed dependencies in `atom`, not generators.
    inherited_loop_edges: Vec<EdgeIndex>,
    provenance_components: Vec<Local4dComponentProvenance>,
    provenance_branches: Vec<Local4dBranchProvenance>,
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
    has_soft_ancestry: bool,
    inherited_loop_edges: Vec<EdgeIndex>,
    provenance_components: Vec<Local4dComponentProvenance>,
    provenance_branches: Vec<Local4dBranchProvenance>,
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

/// Representation-neutral algebra of a completed sector. Raw Taylor sectors
/// remain authoritative for subsequent outer Taylor operations.
#[derive(Clone, Debug)]
pub(crate) struct CanonicalUvSector {
    pub(crate) terms: Vec<CanonicalUvTerm>,
    pub(crate) classes: Vec<CanonicalUvDenominatorClass>,
    pub(crate) active_components: Vec<(SuBitGraph, SuBitGraph, LoopMomentumBasis)>,
    pub(crate) frozen_lmbs: Vec<LoopMomentumBasis>,
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) struct CanonicalUvDenominatorClass {
    pub(crate) id: UvDenominatorClassId,
    pub(crate) component: usize,
    pub(crate) momentum: Atom,
    pub(crate) signature: MomentumSignature,
    pub(crate) mass_squared: Atom,
    /// The exact denominator polynomial also distinguishes any prescription
    /// carried by the input. The current UV source uses the shared Feynman
    /// prescription; canonicalization introduces no new prescription convention.
    pub(crate) full_expr: Atom,
    pub(crate) members: BTreeSet<EdgeIndex>,
}

#[derive(Clone, Debug)]
pub(crate) struct CanonicalUvTerm {
    pub(crate) numerator: Atom,
    /// Factorized raw numerator witness, used only for the exact projection
    /// certificate and physical-parent diagnostics. Its provenance does not
    /// participate in canonical algebra or occurrence allocation.
    pub(crate) source_numerator: Atom,
    pub(crate) powers: BTreeMap<UvDenominatorClassId, usize>,
    /// Physical incidence is independent of algebraic class identity. Equal
    /// buckets retain one witness, never concatenate their denominator lists.
    pub(crate) source_witness: Vec<FourDDenominator>,
    /// Aligned with `source_witness`; sign means raw momentum = sign * class
    /// momentum. None denotes a residual or component-loop-independent factor.
    pub(crate) source_classes: Vec<Option<(UvDenominatorClassId, i32)>>,
}

type CanonicalUvAlgebraKey = (
    BTreeMap<UvDenominatorClassId, usize>,
    Vec<(EdgeIndex, Atom, Atom, Atom)>,
);

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct UvDenominatorClassKey {
    component: usize,
    signature: MomentumSignature,
    mass_squared: Atom,
    full_expr: Atom,
}

impl CanonicalUvDenominatorClass {
    pub(crate) fn accounted_bytes(&self) -> usize {
        std::mem::size_of::<Self>()
            + [&self.momentum, &self.mass_squared, &self.full_expr]
                .into_iter()
                .map(|atom| atom.as_view().get_byte_size())
                .sum::<usize>()
            + (self.signature.loop_signature.capacity()
                + self.signature.external_signature.capacity())
                * std::mem::size_of::<i32>()
            + self.members.len() * 48
    }

    pub(crate) fn momentum(signature: &MomentumSignature, lmb: &LoopMomentumBasis) -> Atom {
        signature
            .loop_signature
            .iter()
            .zip(&lmb.loop_edges)
            .chain(signature.external_signature.iter().zip(&lmb.ext_edges))
            .filter(|(coefficient, _)| **coefficient != 0)
            .fold(Atom::Zero, |sum, (coefficient, edge)| {
                sum + Atom::num(*coefficient) * function!(GS.emr_mom, usize::from(*edge))
            })
    }

    pub(crate) fn momentum_with_indices(&self, indices: &[Atom]) -> Atom {
        GS.indexed_momentum(&self.momentum, indices)
    }

    fn polynomial_in_class(&self, lmb: &LoopMomentumBasis) -> Result<Atom> {
        let (axis, coefficient) = self
            .signature
            .loop_signature
            .iter()
            .enumerate()
            .find(|(_, coefficient)| **coefficient != 0)
            .ok_or_else(|| eyre!("UV class {} has no active loop coordinate", self.id.0))?;
        let carrier = usize::from(lmb.loop_edges[crate::momentum::sample::LoopIndex(axis)]);
        let remainder = &self.momentum - Atom::num(*coefficient) * function!(GS.emr_mom, carrier);
        // This invertible coordinate change acts only on the denominator
        // polynomial. It recovers one formal H even when D(H) was expanded in
        // physical loop coordinates, without splitting a numerator factor.
        let polynomial = self
            .full_expr
            .replace_map(|view, _, output| {
                let AtomView::Fun(momentum) = view else {
                    return;
                };
                if momentum.get_symbol() != GS.emr_mom
                    || momentum.get_nargs() == 0
                    || usize::try_from(momentum.get(0)).ok() != Some(carrier)
                {
                    return;
                }
                let indices = momentum
                    .iter()
                    .skip(1)
                    .map(|index| index.to_owned())
                    .collect::<Vec<_>>();
                let mut reference =
                    FunctionBuilder::new(GS.emr_mom).add_arg(GS.uv_class_ref(self.id));
                for index in &indices {
                    reference = reference.add_arg(index);
                }
                **output = (reference.finish() - GS.indexed_momentum(&remainder, &indices))
                    / Atom::num(*coefficient);
            })
            .normalize_dots()
            .expand();
        let mut residual = false;
        let collapsed = polynomial
            .replace_map(|view, _, output| {
                let AtomView::Fun(momentum) = view else {
                    return;
                };
                if momentum.get_symbol() != GS.emr_mom || momentum.get_nargs() == 0 {
                    return;
                }
                if GS.uv_class_data(momentum.get(0)) != Some(self.id) {
                    residual = true;
                    return;
                }
                let indices = momentum
                    .iter()
                    .skip(1)
                    .map(|index| index.to_owned())
                    .collect::<Vec<_>>();
                **output = self.momentum_with_indices(&indices);
            })
            .normalize_dots()
            .expand();
        if residual || collapsed != self.full_expr {
            return Err(eyre!(
                "positive denominator for UV class {} is not a polynomial solely in its certified momentum channel",
                self.id.0
            ));
        }
        Ok(polynomial)
    }
}

impl FourDSector {
    pub(crate) fn accounted_bytes(&self) -> usize {
        std::mem::size_of::<Self>()
            + self.atom.as_view().get_byte_size()
            + self.inherited_loop_edges.capacity() * std::mem::size_of::<EdgeIndex>()
            + self.provenance_components.capacity()
                * std::mem::size_of::<Local4dComponentProvenance>()
            + self
                .provenance_components
                .iter()
                .map(|component| {
                    component.component.capacity()
                        + (component.route_loop_edges.capacity()
                            + component.route_external_edges.capacity())
                            * std::mem::size_of::<usize>()
                        + component.route_signatures.capacity()
                            * std::mem::size_of::<Local4dRouteSignature>()
                        + component
                            .route_signatures
                            .iter()
                            .map(|route| {
                                route.loop_signature.capacity()
                                    + route.external_signature.capacity()
                            })
                            .sum::<usize>()
                })
                .sum::<usize>()
            + self.provenance_branches.capacity() * std::mem::size_of::<Local4dBranchProvenance>()
            + self.active_components.capacity()
                * std::mem::size_of::<(SuBitGraph, SuBitGraph, LoopMomentumBasis)>()
            + self
                .active_components
                .iter()
                .map(|(owners, scope, lmb)| {
                    owners.size().div_ceil(8)
                        + scope.size().div_ceil(8)
                        + 64
                        + lmb.accounted_bytes()
                })
                .sum::<usize>()
            + self.frozen_lmbs.capacity() * std::mem::size_of::<LoopMomentumBasis>()
            + self
                .frozen_lmbs
                .iter()
                .map(LoopMomentumBasis::accounted_bytes)
                .sum::<usize>()
    }

    pub(crate) fn new(
        atom: Atom,
        active_components: Vec<(SuBitGraph, SuBitGraph, LoopMomentumBasis)>,
        frozen_lmbs: Vec<LoopMomentumBasis>,
    ) -> Self {
        let mut inherited_loop_edges = active_components
            .iter()
            .flat_map(|(_, _, lmb)| lmb.loop_edges.iter().copied())
            .chain(
                frozen_lmbs
                    .iter()
                    .flat_map(|lmb| lmb.loop_edges.iter().copied()),
            )
            .collect::<Vec<_>>();
        inherited_loop_edges.sort();
        inherited_loop_edges.dedup();
        Self {
            has_soft_ancestry: false,
            atom,
            active_components,
            frozen_lmbs,
            inherited_loop_edges,
            provenance_components: Vec::new(),
            provenance_branches: Vec::new(),
        }
    }

    #[cfg(test)]
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

    /// Normalize only a projection-owned copy, with class membership resolved
    /// separately in each rational term and independent component frame.
    pub(crate) fn canonical_projection(&self, graph: &Graph) -> Result<CanonicalUvSector> {
        let terms = self.physical_terms()?;
        let mut classes = BTreeMap::<UvDenominatorClassKey, BTreeSet<EdgeIndex>>::new();
        let mut term_keys = Vec::with_capacity(terms.len());
        for term in &terms {
            let mut keys = Vec::with_capacity(term.denominators.len());
            for denominator in &term.denominators {
                let components = self
                    .active_components
                    .iter()
                    .enumerate()
                    .filter(|(_, (owners, _, _))| {
                        owners.includes(&graph[&denominator.source_edge].1)
                    })
                    .map(|(component, _)| component)
                    .collect::<Vec<_>>();
                let component = match components.as_slice() {
                    [] => {
                        keys.push(None);
                        continue;
                    }
                    [component] => *component,
                    _ => {
                        return Err(eyre!(
                            "4D denominator owner {} belongs to overlapping active Taylor components {:?}",
                            usize::from(denominator.source_edge),
                            components
                        ));
                    }
                };
                let signature = denominator
                    .momentum_signature_in_lmb(&self.active_components[component].2, true)?;
                if signature
                    .loop_signature
                    .iter()
                    .all(|coefficient| *coefficient == 0)
                {
                    keys.push(None);
                    continue;
                }
                let (signature, sign) = signature.canonical_up_to_sign();
                let key = UvDenominatorClassKey {
                    component,
                    signature,
                    mass_squared: denominator.mass_squared.clone(),
                    full_expr: denominator
                        .polynomial_in_lmb(&self.active_components[component].2)?,
                };
                classes
                    .entry(key.clone())
                    .or_default()
                    .insert(denominator.source_edge);
                keys.push(Some((key, sign)));
            }
            term_keys.push(keys);
        }
        let ids = classes
            .keys()
            .cloned()
            .enumerate()
            .map(|(index, key)| (key, UvDenominatorClassId(index)))
            .collect::<BTreeMap<_, _>>();
        let classes = classes
            .into_iter()
            .map(|(key, members)| CanonicalUvDenominatorClass {
                id: ids[&key],
                component: key.component,
                momentum: CanonicalUvDenominatorClass::momentum(
                    &key.signature,
                    &self.active_components[key.component].2,
                ),
                signature: key.signature,
                mass_squared: key.mass_squared,
                full_expr: key.full_expr,
                members,
            })
            .collect();
        let mut canonical = CanonicalUvSector {
            terms: Vec::new(),
            classes,
            active_components: self.active_components.clone(),
            frozen_lmbs: self.frozen_lmbs.clone(),
        };
        for (term, keys) in terms.into_iter().zip(term_keys) {
            let source_classes = keys
                .into_iter()
                .map(|key| key.map(|(key, sign)| (ids[&key], sign)))
                .collect::<Vec<_>>();
            let mut powers = BTreeMap::new();
            for (class, _) in source_classes.iter().flatten() {
                *powers.entry(*class).or_default() += 1;
            }
            let mut term = CanonicalUvTerm {
                source_numerator: term.numerator.clone(),
                numerator: term.numerator,
                powers,
                source_witness: term.denominators,
                source_classes,
            };
            let neutral = canonical.neutral_numerator(&term.numerator, graph)?;
            term.numerator = canonical.normalize_numerator(&term.numerator, &term, graph)?;
            let projected_neutral = canonical.neutral_numerator(&term.numerator, graph)?;
            // A signed linear momentum can acquire a factored minus sign at
            // this boundary. Normalize each operand independently only when
            // their structural forms differ. Distribute numeric coefficients
            // only; products and powers of graph numerators remain factorized.
            if projected_neutral != neutral
                && projected_neutral.expand_num().collect_factors()
                    != neutral.expand_num().collect_factors()
            {
                crate::debug_tags!(#generation, #uv, #local, #four_d, #trace;
                    file.source_numerator = %neutral,
                    file.projected_numerator = %projected_neutral,
                    "Canonical UV numerator certificate failed"
                );
                return Err(eyre!(
                    "canonical UV numerator does not reproduce its source in the component frame"
                ));
            }
            canonical.terms.push(term);
        }
        canonical.merge_terms();
        for term in &canonical.terms {
            let projected_neutral = canonical.neutral_numerator(&term.numerator, graph)?;
            let neutral = canonical.neutral_numerator(&term.source_numerator, graph)?;
            if projected_neutral != neutral
                && projected_neutral.expand_num().collect_factors()
                    != neutral.expand_num().collect_factors()
            {
                return Err(eyre!(
                    "merged canonical UV numerator does not reproduce its source bucket"
                ));
            }
        }
        Ok(canonical)
    }
}

impl CanonicalUvSector {
    pub(crate) fn accounted_bytes(&self) -> usize {
        std::mem::size_of::<Self>()
            + self.terms.capacity() * std::mem::size_of::<CanonicalUvTerm>()
            + self
                .terms
                .iter()
                .map(|term| {
                    term.numerator.as_view().get_byte_size()
                        + term.source_numerator.as_view().get_byte_size()
                        + term.powers.len() * 64
                        + term.source_witness.capacity() * std::mem::size_of::<FourDDenominator>()
                        + term
                            .source_witness
                            .iter()
                            .map(FourDDenominator::accounted_bytes)
                            .sum::<usize>()
                        + term.source_classes.capacity()
                            * std::mem::size_of::<Option<(UvDenominatorClassId, i32)>>()
                })
                .sum::<usize>()
            + self.classes.capacity() * std::mem::size_of::<CanonicalUvDenominatorClass>()
            + self
                .classes
                .iter()
                .map(CanonicalUvDenominatorClass::accounted_bytes)
                .sum::<usize>()
            + self.active_components.capacity()
                * std::mem::size_of::<(SuBitGraph, SuBitGraph, LoopMomentumBasis)>()
            + self
                .active_components
                .iter()
                .map(|(owners, scope, lmb)| {
                    owners.size().div_ceil(8)
                        + scope.size().div_ceil(8)
                        + 64
                        + lmb.accounted_bytes()
                })
                .sum::<usize>()
            + self.frozen_lmbs.capacity() * std::mem::size_of::<LoopMomentumBasis>()
            + self
                .frozen_lmbs
                .iter()
                .map(LoopMomentumBasis::accounted_bytes)
                .sum::<usize>()
    }

    pub(crate) fn class(&self, id: UvDenominatorClassId) -> &CanonicalUvDenominatorClass {
        &self.classes[id.0]
    }

    /// Exact diagonal certificate in the retained component frames. Positive
    /// typed denominators unwrap only in this proof copy, never in the source.
    pub(crate) fn neutral_numerator(&self, atom: &Atom, graph: &Graph) -> Result<Atom> {
        let mut error = None;
        let result = atom.replace_map(|view, _, output| {
            let AtomView::Fun(momentum) = view else {
                return;
            };
            if momentum.get_symbol() != GS.emr_mom || momentum.get_nargs() < 2 {
                return;
            }
            let indices = momentum
                .iter()
                .skip(1)
                .map(|index| index.to_owned())
                .collect::<Vec<_>>();
            if let Some(class) = GS.uv_class_data(momentum.get(0)) {
                if let Some(class) = self.classes.get(class.0) {
                    **output = class.momentum_with_indices(&indices);
                } else {
                    error = Some(eyre!(
                        "canonical UV numerator refers to an unknown class {}",
                        class.0
                    ));
                }
                return;
            }
            let Some((owner, role, hard)) = GS.uv_momentum_provenance_data(momentum.get(0)) else {
                return;
            };
            if role == UvMomentumProvenanceRole::DenominatorDerivedSoft {
                return;
            }
            let components = self
                .active_components
                .iter()
                .filter(|(owners, _, _)| owners.includes(&graph[&owner].1))
                .collect::<Vec<_>>();
            let [(_, _, lmb)] = components.as_slice() else {
                return;
            };
            let descriptor = FourDDenominator {
                source_edge: owner,
                momentum: hard,
                mass_squared: Atom::Zero,
                full_expr: Atom::Zero,
            };
            match descriptor.momentum_signature_in_lmb(lmb, true) {
                Ok(signature) => {
                    **output = GS.indexed_momentum(
                        &CanonicalUvDenominatorClass::momentum(&signature, lmb),
                        &indices,
                    )
                }
                Err(err) => error = Some(err.into()),
            }
        });
        if let Some(error) = error {
            return Err(error);
        }
        let result = result.replace_map(|view, _, output| {
            let AtomView::Fun(denominator) = view else {
                return;
            };
            if denominator.get_symbol() != GS.den || denominator.get_nargs() != 4 {
                return;
            }
            let component = if let Some(class) = GS.uv_class_data(denominator.get(0)) {
                self.classes.get(class.0).map(|class| class.component)
            } else if let Ok(owner) = usize::try_from(denominator.get(0)).map(EdgeIndex) {
                self.active_components
                    .iter()
                    .position(|(owners, _, _)| owners.includes(&graph[&owner].1))
            } else {
                None
            };
            let full_expr = denominator.get(3).to_owned();
            **output = if let Some(component) = component {
                let descriptor = FourDDenominator {
                    source_edge: EdgeIndex(0),
                    momentum: Atom::Zero,
                    mass_squared: Atom::Zero,
                    full_expr,
                };
                match descriptor.polynomial_in_lmb(&self.active_components[component].2) {
                    Ok(full_expr) => full_expr,
                    Err(err) => {
                        error = Some(err);
                        return;
                    }
                }
            } else {
                full_expr
            };
        });
        // Signed frame changes can leave opposite additive factors, including
        // bases of integer powers. Choose one literal base without expanding
        // numerator products or powers. Power bases are handled at the power
        // node so fractional powers retain their branch-sensitive original form.
        let mut result = result;
        loop {
            let normalized = result.replace_map_bottom_up(|view, context, output| {
                if matches!(view, AtomView::Add(_))
                    && context.parent_type != Some(symbolica::atom::AtomType::Pow)
                {
                    let base = view.to_owned().expand_num();
                    let opposite = (-&base).expand_num();
                    **output = if opposite < base { -opposite } else { base };
                    return;
                }
                let AtomView::Pow(power) = view else {
                    return;
                };
                let Ok(exponent) = i64::try_from(power.get_exp()) else {
                    return;
                };
                if !matches!(power.get_base(), AtomView::Add(_)) {
                    return;
                }
                let base = power.get_base().to_owned().expand_num();
                let opposite = (-&base).expand_num();
                if opposite < base {
                    **output =
                        Atom::num(if exponent % 2 == 0 { 1 } else { -1 }) * opposite.pow(exponent);
                }
            });
            // An extracted sign can expose numerical content to its parent on
            // the next pass. Stop at structural equality, retaining products
            // and compressed powers throughout normalization.
            if normalized == result {
                break;
            }
            result = normalized;
        }
        match error {
            Some(error) => Err(error),
            None => Ok(result),
        }
    }

    fn normalize_numerator(
        &self,
        atom: &Atom,
        term: &CanonicalUvTerm,
        graph: &Graph,
    ) -> Result<Atom> {
        let normalize = |view: AtomView<'_>| -> Result<Option<Atom>> {
            let AtomView::Fun(function) = view else {
                return Ok(None);
            };
            if function.get_symbol() == GS.den && function.get_nargs() == 4 {
                let Ok(owner) = usize::try_from(function.get(0)).map(EdgeIndex) else {
                    return Ok(None);
                };
                let descriptor = FourDDenominator {
                    source_edge: owner,
                    momentum: GS.erase_uv_momentum_provenance(&function.get(1).to_owned()),
                    mass_squared: function.get(2).to_owned(),
                    full_expr: GS.erase_uv_momentum_provenance(&function.get(3).to_owned()),
                };
                for class in term.powers.keys().map(|class| self.class(*class)) {
                    let (owners, _, lmb) = &self.active_components[class.component];
                    if !owners.includes(&graph[&owner].1)
                        || class.mass_squared != descriptor.mass_squared
                    {
                        continue;
                    }
                    let signature = descriptor
                        .momentum_signature_in_lmb(lmb, true)?
                        .canonical_up_to_sign()
                        .0;
                    if signature == class.signature
                        && descriptor.polynomial_in_lmb(lmb)? == class.full_expr
                    {
                        let full_expr = class.polynomial_in_class(lmb)?;
                        return Ok(Some(GS.den(
                            GS.uv_class_ref(class.id),
                            &class.momentum,
                            &class.mass_squared,
                            full_expr,
                        )));
                    }
                }
                // An unmatched positive block keeps its physical owner and
                // polynomial together, including all existing provenance.
                return Ok(Some(view.to_owned()));
            }
            if function.get_symbol() != GS.emr_mom || function.get_nargs() < 2 {
                return Ok(None);
            }
            let (owner, role, hard) = if let Some(provenance) =
                GS.uv_momentum_provenance_data(function.get(0))
            {
                provenance
            } else {
                if matches!(function.get(0), AtomView::Fun(tag) if tag.get_symbol() == GS.uv_momentum_provenance)
                {
                    return Err(eyre!(
                        "unfinished or malformed UV momentum provenance in completed projection: {}",
                        function.get(0).to_owned()
                    ));
                }
                return Ok(None);
            };
            if !matches!(
                role,
                UvMomentumProvenanceRole::TaylorFixed
                    | UvMomentumProvenanceRole::DenominatorDerived
            ) {
                return Ok(None);
            }
            let components = self
                .active_components
                .iter()
                .enumerate()
                .filter(|(_, (owners, _, _))| owners.includes(&graph[&owner].1))
                .map(|(component, _)| component)
                .collect::<Vec<_>>();
            let [component] = components.as_slice() else {
                return Ok(None);
            };
            let descriptor = FourDDenominator {
                source_edge: owner,
                momentum: hard,
                mass_squared: Atom::Zero,
                full_expr: Atom::Zero,
            };
            let signature = descriptor
                .momentum_signature_in_lmb(&self.active_components[*component].2, true)?;
            let (signature, sign) = signature.canonical_up_to_sign();
            let candidates = term
                .powers
                .keys()
                .copied()
                .filter(|class| {
                    self.class(*class).component == *component
                        && self.class(*class).signature == signature
                })
                .collect::<BTreeSet<_>>();
            let matching = if candidates.len() <= 1 {
                candidates
            } else {
                // Equal channels with different masses remain different pools.
                // Prefer the original surviving owner;
                // an ambiguous pinched factor stays a fixed affine carrier.
                term.source_witness
                    .iter()
                    .zip(&term.source_classes)
                    .filter(|(denominator, _)| denominator.source_edge == owner)
                    .filter_map(|(_, class)| class.map(|(id, _)| id))
                    .filter(|class| candidates.contains(class))
                    .collect::<BTreeSet<_>>()
            };
            let reference = if matching.len() == 1 {
                GS.uv_class_ref(*matching.first().unwrap())
            } else {
                // A pinched hard carrier has no pole family. Keep its exact
                // affine payload fixed for the source-coordinate mapper.
                let hard = CanonicalUvDenominatorClass::momentum(
                    &signature,
                    &self.active_components[*component].2,
                );
                GS.uv_momentum_provenance_tag(
                    usize::from(owner),
                    UvMomentumProvenanceRole::PhysicalSourceFixed,
                    &hard,
                )
            };
            let mut momentum = FunctionBuilder::new(GS.emr_mom).add_arg(reference);
            for index in function.iter().skip(1) {
                momentum = momentum.add_arg(index);
            }
            Ok(Some(Atom::num(sign) * momentum.finish()))
        };
        let mut error = None;
        let normalized = atom.replace_map(|view, _, output| {
            if error.is_none() {
                match normalize(view) {
                    Ok(Some(value)) => **output = value,
                    Ok(None) => {}
                    Err(err) => error = Some(err),
                }
            }
        });
        match error {
            Some(error) => Err(error),
            None => Ok(normalized),
        }
    }

    fn merge_terms(&mut self) {
        let mut buckets = BTreeMap::new();
        for term in std::mem::take(&mut self.terms) {
            let key = term.algebra_key();
            match buckets.entry(key) {
                std::collections::btree_map::Entry::Vacant(entry) => {
                    entry.insert(term);
                }
                std::collections::btree_map::Entry::Occupied(mut entry) => {
                    entry.get_mut().numerator += term.numerator;
                    entry.get_mut().source_numerator += term.source_numerator;
                }
            }
        }
        self.terms = buckets
            .into_values()
            .filter(|term| !term.numerator.is_zero())
            .collect();
    }
}

impl CanonicalUvTerm {
    fn algebra_key(&self) -> CanonicalUvAlgebraKey {
        let mut residual = self
            .source_witness
            .iter()
            .zip(&self.source_classes)
            .filter(|(_, class)| class.is_none())
            .map(|(denominator, _)| {
                (
                    denominator.source_edge,
                    denominator.momentum.clone(),
                    denominator.mass_squared.clone(),
                    denominator.full_expr.clone(),
                )
            })
            .collect::<Vec<_>>();
        residual.sort();
        (self.powers.clone(), residual)
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
        // Preserve contribution history before pruning zero projection payloads.
        let sectors = active.iter().chain(&recursive_completion);
        let has_soft_ancestry = sectors.clone().any(|sector| sector.has_soft_ancestry);
        let mut inherited_loop_edges = Vec::new();
        let mut provenance_components = Vec::new();
        let mut provenance_branches = Vec::new();
        for sector in sectors {
            for edge in &sector.inherited_loop_edges {
                if !inherited_loop_edges.contains(edge) {
                    inherited_loop_edges.push(*edge);
                }
            }
            provenance_components.extend(sector.provenance_components.iter().cloned());
            provenance_branches.extend(sector.provenance_branches.iter().cloned());
        }
        active.retain(|sector| !sector.atom.is_zero());
        recursive_completion.retain(|sector| !sector.atom.is_zero());
        Self {
            active,
            recursive_completion,
            atom,
            has_soft_ancestry,
            inherited_loop_edges,
            provenance_components,
            provenance_branches,
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
    #[cfg(test)]
    pub(crate) fn atom(&self) -> &Atom {
        &self.0.atom
    }

    #[cfg(test)]
    fn root() -> Self {
        Self(FourDSectors::active_atom(Atom::one()))
    }

    /// Pole-part subtraction inserts only completed poles; MUV and soft/IR keep
    /// the completed local counterterm with its integrated finite contribution.
    pub(crate) fn recursion_input(
        local: &Local4dCts,
        integrated: &IntegratedCts,
        scheme: ApproximationType,
        is_root: bool,
        lmb: &LoopMomentumBasis,
    ) -> Result<Self> {
        // Once the current hard integrations are complete, their coefficient
        // has no remaining dependence on those coordinates. Both 3D routes use
        // this common nominal frame for the normalized localizing addback;
        // the analytic integration itself consumes the exact active frames.
        let completed = |atom| FourDSector {
            has_soft_ancestry: local.0.has_soft_ancestry,
            ..FourDSector::new(atom, Vec::new(), vec![lmb.clone()])
        };
        let mut projected = match scheme {
            ApproximationType::MUV | ApproximationType::IR => {
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
            ApproximationType::OS => unimplemented!(
                "local on-shell counterterms are deferred until local counterterms can be derived from the 4D expanded representation"
            ),
            scheme => Err(eyre!("No recursive counterterm projection for {scheme}")),
        }?;
        projected.0.has_soft_ancestry |= local.0.has_soft_ancestry;
        for edge in &local.0.inherited_loop_edges {
            if !projected.0.inherited_loop_edges.contains(edge) {
                projected.0.inherited_loop_edges.push(*edge);
            }
        }
        Ok(projected)
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
            .map(|sector| FourDSector {
                atom: &sector.atom * &cograph,
                ..sector.clone()
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

impl Local4dCts {
    pub(crate) fn atom(&self) -> &Atom {
        &self.0.atom
    }

    pub(crate) fn provenance_components(&self) -> &[Local4dComponentProvenance] {
        &self.0.provenance_components
    }

    pub(crate) fn provenance_branches(&self) -> &[Local4dBranchProvenance] {
        &self.0.provenance_branches
    }

    pub(crate) fn from_full_product(factors: impl IntoIterator<Item = Full4dCts>) -> Self {
        let mut products = vec![(FourDSector::new(Atom::one(), Vec::new(), Vec::new()), false)];
        let mut atom = Atom::one();
        let mut has_soft_ancestry = false;
        let mut inherited_loop_edges = Vec::new();
        for factor in factors {
            atom *= &factor.0.atom;
            has_soft_ancestry |= factor.0.has_soft_ancestry;
            for edge in &factor.0.inherited_loop_edges {
                if !inherited_loop_edges.contains(edge) {
                    inherited_loop_edges.push(*edge);
                }
            }
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
                            let mut product = FourDSector::new(
                                left.atom.clone() * right.atom,
                                active_components,
                                frozen_lmbs,
                            );
                            product.has_soft_ancestry =
                                left.has_soft_ancestry || right.has_soft_ancestry;
                            for edge in left
                                .inherited_loop_edges
                                .iter()
                                .chain(&right.inherited_loop_edges)
                            {
                                if !product.inherited_loop_edges.contains(edge) {
                                    product.inherited_loop_edges.push(*edge);
                                }
                            }
                            product.provenance_components = left.provenance_components.clone();
                            product
                                .provenance_components
                                .extend(right.provenance_components);
                            (product, left_active || right_active)
                        })
                })
                .collect();
        }
        let (active, recursive_completion) = products
            .into_iter()
            .partition::<Vec<_>, _>(|(_, is_active)| *is_active);
        let combined_byte_size = atom.as_view().get_byte_size();
        let mut sectors = FourDSectors::with_atom(
            active.into_iter().map(|(sector, _)| sector).collect(),
            recursive_completion
                .into_iter()
                .map(|(sector, _)| sector)
                .collect(),
            atom,
        );
        sectors.has_soft_ancestry |= has_soft_ancestry;
        for edge in inherited_loop_edges {
            if !sectors.inherited_loop_edges.contains(&edge) {
                sectors.inherited_loop_edges.push(edge);
            }
        }
        sectors.provenance_branches = vec![Local4dBranchProvenance {
            branch: Local4dBranch::Combined,
            byte_size: combined_byte_size,
        }];
        Self(sectors)
    }

    pub(crate) fn active_sectors(&self) -> &[FourDSector] {
        &self.0.active
    }

    /// Combine only identical integration bindings, retaining raw recursive
    /// sectors for the next Taylor operation. First occurrence order is stable.
    pub(crate) fn projection_sectors(&self) -> Vec<FourDSector> {
        let mut sectors: Vec<FourDSector> = Vec::new();
        let mut indices = std::collections::HashMap::<_, usize>::new();
        for sector in &self.0.active {
            let key = (&sector.active_components, &sector.frozen_lmbs);
            if let Some(index) = indices.get(&key).copied() {
                sectors[index].atom += &sector.atom;
            } else {
                indices.insert(key, sectors.len());
                sectors.push(sector.clone());
            }
        }
        sectors.retain(|sector| !sector.atom.is_zero());
        sectors
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

    fn group(terms: impl IntoIterator<Item = Self>) -> Vec<Self> {
        let mut buckets = BTreeMap::new();
        for mut term in terms {
            term.denominators.sort_by_cached_key(|denominator| {
                (
                    denominator.source_edge,
                    denominator.momentum.clone(),
                    denominator.mass_squared.clone(),
                    denominator.full_expr.clone(),
                )
            });
            let key = term
                .denominators
                .iter()
                .map(|denominator| {
                    (
                        denominator.source_edge,
                        denominator.momentum.clone(),
                        denominator.mass_squared.clone(),
                        denominator.full_expr.clone(),
                    )
                })
                .collect::<Vec<_>>();
            match buckets.entry(key) {
                std::collections::btree_map::Entry::Vacant(entry) => {
                    entry.insert(term);
                }
                std::collections::btree_map::Entry::Occupied(mut entry) => {
                    entry.get_mut().numerator += term.numerator;
                }
            }
        }
        buckets
            .into_values()
            .filter(|term| !term.numerator.is_zero())
            .collect()
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
                Ok(Self::group(terms))
            }
            AtomView::Mul(mul) => {
                let mut products = vec![Self::numerator(Atom::one())];
                for factor in mul.iter() {
                    let factor_terms = Self::from_view(factor)?;
                    products = Self::group(products.into_iter().flat_map(|left| {
                        factor_terms
                            .iter()
                            .cloned()
                            .map(move |right| Self::product(left.clone(), right))
                    }));
                }
                Ok(products)
            }
            AtomView::Pow(power) if power.get_base().contains_symbol(GS.den) => {
                let (base, exponent) = power.get_base_exp();
                if !matches!(base, AtomView::Add(_) | AtomView::Mul(_)) {
                    return Ok(vec![Self::from_factorized_term(view)?]);
                }
                let Ok(exponent) = usize::try_from(exponent) else {
                    return Ok(vec![Self::from_factorized_term(view)?]);
                };
                let factors = Self::from_view(base)?;
                if factors.iter().all(|term| term.denominators.is_empty()) {
                    return Ok(vec![Self::numerator(view.to_owned())]);
                }
                // Multiply rational shells only. Positive denominator blocks
                // and denominator-free numerator powers remain opaque above.
                let mut products = vec![Self::numerator(Atom::one())];
                for _ in 0..exponent {
                    products = Self::group(products.into_iter().flat_map(|left| {
                        factors
                            .iter()
                            .cloned()
                            .map(move |right| Self::product(left.clone(), right))
                    }));
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
    pub(crate) fn polynomial_in_lmb(&self, lmb: &LoopMomentumBasis) -> Result<Atom> {
        let mut error = None;
        let normalized = self.full_expr.replace_map(|view, _, output| {
            let AtomView::Fun(momentum) = view else {
                return;
            };
            if momentum.get_symbol() != GS.emr_mom || momentum.get_nargs() < 1 {
                return;
            }
            let Ok(owner) = usize::try_from(momentum.get(0)).map(EdgeIndex) else {
                return;
            };
            let descriptor = FourDDenominator {
                source_edge: owner,
                momentum: function!(GS.emr_mom, usize::from(owner)),
                mass_squared: Atom::Zero,
                full_expr: Atom::Zero,
            };
            match descriptor.momentum_signature_in_lmb(lmb, true) {
                Ok(signature) => {
                    let indices = momentum
                        .iter()
                        .skip(1)
                        .map(|index| index.to_owned())
                        .collect::<Vec<_>>();
                    **output = GS.indexed_momentum(
                        &CanonicalUvDenominatorClass::momentum(&signature, lmb),
                        &indices,
                    );
                }
                Err(err) => error = Some(err.into()),
            }
        });
        if let Some(error) = error {
            return Err(error);
        }
        // Only the small denominator polynomial is expanded. Graph numerator
        // sums and products never pass through this denominator certificate.
        Ok(normalized.normalize_dots().expand())
    }

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
        self.uv_rescaled_with_loops(expansion_subgraph, n_loops, lmb, reference_lmb, false, atom)
    }

    fn uv_rescaled_with_loops(
        &self,
        expansion_subgraph: &SuBitGraph,
        n_loops: usize,
        lmb: &LoopMomentumBasis,
        reference_lmb: &LoopMomentumBasis,
        soft_expansion: bool,
        atom: &Atom,
    ) -> Atom {
        let mut scaled_edges = self
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
        if soft_expansion {
            // Numerator-only boundary carriers also belong to the external
            // Taylor jet; unrelated padded graph externals stay fixed.
            scaled_edges.extend(
                self.dummy_less_full_crown(expansion_subgraph)
                    .included_iter()
                    .map(|hedge| self[&hedge]),
            );
            scaled_edges.sort();
            scaled_edges.dedup();
        }
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
                    routed = if lmb.ext_edges.contains(&owner) {
                        function!(GS.emr_mom, usize::from(owner))
                    } else {
                        lmb.loop_atom::<Atom>(owner, GS.emr_mom, &[], true)
                            + lmb.ext_atom::<Atom>(owner, GS.emr_mom, &[], true)
                    };
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
                    if lmb.ext_edges.contains(&owner) {
                        function!(GS.emr_mom, usize::from(owner))
                    } else {
                        lmb.loop_atom::<Atom>(owner, GS.emr_mom, &[], true)
                            + lmb.ext_atom::<Atom>(owner, GS.emr_mom, &[], true)
                    },
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
            **output = if soft_expansion {
                hard_component + soft * GS.rescale
            } else {
                hard_component * &inverse_rescale + soft
            };
        });

        if soft_expansion {
            // S retains physical masses. Only the denominator polynomial is
            // differentiated; its momentum metadata names the resulting hard
            // pole at zero component-external momentum.
            return atomarg
                .replace(GS.den(W_.a_, W_.mom_, W_.mass_, W_.prop_))
                .with_map(|matched| {
                    GS.den(
                        matched.get(W_.a_).unwrap().to_atom(),
                        matched
                            .get(W_.mom_)
                            .unwrap()
                            .to_atom()
                            .replace(GS.rescale)
                            .with(Atom::Zero),
                        matched.get(W_.mass_).unwrap().to_atom(),
                        matched.get(W_.prop_).unwrap().to_atom(),
                    )
                });
        }

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
#[cfg(test)]
fn grow<S: super::ForestNodeLike>(
    integrand: &Full4dCts,
    ctx: &UVCtx<'_>,
    current: &S,
    given: &S,
) -> Result<Atom> {
    let factor = grow_factor(ctx, current, given)?;
    let integrand = (&integrand.0.atom * factor).simplify_metrics();

    debug_tags!(#uv, #integrated, #algebra, #start; log.integrand = integrand, reduced = %current.reduced_subgraph(given).string_label());

    Ok(integrand)
}

fn grow_factor<S: super::ForestNodeLike>(ctx: &UVCtx<'_>, current: &S, given: &S) -> Result<Atom> {
    let reduced = current.reduced_subgraph(given);
    let graph = ctx.graph;

    let mut t_arg = ctx
        .graph
        .numerator(&reduced, given.subgraph())
        .to_d_dim(GS.dim)
        .color_simplify()
        .get_single_atom()
        .unwrap();

    t_arg /= graph.denominator(&reduced, |_| 1);
    Ok(t_arg)
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
) -> Result<(Atom, LoopMomentumBasis, usize)> {
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
        .series_preserving_factors(GS.rescale, Atom::Zero.as_view(), 0, &[])
        .map_err(|error| {
            eyre!(
                "local 4D Taylor expansion failed for {}: {error}",
                current.subgraph().string_label()
            )
        })?;
    debug_tags!(#uv,#integrated, #series;log.res = series, "Series expanded");

    let evalutated = series.replace(GS.rescale).with(Atom::num(1));
    debug_tags!(#uv,#integrated,#series;log.res = evalutated, "Evaluated at t = 1");

    let byte_size = evalutated.as_view().get_byte_size();
    Ok((finalize(evalutated), lmb.clone(), byte_size))
}

#[debug_instrument(
        current = %current.log_display(),
        given = %given.log_display(),
    )]
fn compatible_lmb<S: ForestNodeLike>(
    ctx: &UVCtx<'_>,
    current: &S,
    given: &S,
    integrand: &Full4dCts,
) -> Result<LoopMomentumBasis> {
    let graph = ctx.graph;

    // Preserve the validated MUV/PolePart routing whenever no completed soft
    // ancestry is present. Besides retaining every child loop coordinate,
    // this requires exactly the child's number of basis representatives to lie
    // physically inside the contracted child.  Vakint then sees all remaining
    // parent coordinates in the reduced co-graph.  A parent acting on H must
    // instead preserve the expansion coordinates recorded by that completed
    // soft atom and therefore continues through the canonical path below.
    if matches!(
        current.renormalization_scheme(),
        ApproximationType::MUV | ApproximationType::PolePart
    ) && !integrand.0.has_soft_ancestry
    {
        let is_ordinary_compatible = |candidate: &LoopMomentumBasis| {
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
        if given.subgraph().is_empty() || is_ordinary_compatible(current.lmb()) {
            return Ok(current.lmb().clone());
        }
        return graph
            .generate_loop_momentum_bases_of(current.subgraph())
            .into_iter()
            .find(is_ordinary_compatible)
            .ok_or_else(|| eyre!("no loop momentum basis compatible with nested UV subgraph"));
    }

    let affine_required_edges: Vec<_> = integrand
        .0
        .inherited_loop_edges
        .iter()
        .copied()
        .filter(|edge| current.subgraph().includes(&graph[edge].1))
        .filter(|edge| {
            graph.loop_momentum_basis.edge_signatures[*edge]
                .external
                .iter()
                .any(|sign| sign.is_sign())
        })
        .collect();
    if !affine_required_edges.is_empty() {
        return Err(eyre!(
            "component {} cannot preserve inherited loop coordinates {:?}: their canonical routes carry graph-external momentum, so redefining them would be an affine loop-momentum shift",
            current.subgraph().string_label(),
            affine_required_edges
        ));
    }

    // Keep the loop coordinates actually used by the completed child as a
    // subset of the outer basis so its Taylor representation is unchanged.
    // A child's boundary momenta are expressions in the enclosing coordinates;
    // they must not constrain the parent to promote them to loop generators.
    let is_homogeneous = |candidate: &LoopMomentumBasis| {
        candidate.loop_edges.iter().all(|edge| {
            graph.loop_momentum_basis.edge_signatures[*edge]
                .external
                .iter()
                .all(|sign| sign.is_zero())
        })
    };
    let is_compatible = |candidate: &LoopMomentumBasis| {
        integrand
            .0
            .inherited_loop_edges
            .iter()
            .filter(|edge| current.subgraph().includes(&graph[*edge].1))
            .all(|edge| candidate.loop_edges.contains(edge))
    };
    // A loop-coordinate representative may lie on an edge of the contracted
    // child even when that coordinate parametrizes the surviving cograph in
    // the full canonical chart.  The completed child records every coordinate
    // on which its atom actually depends, so physical edge membership in
    // `given` is neither necessary nor sufficient as an additional chart
    // constraint.  Requiring such a count made Appendix B.1 spuriously switch
    // the same outer operation from [e2,e3] to [e2,e4].
    // A forest's nominal `current.lmb()` is inherited from its enumeration
    // path and can choose a different edge representative for the same
    // component. Start from the sub-basis induced by the graph's canonical
    // full route instead, extending the actual child loop generators whenever
    // this sub-basis does not already contain them.
    let external: SuBitGraph = graph.external_filter();
    let full_internal = graph.full_filter().subtract(&external);
    if current.subgraph() == &full_internal
        && graph.n_loops(&full_internal) == graph.loop_momentum_basis.loop_edges.len()
    {
        if is_compatible(&graph.loop_momentum_basis) && is_homogeneous(&graph.loop_momentum_basis) {
            return Ok(graph.loop_momentum_basis.clone());
        }
    } else if let Ok(canonical) = graph.try_compatible_sub_lmb(
        current.subgraph(),
        graph.dummy_less_full_crown(current.subgraph()),
        &graph.loop_momentum_basis,
    ) && is_compatible(&canonical)
        && is_homogeneous(&canonical)
    {
        return Ok(canonical);
    }

    // Restrict spanning-forest bases to homogeneous changes that retain the
    // inherited generators, then pick a
    // stable edge-lexicographic minimum rather than depending on traversal
    // order when several bases preserve the completed child coordinates.
    graph
        .generate_loop_momentum_bases_of(current.subgraph())
        .into_iter()
        .filter(|candidate| is_compatible(candidate) && is_homogeneous(candidate))
        .min_by_key(|candidate| {
            candidate
                .loop_edges
                .iter()
                .map(|edge| usize::from(*edge))
                .collect::<Vec<_>>()
        })
        .ok_or_else(|| {
            eyre!(
                "no homogeneous loop momentum basis for component {} preserves nested expansion coordinates {:?}",
                current.subgraph().string_label(),
                integrand.0.inherited_loop_edges
            )
        })
}

#[cfg(test)]
fn component_external_edges<S: ForestNodeLike>(
    ctx: &UVCtx<'_>,
    current: &S,
    lmb: &LoopMomentumBasis,
) -> Vec<EdgeIndex> {
    let component_edges = ctx.graph.paired_edges(current.subgraph());
    let component_carriers: Vec<_> = ctx
        .graph
        .dummy_less_full_crown(current.subgraph())
        .included_iter()
        .map(|hedge| ctx.graph[&hedge])
        .collect();
    lmb.ext_edges
        .iter_enumerated()
        .filter(|(external, edge)| {
            component_carriers.contains(edge)
                || component_edges.iter().any(|component_edge| {
                    lmb.edge_signatures[*component_edge].external[*external].is_sign()
                })
        })
        .map(|(_, edge)| *edge)
        .collect()
}

#[cfg(test)]
fn independent_component_external_edges<S: ForestNodeLike>(
    ctx: &UVCtx<'_>,
    current: &S,
    lmb: &LoopMomentumBasis,
    completed: &Atom,
) -> Result<Vec<EdgeIndex>> {
    let component_edges = ctx.graph.paired_edges(current.subgraph());
    let component_carriers: Vec<_> = ctx
        .graph
        .dummy_less_full_crown(current.subgraph())
        .included_iter()
        .map(|hedge| ctx.graph[&hedge])
        .collect();
    if component_carriers.is_empty() {
        return Ok(Vec::new());
    }
    let contains = |edge| {
        completed
            .replace(GS.emr_mom(edge, W_.x___))
            .match_iter()
            .next()
            .is_some()
    };
    // Preserve only external LMB coordinates which survive in the completed
    // child atom.  The internal signatures identify propagator dependence;
    // checking crown membership as well retains a numerator-only coordinate
    // (for example on a tadpole) without inventing a padded cograph carrier.
    let active = lmb
        .ext_edges
        .iter_enumerated()
        .filter(|(external, edge)| {
            contains(**edge)
                && (component_carriers.contains(edge)
                    || component_edges.iter().any(|component_edge| {
                        lmb.edge_signatures[*component_edge].external[*external].is_sign()
                    }))
        })
        .map(|(_, edge)| *edge)
        .collect::<Vec<_>>();
    let unmatched = component_carriers
        .iter()
        .copied()
        .filter(|edge| contains(*edge) && !lmb.ext_edges.iter().any(|candidate| candidate == edge))
        .collect::<Vec<_>>();
    if !unmatched.is_empty() {
        return Err(eyre!(
            "component {} retains component-external momenta {unmatched:?} which are not coordinates of its canonical LMB, so nested expansion provenance is ambiguous",
            current.subgraph().string_label(),
        ));
    }
    Ok(active)
}

#[cfg(test)]
fn internal_expansion_carriers(
    graph: &Graph,
    carriers: impl IntoIterator<Item = EdgeIndex>,
) -> Vec<EdgeIndex> {
    let mut preferred = Vec::new();
    for carrier in carriers {
        let pair = &graph[&carrier].1;
        if pair.is_paired() && !graph.tree_edges.includes(pair) {
            // The two boundary representatives of a self-energy can carry
            // the same homogeneous momentum with opposite orientation. Use
            // the representative already present in the graph's canonical
            // LMB so pure and nested parent terms share one ambient chart.
            // Equality up to sign excludes an affine external-momentum shift.
            let carrier = graph
                .loop_momentum_basis
                .loop_edges
                .iter()
                .copied()
                .find(|candidate| {
                    graph.loop_momentum_basis.edge_signatures[*candidate]
                        .equality_up_to_sign(&graph.loop_momentum_basis.edge_signatures[carrier])
                })
                .unwrap_or(carrier);
            if !preferred.contains(&carrier) {
                preferred.push(carrier);
            }
        }
    }
    preferred
}

fn t_raw<S: ForestNodeLike>(
    integrand: &Atom,
    ctx: &UVCtx<'_>,
    current: &S,
    _given: &S,
    lmb: &LoopMomentumBasis,
) -> Result<Atom> {
    let graph = ctx.graph;
    let active_loop_edges = lmb
        .loop_edges
        .iter()
        .copied()
        .filter(|edge| current.subgraph().includes(&graph[edge].1))
        .collect::<Vec<_>>();
    let expected_loops = graph.n_loops(current.subgraph());
    if active_loop_edges.len() != expected_loops {
        return Err(eyre!(
            "the canonical route {:?} exposes {} loop carriers for component {}, expected {}",
            lmb.loop_edges,
            active_loop_edges.len(),
            current.subgraph().string_label(),
            expected_loops,
        ));
    }

    // A parent always uses the full component scope on the completed child.
    let rescaled = graph.uv_rescaled(
        current.subgraph(),
        active_loop_edges.len(),
        lmb,
        lmb,
        integrand,
    );
    debug_tags!(#uv,#integrated,#rescaled;
        log.res = rescaled,
        n_loops = active_loop_edges.len(),
        "Rescaled expanded"
    );

    let series = rescaled
        .series_preserving_factors(GS.rescale, Atom::Zero.as_view(), 0, &[])
        .wrap_err_with(|| {
            format!(
                "failed to construct the local UV series for component {}",
                current.subgraph().string_label()
            )
        })?;
    debug_tags!(#uv,#integrated, #series;log.res = series, "Series expanded");

    let evaluated = series.replace(GS.rescale).with(Atom::num(1));
    debug_tags!(#uv,#integrated,#series;log.res = evaluated, "Evaluated at t = 1");

    Ok(evaluated)
}

fn tilde_t_raw<S: ForestNodeLike>(
    integrand: &Atom,
    ctx: &UVCtx<'_>,
    current: &S,
    _given: &S,
    lmb: &LoopMomentumBasis,
) -> Result<Atom> {
    // The soft Taylor operator acts on the component's external momenta only;
    // physical masses and the MUV scales remain untouched.
    // Canonical component bases are padded with unrelated graph externals.
    // Scale only independent slots that flow through this component and all
    // of its dependent/independent crown carriers, retaining source ownership.
    let rescaled =
        ctx.graph
            .uv_rescaled_with_loops(current.subgraph(), 0, lmb, lmb, true, integrand);
    tilde_taylor(rescaled, current.dod())
}

fn tilde_taylor(rescaled: Atom, dod: i32) -> Result<Atom> {
    // For a logarithmic component, \widetilde T_{d-1}=\widetilde T_{-1}=0.
    if dod == 0 {
        return Ok(Atom::Zero);
    }

    let series = rescaled
        .series_preserving_factors(GS.rescale, Atom::Zero.as_view(), i64::from(dod - 1), &[])
        .map_err(|error| eyre!("soft Taylor expansion failed: {error}"))?;
    Ok(series.replace(GS.rescale).with(Atom::one()))
}

#[cfg(test)]
fn hat_t_raw<S: ForestNodeLike>(
    integrand: &Atom,
    ctx: &UVCtx<'_>,
    current: &S,
    given: &S,
    lmb: &LoopMomentumBasis,
) -> Result<Atom> {
    Ok(hat_t_raw_with_provenance(integrand, ctx, current, given, lmb)?.0)
}

fn hat_t_raw_with_provenance<S: ForestNodeLike>(
    completed: &Atom,
    ctx: &UVCtx<'_>,
    current: &S,
    given: &S,
    lmb: &LoopMomentumBasis,
) -> Result<(Atom, Vec<Local4dBranchProvenance>)> {
    // Apply U once to the completed atom. Soft ancestry is explicit metadata,
    // so routing never needs to expand tensor networks for a comparison.
    let t = t_raw(completed, ctx, current, given, lmb)?;
    let t_byte_size = t.as_view().get_byte_size();
    debug_tags!(#generation, #uv, #local, #inspect;
        stage = "local_4d_soft_ordinary_branch",
        log.expr = t,
        byte_size = t_byte_size,
        component = %current.subgraph().string_label(),
        "Constructed the ordinary U branch of the local soft operator"
    );
    let tilde_t = tilde_t_raw(completed, ctx, current, given, lmb)?;
    if tilde_t.is_zero() {
        return Ok((
            t,
            vec![
                Local4dBranchProvenance {
                    branch: Local4dBranch::U,
                    byte_size: t_byte_size,
                },
                Local4dBranchProvenance {
                    branch: Local4dBranch::Combined,
                    byte_size: t_byte_size,
                },
            ],
        ));
    }
    let soft_byte_size = tilde_t.as_view().get_byte_size();
    debug_tags!(#generation, #uv, #local, #inspect;
        stage = "local_4d_soft_branch",
        log.expr = tilde_t,
        byte_size = soft_byte_size,
        component = %current.subgraph().string_label(),
        "Constructed the physical-mass S branch of the local soft operator"
    );
    let t_tilde_t = t_raw(&tilde_t, ctx, current, given, lmb)?;
    let overlap_byte_size = t_tilde_t.as_view().get_byte_size();
    debug_tags!(#generation, #uv, #local, #inspect;
        stage = "local_4d_soft_overlap_branch",
        log.expr = t_tilde_t,
        byte_size = overlap_byte_size,
        component = %current.subgraph().string_label(),
        "Constructed U acting on the completed S branch"
    );
    let composite = t + tilde_t - t_tilde_t;
    let combined_byte_size = composite.as_view().get_byte_size();
    debug_tags!(#generation, #uv, #local, #inspect;
        stage = "local_4d_soft_composite",
        log.expr = composite,
        byte_size = combined_byte_size,
        component = %current.subgraph().string_label(),
        "Combined the local soft operator as U + S - US"
    );
    Ok((
        composite,
        vec![
            Local4dBranchProvenance {
                branch: Local4dBranch::U,
                byte_size: t_byte_size,
            },
            Local4dBranchProvenance {
                branch: Local4dBranch::S,
                byte_size: soft_byte_size,
            },
            Local4dBranchProvenance {
                branch: Local4dBranch::US,
                byte_size: overlap_byte_size,
            },
            Local4dBranchProvenance {
                branch: Local4dBranch::Combined,
                byte_size: combined_byte_size,
            },
        ],
    ))
}

fn finalize(raw: Atom) -> Atom {
    // Keep the local Taylor numerator factorized through exact CFF projection.
    // Analytic integration owns its Dirac algebra in integrated::simplify;
    // numerical evaluation contracts the retained tensors after residue mapping.
    // Keep each Taylor term's propagator powers. Collecting factors across the
    // sum clears denominators and manufactures higher-rank numerator factors.
    raw.simplify_metrics()
}

fn component_provenance<S: ForestNodeLike>(
    ctx: &UVCtx<'_>,
    current: &S,
    lmb: &LoopMomentumBasis,
) -> Local4dComponentProvenance {
    Local4dComponentProvenance {
        component: current.subgraph().string_label(),
        scheme: current.renormalization_scheme(),
        dod: current.dod(),
        route_loop_edges: lmb
            .loop_edges
            .iter()
            .map(|edge| usize::from(*edge))
            .collect(),
        route_external_edges: lmb
            .ext_edges
            .iter()
            .map(|edge| usize::from(*edge))
            .collect(),
        route_signatures: ctx
            .graph
            .paired_edges(current.subgraph())
            .into_iter()
            .map(|edge| {
                let signature = &lmb.edge_signatures[edge];
                Local4dRouteSignature {
                    edge: usize::from(edge),
                    loop_signature: signature.internal.to_string(),
                    external_signature: signature.external.to_string(),
                }
            })
            .collect(),
    }
}

pub(crate) fn uv_limit<S: ForestNodeLike, M: ForestNodeLike>(
    integrand: &Full4dCts,
    ctx: &UVCtx<'_>,
    current: &S,
    given: &S,
    marker_current: &M,
    marker_given: &M,
) -> Result<Local4dCts> {
    let started = Instant::now();
    let mut taylor_time = Duration::ZERO;
    let scheme = current.renormalization_scheme();
    if scheme == ApproximationType::OS {
        unimplemented!(
            "local on-shell counterterms are deferred until local counterterms can be derived from the 4D expanded representation"
        );
    }
    if scheme == ApproximationType::MUV && current.dod() > 0 && integrand.0.has_soft_ancestry {
        return Err(eyre!(
            "positive-degree MUV component {} cannot contain a local soft scheme refinement in phase 1; its reduced-cograph subtraction requires the deferred finite scheme-change/integrated policy",
            current.subgraph().string_label()
        ));
    }
    if scheme == ApproximationType::PolePart && integrand.0.has_soft_ancestry {
        return Err(eyre!(
            "PolePart cannot be nested with a local soft refinement in the same wood before an integrated scheme-change policy exists"
        ));
    }
    match scheme {
        ApproximationType::MUV | ApproximationType::PolePart | ApproximationType::IR => {
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
                let co_graph = grow_factor(ctx, current, given)?;
                let grown = (&sector.atom * &co_graph).simplify_metrics();
                let soft_route = scheme == ApproximationType::IR || sector.has_soft_ancestry;
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
                let canonical_component_lmb = if given.subgraph().is_empty() || soft_route {
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
                if !soft_route && retained_loop_edges.len() != expected_retained_loops {
                    return Err(eyre!(
                        "local 4D Taylor sector retains {} prefix-plus-quotient carriers, expected {expected_retained_loops}",
                        retained_loop_edges.len()
                    ));
                }
                let taylor_started = Instant::now();
                let (result, coordinate_lmb, provenance_branches) = if soft_route {
                    let input = Full4dCts(FourDSectors::new(vec![sector.clone()], Vec::new()));
                    let lmb = compatible_lmb(ctx, current, given, &input)?;
                    // The nested forest identity is K_parent(1-K_child): the
                    // parent acts with full component grading on each completed
                    // typed sector. U and H are linear across these sectors.
                    let (raw, provenance_branches) = if scheme == ApproximationType::IR {
                        hat_t_raw_with_provenance(&grown, ctx, current, given, &lmb)?
                    } else {
                        let raw = t_raw(&grown, ctx, current, given, &lmb)?;
                        let byte_size = raw.as_view().get_byte_size();
                        (
                            raw,
                            vec![
                                Local4dBranchProvenance {
                                    branch: Local4dBranch::U,
                                    byte_size,
                                },
                                Local4dBranchProvenance {
                                    branch: Local4dBranch::Combined,
                                    byte_size,
                                },
                            ],
                        )
                    };
                    (finalize(raw), lmb, provenance_branches)
                } else {
                    let (result, lmb, byte_size) = t(
                        &grown,
                        ctx,
                        current,
                        given,
                        &retained_loop_edges,
                        expected_prefix_loops,
                    )?;
                    (
                        result,
                        lmb,
                        vec![
                            Local4dBranchProvenance {
                                branch: Local4dBranch::U,
                                byte_size,
                            },
                            Local4dBranchProvenance {
                                branch: Local4dBranch::Combined,
                                byte_size,
                            },
                        ],
                    )
                };
                taylor_time += taylor_started.elapsed();
                // Exact residues must use the very coordinates in which T
                // produced their hard momenta. Rebuilding a canonical quotient
                // LMB here can choose an equivalent but differently spelled
                // carrier and thereby reintroduce a physical crown shift.
                // Ordinary T preselects a graphic quotient, so its contraction
                // removes prefix directions without inventing external ones.
                // H instead preserves the completed child's actual coordinates:
                // a surviving parent coordinate can lie on a child-owned edge.
                // Demote only the consumed child coordinates in that chart;
                // their retained shifts are fixed during the quotient contour.
                let component_lmb = if given.subgraph().is_empty() {
                    coordinate_lmb.clone()
                } else if soft_route {
                    let mut component_lmb = coordinate_lmb.clone();
                    for index in (0..component_lmb.loop_edges.len()).rev() {
                        let index = LoopIndex(index);
                        if retained_loop_edges.contains(&component_lmb.loop_edges[index]) {
                            component_lmb.put_loop_to_ext(index);
                        }
                    }
                    ctx.graph
                        .canonicalize_lmb_external_order(&mut component_lmb);
                    component_lmb
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
                let mut projected = FourDSector::new(
                    marker.apply(
                        UvOperation::Approx,
                        marker_current.subgraph(),
                        marker_given.subgraph(),
                        &-result,
                    ),
                    active_components,
                    sector.frozen_lmbs.clone(),
                );
                // Normalize the selected result once. Descendants act on this
                // completed signed sector; ordinary controls are built in tests.
                projected.has_soft_ancestry = sector.has_soft_ancestry
                    || provenance_branches
                        .iter()
                        .any(|branch| branch.branch == Local4dBranch::S);
                projected.inherited_loop_edges = sector.inherited_loop_edges.clone();
                for edge in coordinate_lmb.loop_edges.iter().copied() {
                    if !projected.inherited_loop_edges.contains(&edge) {
                        projected.inherited_loop_edges.push(edge);
                    }
                }
                projected.provenance_components =
                    vec![component_provenance(ctx, current, &coordinate_lmb)];
                projected.provenance_branches = provenance_branches;
                Ok(projected)
            };
            let sectors = integrand
                .sectors()
                .map(project)
                .collect::<Result<Vec<_>>>()?;
            // T is linear. Summing the individually framed sectors preserves the
            // aggregate compatibility atom without inventing one coordinate LMB
            // for a disconnected product.
            let mut projected = FourDSectors::new(sectors, Vec::new());
            projected.has_soft_ancestry |= integrand.0.has_soft_ancestry;
            for edge in &integrand.0.inherited_loop_edges {
                if !projected.inherited_loop_edges.contains(edge) {
                    projected.inherited_loop_edges.push(*edge);
                }
            }
            let local = Local4dCts(projected);
            let elapsed = started.elapsed();
            debug_tags!(#generation, #uv, #local, #four_d, #profile;
                stage = "local_4d_construction",
                graph = %ctx.graph.name,
                current = %current.log_display(),
                given = %given.log_display(),
                elapsed_ms = elapsed.as_secs_f64() * 1000.0,
                taylor_ms = taylor_time.as_secs_f64() * 1000.0,
                construction_overhead_ms = elapsed.saturating_sub(taylor_time).as_secs_f64() * 1000.0,
                input_sectors = integrand.0.active.len() + integrand.0.recursive_completion.len(),
                output_sectors = local.0.active.len() + local.0.recursive_completion.len(),
                "Constructed factorized local-4D Taylor sectors"
            );
            Ok(local)
        }
        // OS remains deferred; see the current-boundaries section in
        // docs/architecture/uv-renormalization.typ.
        atype => Err(eyre!("Not yet implemented {:?}", atype)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        dot,
        graph::{GraphThreeDSource, parse::IntoGraph},
        initialisation::test_initialise,
        numerator::{aind::Aind, energy_degree::EnergyPowerAnalyzer},
        uv::approx::projected_4d::Local4dProjectionContext,
        uv::{Spinney, UVgenerationSettings, hedge_poset::OwnedForestNode, uv_graph::UVE},
    };
    use gammaloop_tracing_filter::LogMessage;
    use idenso::{
        color::ColorSimplifier,
        dirac::GammaSimplifier,
        representations::{Bispinor, ColorAdjoint},
    };
    use linnet::half_edge::{
        involution::Hedge,
        subgraph::{InternalSubGraph, ModifySubSet, SubSetOps},
    };
    use spenso::{
        network::{
            parsing::{AtomStructureExt, StrictTensorFilter},
            tags::SPENSO_TAG,
        },
        shadowing::{IntoAtom, symbolica_utils::LogPrint},
        structure::representation::{LibraryRep, Minkowski, RepName},
    };
    use std::collections::{BTreeMap, BTreeSet};
    use symbolica::{domains::rational::Rational, function, symbol};

    struct TestNode {
        subgraph: SuBitGraph,
        lmb: LoopMomentumBasis,
        dod: i32,
        scheme: ApproximationType,
    }

    impl LogMessage for TestNode {
        fn log_display(&self) -> String {
            "local-4d test node".to_owned()
        }
    }

    impl ForestNodeLike for TestNode {
        fn subgraph(&self) -> &SuBitGraph {
            &self.subgraph
        }

        fn lmb(&self) -> &LoopMomentumBasis {
            &self.lmb
        }

        fn dod(&self) -> i32 {
            self.dod
        }

        fn renormalization_scheme(&self) -> ApproximationType {
            self.scheme
        }

        fn topo_order(&self) -> usize {
            0
        }

        fn reduced_subgraph(&self, given: &Self) -> SuBitGraph {
            self.subgraph.subtract(&given.subgraph)
        }
    }

    fn two_point_node(graph: &Graph, scheme: ApproximationType) -> TestNode {
        let external: SuBitGraph = graph.external_filter();
        let subgraph = graph.full_filter().subtract(&external);
        TestNode {
            lmb: graph.lmb_of(&subgraph),
            subgraph,
            dod: 1,
            scheme,
        }
    }

    fn root_node(graph: &Graph) -> TestNode {
        TestNode {
            subgraph: graph.empty_subgraph(),
            lmb: graph.empty_lmb(),
            dod: 0,
            scheme: ApproximationType::MUV,
        }
    }

    fn scalar_two_point_graph() -> Graph {
        dot!(
            digraph scalar_self_energy {
                edge [particle="H" num=1];
                node [num=1];
                ext [style=invis];
                ext -> A:0 [id=0];
                B:1 -> ext [id=1];
                A -> B [id=2];
                A -> B [id=3];
            }
        )
        .unwrap()
    }

    fn nested_hard_degree(
        graph: &Graph,
        parent_subgraph: &SuBitGraph,
        atom: &Atom,
        route: &LoopMomentumBasis,
        hard_edges: &[EdgeIndex],
        label: &str,
    ) -> Rational {
        let lambda_symbol = symbol!("local_4d_nested_hard_lambda");
        let lambda = Atom::var(lambda_symbol);
        let routed = GS
            .erase_uv_momentum_provenance(atom)
            .replace(function!(GS.ct_marker, W_.a_))
            .with(Atom::one())
            .replace_multiple(graph.uv_wrapped_replacement(parent_subgraph, route, &[W_.x___]));
        // Callers pass already-finalized scalar forest atoms.  Re-running the
        // tensor/gamma finalizer on their large sum is both redundant and can
        // dominate the independent Laurent check.
        let mut scaled = routed
            .replace(function!(GS.den, W_.a_, W_.mom_, W_.mass_, W_.prop_))
            .with(Atom::var(W_.prop_))
            .replace(GS.dim)
            .with(Atom::num(4));
        let mut dot_images: Vec<((Atom, Atom), Atom)> = Vec::new();
        scaled = scaled.replace_map(|view, _, output| {
            let AtomView::Fun(dot) = view else {
                return;
            };
            if dot.get_symbol() != SPENSO_TAG.dot || dot.get_nargs() != 2 {
                return;
            }
            let mut arguments = dot.iter();
            let left = arguments.next().unwrap().into_atom();
            let right = arguments.next().unwrap().into_atom();
            let image = dot_images
                .iter()
                .find(|((known_left, known_right), _)| {
                    (known_left == &left && known_right == &right)
                        || (known_left == &right && known_right == &left)
                })
                .map(|(_, image)| image.clone())
                .unwrap_or_else(|| {
                    // Generic exact integers keep this an independent hard-ray
                    // oracle without making Symbolica normalize a large
                    // multivariate rational function.  Equal scalar products
                    // still receive one stable value.
                    let image = Atom::num(101 + dot_images.len() as i64);
                    dot_images.push(((left.clone(), right.clone()), image.clone()));
                    image
                });
            let hard_power = hard_edges
                .iter()
                .map(|edge| {
                    usize::from(paper_has_pattern(&left, GS.emr_mom(*edge, W_.x___)))
                        + usize::from(paper_has_pattern(&right, GS.emr_mom(*edge, W_.x___)))
                })
                .sum::<usize>();
            **output = image / lambda.pow(hard_power as i64);
        });
        assert!(
            scaled
                .replace(function!(SPENSO_TAG.dot, W_.a_, W_.b_))
                .match_iter()
                .next()
                .is_none(),
            "{label}: an unrecognized routed dot product escaped scalarization"
        );
        for edge in hard_edges {
            scaled = scaled
                .replace(GS.emr_mom(*edge, W_.x___))
                .with(GS.emr_mom(*edge, W_.x___) / &lambda);
        }
        let mut component_images: Vec<(Atom, Atom)> = Vec::new();
        scaled = scaled.replace_map(|view, _, output| {
            let AtomView::Fun(component) = view else {
                return;
            };
            if component.get_symbol() != GS.emr_mom {
                return;
            }
            let key = view.into_atom();
            let image = component_images
                .iter()
                .find(|(known, _)| known == &key)
                .map(|(_, image)| image.clone())
                .unwrap_or_else(|| {
                    let image = Atom::num(211 + component_images.len() as i64);
                    component_images.push((key, image.clone()));
                    image
                });
            **output = image;
        });
        for (_, edge, _) in graph.iter_edges_of(parent_subgraph) {
            let mass = graph[edge].mass_atom();
            if !mass.is_zero() {
                scaled = scaled
                    .replace(mass)
                    .with(Atom::num(307 + usize::from(edge) as i64));
            }
        }
        scaled = scaled
            .replace(GS.m_uv_expansion)
            .with(Atom::num(401))
            .replace(GS.m_uv_vacuum)
            .with(Atom::num(409));
        scaled *= lambda.pow(-4 * hard_edges.len() as i64);
        let scaled = scaled.together().cancel();
        assert!(!scaled.is_zero(), "{label} is identically zero");
        scaled
            .series(lambda_symbol, Atom::Zero, 4)
            .unwrap_or_else(|error| panic!("{label}: hard series failed: {error}"))
            .get_trailing_exponent()
    }

    fn padded_two_point_graph() -> Graph {
        dot!(
            digraph padded_scalar_self_energies {
                edge [particle="H" num=1];
                node [num=1];
                ext [style=invis];
                ext -> A:0 [id=0];
                B:1 -> ext [id=1];
                A -> B [id=2];
                A -> B [id=3];
                ext -> C:8 [id=4];
                D:9 -> ext [id=5];
                C -> D [id=6];
                C -> D [id=7];
            }
        )
        .unwrap()
    }

    fn paper_appendix_b1_graph() -> Graph {
        include_str!(
            "../../../../../tests/resources/graphs/paper_appendix_b1_nested_gluon_self_energy.dot"
        )
        .into_graph(&crate::utils::load_generic_model("sm"))
        .unwrap()
    }

    fn paper_appendix_b1_massless_graph() -> Graph {
        include_str!(
            "../../../../../tests/resources/graphs/paper_appendix_b1_nested_gluon_self_energy.dot"
        )
        // Keep the checked-in GammaLoop route p=GL e2, k=GL e3 used to
        // realize the loop coordinates of B.1 and B.15, changing only the
        // fermion species for B.10-B.16.
        .replace("particle=t", "particle=d")
        .into_graph(&crate::utils::load_generic_model("sm"))
        .unwrap()
    }

    fn paper_appendix_b1_node(
        graph: &Graph,
        edges: impl IntoIterator<Item = usize>,
        dod: i32,
    ) -> TestNode {
        let mut subgraph = graph.empty_subgraph::<SuBitGraph>();
        for edge in edges {
            subgraph.add(graph[&EdgeIndex(edge)].1);
        }
        let external: SuBitGraph = graph.external_filter();
        let internal = graph.full_filter().subtract(&external);
        let lmb = if subgraph == internal {
            graph.loop_momentum_basis.clone()
        } else {
            graph
                .underlying
                .try_compatible_sub_lmb(
                    &subgraph,
                    graph.dummy_less_full_crown(&subgraph).subtract(&external),
                    &graph.loop_momentum_basis,
                )
                .unwrap()
        };
        TestNode {
            // Preserve the checked-in GammaLoop basis p=GL e2 and k=GL e3,
            // including for the nested child whose independent boundary
            // carrier is p.
            lmb,
            subgraph,
            dod,
            scheme: ApproximationType::IR,
        }
    }

    #[test]
    fn paper_appendix_b1_fixture_has_the_1pi_components_consistent_with_b20() {
        test_initialise().unwrap();
        let graph = paper_appendix_b1_graph();
        // arXiv:2203.11038v1 has an internal edge-label inconsistency here.
        // B.3 prints gamma_2={e1,e2,e4,e5}, which is not the 1PI fermion
        // cycle, and B.5 consequently describes a massive tadpole.  The
        // topology and the surviving (k-p)^2 cograph in B.20 instead require
        // gamma_2={e1,e2,e3,e4}.  The checked-in GL-to-paper edge map is
        // e2->e2, e3->e4, e4->e5, e5->e3 and e6->e1, so the physically
        // consistent child is GL [2,3,5,6], leaving GL e4=paper e5.
        let mut labels = graph
            .spinneys(&graph.full_filter())
            .into_iter()
            .map(|spinney| {
                (
                    spinney.filter.string_label(),
                    graph
                        .iter_edges_of(&spinney.filter)
                        .map(|(_, edge, _)| usize::from(edge))
                        .collect::<Vec<_>>(),
                    graph.compute_dod(&spinney.filter),
                    graph.boundary_pdg_set(&spinney.filter),
                    graph.internal_pdg_set(&spinney.filter),
                )
            })
            .collect::<Vec<_>>();
        labels.sort_by(|left, right| left.0.cmp(&right.0));
        assert_eq!(
            labels
                .into_iter()
                .map(|(_, edges, dod, boundary, internal)| (edges, dod, boundary, internal))
                .collect::<Vec<_>>(),
            vec![
                (vec![], 0, Default::default(), Default::default()),
                (
                    vec![2, 3, 5, 6],
                    0,
                    [-21, 21].into_iter().collect(),
                    [6].into_iter().collect(),
                ),
                (
                    vec![2, 3, 4, 5, 6],
                    2,
                    [-21, 21].into_iter().collect(),
                    [6, 21].into_iter().collect(),
                ),
                (
                    vec![3, 4],
                    1,
                    [-6, 6].into_iter().collect(),
                    [6, 21].into_iter().collect(),
                ),
            ]
        );
    }

    #[test]
    fn appendix_b1_compatible_lmb_is_independent_of_current_full_chart() {
        test_initialise().unwrap();
        let graph = paper_appendix_b1_graph();
        let canonical = graph.loop_momentum_basis.clone();
        assert_eq!(
            canonical.loop_edges.raw,
            vec![EdgeIndex(2), EdgeIndex(3)],
            "the checked-in Appendix-B.1 route is p=e2, k=e3"
        );
        let alternate = graph
            .generate_loop_momentum_bases()
            .into_iter()
            .find(|candidate| candidate.loop_edges.raw == vec![EdgeIndex(2), EdgeIndex(4)])
            .expect("Appendix B.1 must admit the alternate homogeneous [e2,e4] chart");

        let mut canonical_current = paper_appendix_b1_node(&graph, 2..=6, 2);
        canonical_current.lmb = canonical.clone();
        let mut alternate_current = paper_appendix_b1_node(&graph, 2..=6, 2);
        alternate_current.lmb = alternate;
        let root = root_node(&graph);
        let settings = UVgenerationSettings {
            generate_integrated: false,
            ..Default::default()
        };
        let ctx = UVCtx::new(&graph, &settings);

        let from_canonical =
            compatible_lmb(&ctx, &canonical_current, &root, &Full4dCts::root()).unwrap();
        let from_alternate =
            compatible_lmb(&ctx, &alternate_current, &root, &Full4dCts::root()).unwrap();
        assert_eq!(from_canonical, canonical);
        assert_eq!(from_alternate, canonical);

        for (edges, dod) in [(vec![3, 4], 1), (vec![2, 3, 5, 6], 0)] {
            let mut child = paper_appendix_b1_node(&graph, edges.clone(), dod);
            child.scheme = ApproximationType::MUV;
            let child_local =
                uv_limit(&Full4dCts::root(), &ctx, &child, &root, &child, &root).unwrap();
            let completed_child = Full4dCts::from_factorized_local(&child_local);
            if edges == vec![2, 3, 5, 6] {
                assert_eq!(
                    completed_child.0.inherited_loop_edges,
                    vec![EdgeIndex(2)],
                    "the degree-zero gamma_2 child retains only its actual loop coordinate"
                );
                assert!(
                    !paper_has_pattern(completed_child.atom(), GS.emr_mom(EdgeIndex(3), W_.x___)),
                    "physical membership of e3 in gamma_2 must not become a false coordinate dependency"
                );
            }
            let nested_route =
                compatible_lmb(&ctx, &canonical_current, &child, &completed_child).unwrap();
            assert_eq!(
                nested_route.loop_edges, canonical.loop_edges,
                "Appendix-B.1 child {edges:?} retained coordinates {:?} and changed the outer route",
                completed_child.0.inherited_loop_edges,
            );
        }
    }

    #[test]
    fn paper_appendix_b1_1pi_component_consistent_with_b20_cograph_pair_cancels() {
        test_initialise().unwrap();
        std::thread::Builder::new()
            .stack_size(64 * 1024 * 1024)
            .spawn(|| {
                let graph = paper_appendix_b1_graph();
                let gamma_two = paper_appendix_b1_node(&graph, [2, 3, 5, 6], 0);
                let gamma = paper_appendix_b1_node(&graph, [2, 3, 4, 5, 6], 2);
                let root = root_node(&graph);
                let settings = UVgenerationSettings {
                    generate_integrated: false,
                    ..Default::default()
                };
                let ctx = UVCtx::new(&graph, &settings);

                let child_ct = uv_limit(
                    &Full4dCts::root(),
                    &ctx,
                    &gamma_two,
                    &root,
                    &gamma_two,
                    &root,
                )
                .unwrap();
                let completed_child = Full4dCts::from_factorized_local(&child_ct);
                let reduced = gamma.reduced_subgraph(&gamma_two);
                let reduced_edges = graph.paired_edges(&reduced);
                assert_eq!(
                    reduced_edges,
                    vec![EdgeIndex(4)],
                    "the 1PI child consistent with B.20 leaves GL e4=paper e5 as its cograph"
                );
                let outer_lmb = compatible_lmb(&ctx, &gamma, &gamma_two, &completed_child).unwrap();
                assert!(
                    outer_lmb.edge_signatures[EdgeIndex(4)]
                        .external
                        .iter()
                        .all(|sign| sign.is_zero()),
                    "the GL e4=paper e5 cograph momentum must have no component-external shift"
                );
                let child_family = grow(&completed_child, &ctx, &gamma, &gamma_two).unwrap();
                let nested_family = uv_limit(
                    &completed_child,
                    &ctx,
                    &gamma,
                    &gamma_two,
                    &gamma,
                    &gamma_two,
                )
                .unwrap();
                assert!(
                    !finalize(child_family.clone()).is_zero(),
                    "the B.20-consistent child forest family must be nonzero before cancellation"
                );
                assert!(
                    !finalize(nested_family.atom().clone()).is_zero(),
                    "the B.20-consistent nested forest family must be nonzero before cancellation"
                );

                // Gamma/gamma_2 is GL e4=paper e5, the massless (k-p) cograph
                // displayed in B.20, with no external shift. In the canonical
                // full chart its momentum is Q(e2)-Q(e3), whereas the grown
                // child family still carries the equivalent edge-local Q(e4)
                // spelling. Route that family before comparing the two signed
                // atoms; choosing e4 as a parent loop coordinate merely to make
                // this cancellation syntactic would make the parent operator
                // depend on its forest history.
                let child_family = child_family.replace_multiple(graph.uv_wrapped_replacement(
                    gamma.subgraph(),
                    &outer_lmb,
                    &[W_.x___],
                ));
                paper_assert_normalized_zero(
                    "Appendix B.20-consistent production cograph forest pair",
                    child_family + nested_family.atom(),
                );
            })
            .unwrap()
            .join()
            .unwrap();
    }

    fn three_level_nested_scalar_graph() -> Graph {
        // The child numerator is linear in its e2 boundary carrier so a real
        // completed H atom, rather than fabricated provenance, fixes the
        // coordinate that the middle and outer operations must preserve.
        dot!(
            digraph three_level_nested_scalar {
                edge [particle=scalar_1 num=1]
                node [num=1]
                ext [style=invis]
                ext -> A:0 [id=0]
                F:1 -> ext [id=1]
                A:2 -> B:3 [id=2 lmb_id=1]
                B:4 -> C:5 [id=3 lmb_id=0 num="Q(2,spenso::cind(0))"]
                B:6 -> C:7 [id=4]
                C:8 -> D:9 [id=5]
                A:10 -> D:11 [id=6]
                D:12 -> E:13 [id=7]
                E:14 -> F:15 [id=8]
                A:16 -> F:17 [id=9 lmb_id=2]
            },
            "scalars"
        )
        .unwrap()
    }

    #[test]
    fn top_bubble_boundary_carrier_uses_the_graph_canonical_parent_route() {
        test_initialise().unwrap();
        let graph: Graph =
            include_str!("../../../../../tests/resources/graphs/gamma_star_ddbar_top_bubble.dot")
                .into_graph(&crate::utils::load_generic_model("sm"))
                .unwrap();
        let component = |edges: std::ops::RangeInclusive<usize>| {
            let mut subgraph = graph.empty_subgraph::<SuBitGraph>();
            for edge in edges {
                subgraph.add(graph[&EdgeIndex(edge)].1);
            }
            subgraph
        };
        let child_subgraph = component(7..=8);
        let child_lmb = graph
            .try_compatible_sub_lmb(
                &child_subgraph,
                graph.dummy_less_full_crown(&child_subgraph),
                &graph.loop_momentum_basis,
            )
            .unwrap();
        let child = TestNode {
            lmb: child_lmb,
            subgraph: child_subgraph,
            dod: 2,
            scheme: ApproximationType::IR,
        };
        let parent_subgraph = component(3..=8);
        let parent_lmb = graph
            .try_compatible_sub_lmb(
                &parent_subgraph,
                graph.dummy_less_full_crown(&parent_subgraph),
                &graph.loop_momentum_basis,
            )
            .unwrap();
        let parent = TestNode {
            lmb: parent_lmb,
            subgraph: parent_subgraph,
            dod: 0,
            scheme: ApproximationType::MUV,
        };
        let expected = vec![EdgeIndex(6), EdgeIndex(8)];
        assert_eq!(graph.loop_momentum_basis.loop_edges.raw, expected);
        let boundary_signature = &graph.loop_momentum_basis.edge_signatures[EdgeIndex(5)];
        let canonical_signature = &graph.loop_momentum_basis.edge_signatures[EdgeIndex(6)];
        assert!(boundary_signature.equality_up_to_sign(canonical_signature));
        assert!(
            boundary_signature
                .external
                .iter()
                .chain(canonical_signature.external.iter())
                .all(|sign| sign.is_zero()),
            "canonicalizing the top-bubble boundary carrier must not hide an affine shift",
        );
        let canonicalized = internal_expansion_carriers(&graph, [EdgeIndex(5)]);
        assert_eq!(
            canonicalized,
            vec![EdgeIndex(6)],
            "the two top-bubble boundary representatives must use q_g=e6 from the graph LMB",
        );

        let settings = UVgenerationSettings {
            generate_integrated: false,
            ..Default::default()
        };
        let ctx = UVCtx::new(&graph, &settings);
        let root = root_node(&graph);
        let direct = compatible_lmb(&ctx, &parent, &root, &Full4dCts::root()).unwrap();
        let mut completed_child = Full4dCts::root();
        completed_child.0.inherited_loop_edges = child
            .lmb
            .loop_edges
            .iter()
            .copied()
            .chain(canonicalized)
            .collect();
        let nested = compatible_lmb(&ctx, &parent, &child, &completed_child).unwrap();
        assert_eq!(direct.loop_edges.raw, expected);
        assert_eq!(nested.loop_edges.raw, expected);
    }

    #[test]
    fn top_bubble_actual_fermion_child_h2_removes_both_soft_jet_coefficients() {
        // Try alternative common-factor groupings: greedily picking one can
        // hide the cancellation. Every rewrite is exact and keeps products
        // factored; function arguments and denominator owners remain opaque.
        fn factorized_is_zero(atom: Atom, visited: &mut std::collections::HashSet<Atom>) -> bool {
            fn collect(mut atom: Atom) -> Atom {
                for _ in 0..32 {
                    let next = atom.collect_factors().expand_num();
                    if next == atom {
                        return next;
                    }
                    atom = next;
                }
                panic!("soft coefficient factor collection did not converge");
            }

            let atom = collect(atom);
            if atom.is_zero() {
                return true;
            }
            if !visited.insert(atom.clone()) {
                return false;
            }
            assert!(
                visited.len() <= 256,
                "soft coefficient factorization search exceeded its bound"
            );
            let mut sums = Vec::new();
            atom.visitor(&mut |view| {
                if matches!(view, AtomView::Fun(_)) {
                    return false;
                }
                if matches!(view, AtomView::Add(_)) {
                    sums.push(view.to_owned());
                }
                true
            });
            for sum in sums {
                let AtomView::Add(sum_view) = sum.as_view() else {
                    unreachable!();
                };
                let terms = sum_view.iter().collect::<Vec<_>>();
                for (i, left) in terms.iter().enumerate() {
                    for (j, right) in terms.iter().enumerate().skip(i + 1) {
                        let original_pair = *left + *right;
                        let pair = collect(original_pair.clone());
                        if pair == original_pair {
                            continue;
                        }
                        let regrouped = Atom::add_many(
                            terms
                                .iter()
                                .enumerate()
                                .filter_map(|(k, term)| (k != i && k != j).then_some(*term))
                                .chain([pair.as_view()]),
                        );
                        let candidate = collect(atom.replace_map(|view, context, output| {
                            if context.function_level == 0 && view == sum.as_view() {
                                **output = regrouped.clone();
                            }
                        }));
                        if candidate.as_view().get_byte_size() < atom.as_view().get_byte_size()
                            && factorized_is_zero(candidate, visited)
                        {
                            return true;
                        }
                    }
                }
            }
            false
        }

        std::thread::Builder::new()
            .name("top-bubble-4d-soft-jet".to_string())
            .stack_size(64 * 1024 * 1024)
            .spawn(|| {
                test_initialise().unwrap();
                let graph: Graph = include_str!(
                    "../../../../../tests/resources/graphs/gamma_star_ddbar_top_bubble.dot"
                )
                .into_graph(&crate::utils::load_generic_model("sm"))
                .unwrap();
                let mut child_subgraph = graph.empty_subgraph::<SuBitGraph>();
                for edge in [EdgeIndex(7), EdgeIndex(8)] {
                    child_subgraph.add(graph[&edge].1);
                }
                let child = TestNode {
                    lmb: graph
                        .try_compatible_sub_lmb(
                            &child_subgraph,
                            graph.dummy_less_full_crown(&child_subgraph),
                            &graph.loop_momentum_basis,
                        )
                        .unwrap(),
                    subgraph: child_subgraph,
                    dod: 2,
                    scheme: ApproximationType::IR,
                };
                let given = root_node(&graph);
                let settings = UVgenerationSettings {
                    generate_integrated: false,
                    ..Default::default()
                };
                let ctx = UVCtx::new(&graph, &settings);
                let lmb = compatible_lmb(&ctx, &child, &given, &Full4dCts::root()).unwrap();
                let bubble = grow(&Full4dCts::root(), &ctx, &child, &given).unwrap();
                let ordinary = t_raw(&bubble, &ctx, &child, &given, &lmb).unwrap();
                let soft = tilde_t_raw(&bubble, &ctx, &child, &given, &lmb).unwrap();
                let overlap = t_raw(&soft, &ctx, &child, &given, &lmb).unwrap();
                let refined = &ordinary + &soft - &overlap;
                paper_assert_normalized_zero(
                    "the actual top-bubble H2 uses the production U+S-US composition",
                    &refined - hat_t_raw(&bubble, &ctx, &child, &given, &lmb).unwrap(),
                );

                // Extract the two coefficients independently from the actual
                // routed fermion numerator and denominators. This deliberately
                // does not call `tilde_t_raw` on either remainder: the oracle
                // exposes the q_g^0 and q_g^1 coefficients separately.
                let soft_coefficients = |integrand: &Atom| {
                    let reduced = child.reduced_subgraph(&given);
                    let replacements =
                        graph.uv_wrapped_replacement(&reduced, &lmb, &[W_.x___]);
                    let mut deformed = GS.erase_uv_momentum_provenance(integrand)
                        .replace_multiple(&replacements);
                    let external_edges = component_external_edges(&ctx, &child, &lmb);
                    assert!(
                        external_edges.contains(&EdgeIndex(6)),
                        "the 4D soft jet must scale the selected q_g=e6 carrier: {external_edges:?}",
                    );
                    for edge in external_edges {
                        deformed = deformed
                            .replace(GS.emr_mom(edge, W_.x___))
                            .with(GS.emr_mom(edge, W_.x___) * GS.rescale);
                    }
                    let series = deformed.series(GS.rescale, Atom::Zero, 1).unwrap();
                    [0, 1].map(|order| {
                        let coefficient = series
                            .coefficient(Rational::from(order))
                            .expect("requested coefficient is within series precision");
                        finalize(coefficient.replace(GS.dim).with(Atom::num(4)))
                            .collect_factors()
                    })
                };

                let ordinary_remainder = soft_coefficients(&(&bubble - &ordinary));
                assert!(
                    !factorized_is_zero(ordinary_remainder[0].clone(), &mut Default::default()),
                    "ordinary U2 must leave a nonzero q_g^0 coefficient in the actual massive-top bubble",
                );
                let refined_remainder = soft_coefficients(&(&bubble - &refined));
                assert!(
                    refined_remainder.iter().all(|coefficient| {
                        factorized_is_zero(coefficient.clone(), &mut Default::default())
                    }),
                    "(1-H2)Pi must have vanishing q_g^0 and q_g^1 coefficients; got {refined_remainder:?}",
                );
            })
            .unwrap()
            .join()
            .unwrap();
    }

    #[test]
    fn nested_soft_boundary_does_not_constrain_the_parent_loop_basis() {
        test_initialise().unwrap();
        let graph: Graph =
            include_str!("../../../../../tests/resources/graphs/dgse_local_ir_offshell.dot")
                .into_graph(&crate::utils::load_generic_model("sm"))
                .unwrap();
        let mut child_subgraph = graph.empty_subgraph::<SuBitGraph>();
        for edge in [EdgeIndex(5), EdgeIndex(6)] {
            child_subgraph.add(graph[&edge].1);
        }
        let child = TestNode {
            lmb: graph.lmb_of(&child_subgraph),
            subgraph: child_subgraph,
            dod: 1,
            scheme: ApproximationType::IR,
        };
        let external: SuBitGraph = graph.external_filter();
        let parent_subgraph = graph.full_filter().subtract(&external);
        let parent = TestNode {
            lmb: graph.lmb_of(&parent_subgraph),
            subgraph: parent_subgraph,
            dod: 0,
            scheme: ApproximationType::MUV,
        };
        let settings = UVgenerationSettings {
            generate_integrated: false,
            ..Default::default()
        };
        let ctx = UVCtx::new(&graph, &settings);
        let given = root_node(&graph);
        let grown = grow(&Full4dCts::root(), &ctx, &child, &given).unwrap();
        let completed = hat_t_raw(&grown, &ctx, &child, &given, child.lmb()).unwrap();
        let preferred_expansion_carriers = internal_expansion_carriers(
            &graph,
            independent_component_external_edges(&ctx, &child, child.lmb(), &completed).unwrap(),
        );
        assert_eq!(preferred_expansion_carriers.len(), 1);
        let carrier = preferred_expansion_carriers[0];
        assert!(!parent.lmb.loop_edges.contains(&carrier));
        let child_local =
            uv_limit(&Full4dCts::root(), &ctx, &child, &given, &child, &given).unwrap();
        assert!(
            !child_local.0.inherited_loop_edges.contains(&carrier),
            "the child boundary is a routed dependency, not an inherited loop generator"
        );
        let completed_child = Full4dCts::from_factorized_local(&child_local);
        let selected = compatible_lmb(&ctx, &parent, &child, &completed_child).unwrap();

        assert!(selected.loop_edges.contains(&carrier));
        assert!(
            child
                .lmb
                .loop_edges
                .iter()
                .all(|edge| selected.loop_edges.contains(edge))
        );

        // The completed soft refinement is child-UV finite. On the canonical
        // homogeneous route a logarithmic containing component has no
        // coefficient left to project. This is the 4D oracle for checking
        // that the 3D ordered replay uses the same containing route and active
        // variables without changing its orientation-resolved CFF.
        let ordinary = t_raw(&grown, &ctx, &child, &given, child.lmb()).unwrap();
        let refinement = &completed - &ordinary;
        assert!(
            !refinement.expand().is_zero(),
            "the DGSE H-U child refinement must be non-vacuous"
        );
        let co_graph = grow_factor(&ctx, &parent, &child).unwrap();
        let parent_image =
            t_raw(&(&refinement * co_graph), &ctx, &parent, &child, &selected).unwrap();
        paper_assert_normalized_zero(
            "DGSE U_parent,0 annihilates the completed H-U child refinement",
            parent_image,
        );
    }

    #[test]
    fn nested_ordinary_u_boundary_does_not_constrain_the_h_parent_loop_basis() {
        test_initialise().unwrap();
        let graph: Graph =
            include_str!("../../../../../tests/resources/graphs/dgse_local_ir_offshell.dot")
                .into_graph(&crate::utils::load_generic_model("sm"))
                .unwrap();
        let mut child_subgraph = graph.empty_subgraph::<SuBitGraph>();
        for edge in [EdgeIndex(5), EdgeIndex(6)] {
            child_subgraph.add(graph[&edge].1);
        }
        let child = TestNode {
            lmb: graph.lmb_of(&child_subgraph),
            subgraph: child_subgraph,
            dod: 1,
            scheme: ApproximationType::MUV,
        };
        let external: SuBitGraph = graph.external_filter();
        let parent_subgraph = graph.full_filter().subtract(&external);
        let parent = TestNode {
            lmb: graph.lmb_of(&parent_subgraph),
            subgraph: parent_subgraph,
            dod: 1,
            scheme: ApproximationType::IR,
        };
        let settings = UVgenerationSettings {
            generate_integrated: false,
            ..Default::default()
        };
        let ctx = UVCtx::new(&graph, &settings);
        let given = root_node(&graph);
        let grown = grow(&Full4dCts::root(), &ctx, &child, &given).unwrap();
        let completed = t_raw(&grown, &ctx, &child, &given, child.lmb()).unwrap();
        let carriers = internal_expansion_carriers(
            &graph,
            independent_component_external_edges(&ctx, &child, child.lmb(), &completed).unwrap(),
        );
        assert_eq!(carriers.len(), 1);
        let carrier = carriers[0];
        assert!(
            !parent.lmb.loop_edges.contains(&carrier),
            "the fixture must force the H parent away from its default route"
        );

        let child_local =
            uv_limit(&Full4dCts::root(), &ctx, &child, &given, &child, &given).unwrap();
        assert!(
            !child_local.0.inherited_loop_edges.contains(&carrier),
            "a positive-degree U child must inherit only its actual loop generators"
        );
        let completed_child = Full4dCts::from_factorized_local(&child_local);
        let selected = compatible_lmb(&ctx, &parent, &child, &completed_child).unwrap();
        let pure = compatible_lmb(&ctx, &parent, &given, &Full4dCts::root()).unwrap();
        assert_eq!(
            selected, pure,
            "the H parent must use the same basis with or without the U child"
        );

        let parent_local =
            uv_limit(&completed_child, &ctx, &parent, &child, &parent, &child).unwrap();
        assert_eq!(
            parent_local.0.provenance_components[0].route_loop_edges,
            selected
                .loop_edges
                .iter()
                .map(|edge| usize::from(*edge))
                .collect::<Vec<_>>(),
            "the production H parent must extend the child's loop basis"
        );
    }

    #[test]
    fn numerator_only_nested_carrier_is_preserved_from_the_completed_atom() {
        test_initialise().unwrap();
        let graph: Graph = dot!(
            digraph nested_numerator_only_tadpole {
                edge [particle=scalar_1 num=1]
                node [num=1]
                ext [style=invis]
                ext -> A:0 [id=0]
                A:1 -> ext [id=1]
                A:2 -> B:3 [id=2]
                B:4 -> B:5 [id=3 lmb_id=0 num="Q(4,spenso::cind(0))"]
                B:6 -> A:7 [id=4 lmb_id=1]
            },
            "scalars"
        )
        .unwrap();
        let mut child_subgraph = graph.empty_subgraph::<SuBitGraph>();
        child_subgraph.add(graph[&EdgeIndex(3)].1);
        let external: SuBitGraph = graph.external_filter();
        let child_lmb = graph
            .underlying
            .try_compatible_sub_lmb(
                &child_subgraph,
                graph
                    .dummy_less_full_crown(&child_subgraph)
                    .subtract(&external),
                &graph.loop_momentum_basis,
            )
            .unwrap();
        let carrier = EdgeIndex(4);
        let carrier_slot = child_lmb
            .ext_edges
            .iter_enumerated()
            .find(|(_, edge)| **edge == carrier)
            .map(|(external, _)| external)
            .expect("the parent-loop edge e4 must be the child boundary coordinate");
        assert!(
            graph
                .paired_edges(&child_subgraph)
                .iter()
                .all(|edge| child_lmb.edge_signatures[*edge].external[carrier_slot].is_zero()),
            "the tadpole propagator must have rank-zero dependence on its boundary coordinate"
        );
        let child = TestNode {
            lmb: child_lmb,
            subgraph: child_subgraph,
            dod: 1,
            scheme: ApproximationType::IR,
        };
        let given = root_node(&graph);
        let settings = UVgenerationSettings {
            generate_integrated: false,
            ..Default::default()
        };
        let ctx = UVCtx::new(&graph, &settings);
        let grown = grow(&Full4dCts::root(), &ctx, &child, &given).unwrap();
        let completed = hat_t_raw(&grown, &ctx, &child, &given, child.lmb()).unwrap();
        assert!(paper_has_pattern(
            &completed,
            GS.emr_mom(carrier, GS.cind(0))
        ));
        assert_eq!(
            independent_component_external_edges(&ctx, &child, child.lmb(), &completed,).unwrap(),
            vec![carrier]
        );
        let local = uv_limit(&Full4dCts::root(), &ctx, &child, &given, &child, &given).unwrap();
        assert!(
            !local.0.inherited_loop_edges.contains(&carrier),
            "a numerator-only boundary dependency must not become a child loop generator"
        );
    }

    #[test]
    fn degree_zero_ir_reuses_ordinary_t_basis_metadata() {
        test_initialise().unwrap();
        let graph: Graph = include_str!("../../../../../tests/resources/graphs/dgse.dot")
            .into_graph(&crate::utils::load_generic_model("sm"))
            .unwrap();
        let mut subgraph = graph.empty_subgraph::<SuBitGraph>();
        for edge in [EdgeIndex(5), EdgeIndex(6)] {
            subgraph.add(graph[&edge].1);
        }
        let ir = TestNode {
            lmb: graph.lmb_of(&subgraph),
            subgraph: subgraph.clone(),
            dod: 0,
            scheme: ApproximationType::IR,
        };
        let ordinary = TestNode {
            lmb: graph.lmb_of(&subgraph),
            subgraph,
            dod: 0,
            scheme: ApproximationType::MUV,
        };
        let given = root_node(&graph);
        let settings = UVgenerationSettings {
            generate_integrated: false,
            ..Default::default()
        };
        let ctx = UVCtx::new(&graph, &settings);

        let ir_ct = uv_limit(&Full4dCts::root(), &ctx, &ir, &given, &ir, &given).unwrap();
        let ordinary_ct = uv_limit(
            &Full4dCts::root(),
            &ctx,
            &ordinary,
            &given,
            &ordinary,
            &given,
        )
        .unwrap();

        assert_eq!(
            ir_ct.0.inherited_loop_edges,
            ordinary_ct.0.inherited_loop_edges
        );
        assert!(!ir_ct.0.has_soft_ancestry);
        assert!(!ordinary_ct.0.has_soft_ancestry);
    }

    #[test]
    fn deep_nesting_propagates_the_basis_selected_by_the_completed_middle_atom() {
        test_initialise().unwrap();
        let graph = three_level_nested_scalar_graph();
        let mut child_subgraph = graph.empty_subgraph::<SuBitGraph>();
        for edge in [EdgeIndex(3), EdgeIndex(4)] {
            child_subgraph.add(graph[&edge].1);
        }
        let mut middle_subgraph = graph.empty_subgraph::<SuBitGraph>();
        for edge in [2, 3, 4, 5, 6].map(EdgeIndex) {
            middle_subgraph.add(graph[&edge].1);
        }
        let external: SuBitGraph = graph.external_filter();
        let outer_subgraph = graph.full_filter().subtract(&external);
        let nested_lmb = |subgraph: &SuBitGraph| {
            graph
                .try_compatible_sub_lmb(
                    subgraph,
                    graph.dummy_less_full_crown(subgraph),
                    &graph.loop_momentum_basis,
                )
                .unwrap()
        };
        let child = TestNode {
            lmb: nested_lmb(&child_subgraph),
            subgraph: child_subgraph,
            dod: 1,
            scheme: ApproximationType::IR,
        };
        let middle = TestNode {
            lmb: nested_lmb(&middle_subgraph),
            subgraph: middle_subgraph,
            dod: 1,
            scheme: ApproximationType::IR,
        };
        let outer = TestNode {
            lmb: nested_lmb(&outer_subgraph),
            subgraph: outer_subgraph,
            dod: 0,
            scheme: ApproximationType::MUV,
        };
        let settings = UVgenerationSettings {
            generate_integrated: false,
            ..Default::default()
        };
        let ctx = UVCtx::new(&graph, &settings);
        let given = root_node(&graph);
        let child_local =
            uv_limit(&Full4dCts::root(), &ctx, &child, &given, &child, &given).unwrap();
        assert_eq!(child_local.0.inherited_loop_edges, child.lmb.loop_edges.raw);
        let completed_child = Full4dCts::from_factorized_local(&child_local);
        let selected_middle = compatible_lmb(&ctx, &middle, &child, &completed_child).unwrap();
        assert_eq!(
            selected_middle.loop_edges, middle.lmb.loop_edges,
            "the middle route must extend the child without changing the pure-parent chart"
        );
        assert!(
            child
                .lmb
                .loop_edges
                .iter()
                .all(|edge| selected_middle.loop_edges.contains(edge))
        );

        let middle_local =
            uv_limit(&completed_child, &ctx, &middle, &child, &middle, &child).unwrap();
        assert!(middle_local.0.has_soft_ancestry);
        assert!(
            selected_middle
                .loop_edges
                .iter()
                .all(|edge| middle_local.0.inherited_loop_edges.contains(edge)),
            "the production H middle atom must record the route on which it was completed"
        );
        let completed_middle = Full4dCts::from_factorized_local(&middle_local);
        assert!(completed_middle.0.has_soft_ancestry);
        let selected_outer = compatible_lmb(&ctx, &outer, &middle, &completed_middle).unwrap();
        assert!(
            selected_middle
                .loop_edges
                .iter()
                .all(|edge| selected_outer.loop_edges.contains(edge)),
            "the grandparent must extend the basis actually selected by the completed middle atom"
        );
    }

    #[test]
    fn nested_expansion_rejects_affine_carrier_promotion() {
        test_initialise().unwrap();

        // If h=k+p is promoted to a loop coordinate, scaling p while holding h
        // fixed induces k -> k+(1-u)p.  For f=k^2, the finite jets U_1 and S_0
        // then depend on their order:
        //
        //   U_1 S_0 f = k^2-p^2,       S_0 U_1 f = k^2.
        //
        // This is the algebraic obstruction behind the production routing
        // rejection below.  Homogeneous routes do not induce the additive
        // shift and are covered by the commuting-jet tests separately.
        let loop_momentum = Atom::var(symbol!("affine_route_loop_momentum"));
        let external_momentum = Atom::var(symbol!("affine_route_external_momentum"));
        let parameter = Atom::var(GS.rescale);
        let bare = loop_momentum.pow(2);
        let soft_bare = paper_series_through(
            bare.replace(external_momentum.clone())
                .with(&external_momentum * &parameter),
            0,
        );
        let shifted_loop = &loop_momentum + (Atom::one() - &parameter) * &external_momentum;
        let uv_after_soft = paper_series_through(
            soft_bare
                .replace(loop_momentum.clone())
                .with(shifted_loop.clone()),
            1,
        );
        let uv_bare =
            paper_series_through(bare.replace(loop_momentum.clone()).with(shifted_loop), 1);
        let soft_after_uv = paper_series_through(
            uv_bare
                .replace(external_momentum.clone())
                .with(&external_momentum * parameter),
            0,
        );
        assert_eq!(
            (&uv_after_soft - (&bare - external_momentum.pow(2))).expand(),
            Atom::Zero,
            "U_1 S_0 at fixed affine carrier h=k+p"
        );
        assert_eq!(
            (&soft_after_uv - &bare).expand(),
            Atom::Zero,
            "S_0 U_1 at fixed affine carrier h=k+p"
        );
        assert_eq!(
            (&uv_after_soft - &soft_after_uv + external_momentum.pow(2)).expand(),
            Atom::Zero,
            "affine finite jets must exhibit the predicted noncommutator"
        );

        let graph: Graph = dot!(
            digraph nested_soft_affine_routing {
                edge [particle=scalar_1 num=1]
                node [num=1]
                ext [style=invis]
                ext -> v1:0 [id=0]
                v4:1 -> ext [id=1]
                v1 -> v2 [id=2]
                v2 -> v3 [id=3 lmb_id=0]
                v2 -> v3 [id=4]
                v3 -> v4 [id=5]
                v1 -> v4 [id=6 lmb_id=1]
            },
            "scalars"
        )
        .unwrap();
        let mut child_subgraph = graph.empty_subgraph::<SuBitGraph>();
        for edge in [EdgeIndex(3), EdgeIndex(4)] {
            child_subgraph.add(graph[&edge].1);
        }
        let external: SuBitGraph = graph.external_filter();
        let parent_subgraph = graph.full_filter().subtract(&external);
        let child = TestNode {
            lmb: graph.lmb_of(&child_subgraph),
            subgraph: child_subgraph,
            dod: 2,
            scheme: ApproximationType::IR,
        };
        let parent = TestNode {
            lmb: graph.lmb_of(&parent_subgraph),
            subgraph: parent_subgraph.clone(),
            dod: 2,
            scheme: ApproximationType::IR,
        };
        let affine_carrier = EdgeIndex(2);
        assert!(
            graph.loop_momentum_basis.edge_signatures[affine_carrier]
                .external
                .iter()
                .any(|sign| sign.is_sign()),
            "the test carrier must have a graph-external offset in the canonical route"
        );
        let affine_lmb = graph
            .generate_loop_momentum_bases_of(&parent_subgraph)
            .into_iter()
            .find(|candidate| candidate.loop_edges.contains(&affine_carrier))
            .expect("the test parent must admit the affine carrier as a nominal loop edge");
        let affine_parent = TestNode {
            lmb: affine_lmb,
            subgraph: parent_subgraph,
            dod: 2,
            scheme: ApproximationType::IR,
        };
        let settings = UVgenerationSettings {
            generate_integrated: false,
            ..Default::default()
        };
        let ctx = UVCtx::new(&graph, &settings);
        let root = root_node(&graph);
        let selected = compatible_lmb(&ctx, &affine_parent, &root, &Full4dCts::root()).unwrap();
        assert!(
            selected.loop_edges.iter().all(|edge| {
                graph.loop_momentum_basis.edge_signatures[*edge]
                    .external
                    .iter()
                    .all(|sign| sign.is_zero())
            }),
            "a primitive IR operation must replace an affine nominal LMB with a homogeneous one"
        );
        let child_local = uv_limit(&Full4dCts::root(), &ctx, &child, &root, &child, &root).unwrap();
        let completed_child = Full4dCts::from_factorized_local(&child_local);
        assert_eq!(completed_child.0.inherited_loop_edges, vec![EdgeIndex(3)]);
        assert!(paper_has_pattern(
            completed_child.atom(),
            GS.emr_mom(affine_carrier, W_.x___)
        ));
        let nested = compatible_lmb(&ctx, &parent, &child, &completed_child).unwrap();
        assert_eq!(
            nested, selected,
            "pure and nested parent terms must use the same expansion coordinates"
        );
        assert!(!nested.loop_edges.contains(&affine_carrier));
        assert!(
            completed_child
                .0
                .inherited_loop_edges
                .iter()
                .all(|edge| nested.loop_edges.contains(edge))
        );

        // q=r+p stays a boundary expression. The parent soft deformation is
        // q -> r+u*p at fixed inherited loops, never q -> u*q.
        let indices = [GS.cind(0)];
        let hard = nested.loop_atom(affine_carrier, GS.emr_mom, &indices, true);
        let soft = nested.ext_atom(affine_carrier, GS.emr_mom, &indices, true);
        assert!(!hard.is_zero() && !soft.is_zero());
        let boundary = GS.emr_mom(affine_carrier, GS.cind(0)).pow(2);
        let deformed =
            graph.uv_rescaled_with_loops(parent.subgraph(), 0, &nested, &nested, true, &boundary);
        let deformed = GS.erase_uv_momentum_provenance(&deformed);
        let parameter = Atom::var(GS.rescale);
        assert_eq!(
            (&deformed - (&hard + &parameter * &soft).pow(2)).expand(),
            Atom::Zero
        );
        assert!(
            !(&deformed - (&parameter * (&hard + &soft)).pow(2))
                .expand()
                .is_zero()
        );
        // The same parent operation must also accept the completed child itself.
        let nested_local =
            uv_limit(&completed_child, &ctx, &parent, &child, &parent, &child).unwrap();
        assert!(!nested_local.atom().is_zero());

        // An actual frozen loop coordinate may still not be redefined by an
        // affine shift. Boundary dependencies no longer enter this metadata.
        let mut completed_child = Full4dCts::root();
        completed_child.0.inherited_loop_edges = vec![affine_carrier];
        let error = compatible_lmb(&ctx, &parent, &child, &completed_child)
            .expect_err("an affine child carrier must not become a fixed parent loop coordinate");
        let message = error.to_string();
        assert!(
            message.contains("affine loop-momentum shift"),
            "unexpected routing error: {message}"
        );
        assert!(message.contains(&parent.subgraph().string_label()));
    }

    #[test]
    fn four_photon_self_energy_inherits_loops_in_the_original_graph_basis() {
        test_initialise().unwrap();
        let model = crate::utils::load_generic_model("sm");
        for source in [
            include_str!(
                "../../../../../examples/cli/aa_aa/3L/graphs/processes/amplitudes/aa_aa/3L/GL262.dot"
            ),
            include_str!(
                "../../../../../examples/cli/aa_aa/3L/graphs/processes/amplitudes/aa_aa/3L/GL256.dot"
            ),
        ] {
            let graph: Graph = source.into_graph(&model).unwrap();
            assert_eq!(
                graph.loop_momentum_basis.loop_edges.raw,
                vec![EdgeIndex(4), EdgeIndex(9), EdgeIndex(12)]
            );
            let node = |edges: &[usize], dod| {
                let mut subgraph = graph.empty_subgraph::<SuBitGraph>();
                for edge in edges {
                    subgraph.add(graph[&EdgeIndex(*edge)].1);
                }
                let lmb = graph
                    .try_compatible_sub_lmb(
                        &subgraph,
                        graph.dummy_less_full_crown(&subgraph),
                        &graph.loop_momentum_basis,
                    )
                    .unwrap();
                TestNode {
                    subgraph,
                    lmb,
                    dod,
                    scheme: ApproximationType::IR,
                }
            };
            let child = node(&[9, 10], 2);
            let parent = node(&[4, 5, 7, 9, 10], 1);
            let outer = node(&(4..=13).collect::<Vec<_>>(), 0);
            let root = root_node(&graph);
            let settings = UVgenerationSettings {
                generate_integrated: false,
                ..Default::default()
            };
            let ctx = UVCtx::new(&graph, &settings);
            let child_local =
                uv_limit(&Full4dCts::root(), &ctx, &child, &root, &child, &root).unwrap();
            let completed_child = Full4dCts::from_factorized_local(&child_local);
            assert_eq!(completed_child.0.inherited_loop_edges, vec![EdgeIndex(9)]);
            let selected = compatible_lmb(&ctx, &parent, &child, &completed_child).unwrap();
            assert_eq!(selected.loop_edges.raw, vec![EdgeIndex(4), EdgeIndex(9)]);
            assert_eq!(
                selected,
                compatible_lmb(&ctx, &parent, &root, &Full4dCts::root()).unwrap()
            );
            for edge in [EdgeIndex(5), EdgeIndex(7)] {
                assert!(!selected.loop_edges.contains(&edge));
                assert!(
                    graph.loop_momentum_basis.edge_signatures[edge]
                        .external
                        .iter()
                        .any(|sign| sign.is_sign())
                );
            }
            let parent_local =
                uv_limit(&completed_child, &ctx, &parent, &child, &parent, &child).unwrap();
            assert!(!parent_local.atom().is_zero());
            let completed_parent = Full4dCts::from_factorized_local(&parent_local);
            let selected_outer = compatible_lmb(&ctx, &outer, &parent, &completed_parent).unwrap();
            assert_eq!(
                selected_outer.loop_edges,
                graph.loop_momentum_basis.loop_edges
            );
            assert_eq!(
                selected_outer,
                compatible_lmb(&ctx, &outer, &root, &Full4dCts::root()).unwrap()
            );
        }
    }

    #[test]
    fn tilde_taylor_has_degree_d_minus_one_and_keeps_physical_mass() {
        test_initialise().unwrap();
        let mass = Atom::var(symbol!("local_4d_test_mass"));
        let linear = Atom::var(symbol!("local_4d_test_linear"));
        let quadratic = Atom::var(symbol!("local_4d_test_quadratic"));
        let rescaled = &mass + &linear * GS.rescale + quadratic * Atom::var(GS.rescale).pow(2);

        let degree_zero = tilde_taylor(rescaled.clone(), 1).unwrap();
        let actual = tilde_taylor(rescaled, 2).unwrap();

        assert_eq!((degree_zero - &mass).expand(), Atom::Zero);
        assert_eq!((actual.clone() - mass - linear).expand(), Atom::Zero);
        assert_eq!(
            actual
                .clone()
                .replace(GS.m_uv_expansion)
                .with(Atom::num(37)),
            actual,
            "the soft operator must not introduce mUV"
        );
    }

    #[test]
    fn soft_projection_does_not_scale_padded_cograph_externals() {
        test_initialise().unwrap();
        let graph = padded_two_point_graph();
        let external: SuBitGraph = graph.external_filter();
        let internal = graph.full_filter().subtract(&external);
        let component = graph
            .connected_components(&internal)
            .into_iter()
            .find(|component| graph.paired_edges(component).contains(&EdgeIndex(2)))
            .unwrap();
        let current = TestNode {
            lmb: graph.lmb_of(&component),
            subgraph: component,
            dod: 1,
            scheme: ApproximationType::IR,
        };
        let given = root_node(&graph);
        let settings = UVgenerationSettings {
            generate_integrated: false,
            ..Default::default()
        };
        let component_edges = graph.paired_edges(current.subgraph());
        let component_carriers: Vec<_> = graph
            .dummy_less_full_crown(current.subgraph())
            .included_iter()
            .map(|hedge| graph[&hedge])
            .collect();
        let active = current
            .lmb
            .ext_edges
            .iter_enumerated()
            .find(|(external, edge)| {
                component_carriers.contains(edge)
                    || component_edges.iter().any(|component_edge| {
                        current.lmb.edge_signatures[*component_edge].external[*external].is_sign()
                    })
            })
            .map(|(_, edge)| *edge)
            .expect("the component must expose an independent external carrier");
        let padded = current
            .lmb
            .ext_edges
            .iter_enumerated()
            .find(|(external, edge)| {
                !component_carriers.contains(edge)
                    && component_edges.iter().all(|component_edge| {
                        !current.lmb.edge_signatures[*component_edge].external[*external].is_sign()
                    })
            })
            .map(|(_, edge)| *edge)
            .expect("the canonical component LMB must retain a padded cograph external");
        let active_momentum = GS.emr_mom(active, Atom::num(0));
        let padded_momentum = GS.emr_mom(padded, Atom::num(0));

        let projected = tilde_t_raw(
            &(&active_momentum + &padded_momentum),
            &UVCtx::new(&graph, &settings),
            &current,
            &given,
            current.lmb(),
        )
        .unwrap();

        assert_eq!((projected - padded_momentum).expand(), Atom::Zero);
    }

    #[test]
    fn logarithmic_tilde_is_zero_so_hat_is_t() {
        assert!(
            tilde_taylor(Atom::var(GS.m_uv_vacuum), 0)
                .unwrap()
                .is_zero()
        );
    }

    #[test]
    fn local_soft_provenance_records_the_selected_route_and_branch_sizes() {
        test_initialise().unwrap();
        let graph = scalar_two_point_graph();
        let given = root_node(&graph);
        let mut current = two_point_node(&graph, ApproximationType::IR);
        current.dod = 2;
        let settings = UVgenerationSettings {
            generate_integrated: false,
            ..Default::default()
        };
        let ctx = UVCtx::new(&graph, &settings);
        let root = Full4dCts::root();
        let selected = compatible_lmb(&ctx, &current, &given, &root).unwrap();
        let local = uv_limit(&root, &ctx, &current, &given, &current, &given).unwrap();

        assert!(!root.0.has_soft_ancestry);
        assert!(local.0.has_soft_ancestry);
        let completed = Full4dCts::from_factorized_local(&local);
        assert!(completed.0.has_soft_ancestry);
        for factors in [[root.clone(), completed.clone()], [completed, root]] {
            let product = Local4dCts::from_full_product(factors);
            assert!(product.0.has_soft_ancestry);
            assert_eq!(product.atom(), local.atom());
        }

        let [component] = local.provenance_components() else {
            panic!("a connected local operation must export exactly one component route")
        };
        assert_eq!(component.component, current.subgraph().string_label());
        assert_eq!(component.scheme, ApproximationType::IR);
        assert_eq!(component.dod, 2);
        assert_eq!(
            component.route_loop_edges,
            selected
                .loop_edges
                .iter()
                .map(|edge| usize::from(*edge))
                .collect::<Vec<_>>()
        );
        assert_eq!(
            component.route_external_edges,
            selected
                .ext_edges
                .iter()
                .map(|edge| usize::from(*edge))
                .collect::<Vec<_>>()
        );
        assert_eq!(
            component
                .route_signatures
                .iter()
                .map(|signature| signature.edge)
                .collect::<Vec<_>>(),
            graph
                .paired_edges(current.subgraph())
                .into_iter()
                .map(usize::from)
                .collect::<Vec<_>>()
        );
        assert_eq!(
            local
                .provenance_branches()
                .iter()
                .map(|branch| branch.branch)
                .collect::<Vec<_>>(),
            vec![
                Local4dBranch::U,
                Local4dBranch::S,
                Local4dBranch::US,
                Local4dBranch::Combined,
            ]
        );
        assert!(
            local
                .provenance_branches()
                .iter()
                .all(|branch| branch.byte_size > 0)
        );
    }

    #[test]
    fn production_hat_operator_satisfies_complement_and_refinement_identities() {
        test_initialise().unwrap();
        let graph = scalar_two_point_graph();
        let given = root_node(&graph);
        let settings = UVgenerationSettings {
            generate_integrated: false,
            ..Default::default()
        };
        let ctx = UVCtx::new(&graph, &settings);
        let physical_mass = graph[EdgeIndex(2)].mass_atom();
        assert!(!physical_mass.is_zero());

        for degree in 0..=2 {
            let mut current = two_point_node(&graph, ApproximationType::IR);
            current.dod = degree;
            let root = Full4dCts::root();
            let lmb = compatible_lmb(&ctx, &current, &given, &root).unwrap();
            let grown = grow(&root, &ctx, &current, &given).unwrap();
            let ordinary = t_raw(&grown, &ctx, &current, &given, &lmb).unwrap();
            let soft = tilde_t_raw(&grown, &ctx, &current, &given, &lmb).unwrap();
            let overlap = t_raw(&soft, &ctx, &current, &given, &lmb).unwrap();
            let soft_after_ordinary = tilde_t_raw(&ordinary, &ctx, &current, &given, &lmb).unwrap();
            let composite = hat_t_raw(&grown, &ctx, &current, &given, &lmb).unwrap();
            let soft_remainder = &grown - &soft;
            let uv_of_soft_remainder =
                t_raw(&soft_remainder, &ctx, &current, &given, &lmb).unwrap();

            paper_assert_normalized_zero(
                &format!("production T-hat branch composition at d={degree}"),
                &composite - &ordinary - &soft + &overlap,
            );
            paper_assert_normalized_zero(
                &format!("1-H=(1-U)(1-S) at d={degree}"),
                (&grown - &composite) - (soft_remainder - uv_of_soft_remainder),
            );
            paper_assert_normalized_zero(
                &format!("H-U=(1-U)S at d={degree}"),
                (&composite - &ordinary) - (&soft - &overlap),
            );
            // In one homogeneous chart U filters total (external-momentum,
            // physical-mass) degree while S filters external-momentum degree
            // only.  Literal/terminal M is degree zero, so both compositions
            // select the same intersection of bidegrees.  The implementation
            // nevertheless retains the paper's explicit U(S) ordering.
            paper_assert_normalized_zero(
                &format!("US=SU for homogeneous soft/UV jets at d={degree}"),
                &overlap - &soft_after_ordinary,
            );

            for (label, branch) in [
                ("ordinary", &ordinary),
                ("soft", &soft),
                ("overlap", &overlap),
                ("composite", &composite),
            ] {
                assert!(
                    !paper_has_pattern(branch, Atom::var(GS.rescale)),
                    "the evaluated {label} branch at d={degree} retains the Taylor parameter"
                );
            }
            if degree == 0 {
                assert!(soft.is_zero());
                paper_assert_normalized_zero("H_0=U_0", composite - ordinary);
            } else {
                assert!(paper_has_pattern(&soft, physical_mass.clone()));
                assert!(!paper_has_muv(&soft));
                assert!(paper_has_muv(&ordinary));
                assert!(paper_has_muv(&overlap));
            }
        }
    }

    #[test]
    fn production_nested_scheme_matrix_applies_supported_local_operator_pairs() {
        let test = std::thread::Builder::new()
            .name("production-nested-scheme-matrix".to_string())
            .stack_size(64 * 1024 * 1024)
            .spawn(|| {
                test_initialise().unwrap();
                let graph: Graph = dot!(
                    digraph nested_scalar_scheme_matrix {
                        edge [particle=scalar_1 num=1]
                        node [num=1]
                        ext [style=invis]
                        ext -> v1:0 [id=0]
                        v4:1 -> ext [id=1]
                        v1 -> v2 [id=2 lmb_id=1 num="Q(2,spenso::cind(0))^2"]
                        v2 -> v3 [id=3 lmb_id=0 num="Q(3,spenso::cind(0))"]
                        v2 -> v3 [id=4]
                        v3 -> v4 [id=5]
                        v1 -> v4 [id=6]
                    },
                    "scalars"
                )
                .unwrap();
                let mut child_subgraph = graph.empty_subgraph::<SuBitGraph>();
                for edge in [EdgeIndex(3), EdgeIndex(4)] {
                    child_subgraph.add(graph[&edge].1);
                }
                let external: SuBitGraph = graph.external_filter();
                let parent_subgraph = graph.full_filter().subtract(&external);
                assert_eq!(graph.compute_dod(&child_subgraph), 1);
                assert_eq!(graph.compute_dod(&parent_subgraph), 1);
                let mut nonnegative_connected = graph
                    .spinneys(&parent_subgraph)
                    .into_iter()
                    .filter(|spinney| !spinney.filter.is_empty())
                    .filter(|spinney| graph.compute_dod(&spinney.filter) >= 0)
                    .filter(|spinney| {
                        graph
                            .underlying
                            .connected_components(&spinney.filter)
                            .len()
                            == 1
                    })
                    .map(|spinney| {
                        let mut edges = graph
                            .iter_edges_of(&spinney.filter)
                            .map(|(_, edge, _)| usize::from(edge))
                            .collect::<Vec<_>>();
                        edges.sort();
                        edges
                    })
                    .collect::<Vec<_>>();
                nonnegative_connected.sort();
                assert_eq!(
                    nonnegative_connected,
                    vec![vec![2, 3, 4, 5, 6], vec![3, 4]],
                    "the four-family fixture must have only the child and parent UV components"
                );
                let child_lmb = graph
                    .underlying
                    .try_compatible_sub_lmb(
                        &child_subgraph,
                        graph
                            .dummy_less_full_crown(&child_subgraph)
                            .subtract(&external),
                        &graph.loop_momentum_basis,
                    )
                    .unwrap();
                assert!(
                    child_lmb.ext_edges.iter().any(|edge| *edge == EdgeIndex(2)),
                    "the canonical child route must use e2 as its boundary carrier"
                );
                let given = root_node(&graph);
                let settings = UVgenerationSettings {
                    generate_integrated: false,
                    ..Default::default()
                };
                let ctx = UVCtx::new(&graph, &settings);
                let seed = Full4dCts::root();
                let zero = Rational::from(0);

                // This is a physical hard-ray oracle, independent of U/H: expose
                // each exact quadratic propagator, route every family into the same
                // parent coordinates, and scale only the selected loop coordinates.
                let hard_degree = |atom: &Atom,
                                   route: &LoopMomentumBasis,
                                   hard_edges: &[EdgeIndex],
                                   label: &str| {
                    nested_hard_degree(
                        &graph,
                        &parent_subgraph,
                        atom,
                        route,
                        hard_edges,
                        label,
                    )
                };

                let mut matrix = Vec::new();
                let mut production_obstructions = Vec::new();
                for (label, child_scheme, parent_scheme) in [
                    ("U/U", ApproximationType::MUV, ApproximationType::MUV),
                    ("H/U", ApproximationType::IR, ApproximationType::MUV),
                    ("U/H", ApproximationType::MUV, ApproximationType::IR),
                    ("H/H", ApproximationType::IR, ApproximationType::IR),
                ] {
                    let child = TestNode {
                        lmb: child_lmb.clone(),
                        subgraph: child_subgraph.clone(),
                        dod: 1,
                        scheme: child_scheme,
                    };
                    let parent = TestNode {
                        lmb: graph.loop_momentum_basis.clone(),
                        subgraph: parent_subgraph.clone(),
                        dod: 1,
                        scheme: parent_scheme,
                    };
                    let child_local =
                        uv_limit(&seed, &ctx, &child, &given, &child, &given).unwrap();
                    let completed_child = Full4dCts::from_factorized_local(&child_local);
                    let ordinary_node = TestNode {
                        subgraph: child.subgraph.clone(),
                        lmb: child.lmb.clone(),
                        dod: child.dod,
                        scheme: ApproximationType::MUV,
                    };
                    let ordinary_local = uv_limit(&seed, &ctx, &ordinary_node, &given, &ordinary_node, &given).unwrap();
                    let ordinary_child = Full4dCts::from_factorized_local(&ordinary_local);
                    assert!(!ordinary_child.0.has_soft_ancestry);
                    assert_eq!(completed_child.0.has_soft_ancestry, child_scheme == ApproximationType::IR);

                    if child_scheme == ApproximationType::IR
                        && parent_scheme == ApproximationType::MUV
                    {
                        let error = uv_limit(
                            &completed_child,
                            &ctx,
                            &parent,
                            &child,
                            &parent,
                            &child,
                        )
                        .expect_err(
                            "a positive-degree U parent over an H child must be rejected",
                        );
                        let message = error.to_string();
                        assert!(message.contains("positive-degree MUV component"), "{message}");
                        assert!(
                            message.contains("deferred finite scheme-change/integrated policy"),
                            "{message}"
                        );
                        continue;
                    }

                    let pure_parent_route = compatible_lmb(&ctx, &parent, &given, &seed)
                        .unwrap_or_else(|error| panic!("{label} pure-parent route: {error}"));
                    let nested_parent_route =
                        compatible_lmb(&ctx, &parent, &child, &completed_child)
                            .unwrap_or_else(|error| panic!("{label} nested-parent route: {error}"));
                    assert_eq!(
                        pure_parent_route, nested_parent_route,
                        "{label}: pure and nested parent families selected different routes"
                    );
                    assert!(
                        pure_parent_route
                            .loop_edges
                            .iter()
                            .any(|edge| *edge == EdgeIndex(2)),
                        "{label}: the parent route changed the child boundary carrier e2"
                    );
                    if parent_scheme == ApproximationType::MUV {
                        let affine: Vec<_> = completed_child
                            .0.inherited_loop_edges
                            .iter()
                            .copied()
                            .filter(|edge| parent.subgraph().includes(&graph[edge].1))
                            .filter(|edge| {
                                graph.loop_momentum_basis.edge_signatures[*edge]
                                    .external
                                    .iter()
                                    .any(|sign| sign.is_sign())
                            })
                            .collect();
                        assert!(
                            affine.is_empty(),
                            "{label}: MUV-parent route would require affine shifts of {affine:?}"
                        );
                    }

                    // Keep exactly the four signed forest families separate:
                    // I, -K_child I, -K_parent I, and +K_parent K_child I.
                    let f_empty = grow(&seed, &ctx, &parent, &given).unwrap();
                    let f_child = grow(&completed_child, &ctx, &parent, &child).unwrap();
                    let f_parent = uv_limit(&seed, &ctx, &parent, &given, &parent, &given).unwrap();
                    let f_nested =
                        uv_limit(&completed_child, &ctx, &parent, &child, &parent, &child).unwrap();
                    let f_nested_reference =
                        uv_limit(&ordinary_child, &ctx, &parent, &child, &parent, &child).unwrap();
                    let route = pure_parent_route;
                    let normalize = |atom: &Atom| {
                        let routed = GS.erase_uv_momentum_provenance(atom)
                            .replace(function!(GS.ct_marker, W_.a_))
                            .with(Atom::one())
                            .replace_multiple(graph.uv_wrapped_replacement(
                                &parent_subgraph,
                                &route,
                                &[W_.x___],
                            ));
                        finalize(routed)
                            .replace(function!(GS.den, W_.a_, W_.mom_, W_.mass_, W_.prop_))
                            .with(Atom::var(W_.prop_))
                            .replace(GS.dim)
                            .with(Atom::num(4))
                            .together()
                            .cancel()
                            .expand()
                    };
                    let families = [
                        normalize(&f_empty),
                        normalize(&f_child),
                        normalize(f_parent.atom()),
                        normalize(f_nested.atom()),
                    ];
                    for (family_label, family) in ["empty", "child", "parent", "child;parent"]
                        .into_iter()
                        .zip(&families)
                    {
                        assert!(
                            !family.is_zero(),
                            "{label}: the production {family_label} family must be nonzero"
                        );
                        assert!(
                            !paper_has_pattern(family, Atom::var(GS.rescale)),
                            "{label}: {family_label} retains the Taylor parameter"
                        );
                    }
                    // Linearity lets us isolate the H-U child refinement, but
                    // the production parent still applies its single full
                    // component operation to that completed contribution.
                    let nested_input = grow(&completed_child, &ctx, &parent, &child).unwrap();
                    let nested_reference_input =
                        grow(&ordinary_child, &ctx, &parent, &child).unwrap();
                    let nested_refinement_input = &nested_input - &nested_reference_input;
                    let direct_refinement = match parent_scheme {
                        ApproximationType::MUV => t_raw(
                            &nested_refinement_input,
                            &ctx,
                            &parent,
                            &child,
                            &route,
                        )
                        .unwrap(),
                        ApproximationType::IR => hat_t_raw_with_provenance(
                            &nested_refinement_input,
                            &ctx,
                            &parent,
                            &child,
                            &route,
                        )
                        .unwrap()
                        .0,
                        _ => unreachable!(),
                    };
                    let direct_refinement = normalize(&-finalize(direct_refinement));
                    let production_refinement =
                        normalize(f_nested.atom()) - normalize(f_nested_reference.atom());
                    let full_parent_difference =
                        (production_refinement - direct_refinement).together().cancel().expand();
                    assert!(
                        full_parent_difference.is_zero(),
                        "{label}: the nested H-U family is not the full parent image of the completed child refinement"
                    );

                    let child_edges: Vec<_> = child.lmb.loop_edges.iter().copied().collect();
                    let reduced_parent_edges: Vec<_> = route
                        .loop_edges
                        .iter()
                        .copied()
                        .filter(|edge| !child_edges.contains(edge))
                        .collect();
                    assert_eq!(child_edges.len(), 1, "{label}: expected one child loop");
                    assert_eq!(
                        reduced_parent_edges.len(),
                        1,
                        "{label}: expected one reduced-parent loop"
                    );
                    let simultaneous_edges: Vec<_> = route.loop_edges.iter().copied().collect();

                    for (pair_label, remainder, ray) in [
                        (
                            "I-K_child I",
                            &families[0] + &families[1],
                            child_edges.as_slice(),
                        ),
                        (
                            "-K_parent(I-K_child)I",
                            &families[2] + &families[3],
                            child_edges.as_slice(),
                        ),
                        (
                            "-(1-K_parent)K_child I",
                            &families[1] + &families[3],
                            reduced_parent_edges.as_slice(),
                        ),
                    ] {
                        let degree = hard_degree(
                            &remainder,
                            &route,
                            ray,
                            &format!("{label}, {pair_label}"),
                        );
                        if degree <= zero {
                            production_obstructions
                                .push(format!("{label}: {pair_label} -> {degree}"));
                        }
                    }
                    let complete = &families[0] + &families[1] + &families[2] + &families[3];
                    for (ray_label, ray) in [
                        ("child", child_edges.as_slice()),
                        ("reduced-parent", reduced_parent_edges.as_slice()),
                        ("simultaneous", simultaneous_edges.as_slice()),
                    ] {
                        let degree = hard_degree(
                            &complete,
                            &route,
                            ray,
                            &format!("{label}, complete forest, {ray_label}-hard"),
                        );
                        if degree <= zero {
                            production_obstructions.push(format!(
                                "{label}: complete {ray_label} ray -> {degree}"
                            ));
                        }
                    }

                    matrix.push((
                        label,
                        child_scheme,
                        parent_scheme,
                        route,
                        child_edges,
                        reduced_parent_edges,
                        simultaneous_edges,
                        families,
                    ));
                }

                assert!(
                    production_obstructions.is_empty(),
                    "the completed-child parent operation left nested UV obstructions: {production_obstructions:?}"
                );
                assert_eq!(
                    matrix.len(),
                    3,
                    "every supported positive-degree local U/H operator pair must run"
                );

                let by_schemes = |child_scheme, parent_scheme| {
                    matrix
                        .iter()
                        .find(|entry| entry.1 == child_scheme && entry.2 == parent_scheme)
                        .unwrap()
                };
                let uu = by_schemes(ApproximationType::MUV, ApproximationType::MUV);
                let uh = by_schemes(ApproximationType::MUV, ApproximationType::IR);
                let hh = by_schemes(ApproximationType::IR, ApproximationType::IR);
                assert!(
                    !(&uh.7[1] - &hh.7[1]).together().cancel().expand().is_zero(),
                    "the H-child row collapsed to the ordinary U child"
                );
                assert!(
                    !(&uu.7[2] - &uh.7[2]).together().cancel().expand().is_zero(),
                    "the H-parent column collapsed to the ordinary U parent"
                );

                let normalized_zero = |atom: Atom| atom.together().cancel().expand().is_zero();
                let refinements = [
                    (
                        "(H/H)-(U/H), child H-U refinement",
                        hh,
                        &hh.7[1] - &uh.7[1] + &hh.7[3] - &uh.7[3],
                        normalized_zero(&hh.7[0] - &uh.7[0])
                            && normalized_zero(&hh.7[2] - &uh.7[2]),
                    ),
                    (
                        "(U/H)-(U/U), parent H-U refinement",
                        uh,
                        &uh.7[2] - &uu.7[2] + &uh.7[3] - &uu.7[3],
                        normalized_zero(&uh.7[0] - &uu.7[0])
                            && normalized_zero(&uh.7[1] - &uu.7[1]),
                    ),
                ];
                for (label, entry, refinement, unchanged_families_cancel) in refinements {
                    assert!(
                        unchanged_families_cancel,
                        "{label}: a scheme-independent forest family changed"
                    );
                    assert!(
                        !refinement.is_zero(),
                        "{label}: the isolated H-U refinement is vacuous"
                    );
                    assert_eq!(
                        entry.3, uu.3,
                        "{label}: compared refinements use different parent routes"
                    );
                    for (ray_label, ray) in [
                        ("child", entry.4.as_slice()),
                        ("reduced-parent", entry.5.as_slice()),
                        ("simultaneous", entry.6.as_slice()),
                    ] {
                        let degree = hard_degree(
                            &refinement,
                            &entry.3,
                            ray,
                            &format!("{label}, {ray_label}-hard"),
                        );
                        assert!(
                            degree > zero,
                            "{label}: the isolated refinement retains a {ray_label} divergence ({degree})"
                        );
                    }
                }

                // Divergent controls prove that the positive Laurent assertions
                // above do not merely reflect a UV-finite fixture.
                assert!(
                    hard_degree(&uu.7[0], &uu.3, &uu.4, "bare child control") <= zero,
                    "the bare child control is unexpectedly UV finite"
                );
                assert!(
                    hard_degree(&uu.7[0], &uu.3, &uu.6, "bare simultaneous control") <= zero,
                    "the bare simultaneous control is unexpectedly UV finite"
                );
                // The bare graph can be finite with only the reduced-parent
                // coordinate hard because an uncontracted child propagator
                // supplies an extra power.  The child-only forest family below
                // is the non-vacuous divergent control for that ray.
                for entry in [uu, uh] {
                    assert!(
                        hard_degree(
                            &entry.7[2],
                            &entry.3,
                            &entry.4,
                            &format!("{}, pure-parent child-hard control", entry.0),
                        ) <= zero,
                        "{}: the pure-parent child-hard control is unexpectedly finite",
                        entry.0
                    );
                }
                for entry in [uu, hh] {
                    assert!(
                        hard_degree(
                            &entry.7[1],
                            &entry.3,
                            &entry.5,
                            &format!("{}, child-only parent-hard control", entry.0),
                        ) <= zero,
                        "{}: the child-only parent-hard control is unexpectedly finite",
                        entry.0
                    );
                }
            })
            .unwrap();
        test.join().unwrap();
    }

    #[test]
    fn full_parent_u_does_not_project_a_uv_finite_completed_child() {
        test_initialise().unwrap();
        let t = Atom::var(GS.rescale);
        let k = Atom::var(symbol!("nested_compatibility_child_radius"));
        let p = Atom::var(symbol!("nested_compatibility_parent_radius"));
        let child_finite = Atom::one() / (Atom::one() + k.pow(2)).pow(3);
        let linearly_divergent_cograph = &p / (Atom::one() + p.pow(2)).pow(2);
        let product = &child_finite * &linearly_divergent_cograph;

        // Include both four-dimensional loop measures.  The finite child
        // supplies two inverse powers, so the simultaneous full-parent ray is
        // convergent even though the reduced-parent ray remains linearly
        // divergent with the child momentum held fixed.
        let simultaneous = product
            .replace(k.clone())
            .with(&k / &t)
            .replace(p.clone())
            .with(&p / &t)
            / t.pow(8);
        let reduced_parent = product.replace(p.clone()).with(&p / &t) / t.pow(4);
        assert_eq!(
            simultaneous
                .series(GS.rescale, Atom::Zero, 1)
                .unwrap()
                .get_trailing_exponent(),
            Rational::from(1),
        );
        let full_projection = simultaneous
            .series(GS.rescale, Atom::Zero, 0)
            .unwrap()
            .to_atom();
        assert!(
            full_projection.is_zero(),
            "ordinary full U has no divergent coefficient to subtract"
        );
        assert_eq!(
            reduced_parent
                .series(GS.rescale, Atom::Zero, 0)
                .unwrap()
                .get_trailing_exponent(),
            Rational::from(-1),
            "the same finite child coefficient multiplies a reduced-parent divergence"
        );
        let quotient_projection = reduced_parent
            .series(GS.rescale, Atom::Zero, 0)
            .unwrap()
            .to_atom()
            .replace(GS.rescale)
            .with(Atom::one());
        assert!(
            !quotient_projection.is_zero(),
            "the quotient hard-ray diagnostic must distinguish the erroneous production scope"
        );
        let quotient_remainder = (&product - quotient_projection)
            .replace(p.clone())
            .with(&p / &t)
            / t.pow(4);
        assert!(
            quotient_remainder
                .series(GS.rescale, Atom::Zero, 1)
                .unwrap()
                .get_trailing_exponent()
                > 0,
            "the quotient diagnostic projection removes the reduced-parent divergence, but is not a parent forest operator"
        );
    }

    #[test]
    fn full_parent_u_and_h_treat_a_uv_finite_child_refinement_consistently() {
        test_initialise().unwrap();
        let t = Atom::var(GS.rescale);
        let k = Atom::var(symbol!("exact_hu_child_radius"));
        let p = Atom::var(symbol!("exact_hu_parent_radius"));
        let q = Atom::var(symbol!("exact_hu_parent_external"));

        // This is the constant soft coefficient S_0 of a logarithmic child,
        // evaluated once with its physical denominator and once after U S_0
        // introduces the UV denominator.  Hence `child_refinement` is exactly
        // the H-U=S-US part of that child.  Its two denominators differ only by
        // their masses, so it is locally finite in the child hard limit.
        let child_soft = Atom::one() / (Atom::one() + k.pow(2)).pow(2);
        let child_uv_soft = Atom::one() / (Atom::num(4) + k.pow(2)).pow(2);
        let child_refinement = &child_soft - &child_uv_soft;
        let child_hard = child_refinement.replace(k.clone()).with(&k / &t) / t.pow(4);
        assert_eq!(
            child_hard
                .series(GS.rescale, Atom::Zero, 2)
                .unwrap()
                .get_trailing_exponent(),
            Rational::from(2),
            "the exact H-U child refinement must be UV finite"
        );

        // A logarithmically divergent reduced cograph distinguishes the two
        // loop gradings. With both loop radii hard, the finite child lowers the
        // full-parent degree and U_parent has no non-positive coefficient.
        // Holding the child fixed exposes a reduced-parent diagnostic limit,
        // but that limit must not replace the full parent forest operator.
        let cograph = Atom::one() / ((Atom::one() + p.pow(2)) * (Atom::one() + (&p + &q).pow(2)));
        let nested_refinement = &child_refinement * &cograph;
        let full_parent_hard = nested_refinement
            .replace(k.clone())
            .with(&k / &t)
            .replace(p.clone())
            .with(&p / &t)
            / t.pow(8);
        assert_eq!(
            full_parent_hard
                .series(GS.rescale, Atom::Zero, 2)
                .unwrap()
                .get_trailing_exponent(),
            Rational::from(2),
        );
        let full_u_projection = full_parent_hard
            .series(GS.rescale, Atom::Zero, 0)
            .unwrap()
            .to_atom();
        assert!(
            full_u_projection.is_zero(),
            "U_parent must not manufacture a subtraction when its full hard ray is finite"
        );

        let reduced_parent_hard = nested_refinement.replace(p.clone()).with(&p / &t) / t.pow(4);
        assert_eq!(
            reduced_parent_hard
                .series(GS.rescale, Atom::Zero, 0)
                .unwrap()
                .get_trailing_exponent(),
            Rational::from(0),
            "the unprojected H-U insertion has a logarithmic reduced-parent ray"
        );

        let quotient_projection = reduced_parent_hard
            .series(GS.rescale, Atom::Zero, 0)
            .unwrap()
            .to_atom()
            .replace(GS.rescale)
            .with(Atom::one());
        assert!(
            !quotient_projection.is_zero(),
            "the quotient diagnostic must expose the coefficient that an erroneous split would add"
        );

        // For a degree-one H parent, S_0 sets its external momentum q to zero.
        // Its U and US branches both use the full grading and therefore vanish
        // on this UV-finite completed child, while the physical soft branch
        // keeps the child mass and improves the reduced-parent limit.
        let parent_soft = nested_refinement.replace(q.clone()).with(Atom::Zero);
        let full_soft_parent_hard = parent_soft
            .replace(k.clone())
            .with(&k / &t)
            .replace(p.clone())
            .with(&p / &t)
            / t.pow(8);
        let full_soft_projection = full_soft_parent_hard
            .series(GS.rescale, Atom::Zero, 0)
            .unwrap()
            .to_atom();
        assert!(
            full_soft_projection.is_zero(),
            "US must not manufacture a full-parent subtraction on the soft-completed finite child"
        );
        let h_parent_remainder = nested_refinement - parent_soft;
        let h_reduced_parent_hard = h_parent_remainder.replace(p.clone()).with(&p / &t) / t.pow(4);
        assert_eq!(
            h_reduced_parent_hard
                .series(GS.rescale, Atom::Zero, 1)
                .unwrap()
                .get_trailing_exponent(),
            Rational::from(1),
            "the physical parent S branch must improve the reduced-parent limit without a quotient U replay"
        );
    }

    #[test]
    fn full_parent_u0_reexpands_s1_child_and_u1_rejects_it() {
        let test = std::thread::Builder::new()
            .name("full-parent-u-s1-child".to_string())
            .stack_size(64 * 1024 * 1024)
            .spawn(|| {
                test_initialise().unwrap();
                let run = |case: &str, graph: Graph, parent_dod: i32| {
                    let mut child_subgraph = graph.empty_subgraph::<SuBitGraph>();
                    for edge in [EdgeIndex(3), EdgeIndex(4)] {
                        child_subgraph.add(graph[&edge].1);
                    }
                    let external: SuBitGraph = graph.external_filter();
                    let parent_subgraph = graph.full_filter().subtract(&external);
                    assert_eq!(graph.compute_dod(&child_subgraph), 2);
                    assert_eq!(graph.compute_dod(&parent_subgraph), parent_dod);
                    let mut nonnegative_connected = graph
                        .spinneys(&parent_subgraph)
                        .into_iter()
                        .filter(|spinney| !spinney.filter.is_empty())
                        .filter(|spinney| graph.compute_dod(&spinney.filter) >= 0)
                        .filter(|spinney| {
                            graph.underlying.connected_components(&spinney.filter).len() == 1
                        })
                        .map(|spinney| {
                            let mut edges = graph
                                .iter_edges_of(&spinney.filter)
                                .map(|(_, edge, _)| usize::from(edge))
                                .collect::<Vec<_>>();
                            edges.sort();
                            edges
                        })
                        .collect::<Vec<_>>();
                    nonnegative_connected.sort();
                    assert_eq!(
                        nonnegative_connected,
                        vec![vec![2, 3, 4, 5, 6], vec![3, 4]],
                        "the S1 fixture must have only the child and parent UV components"
                    );

                    let child_lmb = graph
                        .underlying
                        .try_compatible_sub_lmb(
                            &child_subgraph,
                            graph
                                .dummy_less_full_crown(&child_subgraph)
                                .subtract(&external),
                            &graph.loop_momentum_basis,
                        )
                        .unwrap();
                    assert!(child_lmb.ext_edges.contains(&EdgeIndex(2)));
                    let given = root_node(&graph);
                    let settings = UVgenerationSettings {
                        generate_integrated: false,
                        ..Default::default()
                    };
                    let ctx = UVCtx::new(&graph, &settings);
                    let seed = Full4dCts::root();
                    let child_h = TestNode {
                        subgraph: child_subgraph.clone(),
                        lmb: child_lmb.clone(),
                        dod: 2,
                        scheme: ApproximationType::IR,
                    };
                    let child_u = TestNode {
                        subgraph: child_subgraph.clone(),
                        lmb: child_lmb,
                        dod: 2,
                        scheme: ApproximationType::MUV,
                    };
                    let parent = TestNode {
                        subgraph: parent_subgraph.clone(),
                        lmb: graph.loop_momentum_basis.clone(),
                        dod: parent_dod,
                        scheme: ApproximationType::MUV,
                    };
                    let child_h_local =
                        uv_limit(&seed, &ctx, &child_h, &given, &child_h, &given).unwrap();
                    let child_u_local =
                        uv_limit(&seed, &ctx, &child_u, &given, &child_u, &given).unwrap();
                    let completed_h = Full4dCts::from_factorized_local(&child_h_local);
                    let completed_u = Full4dCts::from_factorized_local(&child_u_local);
                    let route = compatible_lmb(&ctx, &parent, &given, &seed).unwrap();
                    assert_eq!(
                        route,
                        compatible_lmb(&ctx, &parent, &child_h, &completed_h).unwrap(),
                        "the pure and nested S1 parent families must use one route"
                    );
                    assert_eq!(
                        route,
                        compatible_lmb(&ctx, &parent, &child_u, &completed_u).unwrap(),
                        "the U-reference and S1-refined children must use one route"
                    );
                    let child_edges = child_h.lmb.loop_edges.iter().copied().collect::<Vec<_>>();
                    let reduced_parent_edges = route
                        .loop_edges
                        .iter()
                        .copied()
                        .filter(|edge| !child_edges.contains(edge))
                        .collect::<Vec<_>>();
                    assert_eq!(child_edges, vec![EdgeIndex(3)]);
                    assert_eq!(reduced_parent_edges, vec![EdgeIndex(2)]);

                    let normalize = |atom: &Atom| {
                        let routed = GS.erase_uv_momentum_provenance(atom)
                            .replace(function!(GS.ct_marker, W_.a_))
                            .with(Atom::one())
                            .replace_multiple(graph.uv_wrapped_replacement(
                                &parent_subgraph,
                                &route,
                                &[W_.x___],
                            ));
                        finalize(routed)
                            .replace(function!(GS.den, W_.a_, W_.mom_, W_.mass_, W_.prop_))
                            .with(Atom::var(W_.prop_))
                            .replace(GS.dim)
                            .with(Atom::num(4))
                    };

                    let f_empty = normalize(&grow(&seed, &ctx, &parent, &given).unwrap());
                    let f_child_h =
                        normalize(&grow(&completed_h, &ctx, &parent, &child_h).unwrap());
                    let f_child_u =
                        normalize(&grow(&completed_u, &ctx, &parent, &child_u).unwrap());
                    let f_parent = normalize(
                        uv_limit(&seed, &ctx, &parent, &given, &parent, &given)
                            .unwrap()
                            .atom(),
                    );
                    let f_nested_u = normalize(
                        uv_limit(&completed_u, &ctx, &parent, &child_u, &parent, &child_u)
                            .unwrap()
                            .atom(),
                    );
                    let delta_child = completed_h.atom() - completed_u.atom();
                    assert!(
                        paper_has_pattern(&delta_child, GS.emr_mom(EdgeIndex(2), W_.x___)),
                        "the DOD-2 S1 refinement must retain its boundary-momentum polynomial"
                    );
                    if parent_dod > 0 {
                        let error = uv_limit(
                            &completed_h,
                            &ctx,
                            &parent,
                            &child_h,
                            &parent,
                            &child_h,
                        )
                        .expect_err(
                            "a positive-degree U parent over an S1 child must be rejected",
                        );
                        let message = error.to_string();
                        assert!(message.contains("positive-degree MUV component"), "{message}");
                        assert!(
                            message.contains("deferred finite scheme-change/integrated policy"),
                            "{message}"
                        );
                        return;
                    }

                    // At degree zero U_0=H_0, so the ordinary parent can act
                    // once on the exact completed child without a deferred
                    // finite scheme-change term.
                    let f_nested_h_production = normalize(
                        uv_limit(&completed_h, &ctx, &parent, &child_h, &parent, &child_h)
                            .unwrap()
                            .atom(),
                    );
                    let co_graph = grow(&seed, &ctx, &parent, &child_h).unwrap();
                    let production_nested_delta = &f_nested_h_production - &f_nested_u;
                    let direct_nested_delta = normalize(
                        &-finalize(
                            t_raw(
                                &(&delta_child * &co_graph),
                                &ctx,
                                &parent,
                                &child_h,
                                &route,
                            )
                            .unwrap(),
                        ),
                    );
                    paper_assert_normalized_zero(
                        &format!("{case}: U_parent acts fully on the complete (H-U)_child"),
                        &production_nested_delta - &direct_nested_delta,
                    );

                    let production = [
                        f_empty.clone(),
                        f_child_h.clone(),
                        f_parent.clone(),
                        f_nested_h_production.clone(),
                    ];
                    for (family_label, family) in ["empty", "child", "parent", "child;parent"]
                        .into_iter()
                        .zip(&production)
                    {
                        assert!(!family.is_zero(), "the {family_label} family is vacuous");
                    }
                    let complete =
                        &production[0] + &production[1] + &production[2] + &production[3];
                    // Keep the four forest families distinct and verify the
                    // linear decomposition of K_parent K_child with one full
                    // parent action on the completed child.
                    let direct_recurrence = &f_empty
                        + &f_child_u
                        + &f_parent
                        + &f_nested_u
                        + (&f_child_h - &f_child_u)
                        + &direct_nested_delta;
                    paper_assert_normalized_zero(
                        &format!(
                            "{case}: four families equal the completed-child nested R-operation"
                        ),
                        &complete - &direct_recurrence,
                    );
                    let mut obstructions = Vec::new();
                    for (pair_label, remainder, ray) in [
                        (
                            "I-K_child I",
                            &production[0] + &production[1],
                            child_edges.as_slice(),
                        ),
                        (
                            "-K_parent(I-K_child)I",
                            &production[2] + &production[3],
                            child_edges.as_slice(),
                        ),
                        (
                            "-(1-K_parent)K_child I",
                            &production[1] + &production[3],
                            reduced_parent_edges.as_slice(),
                        ),
                    ] {
                        let degree = nested_hard_degree(
                            &graph,
                            &parent_subgraph,
                            &remainder,
                            &route,
                            ray,
                            pair_label,
                        );
                        if degree <= 0 {
                            obstructions.push(format!("{pair_label}: {degree}"));
                        }
                    }
                    let simultaneous_edges = route.loop_edges.iter().copied().collect::<Vec<_>>();
                    for (ray_label, ray) in [
                        ("child", child_edges.as_slice()),
                        ("reduced-parent", reduced_parent_edges.as_slice()),
                        ("simultaneous", simultaneous_edges.as_slice()),
                    ] {
                        let degree = nested_hard_degree(
                            &graph,
                            &parent_subgraph,
                            &complete,
                            &route,
                            ray,
                            &format!("S1 completed-child forest {ray_label}"),
                        );
                        if degree <= 0 {
                            obstructions.push(format!("complete {ray_label}: {degree}"));
                        }
                    }
                    assert!(
                        obstructions.is_empty(),
                        "the {case} full parent operation left nested UV obstructions: {obstructions:?}"
                    );
                };
                run(
                    "U0",
                    dot!(
                        digraph nested_dod2_child_dod0_parent {
                            edge [particle=scalar_1 num=1]
                            node [num=1]
                            ext [style=invis]
                            ext -> v1:0 [id=0]
                            v4:1 -> ext [id=1]
                            v1 -> v2 [id=2 lmb_id=1]
                            v2 -> v3 [id=3 lmb_id=0 num="Q(3,spenso::cind(0))^2"]
                            v2 -> v3 [id=4]
                            v3 -> v4 [id=5]
                            v1 -> v4 [id=6]
                        },
                        "scalars"
                    )
                    .unwrap(),
                    0,
                );
                run(
                    "U1",
                    dot!(
                        digraph nested_dod2_child_dod1_parent {
                            edge [particle=scalar_1 num=1]
                            node [num=1]
                            ext [style=invis]
                            ext -> v1:0 [id=0]
                            v4:1 -> ext [id=1]
                            v1 -> v2 [id=2 lmb_id=1 num="Q(2,spenso::cind(0))"]
                            v2 -> v3 [id=3 lmb_id=0 num="Q(3,spenso::cind(0))^2"]
                            v2 -> v3 [id=4]
                            v3 -> v4 [id=5]
                            v1 -> v4 [id=6]
                        },
                        "scalars"
                    )
                    .unwrap(),
                    1,
                );
            })
            .unwrap();
        test.join().unwrap();
    }

    #[test]
    fn uv_limit_negates_and_marks_the_completed_hat_atom_once() {
        test_initialise().unwrap();
        let graph = scalar_two_point_graph();
        let given = root_node(&graph);
        let mut current = two_point_node(&graph, ApproximationType::IR);
        current.dod = 2;
        let settings = UVgenerationSettings {
            generate_integrated: false,
            add_marker: true,
            keep_marker: true,
            ..Default::default()
        };
        let ctx = UVCtx::new(&graph, &settings);
        let root = Full4dCts::root();
        let lmb = compatible_lmb(&ctx, &current, &given, &root).unwrap();
        let grown = grow(&root, &ctx, &current, &given).unwrap();
        let raw_hat = hat_t_raw(&grown, &ctx, &current, &given, &lmb).unwrap();
        let expected_unmarked = -finalize(raw_hat);
        let actual = uv_limit(&root, &ctx, &current, &given, &current, &given).unwrap();
        let stripped = actual
            .atom()
            .replace(function!(GS.ct_marker, W_.a_))
            .with(Atom::one());

        paper_assert_normalized_zero(
            "the forest sign is applied to the completed H atom",
            stripped - expected_unmarked,
        );
        let histories: Vec<_> = actual
            .atom()
            .replace(function!(GS.ct_marker, W_.a_))
            .match_iter()
            .map(|matched| matched.get(&W_.a_).unwrap().clone())
            .collect();
        assert!(!histories.is_empty(), "the enabled marker must be present");
        for history in histories {
            assert_eq!(
                history
                    .replace(function!(GS.uv_approx, W_.a_))
                    .match_iter()
                    .count(),
                1,
                "each completed H term must receive exactly one approximation marker: {history}"
            );
        }
    }

    #[test]
    fn ir_local_operator_is_independent_of_integrated_generation() {
        test_initialise().unwrap();
        let graph = scalar_two_point_graph();
        let given = root_node(&graph);
        let current = two_point_node(&graph, ApproximationType::IR);
        let local = [false, true].map(|generate_integrated| {
            let settings = UVgenerationSettings {
                generate_integrated,
                ..Default::default()
            };
            uv_limit(
                &Full4dCts::root(),
                &UVCtx::new(&graph, &settings),
                &current,
                &given,
                &current,
                &given,
            )
            .unwrap()
        });
        assert!(!local[0].atom().is_zero());
        assert!(local[0].0.has_soft_ancestry);
        assert_eq!(local[0], local[1]);
    }

    #[test]
    fn integrated_massless_gauge_soft_jet_vanishes_but_completed_hat_does_not() -> Result<()> {
        test_initialise()?;
        let graph: Graph =
            include_str!("../../../../../tests/resources/graphs/local_ir_gluon_self_energy.dot")
                .into_graph(&crate::utils::load_generic_model("sm"))?;
        assert!(graph[EdgeIndex(2)].mass_atom().is_zero());
        let mut current = two_point_node(&graph, ApproximationType::IR);
        current.dod = 2;
        let given = root_node(&graph);
        let settings = UVgenerationSettings {
            generate_integrated: true,
            ..Default::default()
        };
        let ctx = UVCtx::new(&graph, &settings);
        let root = Full4dCts::root();
        let lmb = compatible_lmb(&ctx, &current, &given, &root)?;
        let grown = grow(&root, &ctx, &current, &given)?;
        let soft = tilde_t_raw(&grown, &ctx, &current, &given, &lmb)?;
        assert!(!soft.is_zero(), "the local soft jet must be nontrivial");
        let ordinary = t_raw(&grown, &ctx, &current, &given, &lmb)?;
        let overlap = t_raw(&soft, &ctx, &current, &given, &lmb)?;
        let local = uv_limit(&root, &ctx, &current, &given, &current, &given)?;
        assert_eq!(local.active_sectors().len(), 1);
        let mut vakint_settings = settings.vakint.true_settings();
        vakint_settings.number_of_terms_in_epsilon_expansion = 3;
        let integrator = crate::uv::approx::integrated::Integrated::new(
            crate::utils::vakint()?,
            &vakint_settings,
        );
        let integrate_branch = |raw: Atom| -> Result<IntegratedCts> {
            let mut sector = local.active_sectors()[0].clone();
            sector.atom = -finalize(raw);
            integrator.run(
                &Local4dCts(FourDSectors::new(vec![sector], Vec::new())),
                &ctx,
                &current,
                &given,
                &current,
                &given,
            )
        };
        let soft = integrate_branch(soft)?;
        let ordinary = integrate_branch(ordinary)?;
        let overlap = integrate_branch(overlap)?;
        let completed = integrator.run(&local, &ctx, &current, &given, &current, &given)?;
        assert!(soft.physical_pole_atom().is_zero());
        assert!(soft.physical_finite_counterterm_atom().is_zero());
        assert!(!completed.physical_pole_atom().is_zero());
        assert!(
            !overlap.physical_pole_atom().is_zero()
                || !overlap.physical_finite_counterterm_atom().is_zero(),
            "the integrated US overlap must contribute to the completed counterterm"
        );
        // Eq. (5.6) removes the integrated S contribution, not H or its US
        // overlap. All branches here carry the same production forest sign.
        for difference in [
            completed.physical_pole_atom() - ordinary.physical_pole_atom()
                + overlap.physical_pole_atom(),
            completed.physical_finite_counterterm_atom()
                - ordinary.physical_finite_counterterm_atom()
                + overlap.physical_finite_counterterm_atom(),
        ] {
            assert!(difference.together().cancel().is_zero());
        }
        let recursive = Full4dCts::recursion_input(
            &local,
            &completed,
            ApproximationType::IR,
            false,
            current.lmb(),
        )?;
        assert_eq!(
            recursive.atom(),
            &(local.atom() + completed.finite_counterterm_atom())
        );
        assert_eq!(recursive.0.active, local.0.active);
        assert!(recursive.0.has_soft_ancestry);
        let [finite] = recursive.0.recursive_completion.as_slice() else {
            panic!("the nonzero integrated soft completion must occur exactly once");
        };
        assert_eq!(finite.atom, completed.finite_counterterm_atom());
        assert_eq!(finite.frozen_lmbs, vec![current.lmb.clone()]);
        Ok(())
    }

    #[test]
    fn integrated_massive_scalar_soft_jet_retains_nonzero_physical_mass() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph massive_scalar_soft_tadpole {
            edge [particle=H num=1]
            node [num=1]
            ext [style=invis]
            ext -> a [id=0]
            a -> a [id=1]
            a -> ext [id=2]
        })?;
        let physical_mass = graph[EdgeIndex(1)].mass_atom();
        assert!(!physical_mass.is_zero());
        assert_ne!(physical_mass, Atom::var(GS.m_uv_vacuum));
        let mut current = two_point_node(&graph, ApproximationType::IR);
        current.dod = 2;
        let given = root_node(&graph);
        let settings = UVgenerationSettings {
            generate_integrated: true,
            ..Default::default()
        };
        let ctx = UVCtx::new(&graph, &settings);
        let root = Full4dCts::root();
        let lmb = compatible_lmb(&ctx, &current, &given, &root)?;
        let grown = grow(&root, &ctx, &current, &given)?;
        let soft = tilde_t_raw(&grown, &ctx, &current, &given, &lmb)?;
        let soft = Local4dCts(FourDSectors::new(
            vec![FourDSector::new(
                -finalize(soft),
                vec![(current.subgraph.clone(), current.subgraph.clone(), lmb)],
                Vec::new(),
            )],
            Vec::new(),
        ));
        let mut vakint_settings = settings.vakint.true_settings();
        vakint_settings.number_of_terms_in_epsilon_expansion = 3;
        let integrated = crate::uv::approx::integrated::Integrated::new(
            crate::utils::vakint()?,
            &vakint_settings,
        )
        .run(&soft, &ctx, &current, &given, &current, &given)?;
        let pole = integrated.physical_pole_atom();
        assert!(
            !pole.is_zero(),
            "a massive scalar soft tadpole is not scaleless"
        );
        assert!(paper_has_pattern(&pole, physical_mass.clone()));
        assert!(!pole.contains_symbol(GS.m_uv_vacuum));
        // The one-loop tadpole pole has mass dimension two. Evaluating at
        // m_phys=2*mUV also distinguishes retained m_phys² from a second square.
        let at_uv = pole.replace(physical_mass.clone()).with(GS.m_uv_vacuum);
        let at_twice_uv = pole
            .replace(physical_mass)
            .with(Atom::num(2) * GS.m_uv_vacuum);
        assert!(
            (at_twice_uv - Atom::num(4) * at_uv)
                .together()
                .cancel()
                .is_zero()
        );
        Ok(())
    }

    #[test]
    fn integrated_mixed_mass_soft_bubble_matches_massive_tadpoles() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph mixed_mass_scalar_soft_bubble {
            edge [particle=H num=1]
            node [num=1]
            ext [style=invis]
            ext -> a:0 [id=0]
            b:1 -> ext [id=1]
            a -> b [id=2 mass=0]
            a -> b [id=3]
        })?;
        let physical_mass = graph[EdgeIndex(3)].mass_atom();
        assert!(graph[EdgeIndex(2)].mass_atom().is_zero());
        assert!(!physical_mass.is_zero());
        assert_ne!(physical_mass, Atom::var(GS.m_uv_vacuum));
        let mass_squared = physical_mass.pow(2);
        let current = two_point_node(&graph, ApproximationType::IR);
        let given = root_node(&graph);
        let settings = UVgenerationSettings {
            generate_integrated: true,
            ..Default::default()
        };
        let ctx = UVCtx::new(&graph, &settings);
        let momentum = function!(
            GS.emr_mom,
            usize::from(current.lmb.loop_edges[LoopIndex(0)])
        );
        let massless = function!(GS.den, 2, &momentum, Atom::Zero);
        let massive = function!(GS.den, 3, -&momentum, &mass_squared);
        let mut vakint_settings = settings.vakint.true_settings();
        vakint_settings.number_of_terms_in_epsilon_expansion = 3;
        let integrator = crate::uv::approx::integrated::Integrated::new(
            crate::utils::vakint()?,
            &vakint_settings,
        );
        let mut tadpoles = Vec::new();
        let mut bubbles = Vec::new();
        for power in [1, 2] {
            for (owners, denominator, results) in [
                (
                    graph.get_edge_subgraph(EdgeIndex(3)),
                    massive.pow(power),
                    &mut tadpoles,
                ),
                (
                    current.subgraph.clone(),
                    &massless * massive.pow(power),
                    &mut bubbles,
                ),
            ] {
                // Both branches carry the production forest minus. The massive
                // tadpole contracts the omitted massless line, while the mixed
                // bubble retains both original owners and the same loop measure.
                let local = Local4dCts(FourDSectors::new(
                    vec![FourDSector::new(
                        -denominator.pow(-1),
                        vec![(owners, current.subgraph.clone(), current.lmb.clone())],
                        Vec::new(),
                    )],
                    Vec::new(),
                ));
                results.push(integrator.run(&local, &ctx, &current, &given, &current, &given)?);
            }
        }
        // In dimensional regularization the massless tadpole vanishes, so
        // [1/(q²(q²-M²))] = T_1/M² and its derivative in M² is
        // [1/(q²(q²-M²)²)] = T_2/M² - T_1/M⁴. Check both Laurent projections;
        // finite_counterterm_atom includes its own additive forest sign.
        for finite in [false, true] {
            let tadpoles = tadpoles
                .iter()
                .map(|integrated| {
                    if finite {
                        integrated.physical_finite_counterterm_atom()
                    } else {
                        integrated.physical_pole_atom()
                    }
                })
                .collect::<Vec<_>>();
            let bubbles = bubbles
                .iter()
                .map(|integrated| {
                    if finite {
                        integrated.physical_finite_counterterm_atom()
                    } else {
                        integrated.physical_pole_atom()
                    }
                })
                .collect::<Vec<_>>();
            for difference in [
                &bubbles[0] - &tadpoles[0] / &mass_squared,
                &bubbles[1] - &tadpoles[1] / &mass_squared + &tadpoles[0] / mass_squared.pow(2),
            ] {
                assert!(difference.together().cancel().is_zero());
            }
            assert!(
                bubbles
                    .iter()
                    .all(|atom| !atom.contains_symbol(GS.m_uv_vacuum))
            );
        }
        assert!(!bubbles[0].physical_pole_atom().is_zero());
        assert!(bubbles[1].physical_pole_atom().is_zero());
        let finite = bubbles[1].physical_finite_counterterm_atom();
        assert!(!finite.is_zero());
        assert!(paper_has_pattern(&finite, physical_mass));
        Ok(())
    }

    #[test]
    fn integrated_serial_mixed_mass_sunset_matches_contracted_topologies() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph serial_mixed_mass_sunset {
            edge [particle=H num=1]
            node [num=1]
            ext [style=invis]
            ext -> a:0 [id=0]
            a:1 -> ext [id=1]
            a -> c [id=2 mass=0]
            c -> b [id=3 lmb_id=0]
            b -> a [id=4 lmb_id=1]
            b -> a [id=5 mass=0]
        })?;
        let physical_mass = graph[EdgeIndex(3)].mass_atom();
        assert_eq!(physical_mass, graph[EdgeIndex(4)].mass_atom());
        assert!(!physical_mass.is_zero());
        assert_ne!(physical_mass, Atom::var(GS.m_uv_vacuum));
        let mass_squared = physical_mass.pow(2);
        let current = two_point_node(&graph, ApproximationType::IR);
        assert_eq!(current.lmb.loop_edges.len(), 2);
        let given = root_node(&graph);
        let settings = UVgenerationSettings {
            generate_integrated: true,
            ..Default::default()
        };
        let ctx = UVCtx::new(&graph, &settings);
        let mut vakint_settings = settings.vakint.true_settings();
        vakint_settings.number_of_terms_in_epsilon_expansion = 3;
        let integrator = crate::uv::approx::integrated::Integrated::new(
            crate::utils::vakint()?,
            &vakint_settings,
        );
        let [q, l] = [3, 4].map(|edge| function!(GS.emr_mom, edge));
        let massless_q = function!(GS.den, 2, &q, Atom::Zero);
        let massive_q = function!(GS.den, 3, &q, &mass_squared);
        let common =
            function!(GS.den, 4, &l, &mass_squared) * function!(GS.den, 5, &q - &l, Atom::Zero);
        let coefficient = parse!("(a+b)*(c+d)");
        let mut integrals = Vec::new();
        for (owners, denominator) in [
            (
                current.subgraph.clone(),
                &massless_q * massive_q.pow(2) * &common,
            ),
            (
                current
                    .subgraph
                    .subtract(&graph.get_edge_subgraph(EdgeIndex(2))),
                massive_q.pow(2) * &common,
            ),
            (
                current
                    .subgraph
                    .subtract(&graph.get_edge_subgraph(EdgeIndex(2))),
                &massive_q * &common,
            ),
            (
                current
                    .subgraph
                    .subtract(&graph.get_edge_subgraph(EdgeIndex(3))),
                &massless_q * &common,
            ),
        ] {
            // The numerator stays factorized. The reference sunsets contract
            // one serial owner while retaining the full source scope and both
            // loop coordinates; every term carries the same forest minus.
            let local = Local4dCts(FourDSectors::new(
                vec![FourDSector::new(
                    -&coefficient / denominator,
                    vec![(owners, current.subgraph.clone(), current.lmb.clone())],
                    Vec::new(),
                )],
                Vec::new(),
            ));
            integrals.push(integrator.run(&local, &ctx, &current, &given, &current, &given)?);
        }
        // The repeated q denominator has the exact decomposition
        // 1/[q²(q²-M²)²] = 1/[M²(q²-M²)²] - 1/[M⁴(q²-M²)] + 1/[M⁴q²].
        // The remaining l and q-l denominators are unchanged, producing two
        // MM0 sunsets (one dotted) and one M00 sunset with a common measure.
        for finite in [false, true] {
            let values = integrals
                .iter()
                .map(|integrated| {
                    if finite {
                        integrated.physical_finite_counterterm_atom()
                    } else {
                        integrated.physical_pole_atom()
                    }
                })
                .collect::<Vec<_>>();
            let expected = &values[1] / &mass_squared - &values[2] / mass_squared.pow(2)
                + &values[3] / mass_squared.pow(2);
            let difference = &values[0] - expected;
            assert!(
                difference.together().cancel().is_zero(),
                "serial mixed-mass sunset disagrees with contracted topologies for finite={finite}: {difference}",
            );
            assert!(!values[0].is_zero());
            assert!(paper_has_pattern(&values[0], physical_mass.clone()));
            assert!(
                values
                    .iter()
                    .all(|value| !value.contains_symbol(GS.m_uv_vacuum))
            );
        }
        Ok(())
    }

    #[test]
    fn integrated_independent_mixed_mass_tensors_keep_laurent_cross_terms() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph independent_mixed_mass_tensors {
            edge [particle=H num=1]
            node [num=1]
            ext [style=invis]
            ext -> a:0 [id=0]
            a:1 -> ext [id=1]
            a -> a [id=2 lmb_id=0 mass=2]
            a -> a [id=3 lmb_id=1 mass=3]
        })?;
        let current = two_point_node(&graph, ApproximationType::IR);
        assert_eq!(current.lmb.loop_edges.len(), 2);
        let given = root_node(&graph);
        let momenta = [2, 3].map(|edge| {
            function!(
                GS.emr_mom,
                edge,
                Minkowski {}.new_rep(GS.dim).to_symbolic([])
            )
        });
        let denominators = [2, 3].map(|edge| {
            function!(
                GS.den,
                edge,
                function!(GS.emr_mom, edge),
                graph[EdgeIndex(edge)].mass_atom().pow(2)
            )
        });
        // This diagnostic integration input couples loop vectors directly;
        // it tests tensor projection, with no claim about forest-local owners.
        let coefficient = parse!("(a+b)*(c+d)");
        let coupled_numerator =
            &coefficient * function!(SPENSO_TAG.dot, &momenta[0], &momenta[1]).pow(2);
        let mut mode_results = Vec::new();
        for project_onto_tensor_integrals in [true, false] {
            let settings = UVgenerationSettings {
                generate_integrated: true,
                project_integrated_uv_cts_onto_tensor_integrals: project_onto_tensor_integrals,
                ..Default::default()
            };
            let ctx = UVCtx::new(&graph, &settings);
            let mut vakint_settings = settings.vakint.true_settings();
            vakint_settings.project_onto_tensor_integrals = project_onto_tensor_integrals;
            // Three terms retain the two-loop finite coefficient and the one-loop
            // O(epsilon) coefficients that multiply the other factor's UV pole.
            vakint_settings.number_of_terms_in_epsilon_expansion = 3;
            let integrator = crate::uv::approx::integrated::Integrated::new(
                crate::utils::vakint()?,
                &vakint_settings,
            );
            let mut expansions = Vec::new();
            for (subgraph, atom) in [
                (
                    current.subgraph.clone(),
                    -&coupled_numerator / (denominators[0].pow(3) * denominators[1].pow(3)),
                ),
                (
                    graph.get_edge_subgraph(EdgeIndex(2)),
                    -function!(SPENSO_TAG.dot, &momenta[0], &momenta[0]) / denominators[0].pow(3),
                ),
                (
                    graph.get_edge_subgraph(EdgeIndex(3)),
                    -function!(SPENSO_TAG.dot, &momenta[1], &momenta[1]) / denominators[1].pow(3),
                ),
            ] {
                let current = TestNode {
                    lmb: graph.lmb_of(&subgraph),
                    subgraph,
                    dod: 0,
                    scheme: ApproximationType::IR,
                };
                // Preserve the coupled numerator as one factor while the vacuum
                // denominator separates into independent physical-mass factors.
                let local = Local4dCts(FourDSectors::new(
                    vec![FourDSector::new(
                        atom,
                        vec![(
                            current.subgraph.clone(),
                            current.subgraph.clone(),
                            current.lmb.clone(),
                        )],
                        Vec::new(),
                    )],
                    Vec::new(),
                ));
                let integrated =
                    integrator.run(&local, &ctx, &current, &given, &current, &given)?;
                expansions.push(
                    (integrated.physical_pole_atom()
                        - integrated.physical_finite_counterterm_atom())
                    .series(GS.dim_epsilon, Atom::Zero, 3)?,
                );
            }
            let [actual, first, second] = expansions.as_slice() else {
                unreachable!()
            };
            let [a1, b1, c1] = [-1_i64, 0, 1].map(|power| {
                first
                    .coefficient(power.into())
                    .expect("requested coefficient is within series precision")
            });
            let [a2, b2, c2] = [-1_i64, 0, 1].map(|power| {
                second
                    .coefficient(power.into())
                    .expect("requested coefficient is within series precision")
            });
            assert!(!a1.is_zero() && !a2.is_zero() && !c1.is_zero() && !c2.is_zero());

            // I[k^mu k^nu] = g^{mu nu} I[k^2]/D implies
            // I[(k1.k2)^2] = I[k1^2] I[k2^2]/D. Each sector above has the forest
            // minus, so the full integral is minus the product of the radial ones.
            // For D=4-2*epsilon, 1/D = 1/4 + epsilon/8 + epsilon^2/16 + ... .
            let double_pole = &a1 * &a2;
            let simple_pole = &a1 * &b2 + &b1 * &a2;
            let finite = &a1 * &c2 + &b1 * &b2 + &c1 * &a2;
            let expected = [
                -&coefficient * &double_pole / Atom::num(4),
                -&coefficient * (&simple_pole / Atom::num(4) + &double_pole / Atom::num(8)),
                -&coefficient
                    * (finite / Atom::num(4)
                        + &simple_pole / Atom::num(8)
                        + &double_pole / Atom::num(16)),
            ];
            assert!(
                !(&expected[2] + &coefficient * &b1 * &b2 / Atom::num(4))
                    .together()
                    .cancel()
                    .is_zero(),
                "the finite oracle must detect premature componentwise finite projection",
            );
            let actual_coefficients = (-2_i64..=0)
                .map(|power| {
                    actual
                        .coefficient(power.into())
                        .expect("requested coefficient is within series precision")
                })
                .collect::<Vec<_>>();
            for ((power, actual), expected) in (-2_i64..=0).zip(&actual_coefficients).zip(expected)
            {
                let difference = actual - expected;
                assert!(
                    difference.together().cancel().is_zero(),
                    "project={project_onto_tensor_integrals}: mixed-mass tensor integration disagrees at epsilon^{power}: {difference}",
                );
            }
            mode_results.push(actual_coefficients);
        }
        let [projected, complete] = mode_results.as_slice() else {
            unreachable!()
        };
        for (power, (projected, complete)) in (-2_i64..=0).zip(projected.iter().zip(complete)) {
            let difference = (projected - complete).expand();
            assert!(
                difference.is_zero(),
                "Vakint input modes disagree at epsilon^{power}: {difference}",
            );
        }
        Ok(())
    }

    #[test]
    #[should_panic(
        expected = "local on-shell counterterms are deferred until local counterterms can be derived from the 4D expanded representation"
    )]
    fn os_dispatch_is_unconditionally_deferred() {
        test_initialise().unwrap();
        let graph = scalar_two_point_graph();
        let given = root_node(&graph);
        let current = two_point_node(&graph, ApproximationType::OS);
        let settings = UVgenerationSettings {
            generate_integrated: false,
            ..Default::default()
        };
        let ctx = UVCtx::new(&graph, &settings);
        let _ = uv_limit(&Full4dCts::root(), &ctx, &current, &given, &current, &given);
    }

    fn paper_series_through(rescaled: Atom, degree: i64) -> Atom {
        rescaled
            .series(GS.rescale, Atom::Zero, degree)
            .unwrap()
            .to_atom()
            .replace(GS.rescale)
            .with(Atom::one())
            .collect_factors()
    }

    fn paper_has_pattern(atom: &Atom, pattern: Atom) -> bool {
        GS.erase_uv_momentum_provenance(atom)
            .replace(pattern)
            .match_iter()
            .next()
            .is_some()
    }

    fn paper_has_muv(atom: &Atom) -> bool {
        paper_has_pattern(atom, Atom::var(GS.m_uv_expansion))
            || paper_has_pattern(atom, Atom::var(GS.m_uv_vacuum))
    }

    fn paper_assert_normalized_zero(label: &str, atom: Atom) {
        let normalized = finalize(
            GS.erase_uv_momentum_provenance(&atom)
                .replace(GS.dim)
                .with(Atom::num(4)),
        )
        .replace_map(|view, _, output| {
            if let AtomView::Pow(power) = view {
                let (base, exponent) = power.get_base_exp();
                if i64::try_from(exponent).is_ok_and(|power| power < 0) {
                    // Canonicalize only denominator polynomials; graph
                    // numerator factors stay undistributed in the oracle.
                    **output = base.factor().pow(exponent);
                }
            }
        })
        .collect_factors()
        // Resolve numeric signs such as A-B-(A-B) and cancel exact rational
        // factors without distributing products of graph numerator factors.
        .expand_num()
        .cancel()
        // Scalar regrouping can expose metric contractions that were hidden
        // inside separate additive factors at the first normalization pass.
        .simplify_metrics();
        assert!(normalized.is_zero(), "{label}: {normalized}");
    }

    fn paper_dot(left: Atom, right: Atom) -> Atom {
        FunctionBuilder::new(SPENSO_TAG.dot)
            .add_arg(left)
            .add_arg(right)
            .finish()
    }

    fn paper_routed_momentum(edge: EdgeIndex) -> Atom {
        FunctionBuilder::new(GS.emr_mom)
            .add_arg(usize::from(edge))
            .add_arg(Minkowski {}.new_rep(4).to_symbolic([]))
            .finish()
    }

    /// Translate a bare graph atom to a selected homogeneous route and expose
    /// its quadratic propagators.  In the UV version all physical masses are
    /// set to zero and every propagator receives the common terminal mass M,
    /// exactly as in the explicitly routed right-hand sides of Appendix B.
    fn paper_routed_rational(
        graph: &Graph,
        subgraph: &SuBitGraph,
        lmb: &LoopMomentumBasis,
        atom: &Atom,
        physical_masses: &[Atom],
        uv_denominators: bool,
    ) -> Atom {
        let uv_mass_squared = Atom::var(GS.m_uv_vacuum).pow(2);
        let mut routed = GS
            .erase_uv_momentum_provenance(atom)
            .replace_multiple(graph.uv_wrapped_replacement(subgraph, lmb, &[W_.x___]))
            .replace(function!(GS.den, W_.a_, W_.mom_, W_.mass_, W_.prop_))
            .with_map(move |matched| {
                let propagator = matched.get(W_.prop_).unwrap().to_atom();
                if uv_denominators {
                    propagator - &uv_mass_squared
                } else {
                    propagator
                }
            });
        if uv_denominators {
            for mass in physical_masses {
                routed = routed.replace(mass.clone()).with(Atom::Zero);
            }
        }
        routed
    }

    fn paper_routed_propagator(
        graph: &Graph,
        routing_subgraph: &SuBitGraph,
        lmb: &LoopMomentumBasis,
        edge: EdgeIndex,
        physical_masses: &[Atom],
        uv_denominator: bool,
    ) -> Atom {
        let mut single_edge = graph.empty_subgraph::<SuBitGraph>();
        single_edge.add(graph[&edge].1);
        paper_routed_rational(
            graph,
            routing_subgraph,
            lmb,
            &graph.denominator(&single_edge, |_| 1),
            physical_masses,
            uv_denominator,
        )
    }

    fn paper_scale_external_momenta(
        mut atom: Atom,
        external_edges: &[EdgeIndex],
        parameter: Symbol,
    ) -> Atom {
        let scale = Atom::var(parameter);
        for edge in external_edges {
            atom = atom
                .replace(GS.emr_mom(*edge, W_.x___))
                .with(GS.emr_mom(*edge, W_.x___) * &scale);
        }

        // Make bilinearity explicit before differentiating.  Symbolica must
        // otherwise regard `dot` as an opaque scalar function and leaves a
        // formal `der(dot,...)`, which is not the tensor derivative used in
        // Appendix B.  Every routed momentum is affine in the external scale,
        // so splitting both arguments into their value and first derivative at
        // zero is exact here.
        atom.replace_map(|view, _, output| {
            if let AtomView::Fun(dot) = view {
                if dot.get_symbol() != SPENSO_TAG.dot || dot.get_nargs() != 2 {
                    return;
                }
                let mut arguments = dot.iter();
                let left = arguments.next().unwrap().to_owned();
                let right = arguments.next().unwrap().to_owned();
                let left_zero = left.replace(parameter).with(Atom::Zero);
                let right_zero = right.replace(parameter).with(Atom::Zero);
                let left_linear = left
                    .derivative(parameter)
                    .replace(parameter)
                    .with(Atom::Zero);
                let right_linear = right
                    .derivative(parameter)
                    .replace(parameter)
                    .with(Atom::Zero);
                let dot_product = |left: Atom, right: Atom| {
                    if left.is_zero() || right.is_zero() {
                        Atom::Zero
                    } else {
                        FunctionBuilder::new(SPENSO_TAG.dot)
                            .add_arg(left)
                            .add_arg(right)
                            .finish()
                    }
                };
                **output = dot_product(left_zero.clone(), right_zero.clone())
                    + &scale
                        * (dot_product(left_linear.clone(), right_zero)
                            + dot_product(left_zero, right_linear.clone()))
                    + scale.pow(2) * dot_product(left_linear, right_linear);
            }
        })
        .collect_factors()
    }

    /// Extract a Taylor coefficient by literal differentiation.  This is the
    /// independent Appendix-B oracle: it neither calls the production series
    /// projectors nor reconstructs H from their U/S/US branch decomposition.
    fn paper_direct_coefficient(atom: &Atom, parameter: Symbol, degree: usize) -> Atom {
        let derivative = match degree {
            0 => atom.clone(),
            1 => atom.derivative(parameter),
            2 => atom.derivative(parameter).derivative(parameter) / 2,
            _ => panic!("the Appendix-B oracle only needs degrees zero through two"),
        };
        derivative
            .replace(parameter)
            .with(Atom::Zero)
            .collect_factors()
    }

    #[derive(Clone)]
    struct PaperDenominatorMapping {
        edge: EdgeIndex,
        physical_symbol: Symbol,
        uv_symbol: Symbol,
        physical_propagator: Atom,
        uv_propagator: Atom,
    }

    fn paper_is_uv_denominator(mass: &Atom, propagator: &Atom) -> bool {
        paper_has_pattern(mass, Atom::var(GS.m_uv_expansion))
            || paper_has_pattern(mass, Atom::var(GS.m_uv_vacuum))
            || paper_has_pattern(propagator, Atom::var(GS.m_uv_expansion))
            || paper_has_pattern(propagator, Atom::var(GS.m_uv_vacuum))
    }

    /// Before replacing a denominator by an atomic paper symbol, verify that
    /// its stored fourth argument is the routed physical or terminal-M
    /// propagator assigned to that edge.  The exact mapped equality below
    /// retains every factor and exponent, and therefore separately verifies
    /// the expected denominator powers/multiplicities in each paper term.
    fn paper_assert_denominator_provenance(
        label: &str,
        atom: &Atom,
        mappings: &[PaperDenominatorMapping],
    ) {
        let mut checked = Vec::<(EdgeIndex, bool, Atom)>::new();
        for matched in atom
            .replace(function!(GS.den, W_.a_, W_.mom_, W_.mass_, W_.prop_))
            .match_iter()
        {
            let edge_atom = matched.get(&W_.a_).unwrap().clone();
            let edge =
                EdgeIndex(usize::try_from(i64::try_from(edge_atom.as_view()).unwrap()).unwrap());
            let mass = matched.get(&W_.mass_).unwrap().clone();
            let propagator = matched.get(&W_.prop_).unwrap().clone();
            let uv = paper_is_uv_denominator(&mass, &propagator);
            if checked.iter().any(|(seen_edge, seen_uv, seen_propagator)| {
                *seen_edge == edge && *seen_uv == uv && seen_propagator == &propagator
            }) {
                continue;
            }
            let mapping = mappings
                .iter()
                .find(|mapping| mapping.edge == edge)
                .unwrap_or_else(|| {
                    panic!(
                        "{label}: no Appendix denominator mapping for e{}",
                        usize::from(edge)
                    )
                });
            let expected = if uv {
                &mapping.uv_propagator
            } else {
                &mapping.physical_propagator
            };
            paper_assert_normalized_zero(
                &format!(
                    "{label}: e{} {} denominator provenance",
                    usize::from(edge),
                    if uv { "terminal-M" } else { "physical" }
                ),
                propagator.clone() - expected,
            );
            checked.push((edge, uv, propagator));
        }
        assert!(
            !checked.is_empty(),
            "{label}: no production denominators were checked"
        );
    }

    /// Replace every production propagator by the paper denominator atom for
    /// its routed edge and branch. Keeping these denominators atomic lets the
    /// tensor/gamma comparison remain exact without factoring polynomials over
    /// the physical mass.
    fn paper_atomic_denominators(atom: Atom, mappings: &[PaperDenominatorMapping]) -> Atom {
        let mappings = mappings.to_vec();
        atom.replace(function!(GS.den, W_.a_, W_.mom_, W_.mass_, W_.prop_))
            .with_map(move |matched| {
                let edge = matched.get(W_.a_).unwrap().to_atom();
                let edge = i64::try_from(edge.as_view()).unwrap();
                let mass = matched.get(W_.mass_).unwrap().to_atom();
                let propagator = matched.get(W_.prop_).unwrap().to_atom();
                let uv = paper_is_uv_denominator(&mass, &propagator);
                let mapping = mappings
                    .iter()
                    .find(|mapping| usize::from(mapping.edge) as i64 == edge)
                    .unwrap_or_else(|| panic!("no Appendix denominator mapping for edge e{edge}"));
                Atom::var(if uv {
                    mapping.uv_symbol
                } else {
                    mapping.physical_symbol
                })
            })
    }

    fn paper_assert_atomic_rhs(
        label: &str,
        actual: Atom,
        expected: Atom,
        mappings: &[PaperDenominatorMapping],
    ) {
        paper_assert_denominator_provenance(label, &actual, mappings);
        let difference = finalize(
            GS.erase_uv_momentum_provenance(&paper_atomic_denominators(actual, mappings))
                .replace(function!(GS.ct_marker, W_.a_))
                .with(Atom::one())
                .replace(GS.dim)
                .with(Atom::num(4))
                - GS.erase_uv_momentum_provenance(&paper_atomic_denominators(expected, mappings))
                    .replace(GS.dim)
                    .with(Atom::num(4)),
        )
        .collect_factors()
        .expand_num()
        .cancel();
        assert!(
            difference.is_zero(),
            "{label}: {}",
            difference.log_print(Some(400))
        );
    }

    #[test]
    fn paper_eq_4_15_normalized_d0_d1_d2() {
        test_initialise().unwrap();
        let q = Atom::var(symbol!("paper_eq_4_15_q"));
        let mass = Atom::var(symbol!("paper_eq_4_15_mass"));
        let x = &mass
            + (&mass + Atom::one()) * &q
            + (&mass + Atom::num(2)) * q.pow(2)
            + (&mass + Atom::num(3)) * q.pow(3);
        let soft_deformation = x.replace(q.clone()).with(&q * GS.rescale);
        let uv_deformation = soft_deformation
            .replace(mass.clone())
            .with(&mass * GS.rescale);

        for degree in 0..=2 {
            let soft = tilde_taylor(soft_deformation.clone(), degree).unwrap();
            let uv = paper_series_through(uv_deformation.clone(), i64::from(degree));
            let uv_soft_deformation = soft
                .replace(q.clone())
                .with(&q * GS.rescale)
                .replace(mass.clone())
                .with(&mass * GS.rescale);
            let uv_soft = paper_series_through(uv_soft_deformation, i64::from(degree));
            let hat = (&uv + &soft - &uv_soft).expand();
            let expected = match degree {
                0 => Atom::Zero,
                1 => &mass + &q,
                2 => &mass + &q + &mass * &q + Atom::num(2) * q.pow(2),
                _ => unreachable!(),
            };

            assert!(
                !paper_has_muv(&hat),
                "the completed free-mass operator at d={degree} must not introduce M"
            );
            assert_eq!((hat - expected).expand(), Atom::Zero, "d={degree}");
            for (label, branch) in [("ordinary", &uv), ("soft", &soft), ("overlap", &uv_soft)] {
                assert!(
                    !paper_has_muv(branch),
                    "the free-mass {label} branch at d={degree} must not introduce M"
                );
            }
            if degree == 0 {
                assert!(soft.is_zero(), "Eq. 4.15 requires T-hat_0 = T_0");
            } else {
                assert!(paper_has_pattern(&soft, mass.clone()));
            }
        }
    }

    #[test]
    fn paper_eq_4_16_quadratic_kernel_has_the_physical_soft_mass_split() {
        test_initialise().unwrap();
        let lambda = Atom::var(GS.rescale);
        let k_squared = Atom::var(symbol!("paper_eq_4_16_k_squared"));
        let mass = Atom::var(symbol!("paper_eq_4_16_mass"));
        let m_uv = Atom::var(GS.m_uv_expansion);
        let linear = Atom::var(symbol!("paper_eq_4_16_linear"));
        let quadratic = Atom::var(symbol!("paper_eq_4_16_quadratic"));
        let physical_denominator = &k_squared - mass.pow(2);
        let uv_denominator = &k_squared - m_uv.pow(2);

        // T_2 scales the physical mass together with the external momentum,
        // whereas tilde-T_1 leaves it in the unexpanded soft denominator.
        let ordinary_kernel = Atom::one()
            / (&uv_denominator
                + &lambda * &linear
                + lambda.pow(2) * (&quadratic - mass.pow(2) + m_uv.pow(2)));
        let ordinary = paper_series_through(ordinary_kernel.clone(), 2);
        let soft = paper_series_through(
            Atom::one() / (&physical_denominator + &lambda * &linear + lambda.pow(2) * &quadratic),
            1,
        );

        // Applying T_2 to the already completed soft branch produces the
        // overlap, including the rescaled difference between m and m_UV.
        let overlap_denominator = &uv_denominator + lambda.pow(2) * (m_uv.pow(2) - mass.pow(2));
        let overlap = paper_series_through(
            Atom::one() / &overlap_denominator - &lambda * &linear / overlap_denominator.pow(2),
            2,
        );
        let actual = ordinary + &soft - overlap;
        let soft_part = Atom::one() / &physical_denominator - &linear / physical_denominator.pow(2);
        let quadratic_part =
            linear.pow(2) / uv_denominator.pow(3) - &quadratic / uv_denominator.pow(2);
        let expected = &soft_part + &quadratic_part;

        paper_assert_normalized_zero(
            "Eq. 4.16 quadratic mass split",
            (actual - expected).cancel(),
        );
        assert!(paper_has_pattern(&soft_part, mass.clone()));
        assert!(!paper_has_muv(&soft_part));
        assert!(paper_has_muv(&quadratic_part));
        assert!(!paper_has_pattern(&quadratic_part, mass.clone()));

        // The highest Taylor coefficient is exactly one half of the second
        // derivative at lambda=0, as displayed explicitly in Eq. 4.16.
        let half_second_derivative = ordinary_kernel
            .derivative(GS.rescale)
            .derivative(GS.rescale)
            .replace(GS.rescale)
            .with(Atom::Zero)
            / 2;
        let expected_second_order = linear.pow(2) / uv_denominator.pow(3)
            - (&quadratic - mass.pow(2) + m_uv.pow(2)) / uv_denominator.pow(2);
        paper_assert_normalized_zero(
            "Eq. 4.16 includes the Taylor factor 1/2",
            (half_second_derivative - expected_second_order).cancel(),
        );
    }

    #[test]
    fn paper_appendix_b_7_gluon_t2_includes_the_linear_numerator_cross_term() {
        test_initialise().unwrap();
        let lambda = Atom::var(GS.rescale);
        let p_squared = Atom::var(symbol!("paper_appendix_b_7_p_squared"));
        let p_dot_q = Atom::var(symbol!("paper_appendix_b_7_p_dot_q"));
        let q_squared = Atom::var(symbol!("paper_appendix_b_7_q_squared"));
        let mass = Atom::var(symbol!("paper_appendix_b_7_mass"));
        let m_uv = Atom::var(GS.m_uv_expansion);
        let physical_numerator = Atom::var(symbol!("paper_appendix_b_7_physical_n0"));
        let physical_linear = Atom::var(symbol!("paper_appendix_b_7_physical_n1"));
        let uv_numerator = Atom::var(symbol!("paper_appendix_b_7_uv_n0"));
        let uv_linear = Atom::var(symbol!("paper_appendix_b_7_uv_n1"));
        let physical_denominator = &p_squared - mass.pow(2);
        let uv_denominator = &p_squared - m_uv.pow(2);

        // The fixed k-dependent denominators and gamma traces in Eq. B.7 are
        // represented by independent numerator atoms.  The displayed
        // p-dependent denominator is D^2(D-2 lambda p.q+lambda^2 q^2).
        let physical_kernel = (&physical_numerator - &lambda * &physical_linear)
            / (physical_denominator.pow(2)
                * (&physical_denominator - Atom::num(2) * &lambda * &p_dot_q
                    + lambda.pow(2) * &q_squared));
        let uv_kernel = (&uv_numerator - &lambda * &uv_linear)
            / (uv_denominator.pow(2)
                * (&uv_denominator - Atom::num(2) * &lambda * &p_dot_q
                    + lambda.pow(2) * &q_squared));

        let physical_low_orders = paper_series_through(physical_kernel, 1);
        let uv_quadratic =
            paper_series_through(uv_kernel.clone(), 2) - paper_series_through(uv_kernel, 1);
        let actual = &physical_low_orders + &uv_quadratic;
        let expected_physical = (&physical_numerator - &physical_linear)
            / physical_denominator.pow(3)
            + Atom::num(2) * &p_dot_q * &physical_numerator / physical_denominator.pow(4);
        let expected_quadratic = (-Atom::num(2) * &p_dot_q * &uv_linear
            - &q_squared * &uv_numerator)
            / uv_denominator.pow(4)
            + Atom::num(4) * p_dot_q.pow(2) * &uv_numerator / uv_denominator.pow(5);

        paper_assert_normalized_zero(
            "Appendix B.7 physical low orders and UV quadratic coefficient",
            (actual - &expected_physical - &expected_quadratic).cancel(),
        );
        assert!(paper_has_pattern(&physical_low_orders, mass));
        assert!(!paper_has_muv(&physical_low_orders));
        assert!(paper_has_muv(&uv_quadratic));
    }

    #[test]
    fn paper_appendix_b_6_b_7_production_outer_h2_matches_the_routed_rhs() {
        test_initialise().unwrap();
        std::thread::Builder::new()
            .stack_size(64 * 1024 * 1024)
            .spawn(|| {
                let graph = paper_appendix_b1_graph();
                let gamma = paper_appendix_b1_node(&graph, [2, 3, 4, 5, 6], 2);
                let root = root_node(&graph);
                let settings = UVgenerationSettings {
                    generate_integrated: false,
                    ..Default::default()
                };
                let ctx = UVCtx::new(&graph, &settings);
                let seed = Full4dCts::root();
                let grown = grow(&seed, &ctx, &gamma, &root).unwrap();
                let lmb = compatible_lmb(&ctx, &gamma, &root, &seed).unwrap();

                let ordinary = t_raw(&grown, &ctx, &gamma, &root, &lmb).unwrap();
                let physical_low_orders = tilde_t_raw(&grown, &ctx, &gamma, &root, &lmb).unwrap();
                let overlap = t_raw(&physical_low_orders, &ctx, &gamma, &root, &lmb).unwrap();
                let uv_quadratic = (&ordinary - &overlap).expand();
                let actual = hat_t_raw(&grown, &ctx, &gamma, &root, &lmb).unwrap();

                // Independently encode Eqs. (B.6)-(B.7) in GammaLoop's
                // homogeneous route.  The first two displayed lines of B.7
                // are the constant and linear q coefficients with the top
                // mass kept physical.  Its final two lines are the quadratic
                // q coefficient after setting the numerator masses to zero
                // and introducing M in every quadratic denominator.  Direct
                // differentiation is intentionally used here instead of any
                // production Taylor helper.
                let mass = graph[EdgeIndex(2)].mass_atom();
                let external_edges = component_external_edges(&ctx, &gamma, &lmb);
                assert_eq!(
                    external_edges,
                    vec![EdgeIndex(0), EdgeIndex(1)],
                    "the routed B.7 oracle must scale the two oriented representatives of q"
                );
                let paper_parameter = symbol!("paper_appendix_b_7_direct_q_scale");
                let routed_numerator = graph
                    .numerator(gamma.subgraph(), root.subgraph())
                    .to_d_dim(GS.dim)
                    .color_simplify()
                    .get_single_atom()
                    .unwrap()
                    .replace_multiple(graph.uv_wrapped_replacement(
                        gamma.subgraph(),
                        &lmb,
                        &[W_.x___],
                    ));
                let physical_numerator_scaled = paper_scale_external_momenta(
                    routed_numerator.clone(),
                    &external_edges,
                    paper_parameter,
                );
                let physical_numerator = [0, 1].map(|degree| {
                    paper_direct_coefficient(&physical_numerator_scaled, paper_parameter, degree)
                });
                let uv_numerator_scaled = paper_scale_external_momenta(
                    routed_numerator.replace(mass.clone()).with(Atom::Zero),
                    &external_edges,
                    paper_parameter,
                );
                let uv_numerator = [0, 1, 2].map(|degree| {
                    paper_direct_coefficient(&uv_numerator_scaled, paper_parameter, degree)
                });
                let physical_p = paper_routed_propagator(
                    &graph,
                    gamma.subgraph(),
                    &lmb,
                    EdgeIndex(2),
                    std::slice::from_ref(&mass),
                    false,
                );
                let physical_k = paper_routed_propagator(
                    &graph,
                    gamma.subgraph(),
                    &lmb,
                    EdgeIndex(3),
                    std::slice::from_ref(&mass),
                    false,
                );
                let physical_k_minus_p = paper_routed_propagator(
                    &graph,
                    gamma.subgraph(),
                    &lmb,
                    EdgeIndex(4),
                    std::slice::from_ref(&mass),
                    false,
                );
                let physical_p_copy = paper_routed_propagator(
                    &graph,
                    gamma.subgraph(),
                    &lmb,
                    EdgeIndex(5),
                    std::slice::from_ref(&mass),
                    false,
                );
                let physical_shifted_p = paper_routed_propagator(
                    &graph,
                    gamma.subgraph(),
                    &lmb,
                    EdgeIndex(6),
                    std::slice::from_ref(&mass),
                    false,
                );
                let physical_shifted_p = paper_scale_external_momenta(
                    physical_shifted_p,
                    &external_edges,
                    paper_parameter,
                );
                let physical_shift = [0, 1, 2].map(|degree| {
                    paper_direct_coefficient(&physical_shifted_p, paper_parameter, degree)
                });
                let uv_p = paper_routed_propagator(
                    &graph,
                    gamma.subgraph(),
                    &lmb,
                    EdgeIndex(2),
                    std::slice::from_ref(&mass),
                    true,
                );
                let uv_k = paper_routed_propagator(
                    &graph,
                    gamma.subgraph(),
                    &lmb,
                    EdgeIndex(3),
                    std::slice::from_ref(&mass),
                    true,
                );
                let uv_k_minus_p = paper_routed_propagator(
                    &graph,
                    gamma.subgraph(),
                    &lmb,
                    EdgeIndex(4),
                    std::slice::from_ref(&mass),
                    true,
                );
                let uv_shifted_p = paper_routed_propagator(
                    &graph,
                    gamma.subgraph(),
                    &lmb,
                    EdgeIndex(6),
                    std::slice::from_ref(&mass),
                    true,
                );
                let uv_shifted_p =
                    paper_scale_external_momenta(uv_shifted_p, &external_edges, paper_parameter);
                let uv_shift = [0, 1, 2]
                    .map(|degree| paper_direct_coefficient(&uv_shifted_p, paper_parameter, degree));

                // Verify independently that the denominator atoms used below
                // are exactly the routed D_p, D_k and D_{k-p} of B.7.
                let p = paper_routed_momentum(EdgeIndex(2));
                let k = paper_routed_momentum(EdgeIndex(3));
                let k_minus_p = &k - &p;
                let uv_mass_squared = Atom::var(GS.m_uv_vacuum).pow(2);
                paper_assert_normalized_zero(
                    "Appendix B.7 physical D_p",
                    &physical_p - paper_dot(p.clone(), p.clone()) + mass.pow(2),
                );
                paper_assert_normalized_zero(
                    "Appendix B.7 physical D_k",
                    &physical_k - paper_dot(k.clone(), k.clone()) + mass.pow(2),
                );
                paper_assert_normalized_zero(
                    "Appendix B.7 physical D_{k-p}",
                    &physical_k_minus_p - paper_dot(k_minus_p.clone(), k_minus_p.clone()),
                );
                paper_assert_normalized_zero(
                    "Appendix B.7 repeated physical D_p",
                    physical_p_copy - &physical_p,
                );
                paper_assert_normalized_zero(
                    "Appendix B.7 shifted physical denominator at q=0",
                    &physical_shift[0] - &physical_p,
                );
                paper_assert_normalized_zero(
                    "Appendix B.7 UV D_p",
                    &uv_p - paper_dot(p.clone(), p) + &uv_mass_squared,
                );
                paper_assert_normalized_zero(
                    "Appendix B.7 UV D_k",
                    &uv_k - paper_dot(k.clone(), k) + &uv_mass_squared,
                );
                paper_assert_normalized_zero(
                    "Appendix B.7 UV D_{k-p}",
                    &uv_k_minus_p - paper_dot(k_minus_p.clone(), k_minus_p) + &uv_mass_squared,
                );
                paper_assert_normalized_zero(
                    "Appendix B.7 shifted UV denominator at q=0",
                    &uv_shift[0] - &uv_p,
                );
                paper_assert_normalized_zero(
                    "Appendix B.7 physical and UV shifts have the same q coefficients",
                    &physical_shift[1] - &uv_shift[1],
                );
                paper_assert_normalized_zero(
                    "Appendix B.7 physical and UV shifts have the same q-squared coefficient",
                    &physical_shift[2] - &uv_shift[2],
                );

                let dp = symbol!("paper_appendix_b_7_Dp_phys");
                let dk = symbol!("paper_appendix_b_7_Dk_phys");
                let dkp = symbol!("paper_appendix_b_7_Dkp_phys");
                let up = symbol!("paper_appendix_b_7_Dp_uv");
                let uk = symbol!("paper_appendix_b_7_Dk_uv");
                let ukp = symbol!("paper_appendix_b_7_Dkp_uv");
                let denominator_mappings = [
                    PaperDenominatorMapping {
                        edge: EdgeIndex(2),
                        physical_symbol: dp,
                        uv_symbol: up,
                        physical_propagator: physical_p.clone(),
                        uv_propagator: uv_p.clone(),
                    },
                    PaperDenominatorMapping {
                        edge: EdgeIndex(5),
                        physical_symbol: dp,
                        uv_symbol: up,
                        physical_propagator: physical_p.clone(),
                        uv_propagator: uv_p.clone(),
                    },
                    PaperDenominatorMapping {
                        edge: EdgeIndex(6),
                        physical_symbol: dp,
                        uv_symbol: up,
                        physical_propagator: physical_shift[0].clone(),
                        uv_propagator: uv_shift[0].clone(),
                    },
                    PaperDenominatorMapping {
                        edge: EdgeIndex(3),
                        physical_symbol: dk,
                        uv_symbol: uk,
                        physical_propagator: physical_k.clone(),
                        uv_propagator: uv_k.clone(),
                    },
                    PaperDenominatorMapping {
                        edge: EdgeIndex(4),
                        physical_symbol: dkp,
                        uv_symbol: ukp,
                        physical_propagator: physical_k_minus_p.clone(),
                        uv_propagator: uv_k_minus_p.clone(),
                    },
                ];
                let dp = Atom::var(dp);
                let dk_dkp = Atom::var(dk) * Atom::var(dkp);
                let up = Atom::var(up);
                let uk_ukp = Atom::var(uk) * Atom::var(ukp);
                let expected_physical = (&physical_numerator[0] + &physical_numerator[1])
                    / (dp.pow(3) * &dk_dkp)
                    - &physical_shift[1] * &physical_numerator[0] / (dp.pow(4) * &dk_dkp);
                let expected_uv_quadratic = &uv_numerator[2] / (up.pow(3) * &uk_ukp)
                    - &uv_shift[1] * &uv_numerator[1] / (up.pow(4) * &uk_ukp)
                    + &uv_numerator[0]
                        * (uv_shift[1].pow(2) / up.pow(5) - &uv_shift[2] / up.pow(4))
                        / &uk_ukp;
                let expected = &expected_physical + &expected_uv_quadratic;

                paper_assert_atomic_rhs(
                    "Appendix B.7 physical tensor/gamma lines",
                    physical_low_orders.clone(),
                    expected_physical.clone(),
                    &denominator_mappings,
                );
                paper_assert_atomic_rhs(
                    "Appendix B.7 UV-quadratic tensor/gamma lines",
                    uv_quadratic.clone(),
                    expected_uv_quadratic.clone(),
                    &denominator_mappings,
                );
                paper_assert_atomic_rhs(
                    "Appendix B.6-B.7 complete routed H2 RHS",
                    actual.clone(),
                    expected.clone(),
                    &denominator_mappings,
                );

                paper_assert_normalized_zero(
                    "Appendix B.6-B.7 production H2 branch split",
                    actual.clone() - &physical_low_orders - &uv_quadratic,
                );
                assert!(
                    !physical_low_orders.is_zero(),
                    "the constant and linear physical-mass part of Eq. B.7 must survive"
                );
                assert!(
                    paper_has_pattern(&physical_low_orders, graph[EdgeIndex(2)].mass_atom()),
                    "the S1 part of Eq. B.7 must retain the physical top mass"
                );
                assert!(
                    !paper_has_muv(&physical_low_orders),
                    "the S1 part of Eq. B.7 must not introduce the UV mass"
                );
                assert!(
                    paper_has_muv(&uv_quadratic),
                    "the homogeneous quadratic part of Eq. B.7 must use UV denominators"
                );
                paper_assert_normalized_zero(
                    "the Eq. B.7 UV remainder has no constant or linear soft jet",
                    tilde_t_raw(&uv_quadratic, &ctx, &gamma, &root, &lmb).unwrap(),
                );

                // This establishes the signed F_Gamma family without mixing
                // the already-finalized production tensor with the raw open
                // edge labels of the paper oracle.  Together with the direct
                // H2-to-B.7 comparison above, it is the transitive comparison
                // of the production forest atom to the independent RHS.
                let signed = uv_limit(&seed, &ctx, &gamma, &root, &gamma, &root).unwrap();
                let stripped = signed
                    .atom()
                    .replace(function!(GS.ct_marker, W_.a_))
                    .with(Atom::one());
                let expected_signed = -finalize(actual);
                paper_assert_normalized_zero(
                    "Appendix B.7 F_Gamma is minus the finalized independent H2 match",
                    stripped - expected_signed,
                );
            })
            .unwrap()
            .join()
            .unwrap();
    }

    #[test]
    fn paper_appendix_b_10_b_12_massless_child_and_local_b_16_nested_branch() {
        test_initialise().unwrap();
        std::thread::Builder::new()
            .stack_size(64 * 1024 * 1024)
            .spawn(|| {
                let graph = paper_appendix_b1_massless_graph();
                let gamma_one = paper_appendix_b1_node(&graph, [3, 4], 1);
                let gamma = paper_appendix_b1_node(&graph, [2, 3, 4, 5, 6], 2);
                let root = root_node(&graph);
                let settings = UVgenerationSettings {
                    generate_integrated: false,
                    ..Default::default()
                };
                let ctx = UVCtx::new(&graph, &settings);
                let seed = Full4dCts::root();
                assert!(
                    graph[EdgeIndex(3)].mass_atom().is_zero(),
                    "the B.10-B.12 IR specialization must use a massless fermion"
                );

                let child_grown = grow(&seed, &ctx, &gamma_one, &root).unwrap();
                let child_lmb = compatible_lmb(&ctx, &gamma_one, &root, &seed).unwrap();
                let child_u = t_raw(
                    &child_grown,
                    &ctx,
                    &gamma_one,
                    &root,
                    &child_lmb,
                )
                .unwrap();
                let child_s =
                    tilde_t_raw(&child_grown, &ctx, &gamma_one, &root, &child_lmb).unwrap();
                let child_us = t_raw(
                    &child_s,
                    &ctx,
                    &gamma_one,
                    &root,
                    &child_lmb,
                )
                .unwrap();
                let child_h = hat_t_raw(&child_grown, &ctx, &gamma_one, &root, &child_lmb).unwrap();
                let child_uv_linear = (&child_u - &child_us).collect_factors();

                // At m=p_os=0, Eq. (4.25) identifies the paper's OS operator
                // with H_1.  Encode the two surviving routed tensor/gamma
                // terms of B.10-B.12 independently: the physical p=0 value
                // and the coefficient linear in p with M in both quadratic
                // denominators.
                let child_external_edges = component_external_edges(&ctx, &gamma_one, &child_lmb);
                assert_eq!(
                    child_external_edges,
                    vec![EdgeIndex(2), EdgeIndex(5)],
                    "the routed B.10-B.12 oracle must scale both oriented representatives of the paper momentum p=e2"
                );
                let child_parameter = symbol!("paper_appendix_b_12_direct_p_scale");
                let child_routed_numerator = graph
                    .numerator(gamma_one.subgraph(), root.subgraph())
                    .to_d_dim(GS.dim)
                    .color_simplify()
                    .get_single_atom()
                    .unwrap()
                    .replace_multiple(graph.uv_wrapped_replacement(
                        gamma_one.subgraph(),
                        &child_lmb,
                        &[W_.x___],
                    ));
                let child_numerator_scaled = paper_scale_external_momenta(
                    child_routed_numerator,
                    &child_external_edges,
                    child_parameter,
                );
                let child_numerator = [0, 1].map(|degree| {
                    paper_direct_coefficient(&child_numerator_scaled, child_parameter, degree)
                });
                let child_physical_k = paper_routed_propagator(
                    &graph,
                    gamma_one.subgraph(),
                    &child_lmb,
                    EdgeIndex(3),
                    &[],
                    false,
                );
                let child_physical_k_minus_p = paper_routed_propagator(
                    &graph,
                    gamma_one.subgraph(),
                    &child_lmb,
                    EdgeIndex(4),
                    &[],
                    false,
                );
                let child_physical_shifted = paper_scale_external_momenta(
                    child_physical_k_minus_p.clone(),
                    &child_external_edges,
                    child_parameter,
                );
                let child_physical_shift = [0, 1].map(|degree| {
                    paper_direct_coefficient(&child_physical_shifted, child_parameter, degree)
                });
                let child_uv_k = paper_routed_propagator(
                    &graph,
                    gamma_one.subgraph(),
                    &child_lmb,
                    EdgeIndex(3),
                    &[],
                    true,
                );
                let child_uv_k_minus_p = paper_routed_propagator(
                    &graph,
                    gamma_one.subgraph(),
                    &child_lmb,
                    EdgeIndex(4),
                    &[],
                    true,
                );
                let child_uv_shifted = paper_scale_external_momenta(
                    child_uv_k_minus_p.clone(),
                    &child_external_edges,
                    child_parameter,
                );
                let child_uv_shift = [0, 1].map(|degree| {
                    paper_direct_coefficient(&child_uv_shifted, child_parameter, degree)
                });

                // Check the independently routed scalar denominators before
                // replacing them by atomic D_k/U_k placeholders.  This is the
                // massless specialization of the two denominators displayed
                // in B.10-B.12.
                let child_p = paper_routed_momentum(EdgeIndex(2));
                let child_k = paper_routed_momentum(EdgeIndex(3));
                let child_k_minus_p = &child_k - &child_p;
                let uv_mass_squared = Atom::var(GS.m_uv_vacuum).pow(2);
                paper_assert_normalized_zero(
                    "Appendix B.10-B.12 physical D_k",
                    &child_physical_k - paper_dot(child_k.clone(), child_k.clone()),
                );
                paper_assert_normalized_zero(
                    "Appendix B.10-B.12 physical D_{k-p}",
                    &child_physical_k_minus_p
                        - paper_dot(child_k_minus_p.clone(), child_k_minus_p.clone()),
                );
                paper_assert_normalized_zero(
                    "Appendix B.10-B.12 physical D_{k-p} at p=0",
                    &child_physical_shift[0] - &child_physical_k,
                );
                paper_assert_normalized_zero(
                    "Appendix B.10-B.12 UV U_k",
                    &child_uv_k - paper_dot(child_k.clone(), child_k) + &uv_mass_squared,
                );
                paper_assert_normalized_zero(
                    "Appendix B.10-B.12 UV U_{k-p}",
                    &child_uv_k_minus_p - paper_dot(child_k_minus_p.clone(), child_k_minus_p)
                        + &uv_mass_squared,
                );
                paper_assert_normalized_zero(
                    "Appendix B.10-B.12 UV U_{k-p} at p=0",
                    &child_uv_shift[0] - &child_uv_k,
                );
                paper_assert_normalized_zero(
                    "Appendix B.10-B.12 physical and UV p-linear shifts agree",
                    &child_physical_shift[1] - &child_uv_shift[1],
                );

                let dk_symbol = symbol!("paper_appendix_b_12_Dk_phys");
                let uk_symbol = symbol!("paper_appendix_b_12_Dk_uv");
                let child_denominator_mappings = [
                    PaperDenominatorMapping {
                        edge: EdgeIndex(3),
                        physical_symbol: dk_symbol,
                        uv_symbol: uk_symbol,
                        physical_propagator: child_physical_k.clone(),
                        uv_propagator: child_uv_k.clone(),
                    },
                    PaperDenominatorMapping {
                        edge: EdgeIndex(4),
                        physical_symbol: dk_symbol,
                        uv_symbol: uk_symbol,
                        physical_propagator: child_physical_shift[0].clone(),
                        uv_propagator: child_uv_shift[0].clone(),
                    },
                ];
                let dk = Atom::var(dk_symbol);
                let uk = Atom::var(uk_symbol);
                let expected_child_s = &child_numerator[0] / dk.pow(2);
                let expected_child_uv_linear = &child_numerator[1] / uk.pow(2)
                    - &child_uv_shift[1] * &child_numerator[0] / uk.pow(3);
                let expected_child_u = &child_numerator[0] / uk.pow(2)
                    + &expected_child_uv_linear;
                let expected_child_h = &expected_child_s + &expected_child_uv_linear;

                paper_assert_atomic_rhs(
                    "Appendix B.11 massless physical tensor/gamma term",
                    child_s.clone(),
                    expected_child_s.clone(),
                    &child_denominator_mappings,
                );
                paper_assert_atomic_rhs(
                    "Appendix B.12 massless UV-linear tensor/gamma term",
                    child_uv_linear.clone(),
                    expected_child_uv_linear.clone(),
                    &child_denominator_mappings,
                );
                paper_assert_atomic_rhs(
                    "Appendix B.10-B.12 complete massless child RHS",
                    child_h.clone(),
                    expected_child_h.clone(),
                    &child_denominator_mappings,
                );
                paper_assert_atomic_rhs(
                    "Appendix B.12 complete ordinary U1 child bracket used by the outer U",
                    child_u.clone(),
                    expected_child_u,
                    &child_denominator_mappings,
                );

                // As a paper-side algebraic identity, at m=p_os=0 the two
                // P_+ and P_- terms in B.11 add to the identity. They are
                // precisely the physical S_0 term in B.11, while B.12
                // continues the same local T_1(gamma_1) formula with its
                // ordinary-UV linear coefficient. Production OS dispatch is
                // nevertheless unconditionally deferred in Phase 1: this is
                // an H_1 oracle, not a massless-OS shortcut. The integrated
                // bar-K construction starts only at B.13.
                paper_assert_normalized_zero(
                    "massless B.10-B.12 specialization is S0 plus the UV-linear remainder",
                    &child_h - &child_s - &child_uv_linear,
                );
                assert!(!child_s.is_zero(), "the massless B.11 S0 term must survive");
                assert!(
                    !paper_has_muv(&child_s),
                    "the massless B.11 S0 term must retain physical denominators"
                );
                assert!(
                    paper_has_muv(&child_uv_linear),
                    "the massless B.12 UV-linear term must use UV denominators"
                );
                paper_assert_normalized_zero(
                    "the massless B.12 UV-linear remainder has no constant soft jet",
                    tilde_t_raw(&child_uv_linear, &ctx, &gamma_one, &root, &child_lmb).unwrap(),
                );

                let child_local =
                    uv_limit(&seed, &ctx, &gamma_one, &root, &gamma_one, &root).unwrap();
                let child_stripped = child_local
                    .atom()
                    .replace(function!(GS.ct_marker, W_.a_))
                    .with(Atom::one());
                paper_assert_normalized_zero(
                    "Appendix B.10-B.12 F_gamma is minus the finalized independent H1 match",
                    child_stripped + finalize(child_h.clone()),
                );
                let completed_child = Full4dCts::from_factorized_local(&child_local);
                let nested_grown = grow(&completed_child, &ctx, &gamma, &gamma_one).unwrap();
                let outer_lmb = compatible_lmb(&ctx, &gamma, &gamma_one, &completed_child).unwrap();
                let outer_u = t_raw(
                    &nested_grown,
                    &ctx,
                    &gamma,
                    &gamma_one,
                    &outer_lmb,
                )
                .unwrap();
                let outer_s = tilde_t_raw(
                    &nested_grown,
                    &ctx,
                    &gamma,
                    &gamma_one,
                    &outer_lmb,
                )
                .unwrap();
                let outer_us = t_raw(
                    &outer_s,
                    &ctx,
                    &gamma,
                    &gamma_one,
                    &outer_lmb,
                )
                .unwrap();
                let outer_h = hat_t_raw(
                    &nested_grown,
                    &ctx,
                    &gamma,
                    &gamma_one,
                    &outer_lmb,
                )
                .unwrap();

                // B.16 has two products.  Its first is the complete massless
                // B.10-B.12 child multiplied by the physical constant/linear
                // q jet of the reduced outer graph.  Its second uses the
                // independently re-expanded UV constant/linear child bracket
                // and the quadratic-q UV bracket of that co-graph.  Construct
                // all four factors directly from the routed bare networks.
                let co_graph = grow_factor(&ctx, &gamma, &gamma_one).unwrap();
                let reduced = gamma.reduced_subgraph(&gamma_one);
                let outer_external_edges = component_external_edges(&ctx, &gamma, &outer_lmb);
                assert_eq!(
                    outer_external_edges,
                    vec![EdgeIndex(0), EdgeIndex(1)],
                    "the routed B.16 oracle must scale the two oriented representatives of q"
                );
                let outer_parameter = symbol!("paper_appendix_b_16_direct_q_scale");
                let co_graph_routed_numerator = graph
                    .numerator(&reduced, gamma_one.subgraph())
                    .to_d_dim(GS.dim)
                    .color_simplify()
                    .get_single_atom()
                    .unwrap()
                    .replace_multiple(graph.uv_wrapped_replacement(
                        &reduced,
                        &outer_lmb,
                        &[W_.x___],
                    ));
                let co_graph_numerator_scaled = paper_scale_external_momenta(
                    co_graph_routed_numerator.clone(),
                    &outer_external_edges,
                    outer_parameter,
                );
                let co_graph_numerator = [0, 1, 2].map(|degree| {
                    paper_direct_coefficient(&co_graph_numerator_scaled, outer_parameter, degree)
                });
                let outer_physical_p =
                    paper_routed_propagator(&graph, &reduced, &outer_lmb, EdgeIndex(2), &[], false);
                let outer_physical_p_copy =
                    paper_routed_propagator(&graph, &reduced, &outer_lmb, EdgeIndex(5), &[], false);
                let outer_physical_p_minus_q =
                    paper_routed_propagator(&graph, &reduced, &outer_lmb, EdgeIndex(6), &[], false);
                let outer_physical_shifted = paper_scale_external_momenta(
                    outer_physical_p_minus_q.clone(),
                    &outer_external_edges,
                    outer_parameter,
                );
                let outer_physical_shift = [0, 1, 2].map(|degree| {
                    paper_direct_coefficient(&outer_physical_shifted, outer_parameter, degree)
                });
                let outer_uv_p =
                    paper_routed_propagator(&graph, &reduced, &outer_lmb, EdgeIndex(2), &[], true);
                let outer_uv_p_copy =
                    paper_routed_propagator(&graph, &reduced, &outer_lmb, EdgeIndex(5), &[], true);
                let outer_uv_p_minus_q =
                    paper_routed_propagator(&graph, &reduced, &outer_lmb, EdgeIndex(6), &[], true);
                let outer_uv_shifted = paper_scale_external_momenta(
                    outer_uv_p_minus_q.clone(),
                    &outer_external_edges,
                    outer_parameter,
                );
                let outer_uv_shift = [0, 1, 2].map(|degree| {
                    paper_direct_coefficient(&outer_uv_shifted, outer_parameter, degree)
                });

                // The reduced graph carries p on e2/e5 and p-q on e6.
                // Assert that identification before using the atomic D_p/U_p
                // placeholders in the two products displayed in B.16.
                let outer_p = paper_routed_momentum(EdgeIndex(2));
                let outer_q = paper_routed_momentum(EdgeIndex(0));
                let outer_p_minus_q = &outer_p - &outer_q;
                paper_assert_normalized_zero(
                    "Appendix B.16 physical D_p",
                    &outer_physical_p - paper_dot(outer_p.clone(), outer_p.clone()),
                );
                paper_assert_normalized_zero(
                    "Appendix B.16 repeated physical D_p",
                    &outer_physical_p_copy - &outer_physical_p,
                );
                paper_assert_normalized_zero(
                    "Appendix B.16 physical D_{p-q}",
                    &outer_physical_p_minus_q
                        - paper_dot(outer_p_minus_q.clone(), outer_p_minus_q.clone()),
                );
                paper_assert_normalized_zero(
                    "Appendix B.16 physical D_{p-q} at q=0",
                    &outer_physical_shift[0] - &outer_physical_p,
                );
                paper_assert_normalized_zero(
                    "Appendix B.16 UV U_p",
                    &outer_uv_p - paper_dot(outer_p.clone(), outer_p) + &uv_mass_squared,
                );
                paper_assert_normalized_zero(
                    "Appendix B.16 repeated UV U_p",
                    &outer_uv_p_copy - &outer_uv_p,
                );
                paper_assert_normalized_zero(
                    "Appendix B.16 UV U_{p-q}",
                    &outer_uv_p_minus_q - paper_dot(outer_p_minus_q.clone(), outer_p_minus_q)
                        + &uv_mass_squared,
                );
                paper_assert_normalized_zero(
                    "Appendix B.16 UV U_{p-q} at q=0",
                    &outer_uv_shift[0] - &outer_uv_p,
                );
                paper_assert_normalized_zero(
                    "Appendix B.16 physical and UV q-linear shifts agree",
                    &outer_physical_shift[1] - &outer_uv_shift[1],
                );
                paper_assert_normalized_zero(
                    "Appendix B.16 physical and UV q-quadratic shifts agree",
                    &outer_physical_shift[2] - &outer_uv_shift[2],
                );

                let dp_symbol = symbol!("paper_appendix_b_16_Dp_phys");
                let up_symbol = symbol!("paper_appendix_b_16_Dp_uv");
                let denominator_mappings = [
                    PaperDenominatorMapping {
                        edge: EdgeIndex(2),
                        physical_symbol: dp_symbol,
                        uv_symbol: up_symbol,
                        physical_propagator: outer_physical_p.clone(),
                        uv_propagator: outer_uv_p.clone(),
                    },
                    PaperDenominatorMapping {
                        edge: EdgeIndex(5),
                        physical_symbol: dp_symbol,
                        uv_symbol: up_symbol,
                        physical_propagator: outer_physical_p.clone(),
                        uv_propagator: outer_uv_p.clone(),
                    },
                    PaperDenominatorMapping {
                        edge: EdgeIndex(6),
                        physical_symbol: dp_symbol,
                        uv_symbol: up_symbol,
                        physical_propagator: outer_physical_shift[0].clone(),
                        uv_propagator: outer_uv_shift[0].clone(),
                    },
                    PaperDenominatorMapping {
                        edge: EdgeIndex(3),
                        physical_symbol: dk_symbol,
                        uv_symbol: uk_symbol,
                        physical_propagator: child_physical_k.clone(),
                        uv_propagator: child_uv_k.clone(),
                    },
                    PaperDenominatorMapping {
                        edge: EdgeIndex(4),
                        physical_symbol: dk_symbol,
                        uv_symbol: uk_symbol,
                        physical_propagator: child_physical_shift[0].clone(),
                        uv_propagator: child_uv_shift[0].clone(),
                    },
                ];
                let dp = Atom::var(dp_symbol);
                let up = Atom::var(up_symbol);
                let dpq_symbol = symbol!("paper_appendix_b_16_Dpq_phys");
                let upq_symbol = symbol!("paper_appendix_b_16_Dpq_uv");
                let bare_co_graph_mappings = [
                    PaperDenominatorMapping {
                        edge: EdgeIndex(2),
                        physical_symbol: dp_symbol,
                        uv_symbol: up_symbol,
                        physical_propagator: outer_physical_p.clone(),
                        uv_propagator: outer_uv_p.clone(),
                    },
                    PaperDenominatorMapping {
                        edge: EdgeIndex(5),
                        physical_symbol: dp_symbol,
                        uv_symbol: up_symbol,
                        physical_propagator: outer_physical_p.clone(),
                        uv_propagator: outer_uv_p.clone(),
                    },
                    PaperDenominatorMapping {
                        edge: EdgeIndex(6),
                        physical_symbol: dpq_symbol,
                        uv_symbol: upq_symbol,
                        physical_propagator: outer_physical_p_minus_q.clone(),
                        uv_propagator: outer_uv_p_minus_q.clone(),
                    },
                ];
                let expected_co_graph_bare =
                    &co_graph_routed_numerator / (dp.pow(2) * Atom::var(dpq_symbol));
                let co_graph_routed = co_graph.replace_multiple(graph.uv_wrapped_replacement(
                    &reduced,
                    &outer_lmb,
                    &[W_.x___],
                ));
                paper_assert_atomic_rhs(
                    "Appendix B.10-B.12 F_gamma bare co-graph factor",
                    co_graph_routed,
                    expected_co_graph_bare,
                    &bare_co_graph_mappings,
                );
                assert_eq!(
                    nested_grown,
                    (child_local.atom() * &co_graph).simplify_metrics(),
                    "Appendix B.10-B.12 F_gamma must be the signed child forest atom times its co-graph",
                );
                let expected_co_graph_physical = (&co_graph_numerator[0] + &co_graph_numerator[1])
                    / dp.pow(3)
                    - &outer_physical_shift[1] * &co_graph_numerator[0] / dp.pow(4);
                // The TeX of B.16 writes its first quadratic numerator as an
                // overall minus multiplying a bracket which itself contains
                // `-q^2`.  Read literally that conflicts with B.7.  Taking the
                // direct second derivative of the bare B.1 integrand fixes
                // the unambiguous B.7 sign: -q^2 N0/U_p^4.
                let expected_co_graph_uv_quadratic = &co_graph_numerator[2] / up.pow(3)
                    - &outer_uv_shift[1] * &co_graph_numerator[1] / up.pow(4)
                    + &co_graph_numerator[0]
                        * (outer_uv_shift[1].pow(2) / up.pow(5) - &outer_uv_shift[2] / up.pow(4));
                // The child is a completed local atom before the parent acts.
                // Reuse the graph-label-preserving H1 and U1 factors only
                // after their independent manual RHS comparisons above.  The
                // open source/sink labels are graph boundary identities, not
                // dummy indices which an independently finalized expression
                // may rename before its co-graph is attached.
                let verified_child_h = finalize(child_h.clone());
                let verified_child_outer_u = finalize(child_u.clone());
                let expected_b_16_physical =
                    (&verified_child_h * &expected_co_graph_physical).simplify_metrics();
                let expected_b_16_uv =
                    (&verified_child_outer_u * &expected_co_graph_uv_quadratic)
                        .simplify_metrics();
                let expected_b_16 = &expected_b_16_physical + &expected_b_16_uv;

                // B.15 identifies T(T(gamma_1)*(Gamma\\gamma_1)) as the first
                // contribution to the nested K operation; B.16 is its explicit
                // local formula.  Keep its physical and UV products distinct
                // so that no regrouping can hide a wrong nested forest term.
                // B.17, not B.16, starts the corresponding integrated bar-K
                // contribution.
                let outer_uv_remainder = (&outer_u - &outer_us).collect_factors();
                // The child boundary labels must remain fixed while its H1
                // and U1 atoms are checked above.  Once a completed child is
                // attached to its co-graph, however, those paired labels are
                // contracted dummies.  Normalize each fully composed side
                // separately so the tensor/color collector gives them the
                // same names; dangling graph labels remain untouched.
                // Denominators are checked first, then made atomic so their
                // stored momenta cannot be reinterpreted as tensor powers.
                assert_eq!(
                    graph.inv(Hedge(8)),
                    Hedge(9),
                    "the B.16 child/co-graph boundary must be the two half-edges of e5"
                );
                let child_crown = graph.dummy_less_full_crown(gamma_one.subgraph());
                assert!(child_crown.includes(&Hedge(8)));
                assert!(!child_crown.includes(&Hedge(9)));
                // Simplifying the e5 identity may retain either member of
                // this pair as the name of the resulting contraction.  Use
                // the co-graph-side member only after the factors are joined.
                let child_side_boundary = Atom::from(Aind::Hedge(8, 0));
                let co_graph_side_boundary = Atom::from(Aind::Hedge(9, 0));
                let assert_b_16_rhs = |label: &str, actual: Atom, expected: Atom| {
                    paper_assert_denominator_provenance(label, &actual, &denominator_mappings);
                    let prepare = |atom: Atom| {
                        finalize(
                            GS.erase_uv_momentum_provenance(&paper_atomic_denominators(
                                atom,
                                &denominator_mappings,
                            ))
                                .replace(function!(GS.ct_marker, W_.a_))
                                .with(Atom::one())
                                .replace(GS.dim)
                                .with(Atom::num(4)),
                        )
                        .replace(child_side_boundary.clone())
                        .with(co_graph_side_boundary.clone())
                        // Closing the child/co-graph boundary can turn an open
                        // color chain and metric into a trace on either side.
                        .simplify_color()
                        // Finish finite index contractions across the joined
                        // tensor sums; scalar numerator factors stay factored.
                        .schoonschip_with_net_full::<Aind>()
                        .expect("finite tensor contractions succeed")
                        .to_dots()
                        .normalize_dots()
                    };
                    // Compare physical momenta after provenance has served its
                    // routing purpose, preserving every graph numerator factor.
                    let mut difference = prepare(actual) - prepare(expected);
                    for _ in 0..8 {
                        let contracted = difference
                            .collect_factors()
                            .expand_num()
                            .schoonschip_with_net_full::<Aind>()
                            .expect("finite tensor contractions succeed")
                            .to_dots()
                            .normalize_dots();
                        // Scalar regrouping can expose another finite contraction.
                        // Keep complete tensor calls (including their momentum
                        // arguments and common numerator spectators) opaque, and
                        // normalize only their scalar Taylor/denominator coefficients.
                        let mut tensors = Vec::new();
                        contracted.visitor(&mut |view| {
                            if let AtomView::Fun(fun) = view {
                                if fun.get_symbol() != SPENSO_TAG.dot
                                    && view.is_tensorial(StrictTensorFilter::ContainsReps)
                                    && !tensors.iter().any(|atom: &Atom| atom.as_view() == view)
                                {
                                    tensors.push(view.to_owned());
                                }
                                return false;
                            }
                            true
                        });
                        // A finite linear tensor jet may have opaque tensor
                        // spectators, but never multiply two numerator sums.
                        let literal_tensor = |view: AtomView<'_>| {
                            tensors.iter().any(|tensor| tensor.as_view() == view)
                        };
                        contracted.visitor(&mut |view| {
                            if matches!(view, AtomView::Fun(_)) {
                                return false;
                            }
                            if let AtomView::Mul(product) = view {
                                let summed_factors = product
                                    .iter()
                                    .filter(|factor| {
                                        tensors.iter().any(|tensor| factor.contains(tensor))
                                            && !literal_tensor(*factor)
                                            && !matches!(factor, AtomView::Pow(power)
                                                if literal_tensor(power.get_base())
                                                    && i64::try_from(power.get_exp()).is_ok())
                                    })
                                    .count();
                                assert!(
                                    summed_factors <= 1,
                                    "{label}: tensor coefficient collection would distribute numerator products"
                                );
                            }
                            true
                        });
                        let coefficients = contracted
                            .coefficient_list_exact(&tensors)
                            .expect("the finite B.16 tensor polynomial stays within its bound");
                        let next = Atom::add_many(coefficients.into_iter().map(|(tensor, scalar)| {
                            tensor * scalar.together().cancel()
                        }))
                        .collect_factors();
                        if next == difference {
                            break;
                        }
                        difference = next;
                        if difference.is_zero() {
                            break;
                        }
                    }
                    assert!(
                        difference.is_zero(),
                        "{label}: {}",
                        difference.log_print(Some(400))
                    );
                };
                assert_b_16_rhs(
                    "Appendix B.16 physical product before regrouping",
                    outer_s.clone(),
                    -&expected_b_16_physical,
                );
                assert_b_16_rhs(
                    "Appendix B.16 UV product before regrouping",
                    outer_uv_remainder.clone(),
                    -&expected_b_16_uv,
                );
                assert_b_16_rhs(
                    "Appendix B.16 complete local nested RHS before forest signing",
                    outer_h.clone(),
                    -&expected_b_16,
                );
                assert!(
                    paper_has_muv(&outer_uv_remainder),
                    "the local B.16 UV product must retain terminal UV denominators"
                );
                paper_assert_normalized_zero(
                    "the local B.16 UV product has no constant or linear outer soft jet",
                    tilde_t_raw(&outer_uv_remainder, &ctx, &gamma, &gamma_one, &outer_lmb).unwrap(),
                );
                paper_assert_normalized_zero(
                    "massless local B.16 branch comparison",
                    &outer_h - &outer_s - &outer_uv_remainder,
                );
                assert!(
                    !outer_s.is_zero(),
                    "the outer physical low-order product of massless B.16 must survive"
                );
                let nested = uv_limit(
                    &completed_child,
                    &ctx,
                    &gamma,
                    &gamma_one,
                    &gamma,
                    &gamma_one,
                )
                .unwrap();
                let nested_stripped = nested
                    .atom()
                    .replace(function!(GS.ct_marker, W_.a_))
                    .with(Atom::one());
                paper_assert_normalized_zero(
                    "Appendix B.16 F_gamma;Gamma is minus the finalized local B.16 match",
                    nested_stripped + finalize(outer_h),
                );
            })
            .unwrap()
            .join()
            .unwrap();
    }

    #[test]
    fn production_dod2_massive_denominator_has_the_paper_mass_coefficient() {
        test_initialise().unwrap();
        let graph: Graph = dot!(
            digraph massive_scalar_tadpole {
                edge [particle="H" num=1];
                node [num=1];
                ext [style=invis];
                ext -> A:0 [id=0];
                A:1 -> ext [id=1];
                A:2 -> A:3 [id=2 lmb_id=0];
            }
        )
        .unwrap();
        let external: SuBitGraph = graph.external_filter();
        let subgraph = graph.full_filter().subtract(&external);
        let lmb = graph.lmb_of(&subgraph);
        let mass = graph[EdgeIndex(2)].mass_atom();
        let integrand = Atom::one() / graph.denominator(&subgraph, |_| 1);
        let n_loops = graph.n_loops(&subgraph);
        assert_eq!(n_loops, 1);
        let rescaled = graph.uv_rescaled(&subgraph, n_loops, &lmb, &lmb, &integrand);
        let series = rescaled.series(GS.rescale, Atom::Zero, 0).unwrap();

        assert_eq!(series.get_trailing_exponent(), Rational::from(-2));
        assert!(!mass.is_zero());
        let leading = series
            .coefficient(Rational::from(-2))
            .expect("requested coefficient is within series precision");
        let finite = series
            .coefficient(Rational::from(0))
            .expect("requested coefficient is within series precision");
        let expected_finite = (mass.pow(2) - Atom::var(GS.m_uv_expansion).pow(2)) * leading.pow(2);
        paper_assert_normalized_zero(
            "production DOD-2 denominator mass coefficient",
            finite - expected_finite,
        );
    }

    #[test]
    fn production_outer_u_grades_a_free_fermion_mass_but_keeps_muv_fixed() {
        test_initialise().unwrap();
        let graph: Graph =
            include_str!("../../../../../tests/resources/graphs/local_os_top_self_energy.dot")
                .into_graph(&crate::utils::load_generic_model("sm"))
                .unwrap();
        let given = root_node(&graph);
        let current = two_point_node(&graph, ApproximationType::MUV);
        let settings = UVgenerationSettings {
            generate_integrated: false,
            ..Default::default()
        };
        let ctx = UVCtx::new(&graph, &settings);
        let reduced = current.reduced_subgraph(&given);
        let numerator = graph
            .numerator(&reduced, given.subgraph())
            .to_d_dim(GS.dim)
            .color_simplify()
            .get_single_atom()
            .unwrap();
        let denominator = graph.denominator(&reduced, |_| 1);
        let physical_mass = graph[EdgeIndex(2)].mass_atom();
        let massless_numerator = numerator.replace(physical_mass.clone()).with(Atom::Zero);
        let fermion_mass_numerator = &numerator - massless_numerator;
        assert!(
            !fermion_mass_numerator.is_zero(),
            "the production top propagator numerator must contain its free physical mass"
        );

        // Check that the expression under study is exactly the production
        // numerator and denominator assembled by `grow`.
        let grown = grow(&Full4dCts::root(), &ctx, &current, &given).unwrap();
        paper_assert_normalized_zero(
            "production grow uses the numerator whose fermion-mass term is tested",
            &grown - (numerator / &denominator).simplify_metrics(),
        );

        let n_loops = graph.n_loops(current.subgraph());
        assert_eq!(n_loops, 1);
        let rescaled = graph.uv_rescaled(&reduced, n_loops, current.lmb(), current.lmb(), &grown);
        let series = rescaled.series(GS.rescale, Atom::Zero, 0).unwrap();
        assert_eq!(
            series.get_trailing_exponent(),
            Rational::from(-1),
            "the DOD-1 fermion self-energy must start at inverse-hard power -1"
        );
        let degree_zero = series
            .coefficient(Rational::from(-1))
            .expect("requested coefficient is within series precision");
        let degree_one = series
            .coefficient(Rational::from(0))
            .expect("requested coefficient is within series precision");
        assert!(
            !paper_has_pattern(&degree_zero, physical_mass.clone()),
            "a physical fermion mass must not occur in the degree-zero U coefficient"
        );
        assert!(
            paper_has_pattern(&degree_one, physical_mass.clone()),
            "the free physical fermion mass must first occur at direct-U degree one"
        );

        let m_uv_expansion = Atom::var(GS.m_uv_expansion);
        assert!(
            paper_has_pattern(&degree_zero, m_uv_expansion.clone()),
            "fixed mUVexp must already define the degree-zero terminal denominator"
        );
        assert!(
            paper_has_pattern(&degree_one, m_uv_expansion),
            "the degree-one coefficient must retain the same terminal mUVexp denominator"
        );

        // Isolating the actual fermion numerator mass before deformation gives
        // exactly the mass-dependent part of the degree-one coefficient. In
        // the inverse chart k -> k/t, exponent -1 is direct Taylor degree zero
        // for this DOD-1 graph and exponent 0 is degree one. The terminal
        // mUVexp denominator is therefore fixed at degree zero, whereas m_t is
        // graded once, as in the direct deformation m_t -> t m_t.
        let mass_rescaled = graph.uv_rescaled(
            &reduced,
            n_loops,
            current.lmb(),
            current.lmb(),
            &(&fermion_mass_numerator / &denominator),
        );
        let mass_series = mass_rescaled.series(GS.rescale, Atom::Zero, 0).unwrap();
        assert_eq!(mass_series.get_trailing_exponent(), Rational::from(0));
        let degree_one_mass_part =
            &degree_one - degree_one.replace(physical_mass.clone()).with(Atom::Zero);
        paper_assert_normalized_zero(
            "the production fermion numerator mass is the degree-one U contribution",
            degree_one_mass_part
                - mass_series
                    .coefficient(Rational::from(0))
                    .expect("requested coefficient is within series precision"),
        );

        let projected = t_raw(&grown, &ctx, &current, &given, current.lmb()).unwrap();
        assert!(!projected.is_zero());
        assert!(paper_has_pattern(&projected, physical_mass));
        assert!(paper_has_pattern(&projected, Atom::var(GS.m_uv_expansion)));
        assert!(!paper_has_pattern(&projected, Atom::var(GS.rescale)));
    }

    #[test]
    fn nested_outer_u_grades_physical_child_mass_and_preserves_terminal_m() {
        test_initialise().unwrap();
        let lambda = Atom::var(GS.rescale);
        let physical_mass = Atom::var(symbol!("nested_outer_physical_mass"));
        let free_inner_m = Atom::var(GS.m_uv_expansion);
        let terminal_m = Atom::var(GS.m_uv_vacuum);
        let k_squared = Atom::var(symbol!("nested_outer_child_k_squared"));
        let terminal_denominator = &k_squared - terminal_m.pow(2);

        // This is the parent-U deformation of a completed child term: the
        // child's physical mass and a free numerator mUVexp introduced by its
        // U/US branch are both graded. The mUV embedded in an existing terminal
        // denominator remains fixed in role and propagator shape.
        let deformed = &lambda * (&free_inner_m + &physical_mass) / &terminal_denominator;
        let series = deformed.series(GS.rescale, Atom::Zero, 1).unwrap();
        let constant = series
            .coefficient(Rational::from(0))
            .expect("requested coefficient is within series precision");
        let linear = series
            .coefficient(Rational::from(1))
            .expect("requested coefficient is within series precision");
        paper_assert_normalized_zero(
            "free numerator masses have no outer-U constant coefficient",
            constant,
        );
        paper_assert_normalized_zero(
            "free physical and inner-U numerator masses are graded by the outer U",
            &linear - (&free_inner_m + &physical_mass) / &terminal_denominator,
        );
        assert!(paper_has_pattern(&linear, terminal_m));

        let held_child = (&free_inner_m + &physical_mass) / &terminal_denominator;
        assert!(
            !(&held_child
                - series
                    .coefficient(Rational::from(0))
                    .expect("requested coefficient is within series precision"))
            .together()
            .cancel()
            .is_zero(),
            "holding the complete child fixed would incorrectly put its free numerator masses in the zeroth parent coefficient"
        );
    }

    #[test]
    fn production_outer_u_grades_free_local_muv_but_preserves_its_terminal_denominators() {
        test_initialise().unwrap();
        let graph = scalar_two_point_graph();
        let given = root_node(&graph);
        let current = two_point_node(&graph, ApproximationType::MUV);
        let settings = UVgenerationSettings {
            generate_integrated: false,
            ..Default::default()
        };
        let ctx = UVCtx::new(&graph, &settings);
        let reduced = current.reduced_subgraph(&given);
        let local_muv = Atom::var(GS.m_uv_expansion);
        let grown = grow(&Full4dCts::root(), &ctx, &current, &given).unwrap();
        let inner_projection = t_raw(&grown, &ctx, &current, &given, current.lmb()).unwrap();
        let completed_child = &local_muv * inner_projection;

        // Equation (4.4) requires an outer Taylor operation to grade a free
        // mUV factor produced in an inner numerator.  In the inverse hard
        // chart it stays fixed while the loop measure and two terminal
        // propagators cancel, so it first occurs at direct Taylor degree one
        // (hard-series exponent zero for this superficially linear parent).
        let rescaled = graph.uv_rescaled(
            &reduced,
            graph.n_loops(current.subgraph()),
            current.lmb(),
            current.lmb(),
            &completed_child,
        );
        let series = rescaled.series(GS.rescale, Atom::Zero, 0).unwrap();
        assert_eq!(series.get_trailing_exponent(), Rational::from(0));
        assert!(
            series
                .coefficient(Rational::from(-1))
                .expect("requested coefficient is within series precision")
                .is_zero()
        );

        // The free mUVexp coefficient accompanies existing terminal
        // denominators whose actual vacuum mass is mUV. The parent must retain
        // those denominators instead of applying a second
        // physical-propagator rearrangement.
        let routed = completed_child.replace_multiple(graph.uv_wrapped_replacement(
            &reduced,
            current.lmb(),
            &[W_.x___],
        ));
        let projected = t_raw(&completed_child, &ctx, &current, &given, current.lmb()).unwrap();
        paper_assert_normalized_zero(
            "outer U grades the free local mUV while retaining terminal denominators",
            projected - routed,
        );
    }

    #[test]
    fn paper_eq_4_19_three_point_constant_and_linear_terms_use_distinct_masses() {
        test_initialise().unwrap();
        let lambda = Atom::var(GS.rescale);
        let k_squared = Atom::var(symbol!("paper_eq_4_19_k_squared"));
        let m_uv = Atom::var(GS.m_uv_expansion);
        let numerator_constant = Atom::var(symbol!("paper_eq_4_19_numerator_constant"));
        let p_one = Atom::var(symbol!("paper_eq_4_19_p_one"));
        let p_two = Atom::var(symbol!("paper_eq_4_19_p_two"));
        let numerator_linear = Atom::var(symbol!("paper_eq_4_19_a_one")) * &p_one
            + Atom::var(symbol!("paper_eq_4_19_a_two")) * &p_two;
        let denominator_linear = Atom::var(symbol!("paper_eq_4_19_b_one")) * p_one
            + Atom::var(symbol!("paper_eq_4_19_b_two")) * p_two;
        let physical_denominator = k_squared.clone();
        let uv_denominator = &k_squared - m_uv.pow(2);

        let ordinary = paper_series_through(
            (&numerator_constant + &lambda * &numerator_linear)
                / (&uv_denominator + &lambda * &denominator_linear),
            1,
        );
        let soft = &numerator_constant / &physical_denominator;
        let overlap = &numerator_constant / &uv_denominator;
        let actual = ordinary + &soft - overlap;
        let linear_part = &numerator_linear / &uv_denominator
            - &numerator_constant * &denominator_linear / uv_denominator.pow(2);
        let expected = &soft + &linear_part;

        paper_assert_normalized_zero(
            "Eq. 4.19 constant/linear mass split",
            (actual - expected).cancel(),
        );
        assert!(!paper_has_muv(&soft));
        assert!(paper_has_muv(&linear_part));
    }

    #[test]
    fn production_nested_toy_forest_terms_match_the_explicit_derivation() {
        test_initialise().unwrap();
        let graph: Graph = dot!(
            digraph nested_soft_toy_production {
                edge [particle=scalar_1 num=1]
                node [num=1]
                ext [style=invis]
                ext -> C:0 [id=0]
                A:1 -> ext [id=1]
                C:9 -> A:2 [id=2]
                B:3 -> C:4 [id=3 lmb_id=0 num="Q(3,spenso::cind(0))+UFO::mass_scalar_1"]
                B:5 -> C:6 [id=4]
                A:7 -> B:8 [id=5 lmb_id=1 num="Q(5,spenso::cind(0))"]
            },
            "scalars"
        )
        .unwrap();
        let mut child_subgraph = graph.empty_subgraph::<SuBitGraph>();
        for edge in [EdgeIndex(3), EdgeIndex(4)] {
            child_subgraph.add(graph[&edge].1);
        }
        let external: SuBitGraph = graph.external_filter();
        let parent_subgraph = graph.full_filter().subtract(&external);
        let nested_lmb = |subgraph: &SuBitGraph| {
            graph
                .underlying
                .try_compatible_sub_lmb(
                    subgraph,
                    graph.dummy_less_full_crown(subgraph).subtract(&external),
                    &graph.loop_momentum_basis,
                )
                .unwrap()
        };
        let child = TestNode {
            lmb: nested_lmb(&child_subgraph),
            subgraph: child_subgraph,
            dod: 1,
            scheme: ApproximationType::IR,
        };
        let parent = TestNode {
            lmb: graph.loop_momentum_basis.clone(),
            subgraph: parent_subgraph,
            dod: 2,
            scheme: ApproximationType::IR,
        };
        let root = root_node(&graph);
        let settings = UVgenerationSettings {
            generate_integrated: false,
            add_marker: true,
            keep_marker: true,
            ..Default::default()
        };
        let ctx = UVCtx::new(&graph, &settings);
        let seed = Full4dCts::root();

        let f_empty = grow(&seed, &ctx, &parent, &root).unwrap();
        let child_local = uv_limit(&seed, &ctx, &child, &root, &child, &root).unwrap();
        let completed_child = Full4dCts::from_factorized_local(&child_local);
        let f_child = grow(&completed_child, &ctx, &parent, &child).unwrap();
        let f_parent = uv_limit(&seed, &ctx, &parent, &root, &parent, &root).unwrap();
        let f_nested = uv_limit(&completed_child, &ctx, &parent, &child, &parent, &child).unwrap();

        let child_loop = EdgeIndex(3);
        let child_boundary = EdgeIndex(5);
        let unused_affine_crown = EdgeIndex(2);
        assert!(
            parent.lmb.edge_signatures[unused_affine_crown]
                .external
                .iter()
                .any(|sign| sign.is_sign()),
            "e2=p+q must be affine in the graph route"
        );
        assert_eq!(
            completed_child.0.inherited_loop_edges,
            vec![child_loop],
            "only the child loop generator constrains the parent; its boundary remains routed"
        );
        for edge in &completed_child.0.inherited_loop_edges {
            assert!(
                paper_has_pattern(completed_child.atom(), GS.emr_mom(*edge, W_.x___)),
                "recorded child expansion coordinate {edge:?} must occur in its atom"
            );
        }
        assert!(
            !paper_has_pattern(
                completed_child.atom(),
                GS.emr_mom(unused_affine_crown, W_.x___)
            ),
            "the unused affine crown carrier must neither occur in nor constrain the child atom"
        );

        let families = [
            ("empty", &f_empty),
            ("child", &f_child),
            ("parent", f_parent.atom()),
            ("child;parent", f_nested.atom()),
        ];
        for (label, family) in families {
            assert!(
                !family.is_zero(),
                "the production {label} family must be retained"
            );
        }
        let marker_depths = |atom: &Atom| {
            atom.replace(function!(GS.ct_marker, W_.a_))
                .match_iter()
                .map(|matched| {
                    matched
                        .get(&W_.a_)
                        .unwrap()
                        .replace(function!(GS.uv_approx, W_.a_))
                        .match_iter()
                        .count()
                })
                .collect::<Vec<_>>()
        };
        assert!(marker_depths(&f_empty).is_empty());
        for (label, family, depth) in [
            ("child", &f_child, 1),
            ("parent", f_parent.atom(), 1),
            ("child;parent", f_nested.atom(), 2),
        ] {
            let depths = marker_depths(family);
            assert!(
                !depths.is_empty(),
                "the {label} family must retain its marker"
            );
            assert!(
                depths.iter().all(|actual| *actual == depth),
                "the {label} marker history must have depth {depth}: {depths:?}"
            );
        }

        // Normalize to the independent scalar coordinates in the nested-soft-toy
        // section of docs/architecture/uv-renormalization.typ. This algebraic
        // mapping does not invoke either Taylor projector.
        let k0 = Atom::var(symbol!("production_nested_toy_k0"));
        let p0 = Atom::var(symbol!("production_nested_toy_p0"));
        let k_vec = Atom::var(symbol!("production_nested_toy_k_vec"));
        let p_vec = Atom::var(symbol!("production_nested_toy_p_vec"));
        let q_vec = Atom::var(symbol!("production_nested_toy_q_vec"));
        let k_squared = Atom::var(symbol!("production_nested_toy_k_squared"));
        let p_squared = Atom::var(symbol!("production_nested_toy_p_squared"));
        let q_squared = Atom::var(symbol!("production_nested_toy_q_squared"));
        let k_dot_p = Atom::var(symbol!("production_nested_toy_k_dot_p"));
        let p_dot_q = Atom::var(symbol!("production_nested_toy_p_dot_q"));
        let k_dot_q = Atom::var(symbol!("production_nested_toy_k_dot_q"));
        let mass = Atom::var(symbol!("production_nested_toy_mass"));
        let m_uv = Atom::var(symbol!("production_nested_toy_m_uv"));
        let parent_route =
            graph.uv_wrapped_replacement(parent.subgraph(), parent.lmb(), &[W_.x___]);
        let normalize = |atom: &Atom| {
            let routed = GS
                .erase_uv_momentum_provenance(atom)
                .replace(function!(GS.ct_marker, W_.a_))
                .with(Atom::one())
                .replace_multiple(&parent_route);
            let indexed = finalize(routed)
                .replace(function!(GS.den, W_.a_, W_.mom_, W_.mass_, W_.prop_))
                .with(Atom::var(W_.prop_))
                .replace(graph[child_loop].mass_atom())
                .with(mass.clone())
                .replace(GS.m_uv_expansion)
                .with(m_uv.clone())
                .replace(GS.m_uv_vacuum)
                .with(m_uv.clone())
                .replace(GS.emr_mom(child_loop, GS.cind(0)))
                .with(k0.clone())
                .replace(GS.emr_mom(child_boundary, GS.cind(0)))
                .with(p0.clone())
                .replace(GS.emr_mom(child_loop, W_.x___))
                .with(k_vec.clone())
                .replace(GS.emr_mom(child_boundary, W_.x___))
                .with(p_vec.clone())
                .replace(GS.emr_mom(EdgeIndex(0), W_.x___))
                .with(q_vec.clone());
            indexed
                .replace_map(|view, _, output| {
                    let AtomView::Fun(dot) = view else {
                        return;
                    };
                    if dot.get_symbol() != SPENSO_TAG.dot || dot.get_nargs() != 2 {
                        return;
                    }
                    let mut arguments = dot.iter();
                    let left = arguments.next().unwrap().into_atom();
                    let right = arguments.next().unwrap().into_atom();
                    let replacement = if left == k_vec && right == k_vec {
                        &k_squared
                    } else if left == p_vec && right == p_vec {
                        &p_squared
                    } else if left == q_vec && right == q_vec {
                        &q_squared
                    } else if (left == k_vec && right == p_vec) || (left == p_vec && right == k_vec)
                    {
                        &k_dot_p
                    } else if (left == p_vec && right == q_vec) || (left == q_vec && right == p_vec)
                    {
                        &p_dot_q
                    } else if (left == k_vec && right == q_vec) || (left == q_vec && right == k_vec)
                    {
                        &k_dot_q
                    } else {
                        panic!("unexpected routed dot product dot({left},{right})")
                    };
                    **output = replacement.clone();
                })
                .together()
                .cancel()
                .expand()
        };

        let actual_empty = normalize(&f_empty);
        let actual_child = normalize(&f_child);
        let actual_parent = normalize(f_parent.atom());
        let actual_nested = normalize(f_nested.atom());

        let v_k = &k_squared - mass.pow(2);
        let v_r = &k_squared - Atom::num(2) * &k_dot_p + &p_squared - mass.pow(2);
        let v_p = &p_squared - mass.pow(2);
        let v_p_plus_q = &p_squared + Atom::num(2) * &p_dot_q + &q_squared - mass.pow(2);
        let u_k = &k_squared - m_uv.pow(2);
        let u_r = &k_squared - Atom::num(2) * &k_dot_p + &p_squared - m_uv.pow(2);
        let u_p = &p_squared - m_uv.pow(2);
        let child_bare = (&k0 + &mass) / (&v_k * &v_r);
        let child_hat = (&k0 + &mass) / v_k.pow(2) + Atom::num(2) * &k0 * &k_dot_p / u_k.pow(3);
        let child_bare_uv = &k0 / (&u_k * &u_r);
        let child_hat_uv = &k0 / u_k.pow(2) + Atom::num(2) * &k0 * &k_dot_p / u_k.pow(3);
        let l_m = Atom::one() / v_p.pow(2) - Atom::num(2) * &p_dot_q / v_p.pow(3);
        let q_u = Atom::num(4) * p_dot_q.pow(2) / u_p.pow(4) - &q_squared / u_p.pow(3);
        let expected_empty = &p0 * &child_bare / (&v_p * &v_p_plus_q);
        let expected_child = -&p0 * &child_hat / (&v_p * &v_p_plus_q);
        let expected_parent = -&p0 * &child_bare * &l_m - &p0 * &child_bare_uv * &q_u;
        let expected_nested = &p0 * &child_hat * &l_m + &p0 * &child_hat_uv * &q_u;

        for (label, actual, expected) in [
            ("F_empty", &actual_empty, &expected_empty),
            ("F_gamma", &actual_child, &expected_child),
            ("F_Gamma", &actual_parent, &expected_parent),
            ("F_gamma;Gamma", &actual_nested, &expected_nested),
        ] {
            let difference = (actual - expected).together().cancel().expand();
            assert!(difference.is_zero(), "{label}: {difference}");
        }

        let lambda_symbol = symbol!("production_nested_toy_lambda");
        let lambda = Atom::var(lambda_symbol);
        let child_remainder = &actual_empty + &actual_child;
        let child_hard = child_remainder
            .replace(k0.clone())
            .with(&k0 / &lambda)
            .replace(k_squared.clone())
            .with(&k_squared / lambda.pow(2))
            .replace(k_dot_p.clone())
            .with(&k_dot_p / &lambda)
            .together()
            .cancel()
            .series(lambda_symbol, Atom::Zero, 5)
            .unwrap();
        assert_eq!(child_hard.get_trailing_exponent(), Rational::from(5));

        let complete = actual_empty + actual_child + actual_parent + actual_nested;
        let parent_hard = complete
            .replace(p0.clone())
            .with(&p0 / &lambda)
            .replace(p_squared.clone())
            .with(&p_squared / lambda.pow(2))
            .replace(k_dot_p.clone())
            .with(&k_dot_p / &lambda)
            .replace(p_dot_q.clone())
            .with(&p_dot_q / &lambda)
            .together()
            .cancel()
            .series(lambda_symbol, Atom::Zero, 5)
            .unwrap();
        assert_eq!(parent_hard.get_trailing_exponent(), Rational::from(5));

        let soft = complete
            .replace(p_dot_q.clone())
            .with(&p_dot_q * &lambda)
            .replace(q_squared.clone())
            .with(&q_squared * lambda.pow(2))
            .together()
            .cancel()
            .series(lambda_symbol, Atom::Zero, 2)
            .unwrap();
        assert_eq!(soft.get_trailing_exponent(), Rational::from(2));

        let simultaneous = complete
            .replace(k0.clone())
            .with(&k0 / &lambda)
            .replace(p0.clone())
            .with(&p0 / &lambda)
            .replace(k_squared.clone())
            .with(&k_squared / lambda.pow(2))
            .replace(p_squared.clone())
            .with(&p_squared / lambda.pow(2))
            .replace(k_dot_p.clone())
            .with(&k_dot_p / lambda.pow(2))
            .replace(p_dot_q.clone())
            .with(&p_dot_q / &lambda)
            .together()
            .cancel()
            .series(lambda_symbol, Atom::Zero, 9)
            .unwrap();
        assert_eq!(simultaneous.get_trailing_exponent(), Rational::from(9));
    }

    #[test]
    fn nested_toy_forest_terms_reexpand_the_completed_soft_child() {
        test_initialise().unwrap();
        let k0 = Atom::var(symbol!("nested_soft_toy_k0"));
        let k1 = Atom::var(symbol!("nested_soft_toy_k1"));
        let p0 = Atom::var(symbol!("nested_soft_toy_p0"));
        let p1 = Atom::var(symbol!("nested_soft_toy_p1"));
        let q0 = Atom::var(symbol!("nested_soft_toy_q0"));
        let q1 = Atom::var(symbol!("nested_soft_toy_q1"));
        let mass = Atom::var(symbol!("nested_soft_toy_mass"));
        let m_uv = Atom::var(symbol!("nested_soft_toy_m_uv"));

        let square = |time: &Atom, space: &Atom| time.pow(2) - space.pow(2);
        let k_squared = square(&k0, &k1);
        let p_squared = square(&p0, &p1);
        let q_squared = square(&q0, &q1);
        let r0 = &k0 - &p0;
        let r1 = &k1 - &p1;
        let r_squared = square(&r0, &r1);
        let p_plus_q_squared = square(&(&p0 + &q0), &(&p1 + &q1));
        let k_dot_p = &k0 * &p0 - &k1 * &p1;
        let p_dot_q = &p0 * &q0 - &p1 * &q1;

        let v_k = &k_squared - mass.pow(2);
        let v_r = &r_squared - mass.pow(2);
        let v_p = &p_squared - mass.pow(2);
        let v_p_plus_q = &p_plus_q_squared - mass.pow(2);
        let u_k = &k_squared - m_uv.pow(2);
        let u_r = &r_squared - m_uv.pow(2);
        let u_p = &p_squared - m_uv.pow(2);

        let child = (&k0 + &mass) / (&v_k * &v_r);
        let completed_child =
            (&k0 + &mass) / v_k.pow(2) + Atom::num(2) * &k0 * &k_dot_p / u_k.pow(3);
        let child_uv_leading = &k0 / (&u_k * &u_r);
        let completed_child_uv_leading =
            &k0 / u_k.pow(2) + Atom::num(2) * &k0 * &k_dot_p / u_k.pow(3);
        let physical_parent_kernel =
            Atom::one() / v_p.pow(2) - Atom::num(2) * &p_dot_q / v_p.pow(3);
        let uv_parent_quadratic =
            Atom::num(4) * p_dot_q.pow(2) / u_p.pow(4) - &q_squared / u_p.pow(3);

        // These are the four forest families.  In particular the nested term
        // applies the parent operator to the complete H_child atom and to its
        // freshly re-expanded ordinary-UV leading coefficient.
        let f_empty = &p0 * &child / (&v_p * &v_p_plus_q);
        let f_child = -&p0 * &completed_child / (&v_p * &v_p_plus_q);
        let f_parent = -&p0 * &child * &physical_parent_kernel
            - &p0 * &child_uv_leading * &uv_parent_quadratic;
        let f_nested = &p0 * &completed_child * &physical_parent_kernel
            + &p0 * &completed_child_uv_leading * &uv_parent_quadratic;

        for (label, term) in [
            ("empty", &f_empty),
            ("child", &f_child),
            ("parent", &f_parent),
            ("nested", &f_nested),
        ] {
            assert!(!term.is_zero(), "the {label} forest term must be retained");
        }

        let child_remainder = &child - &completed_child;
        let expected =
            &p0 * &child_remainder * (Atom::one() / (&v_p * &v_p_plus_q) - &physical_parent_kernel)
                + &p0 * (&completed_child_uv_leading - &child_uv_leading) * &uv_parent_quadratic;
        paper_assert_normalized_zero(
            "four nested toy forest terms",
            (f_empty + f_child + f_parent + f_nested - &expected).cancel(),
        );

        // The child H remainder and the difference between the two parent-UV
        // leading coefficients both start at k^-5.  This is the termwise
        // subdivergence condition needed before the outer forest pair is summed.
        let lambda_symbol = symbol!("nested_soft_toy_lambda");
        let lambda = Atom::var(lambda_symbol);
        for (mass_label, mass_replacement) in
            [("generic mass", None), ("massless", Some(Atom::Zero))]
        {
            for (label, remainder) in [
                ("physical child", child_remainder.clone()),
                (
                    "parent UV coefficient",
                    &child_uv_leading - &completed_child_uv_leading,
                ),
            ] {
                let remainder = match &mass_replacement {
                    Some(value) => remainder.replace(mass.clone()).with(value.clone()),
                    None => remainder,
                };
                let hard_scaled = remainder
                    .replace(k0.clone())
                    .with(&k0 / &lambda)
                    .replace(k1.clone())
                    .with(&k1 / &lambda)
                    .together()
                    .cancel();
                let hard_series = hard_scaled.series(lambda_symbol, Atom::Zero, 5).unwrap();
                assert_eq!(
                    hard_series.get_trailing_exponent(),
                    Rational::from(5),
                    "{mass_label}, {label}"
                );
            }

            // Each parent remainder carries the two extra powers supplied by
            // the q Taylor subtraction, hence p_0 times either kernel is p^-5.
            let parent_remainders = [
                (
                    "physical parent kernel",
                    &p0 * (Atom::one() / (&v_p * &v_p_plus_q) - &physical_parent_kernel),
                ),
                ("ordinary parent kernel", &p0 * &uv_parent_quadratic),
            ];
            for (label, remainder) in parent_remainders {
                let remainder = match &mass_replacement {
                    Some(value) => remainder.replace(mass.clone()).with(value.clone()),
                    None => remainder,
                };
                let hard_scaled = remainder
                    .replace(p0.clone())
                    .with(&p0 / &lambda)
                    .replace(p1.clone())
                    .with(&p1 / &lambda)
                    .together()
                    .cancel();
                let hard_series = hard_scaled.series(lambda_symbol, Atom::Zero, 5).unwrap();
                assert_eq!(
                    hard_series.get_trailing_exponent(),
                    Rational::from(5),
                    "{mass_label}, {label}"
                );
            }

            // The two parent kernels vanish quadratically at the soft point.
            for (label, kernel) in [
                (
                    "physical parent kernel",
                    Atom::one() / (&v_p * &v_p_plus_q) - &physical_parent_kernel,
                ),
                ("ordinary parent kernel", uv_parent_quadratic.clone()),
            ] {
                let kernel = match &mass_replacement {
                    Some(value) => kernel.replace(mass.clone()).with(value.clone()),
                    None => kernel,
                };
                let soft_scaled = kernel
                    .replace(q0.clone())
                    .with(&q0 * &lambda)
                    .replace(q1.clone())
                    .with(&q1 * &lambda)
                    .together()
                    .cancel();
                let soft_series = soft_scaled.series(lambda_symbol, Atom::Zero, 2).unwrap();
                assert_eq!(
                    soft_series.get_trailing_exponent(),
                    Rational::from(2),
                    "{mass_label}, {label}"
                );
            }

            // On a simultaneous two-loop ray, the apparent Lambda^-8 pieces
            // of the two regrouped terms cancel.  The complete forest starts
            // at Lambda^-9, as required by the outer degree-two projection.
            let complete_forest = match &mass_replacement {
                Some(value) => expected.replace(mass.clone()).with(value.clone()),
                None => expected.clone(),
            };
            let simultaneous_scaled = complete_forest
                .replace(k0.clone())
                .with(&k0 / &lambda)
                .replace(k1.clone())
                .with(&k1 / &lambda)
                .replace(p0.clone())
                .with(&p0 / &lambda)
                .replace(p1.clone())
                .with(&p1 / &lambda)
                .together()
                .cancel();
            let simultaneous_series = simultaneous_scaled
                .series(lambda_symbol, Atom::Zero, 9)
                .unwrap();
            assert_eq!(
                simultaneous_series.get_trailing_exponent(),
                Rational::from(9),
                "{mass_label}, simultaneous two-loop hard ray"
            );
        }
    }

    #[test]
    fn rational_shell_groups_mixed_multisets_without_expanding_numerators() -> Result<()> {
        test_initialise()?;
        let (a, b, c, x) = symbol!(
            "uv_bucket::a",
            "uv_bucket::b",
            "uv_bucket::c",
            "uv_bucket::x"
        );
        let denominator = GS.den(0, function!(GS.emr_mom, 0), 1, Atom::var(x));
        let factorized = (Atom::var(a) + b).pow(5) * (Atom::var(c) + x);
        let expression = (Atom::var(a) / denominator.pow(2) + Atom::var(b) / denominator.pow(3))
            * &factorized
            + Atom::var(c) / denominator.pow(2);
        let square = (Atom::var(a) / denominator.pow(2) + Atom::var(b) / denominator.pow(3)).pow(2);
        let other = GS.den(1, function!(GS.emr_mom, 1), 1, Atom::var(x));
        let mut first =
            FourDTerm::from_view((Atom::var(a) / (&denominator * &other)).as_view())?.remove(0);
        let mut opposite = first.clone();
        opposite.denominators.reverse();
        opposite.numerator = -&opposite.numerator;
        let cancelled = FourDTerm::group([first.clone(), opposite]);
        first.numerator = factorized.clone();
        for (terms, expected) in [
            (
                FourDTerm::from_view(expression.as_view())?,
                (Atom::var(a) * &factorized + c) / denominator.pow(2)
                    + Atom::var(b) * &factorized / denominator.pow(3),
            ),
            (
                FourDTerm::from_view(square.as_view())?,
                Atom::var(a).pow(2) / denominator.pow(4)
                    + Atom::num(2) * a * b / denominator.pow(5)
                    + Atom::var(b).pow(2) / denominator.pow(6),
            ),
            (cancelled, Atom::Zero),
            (
                FourDTerm::group([first]),
                factorized / (&denominator * other),
            ),
        ] {
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
            assert_eq!(reconstructed, expected);
        }

        Ok(())
    }

    #[test]
    fn canonical_uv_classes_certify_signed_aliases_and_preserve_raw_sectors() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph canonical_signed_triangle {
            edge [num=1 mass=1]
            node [num=1]
            a -> b [id=0 lmb_id=0]
            b -> c [id=1]
            c -> a [id=2]
        })?;
        let owners = graph.full_filter();
        let hard = |owner| function!(GS.emr_mom, owner);
        let component = |owner| function!(GS.emr_mom, owner, GS.cind(0));
        let tagged = |owner, sign, role| {
            function!(
                GS.emr_mom,
                GS.uv_momentum_provenance_tag(owner, role, Atom::num(sign) * hard(owner)),
                GS.cind(0)
            )
        };
        let den = |owner, sign, mass| {
            GS.den(
                owner,
                Atom::num(sign) * hard(owner),
                mass,
                component(owner).pow(2) - mass,
            )
        };
        let positive = den(0, 1, 1);
        let negative = den(1, -1, 1);
        let fixed = tagged(0, 1, UvMomentumProvenanceRole::TaylorFixed);
        let derived = tagged(1, -1, UvMomentumProvenanceRole::DenominatorDerived);
        let soft = tagged(2, 1, UvMomentumProvenanceRole::DenominatorDerivedSoft);
        let physical = tagged(2, 1, UvMomentumProvenanceRole::PhysicalSourceFixed);
        let spectator = &soft + &physical;
        let atom = (&fixed + &spectator) / &positive + (-derived + &spectator) / &negative;
        let sector = FourDSector::new(
            atom.clone(),
            vec![(owners.clone(), owners, graph.loop_momentum_basis.clone())],
            vec![],
        );
        let projection = sector.canonical_projection(&graph)?;
        assert_eq!(sector.atom, atom);
        let expected =
            Atom::num(2) * (Atom::num(2) * component(0) + soft) / (component(0).pow(2) - 1);
        let unequal = FourDSector::new(
            fixed / positive + tagged(1, 1, UvMomentumProvenanceRole::TaylorFixed) / den(1, 1, 2),
            sector.active_components.clone(),
            vec![],
        )
        .canonical_projection(&graph)?;
        let unequal_expected =
            component(0) / (component(0).pow(2) - 1) + component(0) / (component(0).pow(2) - 2);
        for (projection, expected) in [(projection, expected), (unequal, unequal_expected)] {
            let reconstructed =
                projection
                    .terms
                    .iter()
                    .try_fold(Atom::Zero, |sum, term| -> Result<Atom> {
                        let denominator =
                            term.powers
                                .iter()
                                .fold(Atom::one(), |product, (id, power)| {
                                    product * projection.class(*id).full_expr.pow(*power)
                                });
                        Ok(sum
                            + projection.neutral_numerator(&term.numerator, &graph)? / denominator)
                    })?;
            assert!(
                (reconstructed - expected)
                    .expand_num()
                    .collect_factors()
                    .is_zero()
            );
        }

        Ok(())
    }

    #[test]
    fn canonical_uv_signed_multiloop_carriers_preserve_factorized_numerators() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph canonical_signed_sunset {
            edge [num=1 mass=1]
            node [num=1]
            a -> b [id=0 lmb_id=0]
            a -> b [id=1 lmb_id=1]
            a -> b [id=2]
        })?;
        let hard = -function!(GS.emr_mom, 0) - function!(GS.emr_mom, 1);
        let component = GS.indexed_momentum(&hard, &[GS.cind(0)]);
        let opaque =
            (Atom::var(symbol!("signed_multiloop::a")) + symbol!("signed_multiloop::b")).pow(5);
        let tagged = |sign| {
            function!(
                GS.emr_mom,
                GS.uv_momentum_provenance_tag(
                    2,
                    UvMomentumProvenanceRole::TaylorFixed,
                    Atom::num(sign) * &hard,
                ),
                GS.cind(0)
            )
        };
        let positive = GS.den(2, &hard, 1, component.pow(2) - 1);
        let negative = GS.den(2, -&hard, 1, (-&component).pow(2) - 1);
        let atom = tagged(1) * &opaque / positive.pow(2) - tagged(-1) * &opaque / negative.pow(2);
        let owners = graph.full_filter();
        let sector = FourDSector::new(
            atom,
            vec![(owners.clone(), owners, graph.loop_momentum_basis.clone())],
            vec![],
        );
        let canonical = sector.canonical_projection(&graph)?;
        let reconstructed =
            canonical
                .terms
                .iter()
                .try_fold(Atom::Zero, |sum, term| -> Result<Atom> {
                    let denominator = term
                        .powers
                        .iter()
                        .fold(Atom::one(), |product, (id, power)| {
                            product * canonical.class(*id).full_expr.pow(*power)
                        });
                    Ok(sum + canonical.neutral_numerator(&term.numerator, &graph)? / denominator)
                })?;
        let expected = Atom::num(2) * &component * opaque / (component.pow(2).expand() - 1).pow(2);
        assert_eq!(
            reconstructed.collect_factors(),
            canonical
                .neutral_numerator(&expected, &graph)?
                .collect_factors()
        );
        let base = Atom::var(symbol!("signed_multiloop::a")) - symbol!("signed_multiloop::b");
        let factor =
            (Atom::var(symbol!("signed_multiloop::c")) + symbol!("signed_multiloop::d")).pow(5);
        let left = Atom::num(2) * &factor * &base + Atom::num(3) * &factor;
        let right = -Atom::num(2) * &factor * (-&base).expand_num() + Atom::num(3) * &factor;
        let normalized = canonical.neutral_numerator(&left, &graph)?;
        assert_eq!(normalized, canonical.neutral_numerator(&right, &graph)?);
        assert_eq!(
            normalized,
            canonical.neutral_numerator(&normalized, &graph)?
        );
        for exponent in [-3_i64, 2, 3, 500_000] {
            let first = canonical.neutral_numerator(&base.pow(exponent), &graph)?;
            let opposite = Atom::num(if exponent % 2 == 0 { 1 } else { -1 })
                * (-&base).expand_num().pow(exponent);
            assert_eq!(first, canonical.neutral_numerator(&opposite, &graph)?);
            assert_eq!(first, canonical.neutral_numerator(&first, &graph)?);
        }
        let fractional = base.pow(Atom::num(1) / Atom::num(2));
        assert_eq!(
            fractional,
            canonical.neutral_numerator(&fractional, &graph)?
        );
        Ok(())
    }

    #[test]
    fn canonical_uv_positive_blocks_and_absent_poles_keep_their_roles() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph canonical_positive_triangle {
            edge [num=1 mass=1]
            node [num=1]
            a -> b [id=0 lmb_id=0]
            b -> c [id=1]
            c -> a [id=2]
        })?;
        let owners = graph.full_filter();
        let den = |owner| {
            GS.den(
                owner,
                function!(GS.emr_mom, owner),
                1,
                function!(GS.emr_mom, owner, GS.cind(0)).pow(2) - 1,
            )
        };
        let tag = GS.uv_momentum_provenance_tag(
            2,
            UvMomentumProvenanceRole::TaylorFixed,
            function!(GS.emr_mom, 2),
        );
        let tagged = function!(GS.emr_mom, tag, GS.cind(0));
        let bindings = vec![(owners.clone(), owners, graph.loop_momentum_basis.clone())];
        let projection = FourDSector::new((den(0) + 1) / den(1).pow(2), bindings.clone(), vec![])
            .canonical_projection(&graph)?;
        let q = function!(GS.emr_mom, 0, GS.cind(0));
        let expected = q.pow(2) / (q.pow(2) - 1).pow(2);
        let no_pole = FourDSector::new(tagged, bindings, vec![]).canonical_projection(&graph)?;
        for (projection, expected) in [(projection, expected), (no_pole, q)] {
            let reconstructed =
                projection
                    .terms
                    .iter()
                    .try_fold(Atom::Zero, |sum, term| -> Result<Atom> {
                        let denominator =
                            term.powers
                                .iter()
                                .fold(Atom::one(), |product, (id, power)| {
                                    product * projection.class(*id).full_expr.pow(*power)
                                });
                        Ok(sum
                            + projection.neutral_numerator(&term.numerator, &graph)? / denominator)
                    })?;
            assert_eq!(reconstructed, expected);
        }

        Ok(())
    }

    #[test]
    fn canonical_positive_denominator_recovers_one_class_from_expanded_coordinates() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph canonical_positive_sunset {
            edge [num=1 mass=1]
            node [num=1]
            a -> b [id=0 lmb_id=0]
            a -> b [id=1 lmb_id=1]
            a -> b [id=2]
        })?;
        let owners = graph.full_filter();
        let bindings = vec![(owners.clone(), owners, graph.loop_momentum_basis.clone())];
        let momentum = function!(GS.emr_mom, 2);
        let q = |owner| function!(GS.emr_mom, owner, GS.cind(0));
        // Expand only this denominator polynomial; the numerator's outer
        // positive block remains opaque throughout projection.
        let expanded = (q(0) + q(1)).pow(2).expand() - 1;
        let positive = GS.den(2, &momentum, 1, &expanded);
        let negative = GS.den(2, &momentum, 1, q(2).pow(2) - 1);
        let projection =
            FourDSector::new((positive + 1) / negative.pow(2), bindings.clone(), vec![])
                .canonical_projection(&graph)?;
        let reconstructed =
            projection
                .terms
                .iter()
                .try_fold(Atom::Zero, |sum, term| -> Result<Atom> {
                    let denominator =
                        term.powers
                            .iter()
                            .fold(Atom::one(), |product, (id, power)| {
                                product * projection.class(*id).full_expr.pow(*power)
                            });
                    Ok(sum + projection.neutral_numerator(&term.numerator, &graph)? / denominator)
                })?;
        // Normalize numeric coefficients while keeping numerator products and powers factorized.
        assert_eq!(
            reconstructed.expand_num(),
            ((&expanded + 1) / expanded.pow(2)).expand_num()
        );

        let inconsistent = GS.den(2, momentum, 1, q(0).pow(2) - 1);
        let error = FourDSector::new((&inconsistent + 1) / inconsistent.pow(2), bindings, vec![])
            .canonical_projection(&graph)
            .unwrap_err();
        assert!(
            error
                .to_string()
                .contains("not a polynomial solely in its certified momentum channel")
        );
        Ok(())
    }

    #[test]
    fn projection_sector_grouping_preserves_frozen_domains_and_recursion() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph canonical_sector_triangle {
            edge [num=1 mass=1]
            node [num=1]
            a -> b [id=0 lmb_id=0]
            b -> c [id=1]
            c -> a [id=2]
        })?;
        let owners = graph.full_filter();
        let bindings = vec![(owners.clone(), owners, graph.loop_momentum_basis.clone())];
        let first = FourDSector::new(Atom::num(2), bindings.clone(), vec![]);
        let second = FourDSector::new(Atom::num(3), bindings.clone(), vec![]);
        let frozen = FourDSector::new(
            Atom::num(7),
            bindings,
            vec![graph.loop_momentum_basis.clone()],
        );
        let completion = FourDSector::new(Atom::num(11), vec![], vec![]);
        let local = Local4dCts(FourDSectors::new(
            vec![first.clone(), second, frozen.clone()],
            vec![completion.clone()],
        ));
        let projection = local.projection_sectors();
        let projected_value = projection.iter().fold(Atom::Zero, |sum, sector| {
            let frozen_factor = sector.frozen_lmbs.iter().fold(Atom::one(), |product, lmb| {
                product * GS.localizing_integrand(lmb)
            });
            sum + &sector.atom * frozen_factor
        });
        assert_eq!(
            projected_value,
            Atom::num(5) + Atom::num(7) * GS.localizing_integrand(&graph.loop_momentum_basis)
        );
        assert_eq!(local.atom(), &Atom::num(23));
        assert_eq!(
            Full4dCts::from_factorized_local(&local)
                .0
                .all()
                .fold(Atom::Zero, |sum, sector| sum + &sector.atom),
            Atom::num(23)
        );
        Ok(())
    }

    #[test]
    fn cancelling_projection_sectors_preserve_final_cut_orders_without_energy_maps() -> Result<()> {
        use crate::{
            cff::esurface::RaisedEsurfaceGroup,
            graph::cuts::CutSet,
            settings::global::OrientationPattern,
            uv::approx::{
                OrientationProjection, final_integrand::FinalIntegrandBuilder,
                integrated::IntegratedCts, local_3d::Localizer,
                projected_4d::Projected4dApproximation,
            },
        };

        test_initialise()?;
        let mut graph: Graph = dot!(digraph cancelling_projection_sectors {
            edge [num=1 mass=1];
            node [num=1];
            a -> b [id=0 lmb_id=0];
            a -> b [id=1];
        })?;
        let owners = graph.full_filter();
        let current = OwnedForestNode {
            spinney: Spinney::new(
                InternalSubGraph::cleaned_filter_optimist(owners.clone(), graph.as_ref()),
                &graph,
                &graph.loop_momentum_basis,
            )
            .expect("UV spinney construction succeeds")
            .expect("the bubble has a compatible UV spinney"),
            topo_order: 1,
        };
        let coefficient = (GS.emr_mom(EdgeIndex(0), GS.cind(0)).pow(2) + Atom::one())
            / graph.denominator(&owners, |_| 1);
        let bindings = vec![(owners.clone(), owners, graph.loop_momentum_basis.clone())];
        let cancelled = Local4dCts(FourDSectors::new(
            vec![
                FourDSector::new(coefficient.clone(), bindings.clone(), Vec::new()),
                FourDSector::new(-coefficient, bindings, Vec::new()),
            ],
            Vec::new(),
        ));
        let zero = Local4dCts(FourDSectors::new(Vec::new(), Vec::new()));
        let mut cutset = CutSet::empty(graph.n_hedges());
        // The complete zero result must retain both allowed cut orders even
        // when cancellation leaves no energy contour to generate or select.
        cutset.residue_selector.left_th_cut = Some(RaisedEsurfaceGroup {
            esurface_ids: Vec::new(),
            max_occurence: 2,
        });
        let production = Default::default();
        let pattern = OrientationPattern::default();
        let options = graph.denominator_only_cff_3d_expression_options();
        let localizer = Localizer::new(
            &cutset,
            OrientationProjection::exact(&production, &options, &pattern, true),
        );
        let settings = UVgenerationSettings {
            local_uv_cts_from_expanded_4d_integrands: true,
            ..Default::default()
        };
        for local in [cancelled, zero] {
            let projected = Projected4dApproximation::new(localizer, &mut graph, &settings)
                .project_local_4d(&local, &mut Local4dProjectionContext::default())?;
            let finalized = FinalIntegrandBuilder::new(localizer, &settings).build_projected(
                &mut graph,
                &current,
                &projected,
                &IntegratedCts::root(),
            )?;
            assert_eq!(
                finalized.into_integrands(),
                cutset
                    .residue_selector
                    .generate_allowed_keys()
                    .into_iter()
                    .map(|index| (index, Atom::Zero))
                    .collect(),
            );
        }
        Ok(())
    }

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
            .expect("the scalar bubble has a compatible UV loop-momentum basis")
            .expect("the scalar bubble is retained"),
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
        for scheme in [
            ApproximationType::MUV,
            ApproximationType::PolePart,
            ApproximationType::IR,
        ] {
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

        current.spinney.renormalization_scheme = ApproximationType::PolePart;
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
            .expect("the self-energy bubble has a compatible child LMB")
            .expect("the self-energy bubble is retained"),
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
                .expect("the GL24 UV triangle has a compatible child LMB")
                .expect("the GL24 UV triangle is retained"),
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
            .expect("requested coefficient is within series precision")
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
        let (production, _, _) = t(
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
        // Compare coefficients of the unchanged quartic numerator: production
        // keeps it outside the Laurent-layer sum, while the scalar oracle
        // repeats it in each layer. Neither numerator needs to be expanded.
        let hard_q1 = GS.uv_momentum_provenance_tag(
            Atom::num(usize::from(owners[0])),
            UvMomentumProvenanceRole::TaylorFixed,
            function!(GS.emr_mom, usize::from(owners[0])),
        );
        let numerator_keys = [1, 2].map(|index| {
            function!(
                GS.emr_mom,
                &hard_q1,
                minkowski.to_symbolic([Atom::num(index)])
            )
        });
        assert_eq!(
            production.coefficient_list::<u32>(&numerator_keys),
            scalar_series_oracle.coefficient_list::<u32>(&numerator_keys),
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
            .expect("the source triangle has a compatible loop-momentum basis")
            .expect("the source triangle is retained"),
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
        let (expanded, _, _) = t(
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
        let mut context = Local4dProjectionContext::default();
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
            let preparation =
                graph.prepare_3d_expression_for_4d_term(source, &options, &term.numerator, &[])?;
            graph.generate_3d_expression_for_4d_term(&preparation, Some(&mut context))?;
        }
        // Undotted and dotted sources have different occurrence counts;
        // compatible dotted owner relabellings may share a canonical CFF.
        // Regardless of reuse, each request must preserve its complete residue.
        let cutset = crate::graph::cuts::CutSet::empty(graph.n_hedges());
        let mut has_nonzero_residue = false;
        for (denominators, term) in active_denominators.iter().zip(&terms) {
            let mut values = Vec::new();
            for generation_cache in [Some(&mut context), None] {
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
                            * term.map_exact_source_numerator(&orientation.orientation, None)?;
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
    #[test]
    fn early_color_simplification_preserves_open_and_nested_numerator_boundaries() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph nested_open_color_bubble {
            edge [num=1 mass=1]
            node [num=1]
            a -> b [id=0 lmb_id=0]
            a -> b [id=1]
            b -> c [id=2 lmb_id=1]
            c -> a [id=3]
        })?;
        let [a, b, c, d, u, v]: [Atom; 6] =
            std::array::from_fn(|i| idenso::coad!(8, Atom::from(Aind::Normal(80 + i))));
        let spectator = (GS.emr_mom(EdgeIndex(0), GS.cind(0)) + Atom::one())
            * (GS.emr_mom(EdgeIndex(1), GS.cind(0)) + Atom::num(2));
        graph.underlying[EdgeIndex(0)].num.value =
            idenso::color_f!(&a, &u, &v) * (GS.emr_mom(EdgeIndex(0), GS.cind(0)) + Atom::one());
        graph.underlying[EdgeIndex(1)].num.value =
            idenso::color_f!(&b, &u, &v) * (GS.emr_mom(EdgeIndex(1), GS.cind(0)) + Atom::num(2));
        graph.underlying[EdgeIndex(2)].num.value = spenso::g!(&a, &c);
        graph.underlying[EdgeIndex(3)].num.value = spenso::g!(&b, &d);
        let bubble = graph
            .get_edge_subgraph(EdgeIndex(0))
            .union(&graph.get_edge_subgraph(EdgeIndex(1)));
        let nodes = [bubble, graph.full_filter()].map(|filter| OwnedForestNode {
            spinney: Spinney::with_scheme(
                InternalSubGraph::cleaned_filter_optimist(filter, graph.as_ref()),
                &graph,
                &graph.loop_momentum_basis,
                ApproximationType::MUV,
                0,
            )
            .expect("the nested bubble has a valid internal UV topology")
            .expect("the nested bubble has a compatible loop-momentum basis"),
            topo_order: 0,
        });
        let [inner, outer] = &nodes;
        let empty = OwnedForestNode {
            spinney: Spinney::empty(&graph),
            topo_order: 0,
        };
        let settings = UVgenerationSettings::default();
        let ctx = UVCtx::new(&graph, &settings);
        let raw_inner = graph
            .numerator(inner.subgraph(), empty.subgraph())
            .get_single_atom()?;
        let casimir = idenso::color_cas!(2, &ColorAdjoint {}.new_rep(8));
        let expected_inner = &casimir * spenso::g!(&a, &b) * &spectator;
        // The same existing numerator operation is used by direct 3D. The
        // complete open metric, including both boundary labels, stays explicit.
        let prepared_inner = graph
            .numerator(inner.subgraph(), empty.subgraph())
            .color_simplify()
            .get_single_atom()?;
        assert_eq!(prepared_inner, expected_inner);
        let grown_inner = grow(&Full4dCts::root(), &ctx, inner, &empty)?;
        assert_eq!(
            &grown_inner * graph.denominator(inner.subgraph(), |_| 1),
            expected_inner
        );

        let raw_remainder = graph
            .numerator(&outer.reduced_subgraph(inner), inner.subgraph())
            .get_single_atom()?;
        let raw_full = graph
            .numerator(outer.subgraph(), empty.subgraph())
            .get_single_atom()?;
        assert_eq!(&raw_inner * &raw_remainder, raw_full);
        let inner_input = Full4dCts(FourDSectors::active_atom(grown_inner));
        let nested =
            grow(&inner_input, &ctx, outer, inner)? * graph.denominator(outer.subgraph(), |_| 1);
        let late = raw_full.simplify_color().simplify_metrics();
        assert_eq!(nested, late);
        assert_eq!(nested, casimir * spenso::g!(&c, &d) * &spectator);
        // Closing the remaining boundary commutes with both preparations;
        // the factorized momentum spectator is never distributed.
        let closure = spenso::g!(&c, &d);
        assert_eq!(
            (&nested * &closure).simplify_color(),
            (raw_full * closure).simplify_color()
        );
        assert_eq!(
            graph
                .numerator(inner.subgraph(), empty.subgraph())
                .get_single_atom()?,
            raw_inner
        );
        Ok(())
    }

    #[test]
    fn local_taylor_retains_dirac_traces_with_or_without_analytic_addbacks() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph factorized_local_spin {
            edge [num=1 mass=1]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]
            incoming -> a [id=0]
            a -> b [id=1 lmb_id=0]
            b -> a [id=2]
            b -> outgoing [id=3]
        })?;
        let filter = graph
            .get_edge_subgraph(EdgeIndex(1))
            .union(&graph.get_edge_subgraph(EdgeIndex(2)));
        let current = OwnedForestNode {
            spinney: Spinney::with_scheme(
                InternalSubGraph::cleaned_filter_optimist(filter, graph.as_ref()),
                &graph,
                &graph.loop_momentum_basis,
                ApproximationType::MUV,
                0,
            )?
            .expect("the bubble has a compatible UV loop-momentum basis"),
            topo_order: 0,
        };
        let given = OwnedForestNode {
            spinney: Spinney::empty(&graph),
            topo_order: 0,
        };
        let index = spenso::mink!(4, Atom::from(Aind::Normal(23)));
        let trace = spenso::trace!(
            Bispinor {}.new_rep(4).to_symbolic([]),
            idenso::gamma!(&index),
            idenso::gamma!(&index),
        );
        let (a, b, c, d) = symbol!(
            "local_spin_a",
            "local_spin_b",
            "local_spin_c",
            "local_spin_d"
        );
        let spectator = (Atom::var(a) + b) * (Atom::var(c) + d);
        graph.underlying[EdgeIndex(1)].num.value = &spectator * trace;
        let input = Full4dCts(FourDSectors::active_atom(Atom::one()));
        let mut previous = None;
        for generate_integrated in [false, true] {
            let settings = UVgenerationSettings {
                generate_integrated,
                ..Default::default()
            };
            let local = uv_limit(
                &input,
                &UVCtx::new(&graph, &settings),
                &current,
                &given,
                &current,
                &given,
            )?;
            assert!(local.atom().contains_symbol(AGS.gamma));
            assert!(
                local
                    .atom()
                    .contains_symbol(spenso::network::tags::SPENSO_TAG.trace)
            );
            assert!(
                local
                    .atom()
                    .pattern_match(&spectator.to_pattern(), None, None)
                    .next()
                    .is_some()
            );
            if let Some(previous) = &previous {
                assert_eq!(local.atom(), previous);
            }
            previous = Some(local.atom().clone());
        }
        Ok(())
    }
}
