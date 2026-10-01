//! Explicit, offline certificates for changes made by the ordinary-U optimizer.
//! The physical FullH generation below is an input inventory, not a proof of
//! the leading soft behavior of the complete orientation/forest sum.

use std::{
    collections::{BTreeMap, BTreeSet},
    sync::{Arc, Mutex},
};

use super::super::branches::DirectResidueKey;
use super::*;
use crate::{
    cff::expression::OrientationID,
    graph::{FeynmanGraph, Graph, parse::IntoGraph},
    initialisation::test_initialise,
    integrands::process::param_builder::FnMapEntry,
    model::InputParamCard,
    numerator::pole_certificate::{ScalarPoleCertificate, ScalarPoleFailure},
    settings::global::{GenerationSettings, ThresholdSubtractionSettings},
    uv::{
        CTIdentifier, CTRenormalizationRule, RenormalizationPrescriptionSettings,
        UVgenerationSettings, approx::CutStructure, hedge_poset::Wood,
    },
};

const MAX_CAPTURE_BYTES: usize = 128 * 1024 * 1024;
const MAX_CAPTURES: usize = 4096;
static RECORDING: Mutex<Option<Recording>> = Mutex::new(None);

struct Capture {
    current: String,
    given: String,
    active: Option<String>,
    loop_edges: String,
    key: DirectResidueKey,
    cuts: Vec<String>,
    numerators: Vec<Arc<FnMapEntry>>,
    before: Atom,
    after: Atom,
}

struct Recording {
    graph: String,
    captures: Vec<Capture>,
    bytes: usize,
    dropped: usize,
}

struct RecordingGuard;

impl RecordingGuard {
    fn start(graph: &str) -> Self {
        let mut recording = RECORDING.lock().unwrap();
        assert!(
            recording.is_none(),
            "only one pole audit may record at a time"
        );
        *recording = Some(Recording {
            graph: graph.to_owned(),
            captures: Vec::new(),
            bytes: 0,
            dropped: 0,
        });
        Self
    }

    fn finish(self) -> Recording {
        let recording = RECORDING.lock().unwrap().take().unwrap();
        drop(self);
        recording
    }
}

impl Drop for RecordingGuard {
    fn drop(&mut self) {
        RECORDING.lock().unwrap().take();
    }
}

// This entire module and its call site are absent from production builds.
// A graph-filtered global recorder also sees forest work on Rayon workers.
#[allow(clippy::too_many_arguments)]
pub(super) fn record<S: ForestNodeLike>(
    ctx: &UVCtx<'_>,
    current: &S,
    given: &S,
    key: &DirectResidueKey,
    series: &DirectResidueBranches,
    before: &Atom,
    after: &Atom,
    active: Option<&SuBitGraph>,
    lmb: &LoopMomentumBasis,
) {
    let mut recording = RECORDING.lock().unwrap();
    let Some(recording) = recording.as_mut().filter(|r| r.graph == ctx.graph.name) else {
        return;
    };
    let integrands = series
        .iter_keys()
        .find(|(candidate, _)| *candidate == key)
        .map(|(_, integrands)| integrands)
        .expect("the recorded residue belongs to this series");
    let bytes = before.as_view().get_byte_size()
        + after.as_view().get_byte_size()
        + integrands
            .numerators()
            .iter()
            .map(|entry| entry.rhs.as_view().get_byte_size())
            .sum::<usize>();
    if recording.captures.len() == MAX_CAPTURES
        || bytes > MAX_CAPTURE_BYTES.saturating_sub(recording.bytes)
    {
        recording.dropped += 1;
        return;
    }
    recording.bytes += bytes;
    recording.captures.push(Capture {
        current: current.subgraph().string_label(),
        given: given.subgraph().string_label(),
        active: active.map(|active| active.string_label()),
        loop_edges: format!("{:?}", lmb.loop_edges),
        key: key.clone(),
        // fallible_map does not expose the cut key. Retain every matching key;
        // this fixture has exactly one empty cut, checked below.
        cuts: integrands
            .iter()
            .filter(|(_, atom)| *atom == before)
            .map(|(cut, _)| format!("{cut:?}"))
            .collect(),
        numerators: integrands.numerators().to_vec(),
        before: before.clone(),
        after: after.clone(),
    });
}

impl Capture {
    /// Mirror the production family-shape gate before exact collection. In
    /// particular, never distribute products of graph-numerator sums merely
    /// to make a diagnostic proof succeed.
    fn groups(atom: &Atom, keys: &[Atom]) -> Option<BTreeMap<Atom, Atom>> {
        if atom.as_view().get_byte_size() > 256 * 1024
            || keys.len() > 128
            || [
                OrientationID::symbol(),
                GS.theta,
                GS.orientation_delta,
                Symbol::IF,
            ]
            .into_iter()
            .any(|symbol| atom.contains_symbol(symbol))
        {
            return None;
        }
        let contains_key = |part: AtomView<'_>| keys.iter().any(|key| part.contains(key));
        let literal = |part: AtomView<'_>| keys.iter().any(|key| key.as_view() == part);
        let monomial_factor = |part: AtomView<'_>| {
            literal(part)
                || matches!(part, AtomView::Pow(power)
                    if literal(power.get_base())
                    && i64::try_from(power.get_exp()).is_ok_and(|n| (0..=64).contains(&n)))
        };
        let mut eligible = keys
            .iter()
            .all(|key| matches!(key.as_view(), AtomView::Var(_) | AtomView::Fun(_)));
        let mut monomials = 1usize;
        atom.visitor(&mut |part| {
            if !eligible || literal(part) || !contains_key(part) {
                return false;
            }
            match part {
                AtomView::Add(sum) => {
                    monomials = monomials.saturating_add(sum.get_nargs().saturating_sub(1));
                    eligible = monomials <= 128;
                }
                AtomView::Mul(product) => {
                    let factors = product
                        .iter()
                        .filter(|part| contains_key(*part))
                        .collect_vec();
                    eligible = factors.len() <= 1 || factors.into_iter().all(monomial_factor);
                }
                AtomView::Pow(_) => eligible = monomial_factor(part),
                _ => eligible = false,
            }
            eligible
        });
        eligible
            .then(|| atom.coefficient_list_exact(keys))
            .flatten()
            .map(|g| g.into_iter().collect())
    }
}

/// Shared physical input for the two explicit offline diagnostics.
/// This follows model import before Feynman-rule and canonical-LMB construction.
pub(crate) fn physical_full_h_fixture() -> Result<(Graph, GenerationSettings)> {
    test_initialise()?;
    // Match `import model sm-default`: apply the restriction before deriving
    // propagator masses and Feynman rules, keeping nonzero parameters symbolic.
    let mut model = crate::utils::load_generic_model("sm");
    model.restriction = Some("default".into());
    let mut parameters = InputParamCard::from_str(
        include_str!("../../../../../../assets/models/json/sm/restrict_default.json").to_owned(),
        "json",
    )?;
    model.simplify(&mut parameters)?;
    let graph: Graph = include_str!(
        "../../../../../../tests/resources/graphs/paper_figure_b1_double_triangle_soft_ir.dot"
    )
    .into_graph(&model)?;
    assert_eq!(graph.name, "paper_figure_b1_double_triangle_soft_ir");
    assert_eq!(
        graph
            .loop_momentum_basis
            .loop_edges
            .iter()
            .map(|edge| usize::from(*edge))
            .collect_vec(),
        [6, 5, 8]
    );
    assert!(!graph.global_prefactor.projector.is_one());
    for (pair, _, edge) in graph.iter_edges_of(&graph.full_filter()) {
        if pair.is_paired() {
            assert!(
                edge.data.mass_atom().is_zero(),
                "every internal quark/gluon must be massless"
            );
        }
    }
    let generation = GenerationSettings {
        explicit_orientation_sum_only: true,
        threshold_subtraction: ThresholdSubtractionSettings {
            enable_thresholds: false,
            ..Default::default()
        },
        uv: UVgenerationSettings {
            generate_integrated: false,
            renormalization_prescription: RenormalizationPrescriptionSettings {
                overrides: vec![
                    CTRenormalizationRule::new(
                        CTIdentifier::new([-21, 21].into(), Some([1].into())),
                        ApproximationType::IR,
                    ),
                    CTRenormalizationRule::new(
                        CTIdentifier::new([-22, 22].into(), Some([1, 21].into())),
                        ApproximationType::IR,
                    ),
                ],
                ..Default::default()
            },
            ..Default::default()
        },
        ..Default::default()
    };
    graph.ensure_energy_convergent_cycles(&graph.no_dummy().subtract(&graph.initial_state_cut))?;
    Ok((graph, generation))
}

#[test]
#[ignore = "offline physical FullH algebra certificates; generates all 102 orientations without evaluators"]
fn physical_full_h_scalar_cancellation_certificates() -> Result<()> {
    let (mut graph, generation) = physical_full_h_fixture()?;
    let options = graph.production_cff_3d_expression_options(&generation)?;
    let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
    let contract_edges = graph
        .iter_edges_of(
            &graph
                .tree_edges
                .subtract(&graph.initial_state_cut)
                .subtract(&graph.external_filter::<SuBitGraph>()),
        )
        .map(|(_, edge, _)| edge)
        .collect_vec();
    let production = graph.generate_3d_expression_for_integrand(
        &contract_edges,
        &canonization,
        &options,
        None,
    )?;
    assert_eq!(production.expression.orientations.len(), 102);
    let orientation = OrientationProjection::exact_expression(
        &production,
        &options,
        &generation.orientation_pattern,
        true,
    );
    let wood = Wood::new(CutStructure::empty(&graph), &graph, &generation.uv)?;
    let soft_components = wood
        .graph
        .iter_nodes()
        .filter(|(_, _, spinney)| spinney.renormalization_scheme == ApproximationType::IR)
        .map(|(_, _, spinney)| {
            (
                spinney.dod,
                graph
                    .iter_edges_of(&spinney.subgraph)
                    .map(|(_, edge, _)| usize::from(edge))
                    .collect::<BTreeSet<_>>(),
            )
        })
        .collect_vec();
    assert!(soft_components.contains(&(2, [8, 9].into())));
    assert!(soft_components.contains(&(2, (2..10).collect())));
    let mut forests = wood.unfold();
    debug_tags!(#uv, #pole_audit;
        stage = "pole_audit_generation_start", graph = %graph.name,
        orientations = production.expression.orientations.len(),
        file.settings = ?generation.uv, model = "sm-default",
        file.fixture = include_str!("../../../../../../tests/resources/graphs/paper_figure_b1_double_triangle_soft_ir.dot"),
        file.projector = %graph.global_prefactor.projector.to_plain_string(),
        loop_edges = ?graph.loop_momentum_basis.loop_edges, soft_components = ?soft_components,
        "Generate physical FullH inputs; no evaluator or certificate work"
    );
    let guard = RecordingGuard::start(&graph.name);
    forests.compute(
        &mut graph,
        crate::utils::vakint()?,
        orientation,
        &generation.uv,
    )?;
    let recording = guard.finish();
    assert!(
        !recording.captures.is_empty(),
        "the physical graph must exercise accepted U transformations"
    );

    assert!(
        recording
            .captures
            .iter()
            .any(|capture| !capture.numerators.is_empty()
                && capture
                    .before
                    .contains_symbol(symbol!("gammalooprs::uv::numerator_family"))),
        "the audit must contain real retained SM numerator families"
    );

    let mut fully_certified_roots = 0usize;
    let mut certified_scalar_pairs = 0usize;
    let mut unsupported_scalar_pairs = 0usize;
    let mut removed_factor_pairs = 0usize;
    let mut zero_pairs = 0usize;
    let mut unchanged_pairs = 0usize;
    let mut definitions = std::collections::HashSet::new();
    for (root_id, capture) in recording.captures.iter().enumerate() {
        assert_eq!(
            capture.cuts.len(),
            1,
            "one exact empty-cut owner per captured root"
        );
        debug_tags!(#uv, #pole_audit;
            stage = "pole_audit_root", root_id, graph = %recording.graph,
            current = %capture.current, given = %capture.given, active = ?capture.active,
            loop_edges = %capture.loop_edges, residue = ?capture.key, cuts = ?capture.cuts,
            file.before = %capture.before.to_plain_string(),
            file.after = %capture.after.to_plain_string(),
            "Accepted ordinary-U root change, before rescale is set to one"
        );
        for entry in &capture.numerators {
            if definitions.insert(entry.clone()) {
                debug_tags!(#uv, #pole_audit;
                    stage = "pole_audit_retained_definition", root_id,
                    file.definition = ?entry,
                    "Retained numerator definition; never resolved for this proof"
                );
            }
        }
        let (Some(before), Some(after)) = (
            capture
                .before
                .coefficient_list_exact(&[Atom::var(GS.rescale)]),
            capture
                .after
                .coefficient_list_exact(&[Atom::var(GS.rescale)]),
        ) else {
            debug_tags!(#uv, #pole_audit; stage = "pole_audit_unproven", root_id,
                reason = "unsupported Laurent syntax or collection budget", "Root remains unproven");
            continue;
        };
        let before = before.into_iter().collect::<BTreeMap<_, _>>();
        let after = after.into_iter().collect::<BTreeMap<_, _>>();
        let powers = before
            .keys()
            .chain(after.keys())
            .cloned()
            .collect::<BTreeSet<_>>();
        let mut root_complete = true;
        for power in powers {
            let zero = Atom::Zero;
            let before = before.get(&power).unwrap_or(&zero);
            let after = after.get(&power).unwrap_or(&zero);
            if before == after {
                unchanged_pairs += 1;
                continue;
            }
            let mut keys = BTreeSet::new();
            for atom in [before, after] {
                atom.visitor(&mut |part| {
                    if let AtomView::Fun(fun) = part
                        && fun.get_symbol() == symbol!("gammalooprs::uv::numerator_family")
                    {
                        keys.insert(part.to_owned());
                        return false;
                    }
                    true
                });
            }
            let keys = keys.into_iter().collect_vec();
            let (Some(before_groups), Some(after_groups)) = (
                Capture::groups(before, &keys),
                Capture::groups(after, &keys),
            ) else {
                root_complete = false;
                debug_tags!(#uv, #pole_audit;
                    stage = "pole_audit_unproven", root_id, power = %power,
                    reason = "numerator shape, guard or grouping budget",
                    "Do not distribute unsupported graph-numerator products"
                );
                continue;
            };
            for monomial in before_groups
                .keys()
                .chain(after_groups.keys())
                .collect::<BTreeSet<_>>()
            {
                let before = before_groups.get(monomial).unwrap_or(&zero);
                let after = after_groups.get(monomial).unwrap_or(&zero);
                if before == after {
                    unchanged_pairs += 1;
                    continue;
                }
                match ScalarPoleCertificate::check(before, after) {
                    Ok(certificate) => {
                        assert!(certificate.identity_verified);
                        certified_scalar_pairs += 1;
                        removed_factor_pairs += usize::from(certificate.removed_factor_degree > 0);
                        zero_pairs += usize::from(certificate.scalar_zero);
                        debug_tags!(#uv, #pole_audit;
                            stage = "pole_audit_certificate", root_id, power = %power,
                            file.numerator_monomial = %monomial.to_plain_string(),
                            file.before = %before.to_plain_string(), file.after = %after.to_plain_string(),
                            file.removed_factor = %certificate.removed_factor.to_plain_string(),
                            identity_verified = certificate.identity_verified,
                            removed_factor_degree = certificate.removed_factor_degree,
                            scalar_zero = certificate.scalar_zero,
                            numerator_divisibility_verified = certificate.numerator_divisibility_verified,
                            before_num_terms = certificate.before_num_terms, before_den_terms = certificate.before_den_terms,
                            after_num_terms = certificate.after_num_terms, after_den_terms = certificate.after_den_terms,
                            before_den_degree = certificate.before_den_degree, after_den_degree = certificate.after_den_degree,
                            cross_product_terms = ?certificate.cross_product_terms,
                            file.opaque_leaves = ?certificate.opaque_leaves,
                            "Offline exact scalar identity and removed-factor certificate"
                        );
                    }
                    Err(error) => {
                        if !matches!(
                            error.downcast_ref::<ScalarPoleFailure>(),
                            Some(ScalarPoleFailure::Unsupported(_) | ScalarPoleFailure::Budget(_))
                        ) {
                            return Err(error);
                        }
                        root_complete = false;
                        unsupported_scalar_pairs += 1;
                        debug_tags!(#uv, #pole_audit;
                            stage = "pole_audit_unproven", root_id, power = %power,
                            file.numerator_monomial = %monomial.to_plain_string(),
                            file.before = %before.to_plain_string(), file.after = %after.to_plain_string(),
                            reason = %error, "Scalar pair remains unproven"
                        );
                    }
                }
            }
        }
        fully_certified_roots += usize::from(root_complete);
    }
    debug_tags!(#uv, #pole_audit;
        stage = "pole_audit_summary", graph = %recording.graph, orientations = 102,
        accepted_roots = recording.captures.len() + recording.dropped,
        captured_roots = recording.captures.len(), dropped_roots = recording.dropped,
        capture_bytes = recording.bytes, fully_certified_roots, certified_scalar_pairs,
        unsupported_roots = recording.captures.len() - fully_certified_roots + recording.dropped,
        unsupported_scalar_pairs, removed_factor_pairs, zero_pairs, unchanged_pairs,
        "Offline coverage of accepted U changes; no claim of a full-H symbolic soft-limit proof"
    );
    assert!(
        certified_scalar_pairs > 0,
        "nonvacuous exact scalar proof coverage is required"
    );
    // No assertion requires finding a removed factor: zero is a valid outcome.
    Ok(())
}
