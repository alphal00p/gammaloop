//! Bare propagators with their original open tensor boundaries.

use eyre::{Result, WrapErr, ensure, eyre};
use linnet::half_edge::subgraph::{SuBitGraph, SubSetLike, SubSetOps};
use symbolica::atom::{Atom, AtomCore};
use vakint::Vakint;

use crate::{
    graph::{FeynmanGraph, Graph},
    model::Model,
    utils::GS,
    uv::{UVgenerationSettings, UltravioletGraph, uv_graph::UVE},
};

use super::{GraphRegion, expression::RStarIntegrand};

/// Integrate a massless, connected 1PI two-point graph without UV subtractions.
///
/// Local vertex and propagator numerators retain their open Lorentz, color and
/// spinor slots. The result uses the original external edge's `Q` momentum,
/// with the two external momenta identified by momentum conservation. Global
/// graph factors and projectors are excluded, as for [`super::RStarPlan`].
///
/// This invokes the physical p-integral backend directly, without building a
/// forest or computing a counterterm. Include `forcer` in the evaluation
/// methods for up to four loops, or RustRed for up to five loops. Propagators must have
/// unit powers. Physical infrared poles, if present, remain in the result.
///
/// `number_of_terms` has Vakint's convention: an L-loop integral is retained
/// through epsilon^(number_of_terms-L-1). Thus L+1 includes the finite term.
/// The normalization, additional per-loop factor and regulator symbols are
/// exactly those of `uv_settings`.
pub fn integrate_propagator(
    graph: &Graph,
    model: &Model,
    vakint: &Vakint,
    uv_settings: &UVgenerationSettings,
    number_of_terms: i64,
) -> Result<Atom> {
    integrate_propagator_with(graph, model, vakint, uv_settings, number_of_terms, None)
}

/// [`integrate_propagator`] with an optional `closure` multiplied into the
/// integrand before any contraction ([`RStarIntegrand::new_with_closure`]).
fn integrate_propagator_with(
    graph: &Graph,
    model: &Model,
    vakint: &Vakint,
    uv_settings: &UVgenerationSettings,
    number_of_terms: i64,
    closure: Option<&Atom>,
) -> Result<Atom> {
    ensure!(
        number_of_terms > 0,
        "a propagator expansion needs at least one term"
    );
    let boundary = graph.external_filter::<SuBitGraph>();
    RStarIntegrand::physical_external_routing(graph, &boundary)?;
    let support = graph.no_dummy().subtract(&boundary);
    let loops = graph.n_loops(&support);
    let max_loops = if matches!(uv_settings.vakint.true_settings().evaluation_order.0.first(), Some(vakint::EvaluationMethod::RustRed(_))) { 5 } else { 4 };
    ensure!((1..=max_loops).contains(&loops),
        "direct propagator integration supports one through {max_loops} loops with this backend");
    ensure!(
        graph.is_connected(&support) && graph.bridges_of(&support).is_empty(),
        "direct propagator integration requires a connected 1PI graph"
    );
    for (_, _, edge) in graph.iter_edges_of(&support) {
        ensure!(
            edge.data.mass_atom(model).is_zero(),
            "direct propagator integration requires massless propagators"
        );
        ensure!(
            edge.data.extra_data.vakint_edge_power.unwrap_or(1) == 1,
            "direct propagator integration currently requires unit propagator powers"
        );
    }
    let marked = graph
        .iter_edges_of(&support)
        .next()
        .map(|(_, edge, _)| edge)
        .ok_or_else(|| eyre!("a propagator graph contains no internal edges"))?;
    let region = GraphRegion {
        remaining: support.clone(),
        support,
        insertions: Vec::new(),
        boundary,
        loops,
    };
    // The marked edge is only a constructor argument here: its auxiliary
    // mass is zero and no Taylor operation is performed.
    let integrand =
        RStarIntegrand::new_with_closure(graph, model, &region, &[], marked, Atom::Zero, closure)?;
    let mut settings = uv_settings.vakint.true_settings();
    settings.number_of_terms_in_epsilon_expansion = number_of_terms;
    integrand
        .integrate_physical(vakint, uv_settings, &settings, 0)
        .wrap_err("while integrating a bare open-tensor propagator")
}

/// The complete form-factor basis of a massless two-point insertion, read
/// off the open slots of its numerator ([`insertion_shape`]).
#[derive(Clone, Debug)]
pub(crate) enum InsertionShape {
    /// One open Dirac line (two bispinor slots, no Lorentz slot) with a
    /// numerator odd in the gamma matrices: `X = A slash(p)` exactly in `d`
    /// dimensions (one vector, no unit-matrix term for an odd line).
    Quark {
        row: Atom,
        column: Atom,
        trace_of_one: Atom,
    },
    /// Two open Lorentz slots of one representation, no spinor slot:
    /// `X^{mu nu} = a g^{mu nu} + b p^mu p^nu`, the complete parity-even
    /// rank-two basis with one vector (no transversality assumed).
    Gluon { slots: [Atom; 2] },
    /// No open Lorentz or spinor slot (a ghost, or a scalar): already closed.
    Scalar,
}

impl InsertionShape {
    pub(crate) fn kind(&self) -> &'static str {
        match self {
            Self::Quark { .. } => "quark",
            Self::Gluon { .. } => "gluon",
            Self::Scalar => "scalar",
        }
    }
}

/// The [`InsertionShape`] of a two-point graph's numerator over its internal
/// edges, or the reason it has none. Slots of any other representation
/// (colour) are left open and untouched by the closure. A parity-odd symbol
/// (`gamma5`, chiral projector, Levi-Civita) declines: the bases above span
/// only the parity-even tensors.
pub(crate) fn insertion_shape(
    graph: &Graph,
) -> Result<std::result::Result<InsertionShape, String>> {
    use idenso::{CookMode, CookSettings, tensor::SymbolicTensor};
    use spenso::structure::{
        abstract_index::AbstractIndex,
        partial::PartialStructureExt,
        representation::{LibraryRep, Minkowski},
        slot::IsAbstractSlot,
    };
    let boundary = graph.external_filter::<SuBitGraph>();
    let support = graph.no_dummy().subtract(&boundary);
    let cooking = CookSettings::indices()
        .with_mode(CookMode::ReversibleEncoding)
        .with_representation_payloads(true, true);
    let numerator = graph
        .numerator(&support, &graph.empty_subgraph())
        .to_d_dim(GS.dim)
        .get_single_atom()?;
    let parity_odd = [
        idenso::dirac::AGS.gamma5,
        idenso::dirac::AGS.projp,
        idenso::dirac::AGS.projm,
        *idenso::epsilon::EPSILON_SYMBOL,
    ];
    if parity_odd
        .iter()
        .any(|symbol| numerator.contains_symbol(*symbol))
    {
        return Ok(Err("parity_odd".into()));
    }
    let tensor = SymbolicTensor::infer(cooking.try_cook(numerator.as_view())?)?;
    let bispinor: LibraryRep = idenso::representations::Bispinor {}.into();
    let minkowski: LibraryRep = Minkowski {}.into();
    let mut spinors = Vec::new();
    let mut lorentz = Vec::new();
    for slot in tensor.structure().logical_slots() {
        let spenso::structure::partial::PartialIndex::Explicit(index) = slot.aind else {
            return Ok(Err("implicit_slot".into()));
        };
        let slot = slot.rep().slot::<AbstractIndex, _>(index);
        if slot.rep().rep == bispinor {
            spinors.push(slot);
        } else if slot.rep().rep == minkowski {
            lorentz.push(slot);
        }
    }
    match (spinors.len(), lorentz.len()) {
        (0, 0) => Ok(Ok(InsertionShape::Scalar)),
        (0, 2) => {
            if lorentz[0].rep() != lorentz[1].rep() {
                return Ok(Err("lorentz_reps".into()));
            }
            Ok(Ok(InsertionShape::Gluon {
                slots: [
                    cooking.uncook(lorentz[0].to_atom().as_view()),
                    cooking.uncook(lorentz[1].to_atom().as_view()),
                ],
            }))
        }
        (2, 0) => {
            if spinors[0].rep() != spinors[1].rep() {
                return Ok(Err("spinor_reps".into()));
            }
            let atoms = spinors
                .iter()
                .map(|slot| cooking.uncook(slot.to_atom().as_view()))
                .collect::<Vec<_>>();
            // As `dirac_boundary_with`: each open spinor slot occurs only in
            // gamma matrices, the one only as their row, the other only as
            // their column.
            let mut rows = [0_usize; 2];
            let mut columns = [0_usize; 2];
            let mut total = [0_usize; 2];
            numerator.visitor(&mut |part| {
                for (k, atom) in atoms.iter().enumerate() {
                    if part == atom.as_view() {
                        total[k] += 1;
                    }
                }
                if let symbolica::atom::AtomView::Fun(fun) = part
                    && fun.get_symbol() == idenso::dirac::AGS.gamma
                    && fun.get_nargs() == 3
                {
                    for (k, atom) in atoms.iter().enumerate() {
                        if fun.get(0) == atom.as_view() {
                            rows[k] += 1;
                        }
                        if fun.get(1) == atom.as_view() {
                            columns[k] += 1;
                        }
                    }
                }
                true
            });
            if (0..2).any(|k| total[k] != rows[k] + columns[k]) {
                return Ok(Err("spinor_outside_gamma".into()));
            }
            let (row, column) = if rows[0] > 0 && columns[0] == 0 && columns[1] > 0 && rows[1] == 0
            {
                (0, 1)
            } else if rows[1] > 0 && columns[1] == 0 && columns[0] > 0 && rows[0] == 0 {
                (1, 0)
            } else {
                return Ok(Err("row_column".into()));
            };
            match super::expression::gamma_parity(numerator.as_view()) {
                Some(1) => {}
                Some(_) => return Ok(Err("gamma_parity:even".into())),
                None => return Ok(Err("gamma_parity:mixed".into())),
            }
            Ok(Ok(InsertionShape::Quark {
                row: atoms[row].clone(),
                column: atoms[column].clone(),
                trace_of_one: cooking.uncook(spinors[0].rep().dim.to_symbolic().as_view()),
            }))
        }
        (spinors, lorentz) => Ok(Err(format!("shape:bis{spinors}_mink{lorentz}"))),
    }
}

/// [`integrate_propagator`] through the closed form factors of the
/// insertion's [`InsertionShape`]: the open tensor `X` is projected onto its
/// complete basis `T_k` before any contraction, the closed scalars are
/// integrated (all channels in one integral, separated by scalar markers),
/// and `X = sum_k F_k T_k` is rebuilt exactly. With `p = Q_E` (the first
/// boundary momentum, [`RStarIntegrand::physical_external_routing`]) and
/// `d = 4 - 2 epsilon`:
/// * quark: `S = Tr[slash(p) X] = Tr[1] A p^2`, `X = S/(Tr[1] p^2) slash(p)`;
/// * gluon: `S1 = g_{mu nu} X^{mu nu} = d a + b p^2`,
///   `S2 = p_mu p_nu X^{mu nu} = (a + b p^2) p^2`, so
///   `a = (S1 - S2/p^2)/(d-1)`, `b = (d S2/p^2 - S1)/((d-1) p^2)`.
///
/// The projectors are regular in epsilon and free of loop momenta, so they
/// commute with the integration and its epsilon expansion: the result equals
/// [`integrate_propagator`]'s through the same depth. Returns the decline
/// reason for a scalar insertion (nothing to close) or an unrecognized shape.
pub(crate) fn integrate_propagator_closed(
    graph: &Graph,
    model: &Model,
    vakint: &Vakint,
    uv_settings: &UVgenerationSettings,
    number_of_terms: i64,
) -> Result<std::result::Result<(Atom, InsertionShape), String>> {
    use spenso::{
        network::{library::symbolic::ETS, tags::SPENSO_TAG},
        structure::representation::{Minkowski, RepName},
    };
    use symbolica::{function, symbol};
    let shape = match insertion_shape(graph)? {
        Ok(InsertionShape::Scalar) => return Ok(Err("scalar".into())),
        Ok(shape) => shape,
        Err(reason) => return Ok(Err(reason)),
    };
    let (external, _) =
        RStarIntegrand::physical_external_routing(graph, &graph.external_filter::<SuBitGraph>())?;
    let external = usize::from(external);
    let minkowski = Minkowski {}.new_rep(GS.dim).to_symbolic([]);
    let p = function!(GS.emr_mom, external, &minkowski);
    let p2 = function!(SPENSO_TAG.dot, &p, &p);
    let inverse_p2 = Atom::one() / &p2;
    let d = Atom::num(4) - Atom::num(2) * Atom::var(GS.dim_epsilon);
    let regularize = crate::uv::approx::integrated::Integrated::dimensionally_regularized;
    let rebuilt = match &shape {
        InsertionShape::Quark {
            row,
            column,
            trace_of_one,
        } => {
            let closure = function!(idenso::dirac::AGS.gamma, column, row, &p);
            let trace = integrate_propagator_with(
                graph,
                model,
                vakint,
                uv_settings,
                number_of_terms,
                Some(&closure),
            )?;
            let coefficient = regularize(&(trace * &inverse_p2).expand()) / trace_of_one;
            coefficient * function!(idenso::dirac::AGS.gamma, row, column, &p)
        }
        InsertionShape::Gluon { slots: [mu, nu] } => {
            let channels = [
                symbol!("gammalooprs::insertion_closure_metric"; Scalar),
                symbol!("gammalooprs::insertion_closure_longitudinal"; Scalar),
            ];
            let metric = ETS.metric(mu, nu);
            let momenta = function!(GS.emr_mom, external, mu) * function!(GS.emr_mom, external, nu);
            let closure = Atom::var(channels[0]) * &metric + Atom::var(channels[1]) * &momenta;
            let integral = integrate_propagator_with(
                graph,
                model,
                vakint,
                uv_settings,
                number_of_terms,
                Some(&closure),
            )?;
            let channel = |k: usize| {
                regularize(
                    &integral
                        .replace(Atom::var(channels[k]))
                        .with(Atom::one())
                        .replace(Atom::var(channels[1 - k]))
                        .with(Atom::Zero),
                )
            };
            let s1 = channel(0);
            let s2_over_p2 = (channel(1) * &inverse_p2).expand();
            let a = (&s1 - &s2_over_p2) / (&d - 1);
            let b = ((&d * &s2_over_p2 - &s1) * &inverse_p2).expand() / (&d - 1);
            a * metric + b * momenta
        }
        InsertionShape::Scalar => unreachable!("declined above"),
    };
    Ok(Ok((rebuilt, shape)))
}

#[cfg(test)]
mod tests {
    use eyre::Result;
    use idenso::representations::Bispinor;
    use linnet::half_edge::involution::EdgeIndex;
    use spenso::{
        network::{library::symbolic::ETS, tags::SPENSO_TAG},
        structure::{
            representation::{Minkowski, RepName},
            slot::{DummyAind, ParseableAind},
        },
    };
    use symbolica::{
        atom::{Atom, Symbol},
        function,
    };
    use vakint::Vakint;

    use crate::{
        finalized_runtime_dot,
        graph::{Autogen, Graph, parse::IntoFinalizedRuntimeGraph},
        initialisation::test_initialise,
        model::Model,
        numerator::aind::Aind,
        utils::{GS, load_generic_model},
        uv::UVgenerationSettings,
    };

    use super::integrate_propagator;
    use crate::uv::rstar::{RStarPlan, tests::assert_same_counterterm};

    fn settings() -> UVgenerationSettings {
        let mut settings = crate::uv::rstar::tests::settings();
        settings.vakint.normalization = "GScheme".into();
        settings
    }

    fn make_massless(graph: &mut Graph) {
        let edges = graph
            .iter_edges()
            .map(|(_, edge, _)| edge)
            .collect::<Vec<_>>();
        for edge in edges {
            let data = &mut graph.underlying[edge];
            data.particle = data.particle.clone().override_mass(Some(Atom::Zero));
        }
    }

    fn bubble(model: &Model) -> Result<Graph> {
        let mut graph: Graph = finalized_runtime_dot!(
            digraph G {
                graph [projector=7];
                node [num=1];
                edge [particle="scalar_1" num=1];
                e [style=invis];
                e -> A:0 [id=2 sink="{ufo_order:0}"];
                B:1 -> e [id=3 source="{ufo_order:0}"];
                A -> B [id=0 source="{ufo_order:1}" sink="{ufo_order:1}"];
                A -> B [id=1 lmb_id=0 source="{ufo_order:2}" sink="{ufo_order:2}"];
            }, model
        )?;
        make_massless(&mut graph);
        graph.overall_factor = Atom::num(11);
        graph.global_prefactor.num = Atom::num(13);
        Ok(graph)
    }

    fn external(graph: &Graph) -> EdgeIndex {
        graph
            .iter_edges()
            .find_map(|(pair, edge, _)| (!pair.is_paired()).then_some(edge))
            .unwrap()
    }

    fn external_square(graph: &Graph) -> Atom {
        let momentum = function!(
            GS.emr_mom,
            usize::from(external(graph)),
            Minkowski {}.new_rep(GS.dim).to_symbolic([])
        );
        function!(SPENSO_TAG.dot, &momentum, &momentum)
    }

    fn scalar_bubble_through_finite(graph: &Graph) -> Atom {
        Atom::one() / Atom::var(GS.dim_epsilon)
            - function!(Symbol::LOG, -external_square(graph) / Atom::var(GS.mu_r_sq))
    }

    #[test]
    fn bare_scalar_bubble_keeps_finite_term_and_excludes_global_factors() -> Result<()> {
        test_initialise()?;
        let model = load_generic_model("scalars");
        let graph = bubble(&model)?;
        let actual = integrate_propagator(&graph, &model, &Vakint::new()?, &settings(), 2)?;
        assert_same_counterterm(&actual, &scalar_bubble_through_finite(&graph))
    }

    #[test]
    fn bare_bubble_retains_open_tensor_and_evanescent_finite_terms() -> Result<()> {
        test_initialise()?;
        let model = load_generic_model("scalars");
        let mut graph = bubble(&model)?;
        let minkowski = Minkowski {}.new_rep(GS.dim);
        let mu = minkowski.to_symbolic([Aind::Hedge(40, 0).to_atom()]);
        let nu = minkowski.to_symbolic([Aind::Hedge(41, 0).to_atom()]);
        let loop_edge = *graph.loop_momentum_basis.loop_edges.first().unwrap();
        let numerator = function!(GS.emr_mom, usize::from(loop_edge), &mu)
            * function!(GS.emr_mom, usize::from(loop_edge), &nu);
        graph.underlying[loop_edge].num = Autogen::explicit(numerator);
        let actual = integrate_propagator(&graph, &model, &Vakint::new()?, &settings(), 2)?;
        let scalar = scalar_bubble_through_finite(&graph);
        let q_mu = function!(GS.emr_mom, usize::from(external(&graph)), &mu);
        let q_nu = function!(GS.emr_mom, usize::from(external(&graph)), &nu);
        // Tensor bubble = [d q_mu q_nu - q² g_mu_nu] B / [4(d-1)].
        let expected = q_mu * q_nu * (&scalar / Atom::num(3) + Atom::one() / Atom::num(18))
            - external_square(&graph)
                * ETS.metric(&mu, &nu)
                * (&scalar / Atom::num(12) + Atom::one() / Atom::num(18));
        assert_same_counterterm(&actual, &expected)
    }

    #[test]
    fn bare_bubble_projects_loop_slash_inside_open_clifford_word() -> Result<()> {
        test_initialise()?;
        let model = load_generic_model("scalars");
        let mut graph = bubble(&model)?;
        let minkowski = Minkowski {}.new_rep(GS.dim);
        let spinor = Bispinor {}.new_rep(4);
        let left = spinor.to_symbolic([Aind::Hedge(40, 0).to_atom()]);
        let right = spinor.to_symbolic([Aind::Hedge(41, 0).to_atom()]);
        let mu = minkowski.to_symbolic([Aind::new_dummy().to_atom()]);
        let loop_edge = *graph.loop_momentum_basis.loop_edges.first().unwrap();
        let momentum = function!(
            GS.emr_mom,
            usize::from(loop_edge),
            minkowski.to_symbolic([])
        );
        let gamma = |argument: &Atom| {
            function!(
                idenso::dirac::AGS.gamma,
                SPENSO_TAG.chain_in,
                SPENSO_TAG.chain_out,
                argument
            )
        };
        graph.underlying[loop_edge].num = Autogen::explicit(spenso::chain!(
            &left,
            &right,
            gamma(&mu),
            gamma(&momentum),
            gamma(&mu)
        ));
        let actual = integrate_propagator(&graph, &model, &Vakint::new()?, &settings(), 2)?;
        let external = function!(
            GS.emr_mom,
            usize::from(external(&graph)),
            minkowski.to_symbolic([])
        );
        // gamma_mu slash(k) gamma^mu = (2-d) slash(k), and the vector
        // bubble is q_mu B/2. The evanescent factor supplies a finite +1.
        let expected = function!(idenso::dirac::AGS.gamma, &left, &right, &external)
            * (Atom::one() - scalar_bubble_through_finite(&graph));
        assert_same_counterterm(&actual, &expected)
    }

    /// GAMMALOOP_INSERTION_CLOSURE: a non-transverse rank-two bubble
    /// (`k_mu k_nu + k^2 g_mu_nu`, both form factors nonzero) rebuilt from its
    /// two closed channels equals the open tensor, with its finite term.
    #[test]
    fn closed_gluon_shape_rebuilds_the_open_rank_two_bubble() -> Result<()> {
        test_initialise()?;
        let model = load_generic_model("scalars");
        let mut graph = bubble(&model)?;
        let minkowski = Minkowski {}.new_rep(GS.dim);
        let mu = minkowski.to_symbolic([Aind::Hedge(40, 0).to_atom()]);
        let nu = minkowski.to_symbolic([Aind::Hedge(41, 0).to_atom()]);
        let loop_edge = *graph.loop_momentum_basis.loop_edges.first().unwrap();
        let k = function!(
            GS.emr_mom,
            usize::from(loop_edge),
            minkowski.to_symbolic([])
        );
        let numerator = function!(GS.emr_mom, usize::from(loop_edge), &mu)
            * function!(GS.emr_mom, usize::from(loop_edge), &nu)
            + Atom::num(3) * function!(SPENSO_TAG.dot, &k, &k) * ETS.metric(&mu, &nu);
        graph.underlying[loop_edge].num = Autogen::explicit(numerator);
        let vakint = Vakint::new()?;
        let open = integrate_propagator(&graph, &model, &vakint, &settings(), 3)?;
        let (closed, shape) =
            super::integrate_propagator_closed(&graph, &model, &vakint, &settings(), 3)?
                .map_err(|reason| eyre::eyre!("declined: {reason}"))?;
        assert_eq!(shape.kind(), "gluon");
        assert!(!open.is_zero());
        // The closed form factors carry 1/(d-1) unexpanded: they agree with
        // the open tensor through the computed depth, eps^(3-1-1).
        use symbolica::atom::AtomCore;
        let depth = |value: &Atom| -> Result<Atom> {
            Ok(value.series(GS.dim_epsilon, 0, 1)?.to_atom().expand())
        };
        assert_same_counterterm(&depth(&closed)?, &depth(&open)?)
    }

    /// GAMMALOOP_INSERTION_CLOSURE: an open Dirac line `gamma_mu slash(k)
    /// gamma^mu` rebuilt from the closed trace `Tr[slash(p) X]` equals the
    /// open matrix; a scalar numerator declines (nothing to close).
    #[test]
    fn closed_quark_shape_rebuilds_the_open_dirac_line() -> Result<()> {
        test_initialise()?;
        let model = load_generic_model("scalars");
        let mut graph = bubble(&model)?;
        let minkowski = Minkowski {}.new_rep(GS.dim);
        let spinor = Bispinor {}.new_rep(4);
        let left = spinor.to_symbolic([Aind::Hedge(40, 0).to_atom()]);
        let right = spinor.to_symbolic([Aind::Hedge(41, 0).to_atom()]);
        let a = spinor.to_symbolic([Aind::new_dummy().to_atom()]);
        let b = spinor.to_symbolic([Aind::new_dummy().to_atom()]);
        let mu = minkowski.to_symbolic([Aind::new_dummy().to_atom()]);
        let loop_edge = *graph.loop_momentum_basis.loop_edges.first().unwrap();
        let momentum = function!(
            GS.emr_mom,
            usize::from(loop_edge),
            minkowski.to_symbolic([])
        );
        let gamma = |i: &Atom, j: &Atom, argument: &Atom| {
            function!(idenso::dirac::AGS.gamma, i, j, argument)
        };
        let vakint = Vakint::new()?;
        assert!(
            super::integrate_propagator_closed(&graph, &model, &vakint, &settings(), 2)?.is_err(),
            "a scalar insertion has nothing to close"
        );
        graph.underlying[loop_edge].num = Autogen::explicit(
            gamma(&left, &a, &mu) * gamma(&a, &b, &momentum) * gamma(&b, &right, &mu),
        );
        let open = integrate_propagator(&graph, &model, &vakint, &settings(), 3)?;
        let (closed, shape) =
            super::integrate_propagator_closed(&graph, &model, &vakint, &settings(), 3)?
                .map_err(|reason| eyre::eyre!("declined: {reason}"))?;
        assert_eq!(shape.kind(), "quark");
        assert!(!open.is_zero());
        assert_same_counterterm(&closed, &open)
    }

    #[test]
    fn bare_nested_bubble_poles_match_rstar_reconstruction() -> Result<()> {
        test_initialise()?;
        let model = load_generic_model("scalars");
        let mut graph: Graph = finalized_runtime_dot!(
            digraph G {
                graph [projector=1];
                node [num=1];
                edge [particle="scalar_1" num=1];
                e [style=invis];
                e -> A:0 [id=4 sink="{ufo_order:0}"];
                C:1 -> e [id=5 source="{ufo_order:0}"];
                A -> B [id=0 source="{ufo_order:1}" sink="{ufo_order:0}"];
                A -> B [id=1 lmb_id=0 source="{ufo_order:2}" sink="{ufo_order:1}"];
                A -> C [id=2 source="{ufo_order:3}" sink="{ufo_order:1}"];
                B -> C [id=3 lmb_id=1 source="{ufo_order:2}" sink="{ufo_order:2}"];
            }, &model
        )?;
        make_massless(&mut graph);
        let settings = settings();
        let vakint = Vakint::new()?;
        let actual = integrate_propagator(&graph, &model, &vakint, &settings, 2)?;
        let plan = RStarPlan::new(&graph, &model, EdgeIndex::from(0), &settings)?;
        let expected = plan.pole_part(&graph, &model, &vakint, &settings)?;
        assert_same_counterterm(&actual, &expected)
    }

    #[test]
    fn bare_propagator_rejects_ignored_propagator_powers() -> Result<()> {
        test_initialise()?;
        let model = load_generic_model("scalars");
        let mut graph = bubble(&model)?;
        graph.underlying[EdgeIndex::from(0)]
            .extra_data
            .vakint_edge_power = Some(2);
        let error =
            integrate_propagator(&graph, &model, &Vakint::new()?, &settings(), 2).unwrap_err();
        assert!(error.to_string().contains("unit propagator powers"));
        Ok(())
    }
}
