use super::{OrientationIntegrandBranch, OrientationIntegrands};
use crate::{
    cff::{
        CutCFFIndex,
        expression::{
            GammaLoopOrientationExpression, OrientationData, OrientationExpression, OrientationID,
        },
        orientations::GraphOrientation,
        surface::LinearEnergyExpr,
    },
    dot,
    graph::{FeynmanGraph, FourDDenominator, Graph, LMBext, cuts::CutSet, parse::IntoGraph},
    initialisation::test_initialise,
    settings::global::OrientationPattern,
    utils::GS,
    uv::{
        ApproximationType, Spinney, UltravioletGraph,
        approx::{
            OrientationProjection,
            direct_3d::{Direct3dCts, DirectResidueBranches},
            local_3d::Localizer,
        },
        hedge_poset::OwnedForestNode,
    },
};
use color_eyre::Result;
use eyre::eyre;
use linnet::half_edge::{
    involution::{EdgeIndex, EdgeVec, Orientation},
    subgraph::{InternalSubGraph, SubSetOps},
};
use std::{collections::BTreeSet, sync::OnceLock};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, FunctionBuilder},
    function,
};
use three_dimensional_reps::CffGenerationContext;
use typed_index_collections::TiVec;

static TWO_EDGE_GRAPH: OnceLock<Graph> = OnceLock::new();

fn orientation(value: i8) -> Orientation {
    match value {
        1 => Orientation::Default,
        -1 => Orientation::Reversed,
        0 => Orientation::Undirected,
        _ => panic!("invalid orientation encoding"),
    }
}

fn edgevec(values: impl IntoIterator<Item = i8>) -> EdgeVec<Orientation> {
    EdgeVec::from_iter(values.into_iter().map(orientation))
}

fn edges(values: impl IntoIterator<Item = usize>) -> Vec<EdgeIndex> {
    values.into_iter().map(EdgeIndex).collect()
}

fn energy_map(
    orientation: EdgeVec<Orientation>,
    edge_energy_map: Vec<LinearEnergyExpr>,
) -> OrientationExpression {
    OrientationExpression {
        data: OrientationData::new(orientation),
        loop_energy_map: Vec::new(),
        edge_energy_map,
        variants: Vec::new(),
    }
}

fn two_edge_graph() -> Result<Graph> {
    Ok(TWO_EDGE_GRAPH
        .get_or_init(|| {
            test_initialise().expect("test model initialization succeeds");
            dot!(
                digraph G {
                    edge [particle="scalar_1"];
                    node [num=1];
                    a -> b [id=0];
                    a -> b [id=1];
                },
                "scalars"
            )
            .expect("two-edge projection graph parses")
        })
        .clone())
}

#[test]
fn production_emr_map_cancels_one_powered_denominator() -> Result<()> {
    let mut graph = two_edge_graph()?;
    let edge = EdgeIndex(0);
    for analysis_numerator in [None, Some(GS.emr_mom(edge, GS.cind(0)))] {
        let ordinary_options = graph.denominator_only_cff_3d_expression_options();
        let ordinary_canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
        let ordinary = graph.generate_3d_expression_for_integrand(
            &[],
            &ordinary_canonization,
            &ordinary_options,
            analysis_numerator.as_ref(),
        )?;
        let ordinary_cff = graph.cff_from_production_expression(
            &ordinary,
            &CutSet::empty(graph.n_hedges()),
            &OrientationPattern::default(),
        )?;
        assert_eq!(
            ordinary_cff.production_prefactor_factor(),
            ordinary.core_global_prefactor_sign.factor(),
            "the typed bridge must be a no-op for scalar and linear-key ordinary production CFFs",
        );
    }
    let on_shell_energy_squared = (1..=3).fold(
        graph.underlying[edge].particle.mass_atom().pow(2),
        |norm_squared, spatial_index| {
            norm_squared + GS.emr_mom(edge, GS.cind(spatial_index)).pow(2)
        },
    );
    let numerator = GS.emr_mom(edge, GS.cind(0)).pow(2) - &on_shell_energy_squared;
    let options = graph.denominator_only_cff_3d_expression_options();
    let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
    let production = graph.generate_3d_expression_for_integrand(
        &[],
        &canonization,
        &options,
        Some(&numerator),
    )?;
    let pattern = OrientationPattern::default();
    let cutset = CutSet::empty(graph.n_hedges());
    let localizer = Localizer::new(
        &cutset,
        OrientationProjection::exact_expression(&production, &options, &pattern, false),
    );
    let contract = graph.empty_subgraph();
    let powered = localizer.projected_cff(
        &mut graph,
        &contract,
        &numerator,
        CffGenerationContext::Standalone,
    )?;
    let mut explicit_sum = Atom::Zero;
    for (selector_id, source_map, integrands) in powered.iter_orientations() {
        let mapped = localizer.map_numerator(&graph, selector_id, source_map, &numerator)?;
        explicit_sum += integrands
            .iter()
            .fold(Atom::Zero, |sum, (_, term)| sum + term * &mapped);
    }
    // The two parallel physical edges carry opposite routing signs. Their
    // acyclic production sectors are therefore (+,-) and (-,+); (+,+) is
    // not a production orientation. The cancellation contacts have zero
    // numerator sampling maps, but their opposite loop lifts still host
    // one matching lower residue in each sector.
    let mut selected_production = edgevec([1, -1]).select(&explicit_sum);
    let regenerated_localizer = Localizer::new(
        &cutset,
        OrientationProjection::exact(
            &production.expression.orientations,
            &options,
            &pattern,
            false,
        ),
    );
    let regenerated = regenerated_localizer.projected_cff(
        &mut graph,
        &contract,
        &numerator,
        CffGenerationContext::Standalone,
    )?;
    let mut regenerated_sum = Atom::Zero;
    for (selector_id, source_map, integrands) in regenerated.iter_orientations() {
        let mapped =
            regenerated_localizer.map_numerator(&graph, selector_id, source_map, &numerator)?;
        regenerated_sum += integrands
            .iter()
            .fold(Atom::Zero, |sum, (_, term)| sum + term * &mapped);
    }
    let mut selected_regenerated = edgevec([1, -1]).select(&regenerated_sum);
    let momentum = FunctionBuilder::new(GS.emr_mom)
        .add_arg(usize::from(edge))
        .finish();
    let mass_squared = graph.underlying[edge].particle.mass_atom().pow(2);
    let denominators = [
        FourDDenominator {
            source_edge: edge,
            momentum: momentum.clone(),
            mass_squared: mass_squared.clone(),
            full_expr: Atom::var(symbolica::symbol!("local_3d_test::first")),
        },
        FourDDenominator {
            source_edge: EdgeIndex(1),
            momentum: -momentum,
            mass_squared: graph.underlying[EdgeIndex(1)].particle.mass_atom().pow(2),
            full_expr: Atom::var(symbolica::symbol!("local_3d_test::second")),
        },
    ];
    let (exact_powered, _) =
        graph.cff_from_4d_denominators(&denominators, &cutset, &options, &numerator)?;
    // A projected post-4D source owns its complete orientation sum. Once
    // the numerator cancels an occurrence, that provenance slot can be
    // undirected and must not be used to filter the surviving contact.
    let mut exact_powered_sum = exact_powered
        .terms
        .values()
        .flat_map(|term| {
            term.orientations.iter().map(|orientation| {
                Ok::<_, eyre::Error>(
                    &orientation.expression
                        * term.map_exact_source_numerator(&orientation.orientation)?,
                )
            })
        })
        .collect::<Result<Vec<_>>>()?
        .into_iter()
        .sum::<Atom>()
        * Atom::num(exact_powered.production_prefactor_factor());
    let (exact_lower, _) =
        graph.cff_from_4d_denominators(&denominators[1..], &cutset, &options, &Atom::one())?;
    // The exact one-denominator source enumerates both contour directions.
    // Compare the remaining-edge direction selected by the (+,-)
    // production sector. Denominator-routing signs and physical selector
    // signs are separate conventions.
    let remaining_lower_occurrence = EdgeIndex(graph.underlying.n_edges());
    let mut exact_lower_sum = exact_lower
        .terms
        .values()
        .flat_map(|term| &term.orientations)
        .find(|orientation| {
            orientation.orientation.data.orientation[remaining_lower_occurrence]
                == Orientation::Default
        })
        .expect("the lower source has the selected remaining-edge direction")
        .expression
        .clone()
        * Atom::num(exact_lower.production_prefactor_factor());
    let mass = graph.underlying[edge].particle.mass_atom();
    for expression in [
        &mut selected_production,
        &mut selected_regenerated,
        &mut exact_powered_sum,
        &mut exact_lower_sum,
    ] {
        *expression = expression.replace(mass.clone()).with(Atom::one());
        for spatial_index in 1..=3 {
            *expression = expression
                .replace(GS.emr_mom(edge, GS.cind(spatial_index)))
                .with(Atom::Zero);
        }
        *expression = expression
            .replace(GS.ose(EdgeIndex(0)))
            .with(Atom::one())
            .replace(GS.ose(EdgeIndex(1)))
            .with(Atom::one());
    }
    let exact_difference = (&exact_powered_sum - &exact_lower_sum).together();
    assert!(
        exact_difference.is_zero(),
        "summed exact-source powered and selected lower channels differ: powered={exact_powered_sum}, lower={exact_lower_sum}, difference={exact_difference}"
    );
    let production_difference = (&selected_production - &exact_lower_sum).together();
    assert!(
        production_difference.is_zero(),
        "production EMR mapping does not cancel the powered denominator: production={selected_production}, lower={exact_lower_sum}, difference={production_difference}"
    );
    let regenerated_difference = (&selected_regenerated - &exact_lower_sum).together();
    assert!(
        regenerated_difference.is_zero(),
        "regenerated graph CFF does not cancel the powered denominator: regenerated={selected_regenerated}, lower={exact_lower_sum}, difference={regenerated_difference}"
    );
    Ok(())
}

#[test]
fn direct_root_preserves_powered_production_entries() -> Result<()> {
    let mut graph = two_edge_graph()?;
    let edge = EdgeIndex(0);
    let on_shell_energy_squared = (1..=3).fold(
        graph.underlying[edge].particle.mass_atom().pow(2),
        |norm_squared, spatial_index| {
            norm_squared + GS.emr_mom(edge, GS.cind(spatial_index)).pow(2)
        },
    );
    let numerator = GS.emr_mom(edge, GS.cind(0)).pow(2) - on_shell_energy_squared;
    let options = graph.denominator_only_cff_3d_expression_options();
    let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
    let production = graph.generate_3d_expression_for_integrand(
        &[],
        &canonization,
        &options,
        Some(&numerator),
    )?;
    let pattern = OrientationPattern::default();
    let cutset = CutSet::empty(graph.n_hedges());
    let raw_root = graph.cff_from_production_expression(&production, &cutset, &pattern)?;
    let production_prefactor = Atom::num(raw_root.production_prefactor_factor());
    let fourddenoms = GS.wrap_tree_denoms(
        graph.denominator(&graph.tree_edges.subtract(&graph.initial_state_cut), |_| -1),
    );

    // The direct root is the stored production CFF itself. A generalized
    // numerator sample may leave an edge undirected, but that is part of
    // its authoritative contour and cannot trigger representative
    // reconstruction. Both generation modes keep that complete map key as
    // an opaque branch while the direct Taylor forest is replayed.
    for explicit_orientation_sum_only in [false, true] {
        let localizer = Localizer::new(
            &cutset,
            OrientationProjection::exact_expression(
                &production,
                &options,
                &pattern,
                explicit_orientation_sum_only,
            ),
        );
        let localized = Direct3dCts::root(&mut graph, localizer)?.branches()?;
        let expected_keys = raw_root
            .terms
            .values()
            .flat_map(|term| &term.orientations)
            .map(|raw| {
                (
                    raw.production_orientation_id
                        .expect("stored production CFF entries retain their original ID"),
                    raw.orientation.edge_energy_map.clone(),
                )
            })
            .fold(Vec::new(), |mut keys, key| {
                if !keys.contains(&key) {
                    keys.push(key);
                }
                keys
            });
        let actual_keys = localized
            .iter_keys()
            .map(|(key, _)| {
                (
                    key.selector_host,
                    production.expression.orientations[key.selector_host]
                        .edge_energy_map
                        .clone(),
                )
            })
            .collect::<Vec<_>>();
        assert_eq!(actual_keys.len(), expected_keys.len());
        assert!(expected_keys.iter().all(|key| actual_keys.contains(key)));
        for (index, term) in &raw_root.terms {
            for raw in &term.orientations {
                let id = raw
                    .production_orientation_id
                    .expect("stored production CFF entries retain their original ID");
                let original = &production.expression.orientations[id];
                assert_eq!(raw.orientation.loop_energy_map, original.loop_energy_map);
                assert_eq!(raw.orientation.edge_energy_map, original.edge_energy_map);
                let (key, integrands) = localized
                    .iter_keys()
                    .find(|(key, _)| key.selector_host == id)
                    .expect("the direct root preserves each stored ID and exact energy map");
                assert_eq!(key.source_edge_energy_map(), None);
                let actual = integrands
                    .iter()
                    .find_map(|(actual_index, atom)| (actual_index == index).then_some(atom))
                    .expect("the direct root preserves every selected CutCFF index");
                let expected = &raw.expression * &production_prefactor * &fourddenoms;
                assert_eq!(actual, &expected);
            }
        }
    }
    Ok(())
}

#[test]
fn orientation_term_keeps_external_selectors_and_adds_internal_ones() {
    let reduced_expression = function!(GS.ose, 0);
    let reduced_orientation = edgevec([1, 0, -1]);
    let valid = vec![edgevec([1, 1, -1]), edgevec([1, -1, -1])];
    let pat = OrientationPattern::default();
    let cutset = CutSet::empty(1);
    let loc = Localizer::new(&cutset, OrientationProjection::new(&valid, &pat));
    let localized = loc
        .localized_orientation_term(&reduced_expression, &reduced_orientation, &edges([1]))
        .unwrap();

    let expected = reduced_expression
        * GS.sign_theta(GS.sign(EdgeIndex(0)))
        * GS.sign_theta(-GS.sign(EdgeIndex(2)))
        * GS.sign_theta(GS.sign(EdgeIndex(1)));
    assert_eq!(localized, expected);
}

#[test]
fn zero_sampling_maps_keep_loop_lifts_as_provenance_only() -> Result<()> {
    let graph = two_edge_graph()?;
    let production = [[1, -1], [-1, 1]]
        .into_iter()
        .map(|directions| {
            energy_map(
                edgevec(directions),
                directions
                    .into_iter()
                    .enumerate()
                    .map(|(edge, direction)| {
                        LinearEnergyExpr::ose(EdgeIndex(edge), i64::from(direction))
                    })
                    .collect(),
            )
        })
        .collect::<TiVec<OrientationID, _>>();
    let options = graph.denominator_only_cff_3d_expression_options();
    let pattern = OrientationPattern::default();
    let cutset = CutSet::empty(graph.n_hedges());
    let localizer = Localizer::new(
        &cutset,
        OrientationProjection::exact(&production, &options, &pattern, false),
    );
    let surviving_edge = EdgeIndex(1);
    let expected_compatible_hosts = vec![OrientationID(0), OrientationID(1)];
    let mut deterministic_hosts = Vec::new();

    for loop_energy in [
        LinearEnergyExpr::ose(surviving_edge, 1),
        LinearEnergyExpr::ose(surviving_edge, -1),
    ] {
        let mut reduced = energy_map(
            edgevec([0, 0]),
            vec![LinearEnergyExpr::zero(), LinearEnergyExpr::zero()],
        );
        reduced.loop_energy_map = vec![loop_energy];
        let representatives =
            localizer.source_selector_representatives(&graph, &reduced, &graph.empty_subgraph())?;
        assert_eq!(
            representatives, expected_compatible_hosts,
            "a loop lift is residue provenance and must not partition an undirected physical prefix",
        );
        let hosted = localizer.localized_orientation_terms(
            &graph,
            &reduced,
            &Atom::one(),
            &graph.empty_subgraph(),
            &[],
            None,
            None,
            Some(&reduced.edge_energy_map),
        )?;
        assert_eq!(hosted.len(), 1);
        assert_eq!(hosted[0].1, Atom::one());
        deterministic_hosts.push(hosted[0].0);
    }

    assert_eq!(
        deterministic_hosts,
        vec![OrientationID(0), OrientationID(0)],
        "the same physical-prefix score deterministically hosts either provenance lift",
    );
    Ok(())
}

#[test]
fn source_selector_is_invariant_under_source_generation_lmb() -> Result<()> {
    let graph = two_edge_graph()?;
    let source_lmb = graph
        .generate_loop_momentum_bases_of(&graph.full_filter())
        .into_iter()
        .find(|lmb| lmb.loop_edges != graph.loop_momentum_basis.loop_edges)
        .ok_or_else(|| eyre!("the selector-hosting fixture needs a non-global LMB"))?;
    let source_edge = *source_lmb
        .loop_edges
        .first()
        .expect("the two-edge graph has one loop");
    let production = [[1, -1], [-1, 1]]
        .into_iter()
        .map(|directions| {
            energy_map(
                edgevec(directions),
                directions
                    .into_iter()
                    .enumerate()
                    .map(|(edge, direction)| {
                        LinearEnergyExpr::ose(EdgeIndex(edge), i64::from(direction))
                    })
                    .collect(),
            )
        })
        .collect::<TiVec<OrientationID, _>>();
    let mut reduced = energy_map(
        edgevec([0, 0]),
        vec![LinearEnergyExpr::zero(), LinearEnergyExpr::zero()],
    );
    reduced.loop_energy_map = vec![LinearEnergyExpr::ose(source_edge, 1)];
    let interpret_source_edge = |lmb: &crate::graph::LoopMomentumBasis| {
        lmb.edge_signatures[source_edge]
            .try_compute_momentum(&reduced.loop_energy_map, &[])
            .expect("the source loop edge depends on its loop coordinate")
            .canonical()
    };
    assert_eq!(
        interpret_source_edge(&source_lmb),
        LinearEnergyExpr::ose(source_edge, 1),
    );
    assert_ne!(
        interpret_source_edge(&graph.loop_momentum_basis),
        LinearEnergyExpr::ose(source_edge, 1),
        "the fixture must distinguish the child source LMB from the graph LMB",
    );

    let options = graph.denominator_only_cff_3d_expression_options();
    let pattern = OrientationPattern::default();
    let cutset = CutSet::empty(graph.n_hedges());
    let localizer = Localizer::new(
        &cutset,
        OrientationProjection::exact(&production, &options, &pattern, true),
    );
    let valid_hosts = BTreeSet::from([OrientationID(0), OrientationID(1)]);
    let mut hosted_by_lift = Vec::new();
    for loop_sign in [1, -1] {
        reduced.loop_energy_map = vec![LinearEnergyExpr::ose(source_edge, loop_sign)];
        hosted_by_lift.push(localizer.localized_orientation_terms(
            &graph,
            &reduced,
            &Atom::one(),
            &graph.empty_subgraph(),
            &[],
            Some(&valid_hosts),
            None,
            Some(&reduced.edge_energy_map),
        )?);
    }

    assert_eq!(hosted_by_lift[0], vec![(OrientationID(0), Atom::one())]);
    assert_eq!(
        hosted_by_lift[1], hosted_by_lift[0],
        "source-LMB loop coordinates remain provenance and cannot change the deterministic physical-prefix host",
    );
    Ok(())
}

#[test]
fn projected_source_sum_uses_a_host_without_imposing_a_global_orientation() -> Result<()> {
    let graph = two_edge_graph()?;
    let production = [[1, -1], [-1, 1]]
        .into_iter()
        .map(|directions| {
            energy_map(
                edgevec(directions),
                directions
                    .into_iter()
                    .enumerate()
                    .map(|(edge, direction)| {
                        LinearEnergyExpr::ose(EdgeIndex(edge), i64::from(direction))
                    })
                    .collect(),
            )
        })
        .collect::<TiVec<OrientationID, _>>();
    let reduced = energy_map(
        edgevec([1, 1]),
        vec![
            LinearEnergyExpr::ose(EdgeIndex(0), 1),
            LinearEnergyExpr::ose(EdgeIndex(1), 1),
        ],
    );
    let options = graph.denominator_only_cff_3d_expression_options();
    let pattern = OrientationPattern::default();
    let cutset = CutSet::empty(graph.n_hedges());
    let direct = Localizer::new(
        &cutset,
        OrientationProjection::exact(&production, &options, &pattern, true),
    );
    assert!(
        direct
            .source_selector_representatives(&graph, &reduced, &graph.empty_subgraph())
            .is_err(),
        "the reduced directions deliberately have no complete production extension",
    );

    let projected = direct.with_independent_source_sum();
    let valid_host = BTreeSet::from([OrientationID(1)]);
    let hosted = projected.localized_orientation_terms(
        &graph,
        &reduced,
        &Atom::one(),
        &graph.empty_subgraph(),
        &[],
        Some(&valid_host),
        None,
        Some(&reduced.edge_energy_map),
    )?;
    assert_eq!(hosted, vec![(OrientationID(1), Atom::one())]);
    Ok(())
}

#[test]
fn projected_source_sum_prefers_a_compatible_physical_host() -> Result<()> {
    let graph = two_edge_graph()?;
    let production = [[1, -1], [-1, 1]]
        .into_iter()
        .map(|directions| {
            energy_map(
                edgevec(directions),
                directions
                    .into_iter()
                    .enumerate()
                    .map(|(edge, direction)| {
                        LinearEnergyExpr::ose(EdgeIndex(edge), i64::from(direction))
                    })
                    .collect(),
            )
        })
        .collect::<TiVec<OrientationID, _>>();
    let reduced = production[OrientationID(0)].clone();
    let options = graph.denominator_only_cff_3d_expression_options();
    let pattern = OrientationPattern::default();
    let cutset = CutSet::empty(graph.n_hedges());
    let projected = Localizer::new(
        &cutset,
        OrientationProjection::exact(&production, &options, &pattern, true),
    )
    .with_independent_source_sum();
    let hosted = projected.localized_orientation_terms(
        &graph,
        &reduced,
        &Atom::one(),
        &graph.empty_subgraph(),
        &[],
        None,
        None,
        Some(&reduced.edge_energy_map),
    )?;
    assert_eq!(hosted, vec![(OrientationID(0), Atom::one())]);
    Ok(())
}

#[test]
fn projected_source_sum_falls_back_when_the_cut_excludes_compatible_hosts() -> Result<()> {
    let graph = two_edge_graph()?;
    let production = [[1, -1], [-1, 1]]
        .into_iter()
        .map(|directions| {
            energy_map(
                edgevec(directions),
                directions
                    .into_iter()
                    .enumerate()
                    .map(|(edge, direction)| {
                        LinearEnergyExpr::ose(EdgeIndex(edge), i64::from(direction))
                    })
                    .collect(),
            )
        })
        .collect::<TiVec<OrientationID, _>>();
    let reduced = production[OrientationID(0)].clone();
    let options = graph.denominator_only_cff_3d_expression_options();
    let pattern = OrientationPattern::default();
    let cutset = CutSet::empty(graph.n_hedges());
    let direct = Localizer::new(
        &cutset,
        OrientationProjection::exact(&production, &options, &pattern, true),
    );
    let valid_host = BTreeSet::from([OrientationID(1)]);
    let direct_hosted = direct.localized_orientation_terms(
        &graph,
        &reduced,
        &Atom::one(),
        &graph.empty_subgraph(),
        &[],
        Some(&valid_host),
        None,
        Some(&reduced.edge_energy_map),
    )?;
    assert!(direct_hosted.is_empty());

    let projected_hosted = direct
        .with_independent_source_sum()
        .localized_orientation_terms(
            &graph,
            &reduced,
            &Atom::one(),
            &graph.empty_subgraph(),
            &[],
            Some(&valid_host),
            None,
            Some(&reduced.edge_energy_map),
        )?;
    assert_eq!(projected_hosted, vec![(OrientationID(1), Atom::one())]);
    Ok(())
}

#[test]
fn nested_powered_contact_keeps_loop_lift_as_provenance_only() -> Result<()> {
    let mut graph = two_edge_graph()?;
    let powered_edge = EdgeIndex(0);
    let on_shell_energy_squared = (1..=3).fold(
        graph.underlying[powered_edge].particle.mass_atom().pow(2),
        |norm_squared, spatial_index| {
            norm_squared + GS.emr_mom(powered_edge, GS.cind(spatial_index)).pow(2)
        },
    );
    let numerator = GS.emr_mom(powered_edge, GS.cind(0)).pow(2) - on_shell_energy_squared;
    let options = graph.denominator_only_cff_3d_expression_options();
    let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
    let production = graph.generate_3d_expression_for_integrand(
        &[],
        &canonization,
        &options,
        Some(&Atom::one()),
    )?;
    let pattern = OrientationPattern::default();
    let cutset = CutSet::empty(graph.n_hedges());
    let localizer = Localizer::new(
        &cutset,
        OrientationProjection::exact_expression(&production, &options, &pattern, true),
    );
    let source_denominator = || FourDDenominator {
        source_edge: powered_edge,
        momentum: FunctionBuilder::new(GS.emr_mom)
            .add_arg(usize::from(powered_edge))
            .finish(),
        mass_squared: graph.underlying[powered_edge].particle.mass_atom().pow(2),
        full_expr: Atom::one(),
    };
    let denominators = [source_denominator(), source_denominator()];
    let (exact, contract_subgraph) =
        graph.cff_from_4d_denominators(&denominators, &cutset, &options, &numerator)?;
    let internal_edges = graph.paired_edges(&contract_subgraph);
    let term = exact
        .terms
        .values()
        .next()
        .expect("the uncut powered source has one CFF term");
    let mut shared_hosts = None;
    let mut contact_count = 0;
    let mut remainder_count = 0;
    let mut saw_zero_canonical_contact = false;
    for orientation in &term.orientations {
        let is_contact = orientation.orientation.variants.iter().any(|variant| {
            variant
                .origin
                .as_deref()
                .is_some_and(|origin| origin.contains("contact"))
        });
        let is_remainder = orientation.orientation.variants.iter().any(|variant| {
            variant
                .origin
                .as_deref()
                .is_some_and(|origin| origin.contains("remainder"))
        });
        assert_ne!(is_contact, is_remainder);
        let canonical_sample = term.map_exact_source_atom(
            &orientation.orientation,
            &GS.emr_mom(powered_edge, GS.cind(0)),
        )?;
        let physical_loop_lift_energies = term.map_exact_source_physical_loop_lift_energies(
            &orientation.orientation,
            [powered_edge],
        )?;
        let powered_loop_lift = physical_loop_lift_energies
            .iter()
            .find_map(|(edge, energy)| (*edge == powered_edge).then_some(energy));
        let powered_loop_lift = powered_loop_lift
            .expect("the contact owner is reconstructed from the exact source loop lift");
        assert!(
            orientation
                .orientation
                .loop_energy_map
                .iter()
                .any(|energy| energy.clone().canonical() != LinearEnergyExpr::zero()),
            "each powered-source branch must retain a nonzero residue loop lift"
        );
        assert!(!powered_loop_lift.is_zero());
        if is_contact && canonical_sample.is_zero() {
            saw_zero_canonical_contact = true;
            assert_eq!(powered_loop_lift, &GS.ose(powered_edge));
        }
        let hosts = localizer.source_selector_representatives(
            &graph,
            &orientation.orientation,
            &contract_subgraph,
        )?;
        let host_set = hosts.iter().copied().collect::<BTreeSet<_>>();
        if let Some(expected) = &shared_hosts {
            assert_eq!(
                &host_set, expected,
                "synthetic contact and remainder directions must not partition physical production hosts"
            );
        } else {
            shared_hosts = Some(host_set);
        }
        if is_contact {
            contact_count += 1;
        } else {
            remainder_count += 1;
        }
        let selected_host = hosts[0];
        let valid_host = BTreeSet::from([selected_host]);
        let hosted = localizer.localized_orientation_terms(
            &graph,
            &orientation.orientation,
            &Atom::one(),
            &contract_subgraph,
            &internal_edges,
            Some(&valid_host),
            orientation.production_orientation_id,
            Some(&orientation.orientation.edge_energy_map),
        )?;
        assert_eq!(hosted.len(), 1);
        assert_eq!(hosted[0].0, selected_host);
    }
    assert!(
        saw_zero_canonical_contact,
        "at least one recursive contact must expose the stale zero numerator sample"
    );
    assert_eq!(contact_count, 3);
    assert_eq!(remainder_count, 2);
    assert_eq!(
        shared_hosts
            .expect("the exact powered source has contour pieces")
            .len(),
        2,
        "every algebraic contour piece may be placed under either physical production host; the caller selects one deterministically"
    );
    Ok(())
}

#[test]
fn exact_projection_distinguishes_affine_maps_with_the_same_coarse_orientation() -> Result<()> {
    let graph = two_edge_graph()?;
    let production = vec![
        energy_map(
            edgevec([1, 1]),
            vec![
                LinearEnergyExpr::ose(EdgeIndex(0), 1),
                LinearEnergyExpr::ose(EdgeIndex(1), 1),
            ],
        ),
        energy_map(
            edgevec([1, 1]),
            vec![
                LinearEnergyExpr::ose(EdgeIndex(0), 2),
                LinearEnergyExpr::ose(EdgeIndex(1), 1),
            ],
        ),
    ]
    .into_iter()
    .collect::<TiVec<OrientationID, _>>();
    let reduced = energy_map(
        edgevec([1, 0]),
        vec![
            LinearEnergyExpr::ose(EdgeIndex(0), 1),
            LinearEnergyExpr::zero(),
        ],
    );
    let options = graph.denominator_only_cff_3d_expression_options();
    let pattern = OrientationPattern::default();
    let cutset = CutSet::empty(graph.n_hedges());
    let localizer = Localizer::new(
        &cutset,
        OrientationProjection::exact(&production, &options, &pattern, false),
    );

    assert_eq!(
        localizer.exact_representatives(
            &graph,
            &reduced,
            &graph.get_edge_subgraph(EdgeIndex(1)),
        )?,
        vec![OrientationID(0)]
    );
    Ok(())
}

#[test]
fn exact_projection_errors_when_no_affine_map_matches() -> Result<()> {
    let graph = two_edge_graph()?;
    let production = vec![energy_map(
        edgevec([1, 1]),
        vec![
            LinearEnergyExpr::ose(EdgeIndex(0), 1),
            LinearEnergyExpr::ose(EdgeIndex(1), 1),
        ],
    )]
    .into_iter()
    .collect::<TiVec<OrientationID, _>>();
    let reduced = energy_map(
        edgevec([1, 0]),
        vec![
            LinearEnergyExpr::ose(EdgeIndex(0), 2),
            LinearEnergyExpr::zero(),
        ],
    );
    let options = graph.denominator_only_cff_3d_expression_options();
    let pattern = OrientationPattern::default();
    let cutset = CutSet::empty(graph.n_hedges());
    let localizer = Localizer::new(
        &cutset,
        OrientationProjection::exact(&production, &options, &pattern, false),
    );

    let error = localizer
        .exact_representatives(&graph, &reduced, &graph.get_edge_subgraph(EdgeIndex(1)))
        .expect_err("a different surviving affine map must not use a coarse fallback");
    assert!(
        error
            .to_string()
            .contains("no production energy map exactly extends")
    );
    Ok(())
}

#[test]
fn source_energy_map_owns_factorized_contact_branch() -> Result<()> {
    let graph = two_edge_graph()?;
    let production = vec![
        energy_map(
            edgevec([1, 1]),
            vec![
                LinearEnergyExpr::ose(EdgeIndex(0), 1),
                LinearEnergyExpr::ose(EdgeIndex(1), 1),
            ],
        ),
        energy_map(
            edgevec([1, -1]),
            vec![
                LinearEnergyExpr::ose(EdgeIndex(0), 1),
                LinearEnergyExpr::ose(EdgeIndex(1), -1),
            ],
        ),
    ]
    .into_iter()
    .collect::<TiVec<OrientationID, _>>();
    let reduced = energy_map(
        edgevec([1, 0]),
        vec![LinearEnergyExpr::uniform_scale(2), LinearEnergyExpr::zero()],
    );
    let options = graph.denominator_only_cff_3d_expression_options();
    let pattern = OrientationPattern::default();
    let cutset = CutSet::empty(graph.n_hedges());
    let localizer = Localizer::new(
        &cutset,
        OrientationProjection::exact(&production, &options, &pattern, false),
    );
    let contract = graph.get_edge_subgraph(EdgeIndex(1));

    let strict_error = localizer
        .exact_representatives(&graph, &reduced, &contract)
        .expect_err("the source map must not weaken strict production-map matching");
    assert!(
        strict_error
            .to_string()
            .contains("no production energy map exactly extends")
    );
    let localized = localizer.localized_orientation_terms(
        &graph,
        &reduced,
        &Atom::one(),
        &contract,
        &[EdgeIndex(1)],
        None,
        None,
        Some(&reduced.edge_energy_map),
    )?;
    assert_eq!(localized, vec![(OrientationID(0), Atom::one())]);
    let selector_index = CutCFFIndex::new_all_none();
    let transient = OrientationIntegrands(vec![OrientationIntegrandBranch {
        selector_id: localized[0].0,
        source_edge_energy_map: Some(reduced.edge_energy_map.clone()),
        integrands: [(selector_index, localized[0].1.clone())]
            .into_iter()
            .collect(),
    }]);
    let keyed = DirectResidueBranches::from_transient(&transient)?;
    assert_eq!(
        keyed.materialize(false)?.iter().next(),
        Some((&selector_index, &Atom::one())),
        "an explicit orientation sum exposes the unchanged pre-T body",
    );
    let selected = keyed.materialize(true)?;
    let selected_body = selected
        .iter()
        .next()
        .expect("the localized residue retains its cut support")
        .1;
    assert_eq!(
        OrientationID(0).select(selected_body.as_view()),
        Atom::one()
    );
    assert_eq!(
        OrientationID(1).select(selected_body.as_view()),
        Atom::Zero,
        "late residue-key materialization must not leak into another key",
    );

    let factor_a = Atom::var(symbolica::symbol!("source_factor_a"));
    let factor_b = Atom::var(symbolica::symbol!("source_factor_b"));
    let edge_energy = GS.emr_mom(EdgeIndex(0), GS.cind(0));
    let numerator = (edge_energy.clone() + &factor_a) * (edge_energy + &factor_b);
    let mapped = localizer.map_numerator(
        &graph,
        OrientationID(0),
        Some(&reduced.edge_energy_map),
        &numerator,
    )?;
    assert_eq!(
        localizer.map_numerator(
            &graph,
            OrientationID(1),
            Some(&reduced.edge_energy_map),
            &numerator,
        )?,
        mapped,
        "production IDs partition complete residue-map keys but do not remap source-owned numerators"
    );
    let sampled_energy = Atom::num(2) * Atom::var(GS.numerator_sampling_scale);
    assert_eq!(
        mapped,
        (sampled_energy.clone() + &factor_a) * (sampled_energy + &factor_b)
    );
    let AtomView::Mul(factors) = mapped.as_view() else {
        panic!("the source-mapped numerator must remain factorized")
    };
    assert_eq!(
        factors
            .iter()
            .filter(|factor| matches!(factor, AtomView::Add(_)))
            .count(),
        2
    );
    assert_ne!(
        localizer.map_numerator(&graph, OrientationID(0), None, &numerator)?,
        mapped
    );
    let coarse = Localizer::new(&cutset, OrientationProjection::new(&[], &pattern));
    assert!(
        coarse
            .map_numerator(
                &graph,
                OrientationID(0),
                Some(&reduced.edge_energy_map),
                &numerator,
            )
            .expect_err("a coarse projector cannot consume an exact source map")
            .to_string()
            .contains("cannot be carried by a coarse orientation projector")
    );

    let index = CutCFFIndex::new_all_none();
    let source = OrientationIntegrands(vec![OrientationIntegrandBranch {
        selector_id: OrientationID(0),
        source_edge_energy_map: Some(reduced.edge_energy_map.clone()),
        integrands: [(index, Atom::num(3))].into_iter().collect(),
    }]);
    let distinct_map = vec![LinearEnergyExpr::zero(), LinearEnergyExpr::zero()];
    let distinct = OrientationIntegrands(vec![OrientationIntegrandBranch {
        selector_id: OrientationID(0),
        source_edge_energy_map: Some(distinct_map.clone()),
        integrands: [(index, Atom::num(5))].into_iter().collect(),
    }]);
    let production_zero = OrientationIntegrands(vec![OrientationIntegrandBranch {
        selector_id: OrientationID(0),
        source_edge_energy_map: None,
        integrands: [(index, Atom::Zero)].into_iter().collect(),
    }]);
    let union = source.zip_add(&distinct)?.zip_add(&production_zero)?;
    assert_eq!(union.iter_orientations().count(), 2);
    assert!(
        union
            .iter_orientations()
            .all(|(_, source_map, _)| source_map.is_some()),
        "a zero production-map identity branch is pruned beside source-owned terms"
    );
    assert!(
        union
            .iter_orientations()
            .any(|(_, map, _)| map == Some(reduced.edge_energy_map.as_slice()))
    );
    assert!(
        union
            .iter_orientations()
            .any(|(_, map, _)| map == Some(distinct_map.as_slice()))
    );

    let same_source = OrientationIntegrands(vec![OrientationIntegrandBranch {
        selector_id: OrientationID(0),
        source_edge_energy_map: Some(reduced.edge_energy_map.clone()),
        integrands: [(index, Atom::num(7))].into_iter().collect(),
    }]);
    let merged = source.zip_add(&same_source)?;
    assert_eq!(merged.iter_orientations().count(), 1);
    assert_eq!(
        merged.iter().next(),
        Some((&index, &Atom::num(10))),
        "only identical selector/source-map branches are coalesced"
    );
    Ok(())
}

#[test]
fn factorized_products_keep_only_shared_production_hosts() -> Result<()> {
    let index = CutCFFIndex::new_all_none();
    let branch = |selector, value| OrientationIntegrandBranch {
        selector_id: OrientationID(selector),
        source_edge_energy_map: None,
        integrands: [(index, Atom::num(value))].into_iter().collect(),
    };
    let left = OrientationIntegrands(vec![branch(0, 2), branch(1, 3)]);
    let right = OrientationIntegrands(vec![branch(1, 5), branch(2, 7)]);
    let mut mapped_hosts = Vec::new();
    let product = left.zip_mul_mapped_factor(&right, |host, _, factor| {
        mapped_hosts.push(host);
        Ok(factor.clone())
    })?;

    assert_eq!(mapped_hosts, vec![OrientationID(1)]);
    assert_eq!(
        product
            .iter_orientations()
            .map(|(host, _, _)| host)
            .collect::<Vec<_>>(),
        vec![OrientationID(1)],
        "factorized multiplication must neither retain nor synthesize unmatched hosts"
    );
    assert_eq!(product.iter().next(), Some((&index, &Atom::num(15))));
    Ok(())
}

fn energy_bounds(graph: &Graph, atom: &Atom) -> Result<Vec<(usize, usize)>> {
    Ok(
        graph.automatic_numerator_energy_degree_bounds_in_atom_excluding_with_min_degree(
            atom,
            [],
            1,
        )?,
    )
}

#[test]
fn outer_cff_capacity_does_not_cancel_between_selector_branches() -> Result<()> {
    let graph = two_edge_graph()?;
    let energy = GS.emr_mom(EdgeIndex(0), GS.cind(0));
    let factorized = (energy.clone() + Atom::num(1)) * (energy.clone() + Atom::num(2));
    let root = CutCFFIndex::new_all_none();
    let branches = OrientationIntegrands(vec![
        OrientationIntegrandBranch {
            selector_id: OrientationID(0),
            source_edge_energy_map: None,
            integrands: [(root, factorized.clone())].into_iter().collect(),
        },
        OrientationIntegrandBranch {
            selector_id: OrientationID(1),
            source_edge_energy_map: None,
            integrands: [(root, -&factorized + &energy)].into_iter().collect(),
        },
    ]);

    assert_eq!(
        energy_bounds(&graph, &branches.factorized_sum())?,
        vec![(0, 1)]
    );
    assert_eq!(
        energy_bounds(&graph, &branches.factorized_capacity_envelope())?,
        vec![(0, 2)],
        "mutually exclusive selector branches need the maximum of their separate ranks"
    );
    assert_eq!(
        branches
            .iter_orientations()
            .next()
            .and_then(|(_, _, integrands)| integrands.iter().next())
            .map(|(_, atom)| atom),
        Some(&factorized),
        "capacity analysis must not expand or rewrite the stored factorized numerator"
    );
    Ok(())
}

#[test]
fn outer_cff_capacity_does_not_cancel_between_cut_orders() -> Result<()> {
    let graph = two_edge_graph()?;
    let energy = GS.emr_mom(EdgeIndex(0), GS.cind(0));
    let cubic = (energy.clone() + Atom::num(1))
        * (energy.clone() + Atom::num(2))
        * (energy.clone() + Atom::num(3));
    let raised = CutCFFIndex {
        lu_cut_order: Some(1),
        ..CutCFFIndex::new_all_none()
    };
    let branches = OrientationIntegrands(vec![OrientationIntegrandBranch {
        selector_id: OrientationID(0),
        source_edge_energy_map: None,
        integrands: [
            (CutCFFIndex::new_all_none(), cubic.clone()),
            (raised, -&cubic + &energy),
        ]
        .into_iter()
        .collect(),
    }]);

    assert_eq!(
        energy_bounds(&graph, &branches.factorized_sum())?,
        vec![(0, 1)]
    );
    assert_eq!(
        energy_bounds(&graph, &branches.factorized_capacity_envelope())?,
        vec![(0, 3)],
        "separately evaluated CutCFFIndex values need the maximum of their separate ranks"
    );
    Ok(())
}

#[test]
fn exact_projection_skips_extensions_excluded_by_the_full_pattern() -> Result<()> {
    let graph = two_edge_graph()?;
    let production = vec![
        energy_map(
            edgevec([1, 1]),
            vec![
                LinearEnergyExpr::ose(EdgeIndex(0), 1),
                LinearEnergyExpr::ose(EdgeIndex(1), 1),
            ],
        ),
        energy_map(
            edgevec([-1, -1]),
            vec![
                LinearEnergyExpr::ose(EdgeIndex(0), -1),
                LinearEnergyExpr::ose(EdgeIndex(1), -1),
            ],
        ),
    ]
    .into_iter()
    .collect::<TiVec<OrientationID, _>>();
    let reduced = energy_map(
        edgevec([1, 0]),
        vec![
            LinearEnergyExpr::ose(EdgeIndex(0), 1),
            LinearEnergyExpr::zero(),
        ],
    );
    let options = graph.denominator_only_cff_3d_expression_options();
    let pattern = OrientationPattern::from_user_pattern("(-1,-1)")?;
    let cutset = CutSet::empty(graph.n_hedges());
    let localizer = Localizer::new(
        &cutset,
        OrientationProjection::exact(&production, &options, &pattern, false),
    );
    let contract = graph.get_edge_subgraph(EdgeIndex(1));

    assert!(
        localizer
            .exact_representatives(&graph, &reduced, &contract)?
            .is_empty(),
        "an exact extension excluded by the full pattern is a zero contribution, not a missing-map error"
    );
    assert!(
        localizer
            .localized_orientation_terms(
                &graph,
                &reduced,
                &Atom::one(),
                &contract,
                &[EdgeIndex(1)],
                None,
                None,
                None,
            )?
            .is_empty(),
        "an excluded extension contributes zero without normalization"
    );
    Ok(())
}

#[test]
fn exact_cff_defers_contracted_edge_patterns_to_full_projection() -> Result<()> {
    let mut graph = two_edge_graph()?;
    let production = vec![
        energy_map(
            edgevec([1, -1]),
            vec![
                LinearEnergyExpr::ose(EdgeIndex(0), 1),
                LinearEnergyExpr::ose(EdgeIndex(1), -1),
            ],
        ),
        energy_map(
            edgevec([-1, 1]),
            vec![
                LinearEnergyExpr::ose(EdgeIndex(0), -1),
                LinearEnergyExpr::ose(EdgeIndex(1), 1),
            ],
        ),
    ]
    .into_iter()
    .collect::<TiVec<OrientationID, _>>();
    let options = graph.denominator_only_cff_3d_expression_options();
    let pattern = OrientationPattern::from_user_pattern("(1,-1)")?;
    let cutset = CutSet::empty(graph.n_hedges());
    let localizer = Localizer::new(
        &cutset,
        OrientationProjection::exact(&production, &options, &pattern, false),
    );
    let to_contract = graph.get_edge_subgraph(EdgeIndex(1));
    let analysis_numerator = graph.production_numerator_atom_for_full_3d_expression();
    let projected = localizer.projected_cff(
        &mut graph,
        &to_contract,
        &analysis_numerator,
        CffGenerationContext::Standalone,
    )?;

    assert_eq!(
        projected
            .iter_orientations()
            .map(|(id, _, _)| id)
            .collect::<Vec<_>>(),
        vec![OrientationID(0)]
    );
    assert!(projected.iter().any(|(_, atom)| !atom.is_zero()));
    Ok(())
}

#[test]
fn coarse_cff_keeps_denominator_only_capacity() -> Result<()> {
    let mut graph = two_edge_graph()?;
    let pattern = OrientationPattern::default();
    let cutset = CutSet::empty(graph.n_hedges());
    let localizer = Localizer::new(&cutset, OrientationProjection::new(&[], &pattern));
    let to_contract = graph.empty_subgraph();
    let unsupported_numerator = (GS.emr_mom(EdgeIndex(0), GS.cind(0)) + Atom::one()).pow(-1);

    let projected = localizer.projected_cff(
        &mut graph,
        &to_contract,
        &unsupported_numerator,
        CffGenerationContext::Standalone,
    )?;

    assert!(projected.iter().any(|(_, atom)| !atom.is_zero()));
    Ok(())
}

#[test]
fn exact_localization_capacity_excludes_replaced_spinney_numerator() -> Result<()> {
    // Share model initialization with the other graph-backed tests in this module.
    let _ = two_edge_graph()?;
    let mut graph: Graph = dot!(
        digraph replaced_spinney_numerator {
            edge [pdg=1000 num=1 mass=1]
            node [num=1]
            a -> b [id=0 lmb_id=0 num="1/(Q(0,spenso::cind(0))+1)"]
            a -> b [id=1]
        },
        "scalars"
    )?;
    let options = graph.denominator_only_cff_3d_expression_options();
    let contracted = graph.tree_edges.subtract(&graph.initial_state_cut);
    let contract_edges = graph.paired_edges(&contracted);
    let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
    let production = graph
        .generate_3d_expression_for_integrand(&contract_edges, &canonization, &options, None)?
        .expression
        .orientations;
    let pattern = OrientationPattern::default();
    let cutset = CutSet::empty(graph.n_hedges());
    let localizer = Localizer::new(
        &cutset,
        OrientationProjection::exact(&production, &options, &pattern, false),
    );
    let subgraph = InternalSubGraph::cleaned_filter_optimist(graph.full_filter(), graph.as_ref());
    let integrated_node = OwnedForestNode {
        spinney: Spinney::with_scheme(
            subgraph,
            &graph,
            &graph.loop_momentum_basis,
            ApproximationType::MUV,
            0,
        )
        .expect("the complete bubble has a compatible loop basis"),
        topo_order: 0,
    };

    let localized = localizer.localize(&Atom::one(), &mut graph, &integrated_node)?;

    assert!(localized.active.iter().any(|(_, atom)| !atom.is_zero()));
    Ok(())
}

#[test]
fn exact_localization_maps_finite_ct_instead_of_leaving_it_unmapped() -> Result<()> {
    // Share model initialization with the other graph-backed tests in this module.
    let _ = two_edge_graph()?;
    let mut graph: Graph = dot!(
        digraph finite_ct_map {
            edge [particle="scalar_1"];
            node [num=1];
            external [style=invis];
            external -> A:0 [id=3];
            C:1 -> external [id=4];
            A -> B [id=0];
            A -> B [id=1];
            B -> C [id=2];
        },
        "scalars"
    )?;
    let options = graph.denominator_only_cff_3d_expression_options();
    let contracted = graph.tree_edges.subtract(&graph.initial_state_cut);
    let contract_edges = graph.paired_edges(&contracted);
    assert!(contract_edges.contains(&EdgeIndex(2)));
    let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
    let production = graph
        .generate_3d_expression_for_integrand(&contract_edges, &canonization, &options, None)?
        .expression
        .orientations;
    let (orientation_id, finite_ct, mapped_finite_ct) = production
            .iter_enumerated()
            .flat_map(|(orientation_id, orientation)| {
                orientation
                    .edge_energy_map
                    .iter()
                    .enumerate()
                    .filter(|(edge, energy)| {
                        contract_edges.contains(&EdgeIndex(*edge))
                            && !energy.external_terms.is_empty()
                    })
                    .map(move |(edge, _)| (orientation_id, orientation, EdgeIndex(edge)))
            })
            .find_map(|(orientation_id, orientation, edge)| {
                let finite_ct = GS.emr_mom(edge, GS.cind(0));
                let replacements = orientation.energy_replacements_gs(&graph);
                let mapped_finite_ct = finite_ct.replace_multiple(&replacements);
                (mapped_finite_ct != finite_ct).then_some((
                    orientation_id,
                    finite_ct,
                    mapped_finite_ct,
                ))
            })
            .ok_or_else(|| {
                eyre!(
                    "the finite-CT fixture must contain a contracted affine map that changes an edge energy"
                )
            })?;
    let pattern =
        OrientationPattern::from_orientation(&production[orientation_id].data.orientation);
    let cutset = CutSet::empty(graph.n_hedges());
    let localizer = Localizer::new(
        &cutset,
        OrientationProjection::exact(&production, &options, &pattern, false),
    );
    let integrated_node = OwnedForestNode {
        spinney: Spinney::empty(&graph),
        topo_order: 0,
    };

    let baseline = localizer.localize(&Atom::one(), &mut graph, &integrated_node)?;
    let localized = localizer.localize(&finite_ct, &mut graph, &integrated_node)?;
    let baseline = baseline
        .active
        .iter_orientations()
        .find_map(|(id, _, integrands)| (id == orientation_id).then_some(integrands))
        .expect("the selected production map has a baseline branch");
    let localized = localized
        .active
        .iter_orientations()
        .find_map(|(id, _, integrands)| (id == orientation_id).then_some(integrands))
        .expect("the selected production map has a finite-CT branch");

    assert_eq!(localized, &baseline.map(|atom| atom * &mapped_finite_ct));
    assert_ne!(localized, &baseline.map(|atom| atom * &finite_ct));
    Ok(())
}

#[test]
fn contracted_exact_extensions_follow_evaluator_orientation_mode() -> Result<()> {
    let graph = two_edge_graph()?;
    let production = vec![
        energy_map(
            edgevec([1, 1]),
            vec![
                LinearEnergyExpr::ose(EdgeIndex(0), 1),
                LinearEnergyExpr::ose(EdgeIndex(1), 1),
            ],
        ),
        energy_map(
            edgevec([1, -1]),
            vec![
                LinearEnergyExpr::ose(EdgeIndex(0), 1),
                LinearEnergyExpr::ose(EdgeIndex(1), -1),
            ],
        ),
    ]
    .into_iter()
    .collect::<TiVec<OrientationID, _>>();
    let reduced = energy_map(
        edgevec([1, 0]),
        vec![
            LinearEnergyExpr::ose(EdgeIndex(0), 1),
            LinearEnergyExpr::zero(),
        ],
    );
    let options = graph.denominator_only_cff_3d_expression_options();
    let pattern = OrientationPattern::default();
    let cutset = CutSet::empty(graph.n_hedges());
    let localizer = Localizer::new(
        &cutset,
        OrientationProjection::exact(&production, &options, &pattern, false),
    );
    let contract = graph.get_edge_subgraph(EdgeIndex(1));
    let localized = localizer.localized_orientation_terms(
        &graph,
        &reduced,
        &Atom::one(),
        &contract,
        &[EdgeIndex(1)],
        None,
        None,
        None,
    )?;
    assert_eq!(localized, vec![(OrientationID(0), Atom::one())]);
    assert_eq!(
        localizer.residue_map_key_selector(localized[0].0),
        OrientationID(0).atom(),
        "an exact reduced residue is localized by its complete map key, not by a physical-theta product"
    );

    let explicit_localizer = Localizer::new(
        &cutset,
        OrientationProjection::exact(&production, &options, &pattern, true),
    );
    let explicit = explicit_localizer.localized_orientation_terms(
        &graph,
        &reduced,
        &Atom::one(),
        &contract,
        &[EdgeIndex(1)],
        None,
        None,
        None,
    )?;
    assert_eq!(explicit, vec![(OrientationID(0), Atom::one())]);
    Ok(())
}

#[test]
fn cut_valid_ids_host_one_inner_representative_per_outer_sector() -> Result<()> {
    let graph = two_edge_graph()?;
    let production = [[1, 1], [1, -1], [-1, 1], [-1, -1]]
        .into_iter()
        .map(|directions| {
            energy_map(
                edgevec(directions),
                directions
                    .into_iter()
                    .enumerate()
                    .map(|(edge, direction)| {
                        LinearEnergyExpr::ose(EdgeIndex(edge), i64::from(direction))
                    })
                    .collect(),
            )
        })
        .collect::<TiVec<OrientationID, _>>();
    let reduced = [1, -1].map(|outer_direction| {
        energy_map(
            edgevec([outer_direction, 0]),
            vec![
                LinearEnergyExpr::ose(EdgeIndex(0), i64::from(outer_direction)),
                LinearEnergyExpr::zero(),
            ],
        )
    });
    let valid_ids = BTreeSet::from([OrientationID(1), OrientationID(2)]);
    let options = graph.denominator_only_cff_3d_expression_options();
    let pattern = OrientationPattern::default();
    let cutset = CutSet::empty(graph.n_hedges());
    let contract = graph.get_edge_subgraph(EdgeIndex(1));

    for explicit_orientation_sum_only in [false, true] {
        let localizer = Localizer::new(
            &cutset,
            OrientationProjection::exact(
                &production,
                &options,
                &pattern,
                explicit_orientation_sum_only,
            ),
        );
        let localized = reduced
            .iter()
            .enumerate()
            .map(|(outer_sector, reduced)| {
                let body = Atom::var(symbolica::symbol!(format!("cut_valid_body_{outer_sector}")));
                let terms = localizer.localized_orientation_terms(
                    &graph,
                    reduced,
                    &body,
                    &contract,
                    &[EdgeIndex(1)],
                    Some(&valid_ids),
                    None,
                    Some(&reduced.edge_energy_map),
                )?;
                Ok((body, terms))
            })
            .collect::<Result<Vec<_>>>()?;

        assert!(localized.iter().all(|(_, terms)| terms.len() == 1));
        assert_eq!(localized[0].1[0].0, OrientationID(1));
        assert_eq!(localized[1].1[0].0, OrientationID(2));
        assert!(
            localized
                .iter()
                .flat_map(|(_, terms)| terms)
                .all(|(id, _)| {
                    valid_ids.contains(id) && *id != OrientationID(0) && *id != OrientationID(3)
                })
        );

        for (outer_sector, (body, terms)) in localized.iter().enumerate() {
            let (id, expression) = &terms[0];
            assert_eq!(expression, body);
            if explicit_orientation_sum_only {
                assert_eq!(localizer.residue_map_key_selector(*id), Atom::one());
            } else {
                let selector = localizer.residue_map_key_selector(*id);
                assert_eq!(
                    id.select(selector.as_view()),
                    Atom::one(),
                    "the valid representative owns its complete map-key sector exactly once"
                );
                let invalid_id = [OrientationID(0), OrientationID(3)][outer_sector];
                assert_eq!(
                    invalid_id.select(selector.as_view()),
                    Atom::Zero,
                    "a cut-invalid map key must not receive the body"
                );
            }
        }
    }
    Ok(())
}

#[test]
fn stored_root_residue_keeps_its_production_orientation_diagonal() -> Result<()> {
    let graph = two_edge_graph()?;
    let production = vec![
        energy_map(
            edgevec([1, 1]),
            vec![
                LinearEnergyExpr::ose(EdgeIndex(0), 1),
                LinearEnergyExpr::ose(EdgeIndex(1), 1),
            ],
        ),
        energy_map(
            edgevec([1, -1]),
            vec![
                LinearEnergyExpr::ose(EdgeIndex(0), 1),
                LinearEnergyExpr::ose(EdgeIndex(1), -1),
            ],
        ),
    ]
    .into_iter()
    .collect::<TiVec<OrientationID, _>>();
    let reduced = energy_map(
        edgevec([1, 0]),
        vec![
            LinearEnergyExpr::ose(EdgeIndex(0), 1),
            LinearEnergyExpr::zero(),
        ],
    );
    let options = graph.denominator_only_cff_3d_expression_options();
    let pattern = OrientationPattern::default();
    let cutset = CutSet::empty(graph.n_hedges());
    let localizer = Localizer::new(
        &cutset,
        OrientationProjection::exact(&production, &options, &pattern, false),
    );
    let contract = graph.get_edge_subgraph(EdgeIndex(1));

    let localized = localizer.localized_orientation_terms(
        &graph,
        &reduced,
        &Atom::one(),
        &contract,
        &[EdgeIndex(1)],
        None,
        Some(OrientationID(1)),
        None,
    )?;

    assert_eq!(localized.len(), 1);
    assert_eq!(localized[0].0, OrientationID(1));
    assert_eq!(localized[0].1, Atom::one());
    let selected = localizer.residue_map_key_selector(localized[0].0);
    assert_eq!(
        OrientationID(0).select(selected.as_view()),
        Atom::Zero,
        "a stored root residue must not leak into another residue-map key"
    );
    assert_eq!(OrientationID(1).select(selected.as_view()), Atom::one());

    let explicit_localizer = Localizer::new(
        &cutset,
        OrientationProjection::exact(&production, &options, &pattern, true),
    );
    let explicit = explicit_localizer.localized_orientation_terms(
        &graph,
        &reduced,
        &Atom::one(),
        &contract,
        &[EdgeIndex(1)],
        None,
        Some(OrientationID(1)),
        None,
    )?;
    assert_eq!(explicit, vec![(OrientationID(1), Atom::one())]);
    Ok(())
}

#[test]
fn stored_generalized_root_map_preserves_its_own_selector() -> Result<()> {
    let graph = two_edge_graph()?;
    let production = vec![
        energy_map(
            edgevec([1, 1]),
            vec![
                LinearEnergyExpr::ose(EdgeIndex(0), 1),
                LinearEnergyExpr::ose(EdgeIndex(1), 1),
            ],
        ),
        energy_map(
            edgevec([1, -1]),
            vec![
                LinearEnergyExpr::ose(EdgeIndex(0), 1),
                LinearEnergyExpr::ose(EdgeIndex(1), -1),
            ],
        ),
        energy_map(
            edgevec([1, 0]),
            vec![
                LinearEnergyExpr::ose(EdgeIndex(0), 1),
                LinearEnergyExpr::zero(),
            ],
        ),
    ]
    .into_iter()
    .collect::<TiVec<OrientationID, _>>();
    let reduced = production[OrientationID(2)].clone();
    let options = graph.denominator_only_cff_3d_expression_options();
    let pattern = OrientationPattern::default();
    let cutset = CutSet::empty(graph.n_hedges());
    let localizer = Localizer::new(
        &cutset,
        OrientationProjection::exact(&production, &options, &pattern, false),
    );
    let contract = graph.empty_subgraph();

    let localized = localizer.localized_orientation_terms(
        &graph,
        &reduced,
        &Atom::one(),
        &contract,
        &[],
        None,
        Some(OrientationID(2)),
        Some(&reduced.edge_energy_map),
    )?;

    assert_eq!(
        localized,
        vec![(OrientationID(2), Atom::one())],
        "a stored generalized root map keeps its original ID in the sparse selector sidecar",
    );

    let cut_valid_ids = BTreeSet::from([OrientationID(1), OrientationID(2)]);
    let cut_localized = localizer.localized_orientation_terms(
        &graph,
        &reduced,
        &Atom::one(),
        &contract,
        &[],
        Some(&cut_valid_ids),
        Some(OrientationID(2)),
        Some(&reduced.edge_energy_map),
    )?;
    assert_eq!(
        cut_localized,
        vec![(OrientationID(2), Atom::one())],
        "cut filtering must preserve a cut-valid stored root ID instead of re-hosting it",
    );
    assert_eq!(
        OrientationID(2).select(
            localizer
                .residue_map_key_selector(cut_localized[0].0)
                .as_view(),
        ),
        Atom::one(),
        "the branch's complete map key must select it without resolving its undirected edge",
    );

    let explicit_localizer = Localizer::new(
        &cutset,
        OrientationProjection::exact(&production, &options, &pattern, true),
    );
    assert_eq!(
        explicit_localizer.localized_orientation_terms(
            &graph,
            &reduced,
            &Atom::one(),
            &contract,
            &[],
            None,
            Some(OrientationID(2)),
            Some(&reduced.edge_energy_map),
        )?,
        vec![(OrientationID(2), Atom::one())],
        "an explicit orientation sum keeps the stored generalized branch exactly once without a selector",
    );
    assert_eq!(
        explicit_localizer.localized_orientation_terms(
            &graph,
            &reduced,
            &Atom::one(),
            &contract,
            &[],
            Some(&cut_valid_ids),
            Some(OrientationID(2)),
            Some(&reduced.edge_energy_map),
        )?,
        vec![(OrientationID(2), Atom::one())],
        "the explicit sum retains the stored cut-valid branch metadata without a selector",
    );
    Ok(())
}

#[test]
fn contracted_uv_source_directions_use_late_residue_key_selectors() -> Result<()> {
    let graph = two_edge_graph()?;
    let production = vec![
        energy_map(
            edgevec([1, 1]),
            vec![
                LinearEnergyExpr::ose(EdgeIndex(0), 1),
                LinearEnergyExpr::ose(EdgeIndex(1), 1),
            ],
        ),
        energy_map(
            edgevec([-1, -1]),
            vec![
                LinearEnergyExpr::ose(EdgeIndex(0), -1),
                LinearEnergyExpr::ose(EdgeIndex(1), -1),
            ],
        ),
    ]
    .into_iter()
    .collect::<TiVec<OrientationID, _>>();
    let reduced = energy_map(
        edgevec([-1, 1]),
        vec![
            LinearEnergyExpr::ose(EdgeIndex(0), -1),
            LinearEnergyExpr::ose(EdgeIndex(1), 1),
        ],
    );
    let options = graph.denominator_only_cff_3d_expression_options();
    let pattern = OrientationPattern::default();
    let cutset = CutSet::empty(graph.n_hedges());
    let localizer = Localizer::new(
        &cutset,
        OrientationProjection::exact(&production, &options, &pattern, false),
    );
    let contract = graph.full_filter();
    let localized = localizer.localized_orientation_terms(
        &graph,
        &reduced,
        &Atom::one(),
        &contract,
        &[EdgeIndex(0), EdgeIndex(1)],
        None,
        None,
        None,
    )?;

    assert_eq!(localized, vec![(OrientationID(0), Atom::one())]);
    let index = CutCFFIndex::new_all_none();
    let transient = OrientationIntegrands(vec![OrientationIntegrandBranch {
        selector_id: localized[0].0,
        source_edge_energy_map: None,
        integrands: [(index, localized[0].1.clone())].into_iter().collect(),
    }]);
    let keyed = DirectResidueBranches::from_transient(&transient)?;
    assert_eq!(
        keyed.materialize(false)?.iter().next(),
        Some((&index, &Atom::one())),
        "contracted directions do not add selector atoms to the pre-T body",
    );
    let selected = keyed.materialize(true)?;
    let selected_body = selected
        .iter()
        .next()
        .expect("the localized residue retains its cut support")
        .1;
    assert_eq!(
        OrientationID(0).select(selected_body.as_view()),
        Atom::one()
    );
    assert_eq!(
        OrientationID(1).select(selected_body.as_view()),
        Atom::Zero,
        "late residue-key materialization selects exactly one complete production key",
    );
    Ok(())
}
