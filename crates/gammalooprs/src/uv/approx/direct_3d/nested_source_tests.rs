//! Exact source boundary for a completed OS component with a live spectator.

use super::*;
use crate::{
    initialisation::test_initialise,
    uv::approx::local_4d::on_shell_forest_tests::three_loop_quark_fixture,
};

#[test]
fn proper_two_loop_os_source_keeps_contact_rows_and_live_crown() -> Result<()> {
    test_initialise()?;
    let (mut graph, parent, child, root) = three_loop_quark_fixture()?;
    assert_eq!(graph.n_loops(&graph.full_filter()), 3);
    assert_eq!(graph.n_loops(parent.subgraph()), 2);
    assert_eq!(graph.n_loops(child.subgraph()), 1);
    assert_eq!(
        graph
            .loop_momentum_basis
            .loop_edges
            .iter()
            .copied()
            .collect::<Vec<_>>(),
        [EdgeIndex(2), EdgeIndex(3), EdgeIndex(4)],
    );
    assert_eq!(
        graph.loop_momentum_basis.edge_signatures[EdgeIndex(3)],
        graph.loop_momentum_basis.edge_signatures[EdgeIndex(6)],
    );
    let numerator = graph
        .numerator(parent.subgraph(), root.subgraph())?
        .state
        .expr;
    let source = FourDTerm::from_view(
        graph
            .denominator(parent.subgraph(), &graph.model, |_| -1)
            .as_view(),
    )?;
    let [source] = source.as_slice() else {
        return Err(eyre!("the original parent has one propagator product"));
    };
    assert_eq!(source.denominators.len(), 5);
    assert!(source.numerator.expression().is_one());
    for denominator in &source.denominators {
        assert_eq!(
            denominator.mass_squared,
            graph[denominator.source_edge]
                .mass_atom(&graph.model)
                .pow(2),
        );
        let momentum = function!(GS.emr_mom, usize::from(denominator.source_edge));
        assert!(denominator.momentum == momentum || denominator.momentum == -momentum);
    }
    let uv_edges = graph.paired_edges(parent.subgraph());
    let crown = graph.dummy_stripped_external_flows_of(parent.subgraph());
    let mut options = graph.denominator_only_cff_3d_expression_options();
    options.cff_generation_context = CffGenerationContext::EmbeddedCffFactor;
    let cutset = CutSet::empty(graph.n_hedges());
    // This exact source constructor certifies its original denominator
    // momenta, masses, powers and retained physical crown before generation.
    // Its capacity comes from the full physical factorized numerator.
    let (cff, _) = graph.cff_from_4d_denominators_in_uv_sub_lmb(
        &source.denominators,
        uv_edges,
        crown.included_iter(),
        parent.lmb(),
        ExactUvSubLmbFrame::RetainedPhysicalCrown,
        &cutset,
        &options,
        numerator.expression(),
        None,
        &[],
        &BTreeMap::new(),
    )?;
    let term = cff
        .terms
        .get(&CutCFFIndex::new_all_none())
        .expect("the source is an uncut component contour");
    let (body, parameters, rows) =
        term.prepare_exact_source_numerator(&DirectResidueBranches::numerator_scope().1)?;
    let certified_rows = term.certified_original_source_energy_rows(
        numerator.expression(),
        graph.underlying.n_edges(),
    )?;
    assert_eq!(rows.len(), term.orientations.len());
    let first = GS.emr_mom(EdgeIndex(3), GS.cind(0));
    let second = GS.emr_mom(EdgeIndex(6), GS.cind(0));
    let live_crown = GS.emr_mom(EdgeIndex(2), GS.cind(0));
    let mut distinct_contact = false;
    for ((orientation, row), energy_row) in term.orientations.iter().zip(rows).zip(certified_rows) {
        let mapped = body.replace_multiple(parameters.iter().zip(row).map(|(parameter, value)| {
            Replacement::new(parameter.to_pattern(), value.to_pattern())
        }));
        let independent =
            term.map_exact_source_atom(&orientation.orientation, numerator.expression())?;
        assert_eq!(
            mapped, independent,
            "the full factorized numerator must reconstruct from its component source rows"
        );
        assert_eq!(
            term.map_exact_source_atom(&orientation.orientation, &live_crown)?,
            live_crown,
            "the surrounding graph's p0 is outside this component contour",
        );
        let first = term.map_exact_source_atom(&orientation.orientation, &first)?;
        let second = term.map_exact_source_atom(&orientation.orientation, &second)?;
        assert_eq!(energy_row[3], first);
        assert_eq!(energy_row[6], second);
        assert_eq!(energy_row[2], live_crown);
        // Expansion here is limited to one linear energy-map difference;
        // no graph numerator or rational CFF sum is distributed.
        distinct_contact |= !(first - second).expand().is_zero();
    }
    assert!(
        distinct_contact,
        "the equal-momentum e3/e6 occurrences need separate raised/contact samples"
    );
    Ok(())
}
