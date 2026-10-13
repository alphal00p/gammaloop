//! Original numerator occurrences remain distinct under raised-pole sampling.

use super::*;
use crate::{
    graph::FeynmanGraph,
    initialisation::test_initialise,
    settings::global::OrientationPattern,
    uv::approx::{OrientationProjection, local_4d::on_shell_forest_tests::nested_quark_fixture},
};

#[test]
fn complete_os_parent_keeps_raised_quark_numerator_occurrences() -> Result<()> {
    test_initialise()?;
    let (mut graph, parent, child, root) = nested_quark_fixture(false)?;
    assert_eq!(
        graph.loop_momentum_basis.edge_signatures[EdgeIndex(2)],
        graph.loop_momentum_basis.edge_signatures[EdgeIndex(5)],
    );
    assert_eq!(
        graph[EdgeIndex(2)].mass_atom(),
        graph[EdgeIndex(5)].mass_atom()
    );
    let prepared = DirectOnShellPreparation::new(&graph, &parent, &child, parent.lmb())?;
    let owners = |atom: &Atom, tagged: bool| {
        let mut owners = BTreeSet::new();
        let _ = atom.replace_map(|view, _, _| {
            let AtomView::Fun(momentum) = view else {
                return;
            };
            if momentum.get_symbol() != GS.emr_mom || momentum.get_nargs() == 0 {
                return;
            }
            if tagged {
                if let Some((
                    owner,
                    crate::utils::symbols::UvMomentumProvenanceRole::PhysicalSourceFixed,
                    _,
                )) = GS.uv_momentum_provenance_data(momentum.get(0))
                {
                    owners.insert(owner);
                }
            } else if let Ok(owner) = usize::try_from(momentum.get(0)) {
                owners.insert(EdgeIndex(owner));
            }
        });
        owners
    };
    let expected = BTreeSet::from([EdgeIndex(2), EdgeIndex(5)]);
    for branch in &prepared.branches {
        assert!(owners(branch, true).is_superset(&expected));
        let original = prepared.original_source_numerator(branch);
        assert!(owners(&original, false).is_superset(&expected));
        let erased = prepared.scope_masses(&GS.erase_uv_momentum_provenance(branch));
        assert!(
            !owners(&erased, false).contains(&EdgeIndex(5)),
            "the negative control must expose the collapsed raised occurrence"
        );
        assert_eq!(
            original.replace_multiple(graph.uv_wrapped_replacement(
                parent.subgraph(),
                parent.lmb(),
                &[W_.x___],
            )),
            erased,
            "restoring source owners preserves the physical diagonal and projector order",
        );
    }

    let numerator = graph
        .numerator(&graph.full_filter(), root.subgraph())
        .get_single_atom()
        .expect("the full QCD numerator exists");
    let options = graph.denominator_only_cff_3d_expression_options();
    let surfaces = graph.get_esurface_canonization(&graph.loop_momentum_basis);
    let production =
        graph.generate_3d_expression_for_integrand(&[], &surfaces, &options, Some(&numerator))?;
    let pattern = OrientationPattern::default();
    let cut = CutSet::empty(graph.n_hedges());
    let orientation =
        OrientationProjection::exact_expression(&production, &options, &pattern, true);
    let source = Direct3dCts::root(&graph, Localizer::new(&cut, orientation))?.branches()?;
    let witness = GS.emr_mom(EdgeIndex(5), GS.cind(0)) - GS.emr_mom(EdgeIndex(2), GS.cind(0));
    let (body, parameters, rows) = source.prepare_numerator(
        orientation,
        &graph,
        &witness,
        DirectResidueBranches::numerator_scope().0,
    )?;
    let mut distinct_contact = false;
    for ((key, _), (row_key, arguments)) in source.iter_keys().zip(rows) {
        assert_eq!(key, &row_key);
        let mapped = body.replace_multiple(
            parameters
                .iter()
                .zip(arguments)
                .map(|(parameter, argument)| Replacement::new(parameter.to_pattern(), argument)),
        );
        let expected = key.map_numerator(orientation, &graph, &witness)?;
        assert!(
            (mapped - &expected).expand().is_zero(),
            "the shared numerator template must retain the complete source row"
        );
        let diagonal = expected
            .replace(GS.ose(EdgeIndex(5)))
            .with(GS.ose(EdgeIndex(2)));
        distinct_contact |= !diagonal.expand().is_zero();
    }
    assert!(
        distinct_contact,
        "the actual raised channel must sample e2/e5 independently even at equal physical energies"
    );
    Ok(())
}
