//! Spinneys and woods built from user-given selections, compared with classification.

use std::error::Error;

use feynkit_graph::{DiagramError, FeynmanDiagram, LoopMomentumBasis};
use feynkit_model::Model;
use feynkit_uv::{
    Counterterm, CountertermError, CountertermScheme, RenormalizationPrescriptionSettings, Spinney,
    UvGraph, Wood,
};
use linnet::half_edge::involution::{EdgeIndex, HedgePair};
use linnet::half_edge::subgraph::{InternalSubGraph, ModifySubSet, SuBitGraph};

/// The selection of the edges `ids`, with both halves of each.
fn lines(diagram: &FeynmanDiagram, ids: &[usize]) -> InternalSubGraph {
    let graph = diagram.hedge_graph();
    let mut filter: SuBitGraph = graph.empty_subgraph();
    for (pair, EdgeIndex(edge), _) in graph.iter_edges() {
        if let (true, HedgePair::Paired { source, sink }) = (ids.contains(&edge), pair) {
            filter.add(source);
            filter.add(sink);
        }
    }
    InternalSubGraph::cleaned_filter_optimist(filter, graph)
}

fn classified(diagram: &FeynmanDiagram) -> Vec<Spinney<LoopMomentumBasis>> {
    diagram
        .classify_spinneys(
            &diagram.hedge_graph().full_filter(),
            diagram.model(),
            &RenormalizationPrescriptionSettings::default(),
            diagram.loop_momentum_basis(),
        )
        .unwrap()
}

#[test]
fn convergent_factors_are_told_apart_and_never_contracted() -> Result<(), Box<dyn Error>> {
    // Quartic prism: two vertex-disjoint triangles, each convergent (4 - 6),
    // with a logarithmic aggregate degree.
    let prism = FeynmanDiagram::from_dot(
        Model::phi4(),
        r#"digraph prism { edge [particle=phi]; ext [style=invis];
            a -> b [id=0]; b -> c [id=1]; c -> a [id=2];
            d -> e [id=3]; e -> f [id=4]; f -> d [id=5];
            a -> d [id=6]; b -> e [id=7]; c -> f [id=8];
            ext -> a; ext -> b; ext -> c; d -> ext; e -> ext; f -> ext; }"#,
    )?;
    // Cubic convergent triangle with a divergent tadpole hanging off a bridge.
    let kite = FeynmanDiagram::from_dot(
        Model::phi3(),
        r#"digraph kite { edge [particle=phi]; ext [style=invis];
            a -> b [id=0]; b -> c [id=1]; c -> a [id=2]; a -> t [id=3]; t -> t [id=4];
            ext -> b; c -> ext; }"#,
    )?;
    for (diagram, region) in [
        (&prism, [0, 1, 2, 3, 4, 5].as_slice()),
        (&kite, &[0, 1, 2, 4]),
    ] {
        let basis = diagram.loop_momentum_basis();
        // Hosts build collective regions as atomic Taylor sources ...
        let collective = Spinney::new(lines(diagram, region), diagram, basis)?
            .expect("the union is a collective region");
        assert_eq!(collective.n_components(), 2);
        // ... but classification never retains them, and counterterms reject them.
        assert!(!collective.has_divergent_factors(diagram)?);
        assert!(
            classified(diagram)
                .iter()
                .all(|spinney| spinney.has_divergent_factors(diagram).unwrap())
        );
        assert!(matches!(
            Counterterm::new([collective], CountertermScheme::MsBar, diagram, basis),
            Err(CountertermError::<DiagramError>::ConvergentFactor { .. })
        ));
    }
    let kept = classified(&kite)
        .into_iter()
        .filter(|spinney| !spinney.is_empty())
        .collect::<Vec<_>>();
    assert_eq!(kept.len(), 1);
    assert_eq!(kept[0].edges(&kite), [EdgeIndex(4)]);
    Ok(())
}

#[test]
fn disjoint_products_complete_user_given_woods() -> Result<(), Box<dyn Error>> {
    for (model, dot) in [
        // Two tadpoles.
        (
            Model::phi3(),
            "digraph spectacles { edge [particle=phi]; a -> a; a -> b; b -> b; }",
        ),
        // Two quartic bubbles joined by a line: the full graph is convergent.
        (
            Model::phi4(),
            "digraph pair { edge [particle=phi]; ext [style=invis]; ext -> a; ext -> a; \
             a -> b; a -> b; b -> c; b -> ext; c -> d; c -> d; c -> ext; d -> ext; d -> ext; }",
        ),
    ] {
        let diagram = FeynmanDiagram::from_dot(model, dot)?;
        let basis = diagram.loop_momentum_basis();
        let all = classified(&diagram);
        let mut given = all
            .iter()
            .filter(|spinney| spinney.n_components() == 1)
            .cloned()
            .collect::<Vec<_>>();
        // The wood is literal: the given connected spinneys have no product.
        let literal = Wood::from_spinneys(given.clone(), &diagram, basis)?;
        assert_eq!(literal.unfold().n_forests(), 3);
        // Their products, as classification builds them, restore it.
        let products = Wood::disjoint_products(&given, &diagram, basis)?;
        assert_eq!(products.len(), 1);
        assert!(Wood::disjoint_products(&all, &diagram, basis)?.is_empty());
        given.extend(products);
        let closed = Wood::from_spinneys(given, &diagram, basis)?;
        let enumerated = Wood::from_spinneys(all, &diagram, basis)?;
        assert_eq!(closed.n_spinneys(), enumerated.n_spinneys());
        assert_eq!(closed.unfold().n_forests(), 4);
    }
    Ok(())
}
