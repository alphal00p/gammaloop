use super::*;
use crate::network::NetworkState;
use crate::structure::representation::{Lorentz, Minkowski, RepName};
use crate::{chain, slot, trace};
use symbolica::{function, symbol};

fn mink4() -> Representation<Minkowski> {
    Minkowski {}.new_rep(4)
}

fn chain_factor(name: Symbol) -> Atom {
    FunctionBuilder::new(name)
        .add_arg(Atom::var(SPENSO_TAG.chain_in))
        .add_arg(Atom::var(SPENSO_TAG.chain_out))
        .finish()
}

fn chain_factor_with_external(name: Symbol, external: Atom) -> Atom {
    FunctionBuilder::new(name)
        .add_arg(external)
        .add_arg(Atom::var(SPENSO_TAG.chain_in))
        .add_arg(Atom::var(SPENSO_TAG.chain_out))
        .finish()
}

fn opaque_fast_settings() -> ParseSettings {
    ParseSettings {
        shorthand_parsing: ShorthandParsing::Opaque {
            inference: StructureInferenceMode::Fast,
        },
        ..Default::default()
    }
}

fn opaque_expanded_settings() -> ParseSettings {
    ParseSettings {
        shorthand_parsing: ShorthandParsing::Opaque {
            inference: StructureInferenceMode::Expanded,
        },
        ..Default::default()
    }
}

#[test]
fn parse_chain_as_opaque_tensor() {
    let rep = mink4();
    let expr = chain!(
        slot!(rep, i),
        slot!(rep, j),
        chain_factor(symbol!("f")),
        chain_factor(symbol!("g")),
    );

    let parsed = expr
        .parse_to_atom_net::<AbstractIndex>(&opaque_fast_settings())
        .unwrap();

    assert_eq!(parsed.state, NetworkState::SelfDualTensor);
    assert_eq!(parsed.graph.dangling_indices().len(), 2);
}

#[test]
fn parse_chain_as_opaque_tensor_with_expanded_structure() {
    let rep = mink4();
    let expr = chain!(
        slot!(rep, i),
        slot!(rep, j),
        chain_factor(symbol!("f")),
        chain_factor(symbol!("g")),
    );

    let parsed = expr
        .parse_to_atom_net::<AbstractIndex>(&opaque_expanded_settings())
        .unwrap();

    assert_eq!(parsed.state, NetworkState::SelfDualTensor);
    assert_eq!(parsed.graph.n_nodes(), 1);
    assert_eq!(parsed.graph.dangling_indices().len(), 2);
}

#[test]
fn parse_trace_materializes_closed_links() {
    let rep = mink4();
    let expr = trace!(&rep, chain_factor(symbol!("f")), chain_factor(symbol!("g")));

    let parsed = expr
        .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();

    assert!(
        parsed.state.is_scalar(),
        "expected compact Schoonschip vector product `{expr}` to parse as a scalar; got state {:?} with dangling indices {:?}",
        parsed.state,
        parsed.graph.dangling_indices()
    );
    assert!(parsed.graph.dangling_indices().is_empty());
}

#[test]
fn parse_trace_closes_links_and_keeps_external_indices() {
    let trace_rep = Lorentz {}.new_rep(4);
    let external_rep = mink4();
    let expr = trace!(
        &trace_rep,
        chain_factor_with_external(symbol!("f"), slot!(external_rep, a).to_atom()),
        chain_factor_with_external(symbol!("g"), slot!(external_rep, b).to_atom()),
        chain_factor_with_external(symbol!("h"), slot!(external_rep, c).to_atom()),
    );

    let parsed = expr
        .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();

    assert_eq!(parsed.state, NetworkState::SelfDualTensor);
    assert_eq!(parsed.graph.dangling_indices().len(), 3);
}

#[test]
fn parse_trace_as_opaque_tensor_with_fast_structure() {
    let trace_rep = Lorentz {}.new_rep(4);
    let external_rep = mink4();
    let expr = trace!(
        &trace_rep,
        chain_factor_with_external(symbol!("f"), slot!(external_rep, a).to_atom()),
        chain_factor_with_external(symbol!("g"), slot!(external_rep, b).to_atom()),
        chain_factor_with_external(symbol!("h"), slot!(external_rep, c).to_atom()),
    );

    let parsed = expr
        .parse_to_atom_net::<AbstractIndex>(&opaque_fast_settings())
        .unwrap();

    assert_eq!(parsed.state, NetworkState::SelfDualTensor);
    assert_eq!(parsed.graph.dangling_indices().len(), 3);
}

#[test]
fn parse_trace_as_opaque_tensor_with_expanded_structure() {
    let trace_rep = Lorentz {}.new_rep(4);
    let external_rep = mink4();
    let expr = trace!(
        &trace_rep,
        chain_factor_with_external(symbol!("f"), slot!(external_rep, a).to_atom()),
        chain_factor_with_external(symbol!("g"), slot!(external_rep, b).to_atom()),
        chain_factor_with_external(symbol!("h"), slot!(external_rep, c).to_atom()),
    );

    let parsed = expr
        .parse_to_atom_net::<AbstractIndex>(&opaque_expanded_settings())
        .unwrap();

    assert_eq!(parsed.state, NetworkState::SelfDualTensor);
    assert_eq!(parsed.graph.n_nodes(), 1);
    assert_eq!(parsed.graph.dangling_indices().len(), 3);
}

#[test]
fn parse_four_factor_trace_closes_links_and_keeps_external_indices() {
    let trace_rep = Lorentz {}.new_rep(4);
    let external_rep = mink4();
    let expr = trace!(
        &trace_rep,
        chain_factor_with_external(symbol!("f"), slot!(external_rep, a).to_atom()),
        chain_factor_with_external(symbol!("g"), slot!(external_rep, b).to_atom()),
        chain_factor_with_external(symbol!("h"), slot!(external_rep, c).to_atom()),
        chain_factor_with_external(symbol!("q"), slot!(external_rep, d).to_atom()),
    );

    let parsed = expr
        .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();

    assert_eq!(parsed.state, NetworkState::SelfDualTensor);
    assert_eq!(parsed.graph.dangling_indices().len(), 4);
}

#[test]
fn parse_single_factor_trace_closes_dualizable_links() {
    let rep = Lorentz {}.new_rep(4);
    let expr = trace!(&rep, chain_factor(symbol!("f")));

    let parsed = expr
        .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();

    assert!(parsed.state.is_scalar());
    assert!(parsed.graph.dangling_indices().is_empty());
}

#[test]
fn parse_empty_trace_unwraps_representation_dimension() {
    let rep = mink4();
    let expr = trace!(&rep);

    let parsed = expr
        .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();

    assert_eq!(parsed.state, NetworkState::PureScalar);
    assert!(parsed.graph.dangling_indices().is_empty());
}

#[test]
fn parse_empty_chain_as_endpoint_metric() {
    let rep = mink4();
    let expr = chain!(slot!(rep, i), slot!(rep, j));

    let parsed = expr
        .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();

    assert_eq!(parsed.state, NetworkState::SelfDualTensor);
    assert_eq!(parsed.graph.dangling_indices().len(), 2);
}

#[test]
fn opaque_schoonschip_vectors_stay_scalar() {
    let rep = mink4();
    let expr =
        function!(symbol!("p"), rep.to_symbolic([])) * function!(symbol!("q"), rep.to_symbolic([]));

    let parsed = expr
        .parse_to_atom_net::<AbstractIndex>(&opaque_fast_settings())
        .unwrap();

    assert_eq!(parsed.state, NetworkState::PureScalar);
    assert!(parsed.graph.dangling_indices().is_empty());
}

#[test]
fn parse_slot_metric_as_tensor() {
    let rep = mink4();
    let expr = ETS.metric(slot!(rep, i).to_atom(), slot!(rep, j).to_atom());

    let parsed = expr
        .parse_to_atom_net::<AbstractIndex>(&opaque_fast_settings())
        .unwrap();

    assert_eq!(parsed.state, NetworkState::SelfDualTensor);
    assert_eq!(parsed.graph.n_nodes(), 1);
    assert_eq!(parsed.graph.dangling_indices().len(), 2);
}

#[test]
fn three_argument_metric_inner_product_is_not_parser_syntax() {
    let rep = mink4();
    let expr = function!(
        ETS.metric,
        rep.to_symbolic([]),
        Atom::var(symbol!("p")),
        Atom::var(symbol!("q"))
    );

    let parsed = expr
        .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();

    assert_eq!(parsed.state, NetworkState::PureScalar);
    assert!(parsed.graph.dangling_indices().is_empty());
}

#[test]
fn parse_schoonschipped_metric_product() {
    let rep = mink4();
    let expr = function!(
        ETS.metric,
        function!(symbol!("p"), rep.to_symbolic([])),
        function!(symbol!("q"), rep.to_symbolic([]))
    );

    let parsed = expr
        .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();

    assert!(parsed.state.is_scalar());
    assert!(parsed.graph.dangling_indices().is_empty());
}

#[test]
fn parse_schoonschipped_metric_with_open_slot() {
    let rep = mink4();
    let expr = ETS.metric(
        slot!(rep, i).to_atom(),
        function!(symbol!("p"), rep.to_symbolic([])),
    );

    let parsed = expr
        .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();

    assert_eq!(parsed.state, NetworkState::SelfDualTensor);
    assert_eq!(parsed.graph.dangling_indices().len(), 1);
}

#[test]
fn parse_schoonschipped_higher_rank_tensor_keeps_open_slots() {
    let rep = mink4();
    let expr = function!(
        symbol!("F"),
        slot!(rep, i).to_atom(),
        function!(symbol!("p"), rep.to_symbolic([])),
        slot!(rep, j).to_atom()
    );

    let parsed = expr
        .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();

    assert_eq!(parsed.state, NetworkState::SelfDualTensor);
    assert!(parsed.graph.n_nodes() > 1);
    assert_eq!(parsed.graph.dangling_indices().len(), 2);
}

#[test]
fn opaque_schoonschipped_metric_product_stays_scalar() {
    let rep = mink4();
    let expr = function!(
        ETS.metric,
        function!(symbol!("p"), rep.to_symbolic([])),
        function!(symbol!("q"), rep.to_symbolic([]))
    );

    let parsed = expr
        .parse_to_atom_net::<AbstractIndex>(&opaque_fast_settings())
        .unwrap();

    assert_eq!(parsed.state, NetworkState::PureScalar);
    assert!(parsed.graph.dangling_indices().is_empty());
}

#[test]
fn parse_schoonschipped_metric_sum_product() {
    let rep = mink4();
    let k1 = function!(symbol!("k"), Atom::num(1), rep.to_symbolic([]));
    let k2 = function!(symbol!("k"), Atom::num(2), rep.to_symbolic([]));
    let p = function!(symbol!("p"), Atom::num(3), rep.to_symbolic([]));
    let expr = function!(ETS.metric, k1 + k2, p);

    let parsed = expr
        .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();

    assert!(parsed.state.is_scalar());
    assert!(parsed.graph.dangling_indices().is_empty());
}

#[test]
fn parse_schoonschipped_dot_product() {
    let rep = mink4();
    let expr = function!(
        SPENSO_TAG.dot,
        function!(symbol!("p"), rep.to_symbolic([])),
        function!(symbol!("q"), rep.to_symbolic([]))
    );

    let parsed = expr
        .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();

    assert!(parsed.state.is_scalar());
    assert!(parsed.graph.dangling_indices().is_empty());
}

#[test]
fn opaque_schoonschipped_dot_product_stays_scalar() {
    let rep = mink4();
    let expr = function!(
        SPENSO_TAG.dot,
        function!(symbol!("p"), rep.to_symbolic([])),
        function!(symbol!("q"), rep.to_symbolic([]))
    );

    let parsed = expr
        .parse_to_atom_net::<AbstractIndex>(&opaque_fast_settings())
        .unwrap();

    assert_eq!(parsed.state, NetworkState::PureScalar);
    assert!(parsed.graph.dangling_indices().is_empty());
}

#[test]
fn parse_linear_schoonschipped_dot_product() {
    let rep = mink4();
    let p = function!(symbol!("p"), rep.to_symbolic([]));
    let q = function!(symbol!("q"), rep.to_symbolic([]));
    let r = function!(symbol!("r"), rep.to_symbolic([]));
    let expr = function!(SPENSO_TAG.dot, p + q, r);

    let parsed = expr
        .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();

    assert!(parsed.state.is_scalar());
    assert!(parsed.graph.dangling_indices().is_empty());
}

#[test]
fn parse_chain_materializes_schoonschip_factor_argument() {
    let rep = mink4();
    let compact_vector = function!(symbol!("p"), rep.to_symbolic([]));
    let factor = FunctionBuilder::new(symbol!("f"))
        .add_arg(Atom::var(SPENSO_TAG.chain_in))
        .add_arg(Atom::var(SPENSO_TAG.chain_out))
        .add_arg(&compact_vector)
        .finish();
    let expr = chain!(slot!(rep, i), slot!(rep, j), factor);

    let parsed = expr
        .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();

    assert_eq!(parsed.state, NetworkState::SelfDualTensor);
    assert_eq!(parsed.graph.dangling_indices().len(), 2);
    assert_eq!(parsed.store.tensors.len(), 2);
}
