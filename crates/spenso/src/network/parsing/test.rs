use super::*;
use crate::network::NetworkState;
use crate::network::library::symbolic::ETS;
use crate::structure::OrderedStructure;
use crate::structure::representation::{Lorentz, Minkowski, RepName};
use crate::{broadcast_symbol, chain, mink, p, q, slot, tensor, tensor_symbol, trace, vector};
use symbolica::{atom::FunctionBuilder, function, symbol};

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

fn schoonschip_only_settings() -> ParseSettings {
    ParseSettings {
        shorthand_parsing: ShorthandParsing::expand_schoonschip_only(),
        ..Default::default()
    }
}

fn schoonschip_without_inner_products_settings() -> ParseSettings {
    ParseSettings {
        shorthand_parsing: ShorthandParsing::Expand {
            schoonschip: SchoonschipExpansionMode {
                inner_products: false,
                expand_schoonship: true,
                expand_inside_chains: true,
            },
            trace: true,
            chain: true,
        },
        ..Default::default()
    }
}

#[test]
fn nested_chain_and_trace_scopes_are_rejected_before_parsing() {
    let rep = mink4();
    let factor = chain_factor(tensor_symbol!(nested_scope_factor));
    let inner_chain = SPENSO_TAG.chain(
        rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(31))
            .to_atom(),
        rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(37))
            .to_atom(),
        [factor.clone()],
    );
    let inner_trace = SPENSO_TAG.trace(rep.to_symbolic([]), [factor]);
    let outer_chain = |inner: Atom| {
        SPENSO_TAG.chain(
            rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(41))
                .to_atom(),
            rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(43))
                .to_atom(),
            [inner],
        )
    };
    let outer_trace = |inner: Atom| SPENSO_TAG.trace(rep.to_symbolic([]), [inner]);
    let wrapped_inner = FunctionBuilder::new(broadcast_symbol!("nested_scope_wrapper"))
        .add_arg(Atom::num(2) * inner_trace.clone())
        .finish();

    for expression in [
        outer_chain(inner_chain.clone()),
        outer_chain(inner_trace.clone()),
        outer_trace(inner_chain.clone()),
        outer_trace(inner_trace),
        outer_chain(wrapped_inner.clone()),
    ] {
        assert!(matches!(
            expression.parse_to_atom_net::<AbstractIndex>(&ParseSettings::default()),
            Err(TensorNetworkError::ChainNesting(_))
        ));
    }
    assert!(matches!(
        outer_chain(wrapped_inner).parse_to_atom_net::<AbstractIndex>(&opaque_fast_settings()),
        Err(TensorNetworkError::ChainNesting(_))
    ));
}

#[test]
fn sibling_and_singly_wrapped_chain_scopes_remain_valid() {
    let rep = mink4();
    let make_chain = |name, start, end| {
        SPENSO_TAG.chain(
            rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(start))
                .to_atom(),
            rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(end))
                .to_atom(),
            [chain_factor(SPENSO_TAG.tensor_symbol(name))],
        )
    };
    let first = make_chain("sibling_chain_first", 61, 67);
    let second = make_chain("sibling_chain_second", 71, 73);
    let wrapped = FunctionBuilder::new(broadcast_symbol!("single_chain_wrapper"))
        .add_arg(first.clone())
        .finish();

    assert!(
        (first * second)
            .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
            .is_ok()
    );
    assert!(
        wrapped
            .parse_to_atom_net::<AbstractIndex>(&opaque_fast_settings())
            .is_ok()
    );
}

#[test]
fn direct_structure_inference_rejects_nested_chain_scopes() {
    let rep = mink4();
    let inner = SPENSO_TAG.trace(
        rep.to_symbolic([]),
        [chain_factor(tensor_symbol!(nested_inference_factor))],
    );
    let expression = SPENSO_TAG.chain(
        rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(47))
            .to_atom(),
        rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(53))
            .to_atom(),
        [inner],
    );

    for mode in [
        StructureInferenceMode::Fast,
        StructureInferenceMode::Expanded,
    ] {
        assert!(matches!(
            OrderedStructure::<LibraryRep, AbstractIndex>::structure_from_atom(
                expression.as_view(),
                mode,
            ),
            Err(StructureError::ParsingError(message))
                if message.contains("cannot be nested")
        ));
    }
}

#[test]
fn parse_chain_as_opaque_tensor() {
    let rep = mink4();
    let expr = chain!(
        slot!(rep, i),
        slot!(rep, j),
        chain_factor(tensor_symbol!(parse_factor_f)),
        chain_factor(tensor_symbol!(parse_factor_g)),
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
        chain_factor(tensor_symbol!(parse_factor_f)),
        chain_factor(tensor_symbol!(parse_factor_g)),
    );

    let parsed = expr
        .parse_to_atom_net::<AbstractIndex>(&opaque_expanded_settings())
        .unwrap();

    assert_eq!(parsed.state, NetworkState::SelfDualTensor);
    assert_eq!(parsed.graph.n_nodes(), 1);
    assert_eq!(parsed.graph.dangling_indices().len(), 2);
}

#[test]
fn parse_chain_materializes_schoonschip_endpoints() {
    let rep = mink4();
    let expr = chain!(
        vector!(compact_start, Atom::num(1), rep.to_symbolic([])),
        vector!(compact_end, Atom::num(2), rep.to_symbolic([])),
        chain_factor(tensor_symbol!(parse_endpoint_chain_f)),
    );

    let parsed = expr
        .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();

    assert!(parsed.state.is_scalar());
    assert!(parsed.graph.dangling_indices().is_empty());
}

#[test]
fn parse_chain_materializes_mixed_schoonschip_endpoint() {
    let rep = mink4();
    let expr = chain!(
        slot!(rep, i),
        vector!(compact_end, Atom::num(2), rep.to_symbolic([])),
        chain_factor(tensor_symbol!(parse_mixed_endpoint_chain_f)),
    );

    let parsed = expr
        .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();

    assert_eq!(parsed.state, NetworkState::SelfDualTensor);
    assert_eq!(parsed.graph.dangling_indices().len(), 1);
}

#[test]
fn parse_scalar_prefixed_chain_keeps_scalar_outside_chain() {
    let rep = mink4();
    let expr = Atom::num(-2)
        * chain!(
            slot!(rep, i),
            slot!(rep, j),
            chain_factor(tensor_symbol!(parse_prefixed_chain_f)),
        );

    let parsed = expr
        .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();

    assert_eq!(parsed.state, NetworkState::SelfDualTensor);
    assert_eq!(parsed.graph.dangling_indices().len(), 2);
}

#[test]
fn parse_trace_materializes_closed_links() {
    let rep = mink4();
    let expr = trace!(
        &rep,
        chain_factor(tensor_symbol!(parse_factor_f)),
        chain_factor(tensor_symbol!(parse_factor_g))
    );

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
fn parse_scalar_prefixed_trace_keeps_scalar_outside_trace() {
    let rep = mink4();
    let expr = Atom::num(-2) * trace!(&rep, chain_factor(tensor_symbol!(parse_prefixed_trace_f)));

    let parsed = expr
        .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();

    assert!(parsed.state.is_scalar());
    assert!(parsed.graph.dangling_indices().is_empty());
}

#[test]
fn parse_trace_closes_links_and_keeps_external_indices() {
    let trace_rep = Lorentz {}.new_rep(4);
    let external_rep = mink4();
    let expr = trace!(
        &trace_rep,
        chain_factor_with_external(
            tensor_symbol!(parse_factor_f),
            slot!(external_rep, a).to_atom()
        ),
        chain_factor_with_external(
            tensor_symbol!(parse_factor_g),
            slot!(external_rep, b).to_atom()
        ),
        chain_factor_with_external(
            tensor_symbol!(parse_factor_h),
            slot!(external_rep, c).to_atom()
        ),
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
        chain_factor_with_external(
            tensor_symbol!(parse_factor_f),
            slot!(external_rep, a).to_atom()
        ),
        chain_factor_with_external(
            tensor_symbol!(parse_factor_g),
            slot!(external_rep, b).to_atom()
        ),
        chain_factor_with_external(
            tensor_symbol!(parse_factor_h),
            slot!(external_rep, c).to_atom()
        ),
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
        chain_factor_with_external(
            tensor_symbol!(parse_factor_f),
            slot!(external_rep, a).to_atom()
        ),
        chain_factor_with_external(
            tensor_symbol!(parse_factor_g),
            slot!(external_rep, b).to_atom()
        ),
        chain_factor_with_external(
            tensor_symbol!(parse_factor_h),
            slot!(external_rep, c).to_atom()
        ),
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
        chain_factor_with_external(
            tensor_symbol!(parse_factor_f),
            slot!(external_rep, a).to_atom()
        ),
        chain_factor_with_external(
            tensor_symbol!(parse_factor_g),
            slot!(external_rep, b).to_atom()
        ),
        chain_factor_with_external(
            tensor_symbol!(parse_factor_h),
            slot!(external_rep, c).to_atom()
        ),
        chain_factor_with_external(tensor_symbol!(trace_q), slot!(external_rep, d).to_atom()),
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
    let expr = trace!(&rep, chain_factor(tensor_symbol!(parse_factor_f)));

    let mut parsed = expr
        .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();

    assert!(
        parsed.state.is_scalar(),
        "got state {:?} with dangling indices {:?}",
        parsed.state,
        parsed.graph.dangling_indices()
    );
    assert!(parsed.graph.dangling_indices().is_empty());

    parsed.simple_execute();
    assert!(parsed.result_scalar().is_ok());
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
fn opaque_schoonschip_vectors_parse_as_tensor_scalar() {
    let expr = p!(q!(mink!(4)));

    let parsed = expr
        .parse_to_atom_net::<AbstractIndex>(&opaque_fast_settings())
        .unwrap();

    assert_eq!(parsed.state, NetworkState::Scalar);
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
        Atom::var(symbol!("compact_p")),
        Atom::var(symbol!("compact_q"))
    );

    let parsed = expr
        .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();

    assert_eq!(parsed.state, NetworkState::PureScalar);
    assert!(parsed.graph.dangling_indices().is_empty());
}

#[test]
fn three_argument_dot_is_invalid_parser_syntax() {
    let rep = mink4();
    let expr = function!(
        SPENSO_TAG.dot,
        rep.to_symbolic([]),
        Atom::var(symbol!("compact_p")),
        Atom::var(symbol!("compact_q"))
    );

    let err = expr
        .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap_err();

    assert!(matches!(err, TensorNetworkError::InvalidDotFunction(_)));
}

#[test]
fn malformed_broadcast_arities_are_invalid_parser_syntax() {
    let broadcast = symbol!(
        "spenso::parse_malformed_broadcast",
        tag = SPENSO_TAG.broadcast
    );
    let expressions = [
        FunctionBuilder::new(broadcast).finish(),
        FunctionBuilder::new(broadcast)
            .add_arg(Atom::num(1))
            .add_arg(Atom::num(2))
            .finish(),
    ];

    for expression in expressions {
        let err = expression
            .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap_err();

        assert!(matches!(err, TensorNetworkError::TooManyArgsFunction(_)));
    }
}

#[test]
fn parse_schoonschipped_metric_product() {
    let rep = mink4();
    let expr = function!(
        ETS.metric,
        vector!(compact_p, rep.to_symbolic([])),
        vector!(compact_q, rep.to_symbolic([]))
    );

    let parsed = expr
        .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();

    assert!(parsed.state.is_scalar());
    assert!(parsed.graph.dangling_indices().is_empty());
}

#[test]
fn inner_product_schoonschip_can_be_disabled() {
    let rep = mink4();
    let expr = function!(
        ETS.metric,
        vector!(compact_p, rep.to_symbolic([])),
        vector!(compact_q, rep.to_symbolic([]))
    );

    let parsed = expr
        .parse_to_atom_net::<AbstractIndex>(&schoonschip_without_inner_products_settings())
        .unwrap();

    assert!(parsed.state.is_scalar());
    assert!(parsed.graph.dangling_indices().is_empty());
    assert_eq!(parsed.graph.n_nodes(), 1);
}

#[test]
fn parse_schoonschipped_metric_with_open_slot() {
    let rep = mink4();
    let expr = ETS.metric(
        slot!(rep, i).to_atom(),
        vector!(compact_p, rep.to_symbolic([])),
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
    let expr = tensor!(
        F,
        slot!(rep, i).to_atom(),
        vector!(compact_p, rep.to_symbolic([])),
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
fn schoonschip_only_expands_compact_vectors() {
    let rep = mink4();
    let expr = tensor!(
        F,
        slot!(rep, i).to_atom(),
        vector!(compact_p, rep.to_symbolic([])),
        slot!(rep, j).to_atom()
    );

    let parsed = expr
        .parse_to_atom_net::<AbstractIndex>(&schoonschip_only_settings())
        .unwrap();

    assert_eq!(parsed.state, NetworkState::SelfDualTensor);
    assert!(parsed.graph.n_nodes() > 1);
    assert_eq!(parsed.graph.dangling_indices().len(), 2);
}

#[test]
fn representation_valued_scalar_metadata_is_not_a_compact_vector() {
    let rep = mink4();
    let metadata = function!(symbol!("parser_scalar_rep_metadata"), rep.to_symbolic([]));
    let expr = tensor!(parser_metadata_tensor, metadata, slot!(rep, i).to_atom());

    let parsed = expr
        .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();

    assert_eq!(parsed.state, NetworkState::SelfDualTensor);
    assert_eq!(parsed.graph.n_nodes(), 1);
    assert_eq!(parsed.graph.dangling_indices().len(), 1);
}

#[test]
fn opaque_schoonschipped_metric_product_parses_as_tensor_scalar() {
    let rep = mink4();
    let expr = function!(
        ETS.metric,
        vector!(compact_p, rep.to_symbolic([])),
        vector!(compact_q, rep.to_symbolic([]))
    );

    let parsed = expr
        .parse_to_atom_net::<AbstractIndex>(&opaque_fast_settings())
        .unwrap();

    assert_eq!(parsed.state, NetworkState::Scalar);
    assert!(parsed.graph.dangling_indices().is_empty());
}

#[test]
fn parse_schoonschipped_metric_sum_product() {
    let rep = mink4();
    let k1 = vector!(k, Atom::num(1), rep.to_symbolic([]));
    let k2 = vector!(k, Atom::num(2), rep.to_symbolic([]));
    let p = vector!(compact_p, Atom::num(3), rep.to_symbolic([]));
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
        vector!(compact_p, rep.to_symbolic([])),
        vector!(compact_q, rep.to_symbolic([]))
    );

    let parsed = expr
        .parse_to_atom_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();

    assert!(parsed.state.is_scalar());
    assert!(parsed.graph.dangling_indices().is_empty());
}

#[test]
fn opaque_schoonschipped_dot_product_parses_as_tensor_scalar() {
    let rep = mink4();
    let expr = function!(
        SPENSO_TAG.dot,
        vector!(compact_p, rep.to_symbolic([])),
        vector!(compact_q, rep.to_symbolic([]))
    );

    let parsed = expr
        .parse_to_atom_net::<AbstractIndex>(&opaque_fast_settings())
        .unwrap();

    assert_eq!(parsed.state, NetworkState::Scalar);
    assert!(parsed.graph.dangling_indices().is_empty());
}

#[test]
fn parse_linear_schoonschipped_dot_product() {
    let rep = mink4();
    let p = vector!(compact_p, rep.to_symbolic([]));
    let q = vector!(compact_q, rep.to_symbolic([]));
    let r = vector!(r, rep.to_symbolic([]));
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
    let compact_vector = vector!(compact_p, rep.to_symbolic([]));
    let factor = FunctionBuilder::new(tensor_symbol!(parse_factor_f))
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
