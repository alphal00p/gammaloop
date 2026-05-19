use std::{
    collections::{BTreeMap, BTreeSet},
    path::{Path, PathBuf},
    time::Instant,
};

use gammalooprs::{
    initialisation::test_initialise,
    numerator::{ParsingNet, aind::Aind},
    utils::{F, FUN_LIB, TENSORLIB},
};
use idenso::shorthands::{
    metric::MetricSimplifier,
    schoonschip::{Schoonschip, SchoonschipSettings},
};
use linnet::half_edge::{
    NodeIndex,
    involution::{EdgeData, Hedge, HedgePair, Orientation},
};
use spenso::{
    iterators::IteratableTensor,
    network::graph::{NetworkEdge, NetworkLeaf, NetworkNode, NetworkOp},
    network::library::symbolic::{ETS, ExplicitKey},
    network::{
        ExecutionResult, MinResultRank, Network, Sequential, SmallestDegree, Steps,
        parsing::{ParseSettings, ShadowedStructure, ShorthandParsing, StructureInferenceMode},
        store::NetworkStore,
    },
    structure::slot::{DummyAind, IsAbstractSlot},
    structure::{HasName, HasStructure, TensorStructure},
    tensors::{
        data::DataTensor,
        parametric::{AtomViewOrConcrete, MixedTensor, ParamOrConcrete},
    },
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Indeterminate, Symbol},
    parser::ParseSettings as SymbolicaParseSettings,
    wrap_input,
};

type ConcreteParsingNet = Network<
    NetworkStore<
        spenso::tensors::parametric::MixedTensor<
            F<f64>,
            spenso::network::parsing::ShadowedStructure<Aind>,
        >,
        Atom,
    >,
    ExplicitKey<Aind>,
    Symbol,
    Aind,
>;

type ActualTensor = MixedTensor<F<f64>, ShadowedStructure<Aind>>;

#[derive(Clone, Debug)]
struct BoundarySideSummary {
    nodes: usize,
    leaves: usize,
    sums: usize,
    products: usize,
    max_sum_children: usize,
    tensor_leaves: usize,
    scalar_leaves: usize,
    tensor_logical_entries: usize,
    tensor_flat_entries: usize,
    scalar_bytes: usize,
}

impl BoundarySideSummary {
    fn pressure(&self) -> usize {
        let sum_factor = self.max_sum_children.max(self.sums);
        sum_factor
            .saturating_mul(self.tensor_logical_entries.max(1))
            .saturating_add(self.scalar_bytes)
    }
}

#[derive(Clone, Debug)]
struct ProductBoundaryCandidate {
    product: NodeIndex,
    left_child: NodeIndex,
    right_child: NodeIndex,
    left: BoundarySideSummary,
    right: BoundarySideSummary,
    slot: spenso::structure::representation::LibrarySlot<Aind>,
    source: Hedge,
    sink: Hedge,
    rename_hedge: Hedge,
    rename_child: NodeIndex,
    orientation: Orientation,
    score: usize,
}

fn workspace_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .ancestors()
        .nth(2)
        .expect("gammalooprs crate should live two levels below the workspace root")
        .to_path_buf()
}

fn report(label: &str, start: Instant) {
    eprintln!("{label}: elapsed={:.3?}", start.elapsed());
}

fn strip_symbol_tags(input: &str) -> String {
    let mut stripped = String::with_capacity(input.len());
    let mut rest = input;

    while let Some(start) = rest.find("::{") {
        let (prefix, after_prefix) = rest.split_at(start);
        stripped.push_str(prefix);

        let Some(end) = after_prefix.find("}::") else {
            stripped.push_str(after_prefix);
            return stripped;
        };

        stripped.push_str("::");
        rest = &after_prefix[end + "}::".len()..];
    }

    stripped.push_str(rest);
    stripped
}

fn parse_root_input(filename: &str) -> Atom {
    let path = workspace_root().join(filename);
    let start = Instant::now();
    let input = std::fs::read_to_string(&path)
        .unwrap_or_else(|error| panic!("failed to read {}: {error}", path.display()));
    report(
        &format!("read {} bytes={}", path.display(), input.len()),
        start,
    );

    let start = Instant::now();
    let input = strip_symbol_tags(&input);
    report(&format!("strip_symbol_tags bytes={}", input.len()), start);

    let start = Instant::now();
    let expr = Atom::parse(wrap_input!(&input), SymbolicaParseSettings::default())
        .expect("failed to parse Symbolica expression");
    eprintln!(
        "atom stats: terms={} bytes={}",
        expr.nterms(),
        expr.as_view().get_byte_size()
    );
    report("symbolica_parse", start);

    expr
}

fn collect_horner(expr: &Atom) -> Atom {
    let start = Instant::now();
    let collected = expr.collect_horner::<Indeterminate>(None);
    eprintln!(
        "collect_horner stats: terms={} bytes={}",
        collected.nterms(),
        collected.as_view().get_byte_size()
    );
    report("collect_horner", start);
    collected
}

fn apply_sum46_prepass(label: &str, expr: Atom) -> Atom {
    let Ok(prepass) = std::env::var("SPENSO_SUM46_PREPASS") else {
        return expr;
    };
    if prepass.is_empty() || prepass == "none" {
        return expr;
    }

    let start = Instant::now();
    let simplified = match prepass.as_str() {
        "metrics" => expr.simplify_metrics(),
        "schoonschip" | "schoonschip_partial" => {
            expr.schoonschip_with_net::<false, Aind>(&SchoonschipSettings::partial())
        }
        "schoonschip_full" => {
            expr.schoonschip_with_net::<false, Aind>(&SchoonschipSettings::full())
        }
        "schoonschip_expand" => expr.schoonschip_with_net::<true, Aind>(
            &SchoonschipSettings::partial().with_expanded_contracted_sums(),
        ),
        other => panic!("unknown SPENSO_SUM46_PREPASS={other}"),
    };
    eprintln!(
        "{label} prepass={prepass} elapsed={:.3?} before_terms={} before_bytes={} after_terms={} after_bytes={}",
        start.elapsed(),
        expr.nterms(),
        expr.as_view().get_byte_size(),
        simplified.nterms(),
        simplified.as_view().get_byte_size(),
    );
    simplified
}

fn report_actual_net_stats(label: &str, net: &ParsingNet) {
    let mut concrete_tensors = 0usize;
    let mut param_tensors = 0usize;
    let mut concrete_entries = 0usize;
    let mut param_entries = 0usize;
    let mut param_dense = 0usize;
    let mut param_sparse = 0usize;
    let mut logical_entries = 0usize;
    let mut order_counts = BTreeMap::<usize, usize>::new();
    let mut name_counts = BTreeMap::<String, usize>::new();

    for tensor in &net.store.tensors {
        let name = tensor
            .name()
            .map(|name| name.to_string())
            .unwrap_or_else(|| "<anonymous>".to_owned());
        *name_counts.entry(name).or_default() += 1;

        match tensor {
            ParamOrConcrete::Concrete(tensor) => {
                concrete_tensors += 1;
                concrete_entries += tensor.iter_flat().count();
                *order_counts.entry(tensor.structure().order()).or_default() += 1;
                logical_entries += tensor.structure().size().unwrap_or(0);
            }
            ParamOrConcrete::Param(tensor) => {
                param_tensors += 1;
                param_entries += tensor.tensor.actual_size();
                *order_counts
                    .entry(tensor.tensor.structure().order())
                    .or_default() += 1;
                logical_entries += tensor.tensor.structure().size().unwrap_or(0);
                match &tensor.tensor {
                    DataTensor::Dense(_) => param_dense += 1,
                    DataTensor::Sparse(_) => param_sparse += 1,
                }
            }
        }
    }

    let mut top_names = name_counts.into_iter().collect::<Vec<_>>();
    top_names.sort_by(|left, right| right.1.cmp(&left.1).then_with(|| left.0.cmp(&right.0)));
    top_names.truncate(8);

    eprintln!(
        "{label} actual_net stats: graph_nodes={} graph_edges={} tensors={} concrete_tensors={} param_tensors={} scalars={} logical_entries={} concrete_entries={} param_entries={} param_dense={} param_sparse={} orders={:?} top_names={:?}",
        net.graph.graph.n_nodes(),
        net.graph.graph.n_edges(),
        net.store.tensors.len(),
        concrete_tensors,
        param_tensors,
        net.store.scalar.len(),
        logical_entries,
        concrete_entries,
        param_entries,
        param_dense,
        param_sparse,
        order_counts,
        top_names,
    );
}

fn accumulate_actual_tensor_reference(
    tensor: &ActualTensor,
    refs: &mut usize,
    logical_entries: &mut usize,
    flat_entries: &mut usize,
    nonzero_entries: &mut usize,
    orders: &mut BTreeMap<usize, usize>,
    names: &mut BTreeMap<String, usize>,
    nonzeros_by_name: &mut BTreeMap<String, usize>,
) {
    *refs += 1;
    *logical_entries += tensor.structure().size().unwrap_or(0);
    *flat_entries += tensor.iter_flat().count();
    *orders.entry(tensor.structure().order()).or_default() += 1;

    let name = tensor
        .name()
        .map(|name| name.to_string())
        .unwrap_or_else(|| "<anonymous>".to_owned());
    *names.entry(name.clone()).or_default() += 1;

    let tensor_nonzeros = tensor
        .iter_flat()
        .filter(|(_, value)| match value {
            AtomViewOrConcrete::Atom(atom) => !atom.is_zero(),
            AtomViewOrConcrete::Concrete(_) => true,
        })
        .count();
    *nonzero_entries += tensor_nonzeros;
    *nonzeros_by_name.entry(name).or_default() += tensor_nonzeros;
}

fn report_actual_net_leaf_stats(label: &str, net: &ParsingNet) {
    let mut leaf_nodes = 0usize;
    let mut sum_ops = 0usize;
    let mut product_ops = 0usize;
    let mut neg_ops = 0usize;
    let mut power_ops = 0usize;
    let mut function_ops = 0usize;
    let mut sum_children = 0usize;
    let mut product_children = 0usize;
    let mut max_sum_children = 0usize;
    let mut max_product_children = 0usize;

    let mut local_tensor_leaves = 0usize;
    let mut tensor_sum_leaves = 0usize;
    let mut tensor_sum_terms = 0usize;
    let mut tensor_term_leaves = 0usize;
    let mut tensor_term_scaled = 0usize;
    let mut tensor_term_sum_leaves = 0usize;
    let mut tensor_term_sum_terms = 0usize;
    let mut tensor_term_sum_scaled = 0usize;
    let mut library_leaves = 0usize;
    let mut scalar_leaves = 0usize;

    let mut tensor_refs = 0usize;
    let mut tensor_ref_logical_entries = 0usize;
    let mut tensor_ref_flat_entries = 0usize;
    let mut tensor_ref_nonzero_entries = 0usize;
    let mut library_logical_entries = 0usize;
    let mut scalar_bytes = 0usize;
    let mut tensor_ref_orders = BTreeMap::<usize, usize>::new();
    let mut tensor_ref_names = BTreeMap::<String, usize>::new();
    let mut tensor_ref_nonzeros_by_name = BTreeMap::<String, usize>::new();

    for (node_index, _neigh, node) in net.graph.graph.iter_nodes() {
        match node {
            NetworkNode::Leaf(leaf) => {
                leaf_nodes += 1;
                match leaf {
                    NetworkLeaf::LocalTensor(index) => {
                        local_tensor_leaves += 1;
                        accumulate_actual_tensor_reference(
                            &net.store.tensors[*index],
                            &mut tensor_refs,
                            &mut tensor_ref_logical_entries,
                            &mut tensor_ref_flat_entries,
                            &mut tensor_ref_nonzero_entries,
                            &mut tensor_ref_orders,
                            &mut tensor_ref_names,
                            &mut tensor_ref_nonzeros_by_name,
                        );
                    }
                    NetworkLeaf::TensorSum(indices) => {
                        tensor_sum_leaves += 1;
                        tensor_sum_terms += indices.len();
                        for index in indices {
                            accumulate_actual_tensor_reference(
                                &net.store.tensors[*index],
                                &mut tensor_refs,
                                &mut tensor_ref_logical_entries,
                                &mut tensor_ref_flat_entries,
                                &mut tensor_ref_nonzero_entries,
                                &mut tensor_ref_orders,
                                &mut tensor_ref_names,
                                &mut tensor_ref_nonzeros_by_name,
                            );
                        }
                    }
                    NetworkLeaf::TensorTerm(term) => {
                        tensor_term_leaves += 1;
                        if let Some(scalar) = term.scalar {
                            tensor_term_scaled += 1;
                            scalar_bytes += net.store.scalar[scalar].as_view().get_byte_size();
                        }
                        accumulate_actual_tensor_reference(
                            &net.store.tensors[term.tensor],
                            &mut tensor_refs,
                            &mut tensor_ref_logical_entries,
                            &mut tensor_ref_flat_entries,
                            &mut tensor_ref_nonzero_entries,
                            &mut tensor_ref_orders,
                            &mut tensor_ref_names,
                            &mut tensor_ref_nonzeros_by_name,
                        );
                    }
                    NetworkLeaf::TensorTermSum(terms) => {
                        tensor_term_sum_leaves += 1;
                        tensor_term_sum_terms += terms.len();
                        for term in terms {
                            if let Some(scalar) = term.scalar {
                                tensor_term_sum_scaled += 1;
                                scalar_bytes += net.store.scalar[scalar].as_view().get_byte_size();
                            }
                            accumulate_actual_tensor_reference(
                                &net.store.tensors[term.tensor],
                                &mut tensor_refs,
                                &mut tensor_ref_logical_entries,
                                &mut tensor_ref_flat_entries,
                                &mut tensor_ref_nonzero_entries,
                                &mut tensor_ref_orders,
                                &mut tensor_ref_names,
                                &mut tensor_ref_nonzeros_by_name,
                            );
                        }
                    }
                    NetworkLeaf::LibraryKey { key, .. } => {
                        library_leaves += 1;
                        library_logical_entries += key.structure.size().unwrap_or(0);
                    }
                    NetworkLeaf::Scalar(index) => {
                        scalar_leaves += 1;
                        scalar_bytes += net.store.scalar[*index].as_view().get_byte_size();
                    }
                }
            }
            NetworkNode::Op(op) => {
                let child_count = net.graph.cached_expr_children(node_index).len();
                match op {
                    NetworkOp::Sum => {
                        sum_ops += 1;
                        sum_children += child_count;
                        max_sum_children = max_sum_children.max(child_count);
                    }
                    NetworkOp::Product => {
                        product_ops += 1;
                        product_children += child_count;
                        max_product_children = max_product_children.max(child_count);
                    }
                    NetworkOp::Neg => neg_ops += 1,
                    NetworkOp::Power(_) => power_ops += 1,
                    NetworkOp::Function(_) => function_ops += 1,
                }
            }
        }
    }

    let mut top_tensor_ref_names = tensor_ref_names.into_iter().collect::<Vec<_>>();
    top_tensor_ref_names
        .sort_by(|left, right| right.1.cmp(&left.1).then_with(|| left.0.cmp(&right.0)));
    top_tensor_ref_names.truncate(8);
    let mut top_tensor_ref_nonzeros = tensor_ref_nonzeros_by_name.into_iter().collect::<Vec<_>>();
    top_tensor_ref_nonzeros
        .sort_by(|left, right| right.1.cmp(&left.1).then_with(|| left.0.cmp(&right.0)));
    top_tensor_ref_nonzeros.truncate(8);

    eprintln!(
        "{label} actual_net leaves: leaf_nodes={} sum_ops={} product_ops={} neg_ops={} power_ops={} function_ops={} sum_children={} product_children={} max_sum_children={} max_product_children={} local_tensor_leaves={} tensor_sum_leaves={} tensor_sum_terms={} tensor_term_leaves={} tensor_term_scaled={} tensor_term_sum_leaves={} tensor_term_sum_terms={} tensor_term_sum_scaled={} library_leaves={} scalar_leaves={} tensor_refs={} tensor_ref_logical_entries={} tensor_ref_flat_entries={} tensor_ref_nonzero_entries={} library_logical_entries={} scalar_bytes={} tensor_ref_orders={:?} top_tensor_ref_names={:?} top_tensor_ref_nonzeros={:?}",
        leaf_nodes,
        sum_ops,
        product_ops,
        neg_ops,
        power_ops,
        function_ops,
        sum_children,
        product_children,
        max_sum_children,
        max_product_children,
        local_tensor_leaves,
        tensor_sum_leaves,
        tensor_sum_terms,
        tensor_term_leaves,
        tensor_term_scaled,
        tensor_term_sum_leaves,
        tensor_term_sum_terms,
        tensor_term_sum_scaled,
        library_leaves,
        scalar_leaves,
        tensor_refs,
        tensor_ref_logical_entries,
        tensor_ref_flat_entries,
        tensor_ref_nonzero_entries,
        library_logical_entries,
        scalar_bytes,
        tensor_ref_orders,
        top_tensor_ref_names,
        top_tensor_ref_nonzeros,
    );
}

fn leaf_brief(net: &ParsingNet, node: NodeIndex) -> String {
    match &net.graph.graph[node] {
        NetworkNode::Leaf(NetworkLeaf::LocalTensor(index)) => net.store.tensors[*index]
            .name()
            .map(|name| name.to_string())
            .unwrap_or_else(|| format!("tensor#{index}")),
        NetworkNode::Leaf(NetworkLeaf::TensorSum(indices)) => {
            format!("TensorSum({})", indices.len())
        }
        NetworkNode::Leaf(NetworkLeaf::TensorTerm(term)) => match term.scalar {
            Some(_) => net.store.tensors[term.tensor]
                .name()
                .map(|name| format!("{name}*s"))
                .unwrap_or_else(|| format!("tensor#{}*s", term.tensor)),
            None => net.store.tensors[term.tensor]
                .name()
                .map(|name| name.to_string())
                .unwrap_or_else(|| format!("tensor#{}", term.tensor)),
        },
        NetworkNode::Leaf(NetworkLeaf::TensorTermSum(terms)) => {
            format!("TensorTermSum({})", terms.len())
        }
        NetworkNode::Leaf(NetworkLeaf::LibraryKey { key, .. }) => key
            .structure
            .name()
            .map(|name| format!("lib:{name}"))
            .unwrap_or_else(|| "lib:<anonymous>".to_owned()),
        NetworkNode::Leaf(NetworkLeaf::Scalar(_)) => "scalar".to_owned(),
        NetworkNode::Op(NetworkOp::Sum) => "sum".to_owned(),
        NetworkNode::Op(NetworkOp::Product) => "product".to_owned(),
        NetworkNode::Op(NetworkOp::Neg) => "neg".to_owned(),
        NetworkNode::Op(NetworkOp::Power(power)) => format!("pow({power})"),
        NetworkNode::Op(NetworkOp::Function(function)) => format!("fun:{function}"),
    }
}

fn report_product_interleaving(label: &str, net: &ParsingNet, limit: usize) {
    let mut products = net
        .graph
        .cached_expr_preorder_nodes()
        .into_iter()
        .filter_map(|node| match &net.graph.graph[node] {
            NetworkNode::Op(NetworkOp::Product) => {
                let children = net.graph.cached_expr_children(node);
                Some((node, children))
            }
            _ => None,
        })
        .collect::<Vec<_>>();

    products.sort_by(|left, right| {
        right
            .1
            .len()
            .cmp(&left.1.len())
            .then_with(|| left.0.cmp(&right.0))
    });

    for (rank, (node, children)) in products.into_iter().take(limit).enumerate() {
        let sequence = children
            .iter()
            .map(|child| leaf_brief(net, *child))
            .collect::<Vec<_>>();
        let mut counts = BTreeMap::<String, usize>::new();
        for item in &sequence {
            *counts.entry(item.clone()).or_default() += 1;
        }
        eprintln!(
            "{label} product_interleave rank={rank} node={node} children={} counts={:?} sequence={:?}",
            sequence.len(),
            counts,
            sequence,
        );
    }
}

fn report_concrete_net_stats(label: &str, net: &ConcreteParsingNet) {
    let mut concrete_tensors = 0usize;
    let mut param_tensors = 0usize;
    let mut concrete_entries = 0usize;
    let mut param_entries = 0usize;
    let mut order_counts = BTreeMap::<usize, usize>::new();
    let mut name_counts = BTreeMap::<String, usize>::new();

    for tensor in &net.store.tensors {
        let name = tensor
            .name()
            .map(|name| name.to_string())
            .unwrap_or_else(|| "<anonymous>".to_owned());
        *name_counts.entry(name).or_default() += 1;

        match tensor {
            ParamOrConcrete::Concrete(tensor) => {
                concrete_tensors += 1;
                concrete_entries += tensor.iter_flat().count();
                *order_counts.entry(tensor.structure().order()).or_default() += 1;
            }
            ParamOrConcrete::Param(tensor) => {
                param_tensors += 1;
                param_entries += tensor.tensor.actual_size();
                *order_counts
                    .entry(tensor.tensor.structure().order())
                    .or_default() += 1;
            }
        }
    }

    let mut top_names = name_counts.into_iter().collect::<Vec<_>>();
    top_names.sort_by(|left, right| right.1.cmp(&left.1).then_with(|| left.0.cmp(&right.0)));
    top_names.truncate(8);

    eprintln!(
        "{label} concrete_net stats: graph_nodes={} graph_edges={} tensors={} concrete_tensors={} param_tensors={} scalars={} concrete_entries={} param_entries={} orders={:?} top_names={:?}",
        net.graph.graph.n_nodes(),
        net.graph.graph.n_edges(),
        net.store.tensors.len(),
        concrete_tensors,
        param_tensors,
        net.store.scalar.len(),
        concrete_entries,
        param_entries,
        order_counts,
        top_names,
    );
}

fn collect_expr_subtree_nodes(net: &ParsingNet, root: NodeIndex) -> BTreeSet<NodeIndex> {
    let mut seen = BTreeSet::new();
    let mut stack = vec![root];

    while let Some(node) = stack.pop() {
        if !seen.insert(node) {
            continue;
        }
        stack.extend(net.graph.cached_expr_children(node));
    }

    seen
}

fn summarize_expr_subtree(net: &ParsingNet, root: NodeIndex) -> BoundarySideSummary {
    let nodes = collect_expr_subtree_nodes(net, root);
    let mut summary = BoundarySideSummary {
        nodes: nodes.len(),
        leaves: 0,
        sums: 0,
        products: 0,
        max_sum_children: 0,
        tensor_leaves: 0,
        scalar_leaves: 0,
        tensor_logical_entries: 0,
        tensor_flat_entries: 0,
        scalar_bytes: 0,
    };

    for node in nodes {
        match &net.graph.graph[node] {
            NetworkNode::Leaf(leaf) => {
                summary.leaves += 1;
                match leaf {
                    NetworkLeaf::LocalTensor(index) => {
                        summary.tensor_leaves += 1;
                        let tensor = &net.store.tensors[*index];
                        summary.tensor_logical_entries += tensor.structure().size().unwrap_or(0);
                        summary.tensor_flat_entries += tensor.iter_flat().count();
                    }
                    NetworkLeaf::TensorSum(indices) => {
                        summary.tensor_leaves += indices.len();
                        for index in indices {
                            let tensor = &net.store.tensors[*index];
                            summary.tensor_logical_entries +=
                                tensor.structure().size().unwrap_or(0);
                            summary.tensor_flat_entries += tensor.iter_flat().count();
                        }
                    }
                    NetworkLeaf::TensorTerm(term) => {
                        summary.tensor_leaves += 1;
                        let tensor = &net.store.tensors[term.tensor];
                        summary.tensor_logical_entries += tensor.structure().size().unwrap_or(0);
                        summary.tensor_flat_entries += tensor.iter_flat().count();
                        if let Some(scalar) = term.scalar {
                            summary.scalar_bytes +=
                                net.store.scalar[scalar].as_view().get_byte_size();
                        }
                    }
                    NetworkLeaf::TensorTermSum(terms) => {
                        summary.tensor_leaves += terms.len();
                        for term in terms {
                            let tensor = &net.store.tensors[term.tensor];
                            summary.tensor_logical_entries +=
                                tensor.structure().size().unwrap_or(0);
                            summary.tensor_flat_entries += tensor.iter_flat().count();
                            if let Some(scalar) = term.scalar {
                                summary.scalar_bytes +=
                                    net.store.scalar[scalar].as_view().get_byte_size();
                            }
                        }
                    }
                    NetworkLeaf::LibraryKey { key, .. } => {
                        summary.tensor_leaves += 1;
                        summary.tensor_logical_entries += key.structure.size().unwrap_or(0);
                    }
                    NetworkLeaf::Scalar(index) => {
                        summary.scalar_leaves += 1;
                        summary.scalar_bytes += net.store.scalar[*index].as_view().get_byte_size();
                    }
                }
            }
            NetworkNode::Op(op) => match op {
                NetworkOp::Sum => {
                    summary.sums += 1;
                    summary.max_sum_children = summary
                        .max_sum_children
                        .max(net.graph.cached_expr_children(node).len());
                }
                NetworkOp::Product => summary.products += 1,
                NetworkOp::Neg | NetworkOp::Power(_) | NetworkOp::Function(_) => {}
            },
        }
    }

    summary
}

fn graph_product_boundary_candidates(net: &ParsingNet) -> Vec<ProductBoundaryCandidate> {
    let mut candidates = Vec::new();

    for product in net.graph.cached_expr_preorder_nodes() {
        if !matches!(
            &net.graph.graph[product],
            NetworkNode::Op(NetworkOp::Product)
        ) {
            continue;
        }

        let children = net.graph.cached_expr_children(product);
        if children.len() < 2 {
            continue;
        }

        let child_nodes = children
            .iter()
            .map(|child| (*child, collect_expr_subtree_nodes(net, *child)))
            .collect::<Vec<_>>();
        let summaries = children
            .iter()
            .map(|child| (*child, summarize_expr_subtree(net, *child)))
            .collect::<BTreeMap<_, _>>();
        let mut node_to_child = BTreeMap::new();
        for (child_index, (_child, nodes)) in child_nodes.iter().enumerate() {
            for node in nodes {
                node_to_child.insert(*node, child_index);
            }
        }

        for (pair, _edge_index, edge_data) in net.graph.graph.iter_edges() {
            let HedgePair::Paired { source, sink } = pair else {
                continue;
            };
            let NetworkEdge::Slot(slot) = edge_data.data else {
                continue;
            };

            let source_node = net.graph.graph.node_id(source);
            let sink_node = net.graph.graph.node_id(sink);
            let (Some(source_child), Some(sink_child)) = (
                node_to_child.get(&source_node).copied(),
                node_to_child.get(&sink_node).copied(),
            ) else {
                continue;
            };
            if source_child == sink_child {
                continue;
            }

            let left_child = children[source_child];
            let right_child = children[sink_child];
            let left = summaries[&left_child].clone();
            let right = summaries[&right_child].clone();

            let (rename_hedge, rename_child) = if left.pressure() <= right.pressure() {
                (source, left_child)
            } else {
                (sink, right_child)
            };
            let expensive_pressure = left.pressure().max(right.pressure());
            let cheap_nodes = left.nodes.min(right.nodes).max(1);
            let score = expensive_pressure / cheap_nodes
                + left.max_sum_children.max(right.max_sum_children) * 1_000
                + left
                    .tensor_logical_entries
                    .max(right.tensor_logical_entries);

            candidates.push(ProductBoundaryCandidate {
                product,
                left_child,
                right_child,
                left,
                right,
                slot: *slot,
                source,
                sink,
                rename_hedge,
                rename_child,
                orientation: edge_data.orientation,
                score,
            });
        }
    }

    candidates.sort_by(|left, right| {
        right
            .score
            .cmp(&left.score)
            .then_with(|| left.product.cmp(&right.product))
            .then_with(|| left.source.cmp(&right.source))
    });
    candidates
}

fn log_product_boundary_candidates(label: &str, net: &ParsingNet, limit: usize) {
    let candidates = graph_product_boundary_candidates(net);
    eprintln!("{label} boundary_candidates total={}", candidates.len());
    for (rank, candidate) in candidates.iter().take(limit).enumerate() {
        eprintln!(
            "{label} boundary_candidate rank={rank} score={} product={} slot={} edge=({}->{}) rename_child={} rename_hedge={} left_child={} left_nodes={} left_leaves={} left_sums={} left_max_sum_children={} left_tensor_entries={} left_scalar_bytes={} right_child={} right_nodes={} right_leaves={} right_sums={} right_max_sum_children={} right_tensor_entries={} right_scalar_bytes={}",
            candidate.score,
            candidate.product,
            candidate.slot,
            candidate.source,
            candidate.sink,
            candidate.rename_child,
            candidate.rename_hedge,
            candidate.left_child,
            candidate.left.nodes,
            candidate.left.leaves,
            candidate.left.sums,
            candidate.left.max_sum_children,
            candidate.left.tensor_logical_entries,
            candidate.left.scalar_bytes,
            candidate.right_child,
            candidate.right.nodes,
            candidate.right.leaves,
            candidate.right.sums,
            candidate.right.max_sum_children,
            candidate.right.tensor_logical_entries,
            candidate.right.scalar_bytes,
        );
    }
}

fn reindex_structure_slot(
    structure: ShadowedStructure<Aind>,
    original: spenso::structure::representation::LibrarySlot<Aind>,
    fresh: Aind,
) -> ShadowedStructure<Aind> {
    let indices = structure
        .external_structure_iter()
        .map(|slot| {
            if slot.to_lib() == original {
                fresh
            } else {
                slot.aind()
            }
        })
        .collect::<Vec<_>>();

    structure
        .reindex(&indices)
        .expect("graph-derived staged reindex should preserve arity")
        .structure
}

fn apply_graph_boundary_disconnect(
    net: &mut ParsingNet,
    candidate: &ProductBoundaryCandidate,
) -> (
    spenso::structure::representation::LibrarySlot<Aind>,
    spenso::structure::representation::LibrarySlot<Aind>,
) {
    let original = candidate.slot;
    let mut fresh = original;
    fresh.set_aind(Aind::new_dummy());
    let fresh_aind = fresh.aind();
    let rename_nodes = collect_expr_subtree_nodes(net, candidate.rename_child);

    let mut changed_tensors = BTreeSet::new();
    for node in &rename_nodes {
        let NetworkNode::Leaf(leaf) = &net.graph.graph[*node] else {
            continue;
        };

        match leaf {
            NetworkLeaf::LocalTensor(index) => {
                changed_tensors.insert(*index);
            }
            NetworkLeaf::TensorSum(indices) => {
                changed_tensors.extend(indices.iter().copied());
            }
            NetworkLeaf::TensorTerm(term) => {
                changed_tensors.insert(term.tensor);
            }
            NetworkLeaf::TensorTermSum(terms) => {
                changed_tensors.extend(terms.iter().map(|term| term.tensor));
            }
            NetworkLeaf::LibraryKey { .. } | NetworkLeaf::Scalar(_) => {}
        }
    }

    for index in changed_tensors {
        let tensor = net.store.tensors[index].clone();
        net.store.tensors[index] = tensor.map_same_structure(|structure| {
            reindex_structure_slot(structure, original, fresh_aind)
        });
    }

    let edge_updates = net
        .graph
        .graph
        .iter_edges()
        .filter_map(|(pair, edge_index, edge_data)| {
            let NetworkEdge::Slot(slot) = edge_data.data else {
                return None;
            };
            if *slot != original {
                return None;
            }

            let touches_rename_side = match pair {
                HedgePair::Paired { source, sink } => {
                    rename_nodes.contains(&net.graph.graph.node_id(source))
                        && rename_nodes.contains(&net.graph.graph.node_id(sink))
                }
                HedgePair::Unpaired { hedge, .. } => {
                    rename_nodes.contains(&net.graph.graph.node_id(hedge))
                }
                HedgePair::Split { .. } => false,
            };
            touches_rename_side.then_some(edge_index)
        })
        .collect::<Vec<_>>();

    for edge_index in edge_updates {
        net.graph.graph[edge_index] = NetworkEdge::Slot(fresh);
    }

    net.graph
        .graph
        .split_edge(
            candidate.rename_hedge,
            EdgeData::new(NetworkEdge::Slot(fresh), candidate.orientation),
        )
        .expect("selected product boundary slot edge should be paired");

    (original, fresh)
}

fn reconnect_metric_expr(
    original: spenso::structure::representation::LibrarySlot<Aind>,
    fresh: spenso::structure::representation::LibrarySlot<Aind>,
) -> Atom {
    ETS.metric(original.to_atom(), fresh.to_atom())
}

fn execute_actual_net(label: &str, mut net: ParsingNet) -> Atom {
    let lib = TENSORLIB.read().unwrap();
    spenso::network::profile::reset();
    let start = Instant::now();
    net.execute::<Sequential, SmallestDegree, _, _, _>(&*lib, &*FUN_LIB)
        .unwrap_or_else(|error| panic!("{label} actual execution failed: {error}"));
    report(&format!("{label} actual_execute"), start);
    report_actual_net_stats(&format!("{label} after_execute"), &net);
    spenso::network::profile::report(&format!("{label} after_actual_execute"));

    let result = net
        .result_scalar()
        .unwrap_or_else(|error| panic!("{label} actual scalar result failed: {error}"));
    let result = match result {
        ExecutionResult::One => Atom::num(1),
        ExecutionResult::Zero => Atom::Zero,
        ExecutionResult::Val(value) => value.into_owned(),
    };
    eprintln!(
        "{label} actual_result stats: terms={} bytes={}",
        result.nterms(),
        result.as_view().get_byte_size()
    );
    result
}

fn execute_actual_net_min_result_rank(label: &str, mut net: ParsingNet) -> Atom {
    let lib = TENSORLIB.read().unwrap();
    spenso::network::profile::reset();
    let start = Instant::now();
    net.execute::<Sequential, MinResultRank, _, _, _>(&*lib, &*FUN_LIB)
        .unwrap_or_else(|error| panic!("{label} min-result-rank actual execution failed: {error}"));
    report(&format!("{label} min_result_rank_actual_execute"), start);
    report_actual_net_stats(&format!("{label} min_result_rank_after_execute"), &net);
    spenso::network::profile::report(&format!("{label} after_min_result_rank_actual_execute"));

    let result = net
        .result_scalar()
        .unwrap_or_else(|error| panic!("{label} min-result-rank scalar result failed: {error}"));
    let result = match result {
        ExecutionResult::One => Atom::num(1),
        ExecutionResult::Zero => Atom::Zero,
        ExecutionResult::Val(value) => value.into_owned(),
    };
    eprintln!(
        "{label} min_result_rank_actual_result stats: terms={} bytes={}",
        result.nterms(),
        result.as_view().get_byte_size()
    );
    result
}

fn execute_actual_net_min_result_rank_parallel(label: &str, mut net: ParsingNet) -> Atom {
    let lib = TENSORLIB.read().unwrap();
    spenso::network::profile::reset();
    let start = Instant::now();
    net.execute_parallel::<MinResultRank, _, _, _>(&*lib, &*FUN_LIB)
        .unwrap_or_else(|error| {
            panic!("{label} parallel min-result-rank actual execution failed: {error}")
        });
    report(
        &format!("{label} parallel_min_result_rank_actual_execute"),
        start,
    );
    report_actual_net_stats(
        &format!("{label} parallel_min_result_rank_after_execute"),
        &net,
    );
    spenso::network::profile::report(&format!(
        "{label} after_parallel_min_result_rank_actual_execute"
    ));

    let result = net.result_scalar().unwrap_or_else(|error| {
        panic!("{label} parallel min-result-rank scalar result failed: {error}")
    });
    let result = match result {
        ExecutionResult::One => Atom::num(1),
        ExecutionResult::Zero => Atom::Zero,
        ExecutionResult::Val(value) => value.into_owned(),
    };
    eprintln!(
        "{label} parallel_min_result_rank_actual_result stats: terms={} bytes={}",
        result.nterms(),
        result.as_view().get_byte_size()
    );
    result
}

fn execute_actual_net_min_result_rank_parallel_tensor(
    label: &str,
    mut net: ParsingNet,
) -> ActualTensor {
    let lib = TENSORLIB.read().unwrap();
    spenso::network::profile::reset();
    let start = Instant::now();
    net.execute_parallel::<MinResultRank, _, _, _>(&*lib, &*FUN_LIB)
        .unwrap_or_else(|error| {
            panic!("{label} parallel min-result-rank tensor execution failed: {error}")
        });
    report(
        &format!("{label} parallel_min_result_rank_tensor_execute"),
        start,
    );
    report_actual_net_stats(
        &format!("{label} parallel_min_result_rank_tensor_after_execute"),
        &net,
    );
    spenso::network::profile::report(&format!(
        "{label} after_parallel_min_result_rank_tensor_execute"
    ));

    let result = net.result_tensor(&*lib).unwrap_or_else(|error| {
        panic!("{label} parallel min-result-rank tensor result failed: {error}")
    });
    let tensor = match result {
        ExecutionResult::Val(value) => value.into_owned(),
        ExecutionResult::One => panic!("{label} produced tensor result One"),
        ExecutionResult::Zero => panic!("{label} produced tensor result Zero"),
    };
    eprintln!(
        "{label} tensor_result stats: order={} logical_entries={} flat_entries={}",
        tensor.structure().order(),
        tensor.structure().size().unwrap_or(0),
        tensor.iter_flat().count()
    );
    tensor
}

fn execute_actual_net_steps(label: &str, mut net: ParsingNet, steps: usize) {
    let lib = TENSORLIB.read().unwrap();
    for step in 0..steps {
        spenso::network::profile::reset();
        let start = Instant::now();
        net.execute::<Steps<1>, SmallestDegree, _, _, _>(&*lib, &*FUN_LIB)
            .unwrap_or_else(|error| panic!("{label} actual execution step {step} failed: {error}"));
        report(&format!("{label} actual_execute_step_{step}"), start);
        report_actual_net_stats(&format!("{label} after_step_{step}"), &net);
        spenso::network::profile::report(&format!("{label} after_actual_execute_step_{step}"));
    }
}

fn parse_actual_net(label: &str, expr: &Atom) -> ParsingNet {
    parse_actual_net_with_settings(label, expr, actual_parse_settings())
}

fn parse_actual_net_with_boundary_factor(label: &str, expr: &Atom) -> ParsingNet {
    let mut settings = actual_parse_settings();
    settings.factor_add_boundaries = true;
    parse_actual_net_with_settings(label, expr, settings)
}

fn actual_parse_settings() -> ParseSettings {
    ParseSettings {
        shorthand_parsing: ShorthandParsing::Opaque {
            inference: StructureInferenceMode::Fast,
        },
        ..Default::default()
    }
}

fn parse_actual_net_with_settings(label: &str, expr: &Atom, settings: ParseSettings) -> ParsingNet {
    spenso::network::profile::reset();
    let lib = TENSORLIB.read().unwrap();
    let start = Instant::now();
    let net = ParsingNet::try_from_view(expr.as_view(), &*lib, &settings)
        .unwrap_or_else(|error| panic!("{label} actual parse failed: {error}"));
    report(&format!("{label} actual_network_parse"), start);
    report_actual_net_stats(label, &net);
    spenso::network::profile::report(&format!("{label} after_actual_network_parse"));
    net
}

fn parse_concrete_hep_net(label: &str, expr: &Atom) -> ConcreteParsingNet {
    spenso::network::profile::reset();
    let settings = ParseSettings {
        shorthand_parsing: ShorthandParsing::Opaque {
            inference: StructureInferenceMode::Fast,
        },
        ..Default::default()
    };
    let lib = spenso_hep_lib::hep_lib(F(1.0), F(0.0));
    let start = Instant::now();
    let net = ConcreteParsingNet::try_from_view(expr.as_view(), &lib, &settings)
        .unwrap_or_else(|error| panic!("{label} concrete HEP parse failed: {error}"));
    report(&format!("{label} concrete_hep_network_parse"), start);
    report_concrete_net_stats(label, &net);
    spenso::network::profile::report(&format!("{label} after_concrete_hep_network_parse"));
    net
}

fn execute_concrete_hep_net(label: &str, mut net: ConcreteParsingNet) {
    let lib = spenso_hep_lib::hep_lib(F(1.0), F(0.0));
    let start = Instant::now();
    net.execute::<Sequential, SmallestDegree, _, _, _>(&lib, &*FUN_LIB)
        .unwrap_or_else(|error| panic!("{label} concrete HEP execution failed: {error}"));
    report(&format!("{label} concrete_hep_execute"), start);
    report_concrete_net_stats(&format!("{label} after_concrete_hep_execute"), &net);
}

fn parse_inline_expression(input: &str) -> Atom {
    Atom::parse(wrap_input!(input), SymbolicaParseSettings::default())
        .expect("inline diagnostic expression should parse")
}

fn collect_adds_with_arity(view: AtomView<'_>, arity: usize, out: &mut Vec<Atom>) {
    match view {
        AtomView::Add(add) => {
            if add.get_nargs() == arity {
                out.push(add.as_view().to_owned());
            }
            for term in add.iter() {
                collect_adds_with_arity(term, arity, out);
            }
        }
        AtomView::Mul(mul) => {
            for factor in mul.iter() {
                collect_adds_with_arity(factor, arity, out);
            }
        }
        AtomView::Fun(fun) => {
            for arg in fun.iter() {
                collect_adds_with_arity(arg, arity, out);
            }
        }
        AtomView::Pow(pow) => {
            let (base, exp) = pow.get_base_exp();
            collect_adds_with_arity(base, arity, out);
            collect_adds_with_arity(exp, arity, out);
        }
        AtomView::Num(_) | AtomView::Var(_) => {}
    }
}

fn count_function_occurrences(view: AtomView<'_>, name: &str) -> usize {
    match view {
        AtomView::Fun(fun) => {
            usize::from(fun.get_symbol().get_stripped_name() == name)
                + fun
                    .iter()
                    .map(|arg| count_function_occurrences(arg, name))
                    .sum::<usize>()
        }
        AtomView::Add(add) => add
            .iter()
            .map(|term| count_function_occurrences(term, name))
            .sum(),
        AtomView::Mul(mul) => mul
            .iter()
            .map(|factor| count_function_occurrences(factor, name))
            .sum(),
        AtomView::Pow(pow) => {
            let (base, exp) = pow.get_base_exp();
            count_function_occurrences(base, name) + count_function_occurrences(exp, name)
        }
        AtomView::Num(_) | AtomView::Var(_) => 0,
    }
}

fn add_terms(expr: &Atom) -> Vec<Atom> {
    match expr.as_view() {
        AtomView::Add(add) => add.iter().map(|term| term.to_owned()).collect(),
        _ => panic!("expected an add expression"),
    }
}

fn execute_actual_net_min_result_rank_parallel_summary(
    label: &str,
    mut net: ParsingNet,
) -> std::time::Duration {
    let lib = TENSORLIB.read().unwrap();
    spenso::network::profile::reset();
    let start = Instant::now();
    net.execute_parallel::<MinResultRank, _, _, _>(&*lib, &*FUN_LIB)
        .unwrap_or_else(|error| {
            panic!("{label} parallel min-result-rank actual execution failed: {error}")
        });
    let elapsed = start.elapsed();

    let result_start = Instant::now();
    if let Ok(result) = net.result_scalar() {
        let result = match result {
            ExecutionResult::One => Atom::num(1),
            ExecutionResult::Zero => Atom::Zero,
            ExecutionResult::Val(value) => value.into_owned(),
        };
        let result_elapsed = result_start.elapsed();
        eprintln!(
            "{label} summary kind=scalar execute={elapsed:.3?} result={result_elapsed:.3?} result_terms={} result_bytes={}",
            result.nterms(),
            result.as_view().get_byte_size()
        );
        return elapsed;
    }

    let result = net
        .result_tensor(&*lib)
        .unwrap_or_else(|error| panic!("{label} tensor result failed after execution: {error}"));
    let tensor = match result {
        ExecutionResult::Val(value) => value.into_owned(),
        ExecutionResult::One => {
            let result_elapsed = result_start.elapsed();
            eprintln!(
                "{label} summary kind=tensor_one execute={elapsed:.3?} result={result_elapsed:.3?}"
            );
            return elapsed;
        }
        ExecutionResult::Zero => {
            let result_elapsed = result_start.elapsed();
            eprintln!(
                "{label} summary kind=tensor_zero execute={elapsed:.3?} result={result_elapsed:.3?}"
            );
            return elapsed;
        }
    };

    if tensor.structure().order() == 0 {
        let result: Atom = tensor
            .scalar()
            .unwrap_or_else(|| panic!("{label} rank-0 tensor result did not expose a scalar"))
            .into();
        let result_elapsed = result_start.elapsed();
        eprintln!(
            "{label} summary kind=scalar_from_tensor execute={elapsed:.3?} result={result_elapsed:.3?} result_terms={} result_bytes={}",
            result.nterms(),
            result.as_view().get_byte_size()
        );
        return elapsed;
    }

    let result_elapsed = result_start.elapsed();
    eprintln!(
        "{label} summary kind=tensor execute={elapsed:.3?} result={result_elapsed:.3?} order={} logical_entries={} flat_entries={}",
        tensor.structure().order(),
        tensor.structure().size().unwrap_or(0),
        tensor.iter_flat().count()
    );
    elapsed
}

fn parse_actual_net_quiet(label: &str, expr: &Atom) -> ParsingNet {
    let lib = TENSORLIB.read().unwrap();
    let start = Instant::now();
    let net = ParsingNet::try_from_view(expr.as_view(), &*lib, &actual_parse_settings())
        .unwrap_or_else(|error| panic!("{label} actual parse failed: {error}"));
    eprintln!(
        "{label} quiet_parse elapsed={:.3?} graph_nodes={} tensors={} scalars={}",
        start.elapsed(),
        net.graph.graph.n_nodes(),
        net.store.tensors.len(),
        net.store.scalar.len()
    );
    net
}

fn gamma_ladder_order_mwe_expr() -> Atom {
    parse_inline_expression(
        r#"
        spenso::g(spenso::bis(4,4),spenso::bis(4,5))
        *spenso::gamma(spenso::bis(4,4),spenso::bis(4,23),spenso::mink(4,0))
        *spenso::g(spenso::bis(4,22),spenso::bis(4,23))
        *ϵ(0,spenso::mink(4,0))
        *spenso::gamma(spenso::bis(4,5),spenso::bis(4,6),spenso::mink(4,1))
        *spenso::g(spenso::bis(4,6),spenso::bis(4,7))
        *ϵ(1,spenso::mink(4,1))
        *spenso::gamma(spenso::bis(4,7),spenso::bis(4,12),spenso::mink(4,14))
        *spenso::g(spenso::bis(4,12),spenso::bis(4,13))
        *spenso::gamma(spenso::bis(4,13),spenso::bis(4,16),spenso::mink(4,18))
        *spenso::g(spenso::bis(4,16),spenso::bis(4,17))
        *spenso::gamma(spenso::bis(4,8),spenso::bis(4,17),spenso::mink(4,2))
        *spenso::g(spenso::bis(4,8),spenso::bis(4,9))
        *ϵbar(2,spenso::mink(4,2))
        *spenso::gamma(spenso::bis(4,9),spenso::bis(4,10),spenso::mink(4,3))
        *spenso::g(spenso::bis(4,10),spenso::bis(4,11))
        *ϵbar(3,spenso::mink(4,3))
        *spenso::g(spenso::bis(4,20),spenso::bis(4,21))
        *spenso::gamma(spenso::bis(4,11),spenso::bis(4,20),spenso::mink(4,18))
        *spenso::gamma(spenso::bis(4,21),spenso::bis(4,22),spenso::mink(4,14))
        "#,
    )
}

fn scalar_poly(seed: usize, terms: usize) -> String {
    (0..terms)
        .map(|term| format!("x{seed}_{term}*y{seed}_{term}"))
        .collect::<Vec<_>>()
        .join("+")
}

fn late_tensor_sum_mwe_expanded_sum(leaves: usize, scalar_terms: usize) -> String {
    let tensor = "spenso::g(spenso::mink(4,mu),spenso::mink(4,nu))*Q(0,spenso::mink(4,mu))";
    (0..leaves)
        .map(|leaf| format!("({})*({tensor})", scalar_poly(leaf, scalar_terms)))
        .collect::<Vec<_>>()
        .join("+")
}

fn late_tensor_sum_mwe_expr(leaves: usize, scalar_terms: usize) -> Atom {
    let expanded_sum = late_tensor_sum_mwe_expanded_sum(leaves, scalar_terms);

    parse_inline_expression(&format!("({expanded_sum})*Q(99,spenso::mink(4,nu))"))
}

fn late_tensor_sum_mwe_disconnected_expr(leaves: usize, scalar_terms: usize) -> Atom {
    let expanded_sum = late_tensor_sum_mwe_expanded_sum(leaves, scalar_terms);

    parse_inline_expression(&format!("({expanded_sum})*Q(99,spenso::mink(4,rho))"))
}

fn late_tensor_sum_mwe_reconnect_metric_expr() -> Atom {
    parse_inline_expression("spenso::g(spenso::mink(4,nu),spenso::mink(4,rho))")
}

fn factored_tensor_sum_mwe_expr(leaves: usize, scalar_terms: usize) -> Atom {
    let scalar_sum = (0..leaves)
        .map(|leaf| format!("({})", scalar_poly(leaf, scalar_terms)))
        .collect::<Vec<_>>()
        .join("+");

    parse_inline_expression(&format!(
        "(({scalar_sum})
          *spenso::g(spenso::mink(4,mu),spenso::mink(4,nu))
          *Q(0,spenso::mink(4,mu)))
         *Q(99,spenso::mink(4,nu))"
    ))
}

#[test]
#[ignore = "diagnostic timing for GammaLoop's concrete/mixed tensor evaluator path"]
fn symbolica_expression_input_actual_network_execute() {
    test_initialise().expect("GammaLoop initialization should succeed");
    let expr = parse_root_input("symbolica_expression.txt");
    let net = parse_actual_net("raw", &expr);
    let _ = execute_actual_net("raw", net);
}

#[test]
#[ignore = "diagnostic timing for GammaLoop's concrete/mixed tensor evaluator path with MinResultRank"]
fn symbolica_expression_input_actual_network_min_result_rank_execute() {
    test_initialise().expect("GammaLoop initialization should succeed");
    let expr = parse_root_input("symbolica_expression.txt");
    let net = parse_actual_net("raw_min_result_rank", &expr);
    let _ = execute_actual_net_min_result_rank("raw_min_result_rank", net);
}

#[test]
#[ignore = "diagnostic timing for GammaLoop's concrete/mixed tensor evaluator path with parallel MinResultRank"]
fn symbolica_expression_input_actual_network_parallel_min_result_rank_execute() {
    test_initialise().expect("GammaLoop initialization should succeed");
    let expr = parse_root_input("symbolica_expression.txt");
    let net = parse_actual_net("raw_parallel_min_result_rank", &expr);
    let _ = execute_actual_net_min_result_rank_parallel("raw_parallel_min_result_rank", net);
}

#[test]
#[ignore = "diagnostic timing for the larger spenso_eval_input_0.txt parse path"]
fn spenso_eval_input_0_actual_network_parse() {
    test_initialise().expect("GammaLoop initialization should succeed");
    let expr = parse_root_input("spenso_eval_input_0.txt");
    let _net = parse_actual_net("spenso_eval_input_0", &expr);
}

#[test]
#[ignore = "diagnostic timing for the larger spenso_eval_input_0.txt parallel MinResultRank path"]
fn spenso_eval_input_0_actual_network_parallel_min_result_rank_execute() {
    test_initialise().expect("GammaLoop initialization should succeed");
    let expr = parse_root_input("spenso_eval_input_0.txt");
    let net = parse_actual_net("spenso_eval_input_0_parallel_min_result_rank", &expr);
    let _ = execute_actual_net_min_result_rank_parallel(
        "spenso_eval_input_0_parallel_min_result_rank",
        net,
    );
}

#[test]
#[ignore = "diagnostic timing for each branch of the largest 46-way sum in spenso_eval_input_0.txt"]
fn spenso_eval_input_0_sum46_individual_term_execute() {
    test_initialise().expect("GammaLoop initialization should succeed");
    let expr = parse_root_input("spenso_eval_input_0.txt");

    let mut candidates = Vec::new();
    collect_adds_with_arity(expr.as_view(), 46, &mut candidates);
    assert!(
        !candidates.is_empty(),
        "expected at least one 46-way sum in spenso_eval_input_0.txt"
    );

    eprintln!("sum46 candidates={}", candidates.len());
    for (candidate_index, candidate) in candidates.iter().enumerate() {
        let terms = add_terms(candidate);
        let min_bytes = terms
            .iter()
            .map(|term| term.as_view().get_byte_size())
            .min()
            .unwrap_or(0);
        let max_bytes = terms
            .iter()
            .map(|term| term.as_view().get_byte_size())
            .max()
            .unwrap_or(0);
        let total_bytes = terms
            .iter()
            .map(|term| term.as_view().get_byte_size())
            .sum::<usize>();
        eprintln!(
            "sum46 candidate={candidate_index} bytes={} terms={} child_bytes_min={} child_bytes_max={} child_bytes_total={}",
            candidate.as_view().get_byte_size(),
            candidate.nterms(),
            min_bytes,
            max_bytes,
            total_bytes
        );
    }

    let selected_candidate = std::env::var("SPENSO_SUM46_CANDIDATE")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or_else(|| {
            candidates
                .iter()
                .enumerate()
                .max_by_key(|(_, candidate)| candidate.as_view().get_byte_size())
                .map(|(candidate_index, _)| candidate_index)
                .expect("non-empty candidates")
        });
    let add_expr = candidates
        .get(selected_candidate)
        .unwrap_or_else(|| panic!("SPENSO_SUM46_CANDIDATE={selected_candidate} is out of range"));
    let terms = add_terms(add_expr);
    let add_pattern = add_expr.to_pattern();
    eprintln!(
        "sum46 selected_candidate={selected_candidate} selected_bytes={} selected_terms={} selected_children={}",
        add_expr.as_view().get_byte_size(),
        add_expr.nterms(),
        terms.len()
    );

    let start_index = std::env::var("SPENSO_SUM46_START")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(0)
        .min(terms.len());
    let end_index = std::env::var("SPENSO_SUM46_END")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or_else(|| {
            std::env::var("SPENSO_SUM46_LIMIT")
                .ok()
                .and_then(|value| value.parse::<usize>().ok())
                .map(|limit| start_index.saturating_add(limit))
                .unwrap_or(terms.len())
        })
        .min(terms.len())
        .max(start_index);
    eprintln!("sum46 selected_range start={start_index} end={end_index}");

    let mut execute_times = Vec::with_capacity(end_index - start_index);
    for (term_index, term) in terms
        .iter()
        .enumerate()
        .skip(start_index)
        .take(end_index - start_index)
    {
        let label = format!("sum46_term_{term_index:02}");
        let total_start = Instant::now();
        let replace_start = Instant::now();
        let branch_expr = expr.replace(&add_pattern).once().with(term.to_pattern());
        let replace_elapsed = replace_start.elapsed();
        eprintln!(
            "{label} branch_expr terms={} bytes={} selected_term_terms={} selected_term_bytes={} replace={replace_elapsed:.3?}",
            branch_expr.nterms(),
            branch_expr.as_view().get_byte_size(),
            term.nterms(),
            term.as_view().get_byte_size()
        );
        eprintln!(
            "{label} function_counts selected_delta={} branch_delta={} selected_q3={} branch_q3={}",
            count_function_occurrences(term.as_view(), "δ"),
            count_function_occurrences(branch_expr.as_view(), "δ"),
            count_function_occurrences(term.as_view(), "Q3"),
            count_function_occurrences(branch_expr.as_view(), "Q3"),
        );

        let branch_expr = apply_sum46_prepass(&label, branch_expr);
        let net = parse_actual_net_quiet(&label, &branch_expr);
        if std::env::var_os("SPENSO_SUM46_NETWORK_ONLY").is_some() {
            report_actual_net_stats(&label, &net);
            report_actual_net_leaf_stats(&label, &net);
            report_product_interleaving(&label, &net, 5);
            continue;
        }

        let execute_elapsed = execute_actual_net_min_result_rank_parallel_summary(&label, net);
        execute_times.push(execute_elapsed);
        eprintln!("{label} total={:.3?}", total_start.elapsed());
    }

    if execute_times.is_empty() {
        eprintln!("sum46 execution_summary terms=0");
        return;
    }

    execute_times.sort();
    let sum = execute_times
        .iter()
        .copied()
        .fold(std::time::Duration::ZERO, |sum, elapsed| sum + elapsed);
    let max = execute_times.last().copied().unwrap_or_default();
    let median = execute_times[execute_times.len() / 2];
    eprintln!(
        "sum46 execution_summary terms={} sum={sum:.3?} median={median:.3?} max={max:.3?}",
        execute_times.len()
    );
}

#[test]
#[ignore = "diagnostic timing for boundary-factored spenso_eval_input_0.txt parallel MinResultRank"]
fn spenso_eval_input_0_boundary_factor_actual_network_parallel_min_result_rank_execute() {
    test_initialise().expect("GammaLoop initialization should succeed");
    let expr = parse_root_input("spenso_eval_input_0.txt");
    let net = parse_actual_net_with_boundary_factor(
        "spenso_eval_input_0_boundary_factor_parallel_min_result_rank",
        &expr,
    );
    let _ = execute_actual_net_min_result_rank_parallel(
        "spenso_eval_input_0_boundary_factor_parallel_min_result_rank",
        net,
    );
}

#[test]
#[ignore = "diagnostic product-boundary candidate selection for the larger spenso_eval_input_0.txt"]
fn spenso_eval_input_0_graph_boundary_candidate_diagnostics() {
    test_initialise().expect("GammaLoop initialization should succeed");
    let expr = parse_root_input("spenso_eval_input_0.txt");
    let net = parse_actual_net("spenso_eval_input_0_boundary_candidates", &expr);
    log_product_boundary_candidates("spenso_eval_input_0_boundary_candidates", &net, 24);
}

#[test]
#[ignore = "diagnostic timing for the larger spenso_eval_input_0.txt after Hornering"]
fn spenso_eval_input_0_horner_actual_network_parse() {
    test_initialise().expect("GammaLoop initialization should succeed");
    let expr = parse_root_input("spenso_eval_input_0.txt");
    let expr = collect_horner(&expr);
    let _net = parse_actual_net("spenso_eval_input_0_horner", &expr);
}

#[test]
#[ignore = "diagnostic timing for the larger spenso_eval_input_0.txt after Hornering and parallel MinResultRank"]
fn spenso_eval_input_0_horner_actual_network_parallel_min_result_rank_execute() {
    test_initialise().expect("GammaLoop initialization should succeed");
    let expr = parse_root_input("spenso_eval_input_0.txt");
    let expr = collect_horner(&expr);
    let net = parse_actual_net("spenso_eval_input_0_horner_parallel_min_result_rank", &expr);
    let _ = execute_actual_net_min_result_rank_parallel(
        "spenso_eval_input_0_horner_parallel_min_result_rank",
        net,
    );
}

#[test]
#[ignore = "diagnostic timing for the first few GammaLoop concrete/mixed tensor execution steps"]
fn symbolica_expression_input_actual_network_first_steps() {
    test_initialise().expect("GammaLoop initialization should succeed");
    let expr = parse_root_input("symbolica_expression.txt");
    let net = parse_actual_net("raw_steps", &expr);
    execute_actual_net_steps("raw_steps", net, 8);
}

#[test]
#[ignore = "diagnostic timing for GammaLoop's concrete/mixed tensor evaluator path after Hornering"]
fn symbolica_expression_input_horner_actual_network_execute() {
    test_initialise().expect("GammaLoop initialization should succeed");
    let expr = parse_root_input("symbolica_expression.txt");
    let expr = collect_horner(&expr);
    let net = parse_actual_net("horner", &expr);
    let _ = execute_actual_net("horner", net);
}

#[test]
#[ignore = "diagnostic timing for concrete HEP tensor data execution"]
fn concrete_hep_gamma_trace_probe() {
    test_initialise().expect("GammaLoop initialization should succeed");
    let expr = parse_inline_expression(
        "spenso::gamma(spenso::bis(4,1),spenso::bis(4,2),spenso::mink(4,3))*spenso::gamma(spenso::bis(4,2),spenso::bis(4,1),spenso::mink(4,3))",
    );
    let net = parse_concrete_hep_net("gamma_trace_probe", &expr);
    execute_concrete_hep_net("gamma_trace_probe", net);
}

#[test]
#[ignore = "diagnostic reproduction of the current left-deep SmallestDegree gamma ladder"]
fn gamma_ladder_mwe_smallest_degree_execute() {
    test_initialise().expect("GammaLoop initialization should succeed");
    let expr = gamma_ladder_order_mwe_expr();
    let net = parse_actual_net("gamma_ladder_mwe_smallest_degree", &expr);
    let result = execute_actual_net("gamma_ladder_mwe_smallest_degree", net);
    assert!(result.nterms() > 0);
}

#[test]
#[ignore = "diagnostic comparison order that minimizes intermediate result rank"]
fn gamma_ladder_mwe_min_result_rank_execute() {
    test_initialise().expect("GammaLoop initialization should succeed");
    let expr = gamma_ladder_order_mwe_expr();
    let net = parse_actual_net("gamma_ladder_mwe_min_result_rank", &expr);
    let result = execute_actual_net_min_result_rank("gamma_ladder_mwe_min_result_rank", net);
    assert!(result.nterms() > 0);
}

#[test]
#[ignore = "diagnostic equivalence check for contraction-order output forms"]
fn gamma_ladder_mwe_order_results_match_after_expand() {
    test_initialise().expect("GammaLoop initialization should succeed");
    let expr = gamma_ladder_order_mwe_expr();

    let smallest = execute_actual_net(
        "gamma_ladder_mwe_match_smallest_degree",
        parse_actual_net("gamma_ladder_mwe_match_smallest_degree", &expr),
    );
    let min_result_rank = execute_actual_net_min_result_rank(
        "gamma_ladder_mwe_match_min_result_rank",
        parse_actual_net("gamma_ladder_mwe_match_min_result_rank", &expr),
    );

    let start = Instant::now();
    let diff = (smallest - min_result_rank).expand();
    report("gamma_ladder_mwe_order_diff_expand", start);
    eprintln!(
        "gamma_ladder_mwe_order_diff stats: terms={} bytes={}",
        diff.nterms(),
        diff.as_view().get_byte_size()
    );
    assert!(
        diff.is_zero(),
        "contraction orders produced different scalar values: {}",
        diff
    );
}

#[test]
#[ignore = "diagnostic MWE for a late tensor-valued sum that should be factored before execution"]
fn late_tensor_sum_mwe_expanded_factored_and_hornered() {
    test_initialise().expect("GammaLoop initialization should succeed");
    let expanded = late_tensor_sum_mwe_expr(12, 12);
    let factored = factored_tensor_sum_mwe_expr(12, 12);
    let hornered = collect_horner(&expanded);

    let expanded_result = execute_actual_net_min_result_rank_parallel(
        "late_tensor_sum_mwe_expanded",
        parse_actual_net("late_tensor_sum_mwe_expanded", &expanded),
    );
    let boundary_factor_result = execute_actual_net_min_result_rank_parallel(
        "late_tensor_sum_mwe_boundary_factor",
        parse_actual_net_with_boundary_factor("late_tensor_sum_mwe_boundary_factor", &expanded),
    );
    let factored_result = execute_actual_net_min_result_rank_parallel(
        "late_tensor_sum_mwe_factored",
        parse_actual_net("late_tensor_sum_mwe_factored", &factored),
    );
    let hornered_result = execute_actual_net_min_result_rank_parallel(
        "late_tensor_sum_mwe_hornered",
        parse_actual_net("late_tensor_sum_mwe_hornered", &hornered),
    );

    let start = Instant::now();
    let expanded_factored_diff = (expanded_result.clone() - factored_result).expand();
    report("late_tensor_sum_mwe_expanded_factored_diff_expand", start);
    assert!(
        expanded_factored_diff.is_zero(),
        "expanded and factored MWE forms produced different scalar values: {}",
        expanded_factored_diff
    );

    let start = Instant::now();
    let expanded_boundary_factor_diff = (expanded_result.clone() - boundary_factor_result).expand();
    report(
        "late_tensor_sum_mwe_expanded_boundary_factor_diff_expand",
        start,
    );
    assert!(
        expanded_boundary_factor_diff.is_zero(),
        "expanded and boundary-factored MWE forms produced different scalar values: {}",
        expanded_boundary_factor_diff
    );

    let start = Instant::now();
    let expanded_hornered_diff = (expanded_result - hornered_result).expand();
    report("late_tensor_sum_mwe_expanded_hornered_diff_expand", start);
    assert!(
        expanded_hornered_diff.is_zero(),
        "expanded and Hornered MWE forms produced different scalar values: {}",
        expanded_hornered_diff
    );
}

#[test]
#[ignore = "diagnostic staged orchestration for a late tensor-valued sum"]
fn late_tensor_sum_mwe_staged_disconnect_reconnect_metric() {
    test_initialise().expect("GammaLoop initialization should succeed");
    let direct = late_tensor_sum_mwe_expr(12, 12);
    let disconnected = late_tensor_sum_mwe_disconnected_expr(12, 12);
    let reconnect_metric = late_tensor_sum_mwe_reconnect_metric_expr();

    let direct_result = execute_actual_net_min_result_rank_parallel(
        "late_tensor_sum_mwe_staged_direct",
        parse_actual_net("late_tensor_sum_mwe_staged_direct", &direct),
    );
    let intermediate = execute_actual_net_min_result_rank_parallel_tensor(
        "late_tensor_sum_mwe_staged_disconnected",
        parse_actual_net("late_tensor_sum_mwe_staged_disconnected", &disconnected),
    );
    let intermediate_net: ParsingNet = Network::from_tensor(intermediate);
    let reconnect_net = parse_actual_net(
        "late_tensor_sum_mwe_staged_reconnect_metric",
        &reconnect_metric,
    );
    let staged_result = execute_actual_net_min_result_rank_parallel(
        "late_tensor_sum_mwe_staged_reconnected",
        intermediate_net * reconnect_net,
    );

    let start = Instant::now();
    let diff = (direct_result - staged_result).expand();
    report("late_tensor_sum_mwe_staged_diff_expand", start);
    eprintln!(
        "late_tensor_sum_mwe_staged_diff stats: terms={} bytes={}",
        diff.nterms(),
        diff.as_view().get_byte_size()
    );
    assert!(
        diff.is_zero(),
        "direct and staged MWE forms produced different scalar values: {}",
        diff
    );
}

#[test]
#[ignore = "diagnostic graph-derived staged orchestration for a late tensor-valued sum"]
fn late_tensor_sum_mwe_graph_selected_disconnect_reconnect_metric() {
    test_initialise().expect("GammaLoop initialization should succeed");
    let direct = late_tensor_sum_mwe_expr(12, 12);

    let direct_result = execute_actual_net_min_result_rank_parallel(
        "late_tensor_sum_mwe_graph_selected_direct",
        parse_actual_net("late_tensor_sum_mwe_graph_selected_direct", &direct),
    );

    let mut staged_net = parse_actual_net("late_tensor_sum_mwe_graph_selected_staged", &direct);
    log_product_boundary_candidates("late_tensor_sum_mwe_graph_selected", &staged_net, 8);
    let candidates = graph_product_boundary_candidates(&staged_net);
    let candidate = candidates
        .iter()
        .find(|candidate| {
            candidate
                .left
                .max_sum_children
                .max(candidate.right.max_sum_children)
                >= 12
                && candidate.left.leaves.min(candidate.right.leaves) == 1
        })
        .unwrap_or_else(|| {
            candidates
                .first()
                .expect("MWE should expose at least one product-boundary slot candidate")
        })
        .clone();
    eprintln!(
        "late_tensor_sum_mwe_graph_selected chosen score={} product={} slot={} rename_child={} rename_hedge={}",
        candidate.score,
        candidate.product,
        candidate.slot,
        candidate.rename_child,
        candidate.rename_hedge,
    );
    let (original, fresh) = apply_graph_boundary_disconnect(&mut staged_net, &candidate);
    let intermediate = execute_actual_net_min_result_rank_parallel_tensor(
        "late_tensor_sum_mwe_graph_selected_disconnected",
        staged_net,
    );
    let intermediate_net: ParsingNet = Network::from_tensor(intermediate);
    let reconnect_metric = reconnect_metric_expr(original, fresh);
    let reconnect_net = parse_actual_net(
        "late_tensor_sum_mwe_graph_selected_reconnect_metric",
        &reconnect_metric,
    );
    let staged_result = execute_actual_net_min_result_rank_parallel(
        "late_tensor_sum_mwe_graph_selected_reconnected",
        intermediate_net * reconnect_net,
    );

    let start = Instant::now();
    let diff = (direct_result - staged_result).expand();
    report("late_tensor_sum_mwe_graph_selected_diff_expand", start);
    eprintln!(
        "late_tensor_sum_mwe_graph_selected_diff stats: terms={} bytes={}",
        diff.nterms(),
        diff.as_view().get_byte_size()
    );
    assert!(
        diff.is_zero(),
        "direct and graph-selected staged MWE forms produced different scalar values: {}",
        diff
    );
}
