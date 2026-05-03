use std::{
    collections::{BTreeMap, BTreeSet},
    fmt::Write,
};

use serde::{Deserialize, Serialize};
use symbolica::{atom::AtomView, try_parse};
use thiserror::Error;

use crate::graph_io::{ParsedGraph, ParsedGraphExternalEdge, ParsedGraphInternalEdge};

#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub struct MomentumSignature {
    pub loop_signature: Vec<i32>,
    pub external_signature: Vec<i32>,
}

impl MomentumSignature {
    pub fn negated(&self) -> Self {
        Self {
            loop_signature: self.loop_signature.iter().map(|x| -*x).collect(),
            external_signature: self.external_signature.iter().map(|x| -*x).collect(),
        }
    }

    pub fn canonical_up_to_sign(&self) -> (Self, i32) {
        let negated = self.negated();
        if self <= &negated {
            (self.clone(), 1)
        } else {
            (negated, -1)
        }
    }
}

#[derive(Debug, Error)]
pub enum SignatureParseError {
    #[error("empty momentum expression")]
    EmptyExpression,
    #[error("malformed momentum expression `{0}`")]
    MalformedExpression(String),
    #[error("unsupported momentum token `{token}` in `{expression}`")]
    UnsupportedToken { token: String, expression: String },
    #[error("unknown momentum symbol `{symbol}` in `{expression}`")]
    UnknownSymbol { symbol: String, expression: String },
    #[error(
        "momentum symbol `{symbol}` has coefficient {coefficient}, only -1, 0, +1 are supported"
    )]
    UnsupportedCoefficient { symbol: String, coefficient: i32 },
    #[error("failed to parse propagator expression with Symbolica: {0}")]
    SymbolicaParse(String),
    #[error("propagator pattern `{0}` is not supported by the Rust port yet")]
    UnsupportedPropPattern(String),
    #[error("no propagators found in expression `{0}`")]
    NoPropagators(String),
    #[error("propagator call `{call}` has {actual} arguments, expected 2")]
    InvalidPropagatorArity { call: String, actual: usize },
    #[error(
        "momentum symbol `{symbol}` does not match loop prefix `{loop_prefix}` or external prefixes {external_prefixes:?} in `{expression}`"
    )]
    UnknownPrefixedSymbol {
        symbol: String,
        loop_prefix: String,
        external_prefixes: Vec<String>,
        expression: String,
    },
}

pub type Result<T> = std::result::Result<T, SignatureParseError>;

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ExtractedSignatureExpression {
    pub signatures: Vec<MomentumSignature>,
    pub loop_names: Vec<String>,
    pub external_names: Vec<String>,
    pub masses: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ReconstructDotOptions {
    pub num_vertices: Option<usize>,
    pub require_connected: bool,
    pub max_degree: Option<usize>,
    pub format: ReconstructDotFormat,
    pub graph_engine: Option<String>,
    pub minimize_external_legs: bool,
}

impl Default for ReconstructDotOptions {
    fn default() -> Self {
        Self {
            num_vertices: None,
            require_connected: true,
            max_degree: None,
            format: ReconstructDotFormat::Gammaloop,
            graph_engine: None,
            minimize_external_legs: false,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ReconstructDotFormat {
    Gammaloop,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct EdgePlacement {
    tail: usize,
    head: usize,
}

pub fn extract_signatures_and_masses_from_symbolica_expression(
    expression: &str,
    loop_prefix: &str,
    external_prefixes: &[String],
    prop_pattern: &str,
) -> Result<ExtractedSignatureExpression> {
    if prop_pattern != "prop(q_,m_)" {
        return Err(SignatureParseError::UnsupportedPropPattern(
            prop_pattern.to_string(),
        ));
    }

    let atom = try_parse!(expression)
        .map_err(|error| SignatureParseError::SymbolicaParse(error.to_string()))?;
    let mut calls = Vec::new();
    collect_prop_calls(atom.as_view(), "prop", &mut calls)?;
    if calls.is_empty() {
        return Err(SignatureParseError::NoPropagators(expression.to_string()));
    }

    build_signatures_from_qm_pairs(calls, loop_prefix, external_prefixes)
}

pub fn reconstruct_dot_from_expression(
    expression: &str,
    loop_prefix: &str,
    external_prefixes: &[String],
    prop_pattern: &str,
    options: &ReconstructDotOptions,
) -> Result<(String, ExtractedSignatureExpression)> {
    let extracted = extract_signatures_and_masses_from_symbolica_expression(
        expression,
        loop_prefix,
        external_prefixes,
        prop_pattern,
    )?;
    let dot = reconstruct_dot(&extracted, options)?;
    Ok((dot, extracted))
}

pub fn reconstruct_dot(
    extracted: &ExtractedSignatureExpression,
    options: &ReconstructDotOptions,
) -> Result<String> {
    let signatures = &extracted.signatures;
    if signatures.is_empty() {
        return Err(SignatureParseError::EmptyExpression);
    }

    let n_edges = signatures.len();
    let n_loops = signatures[0].loop_signature.len();
    let n_external = signatures[0].external_signature.len();
    for signature in signatures {
        if signature.loop_signature.len() != n_loops
            || signature.external_signature.len() != n_external
        {
            return Err(SignatureParseError::MalformedExpression(
                "inconsistent signature dimensions".to_string(),
            ));
        }
    }

    let num_vertices = resolve_num_vertices(n_edges, n_loops, options)?;

    let placements = find_edge_placements(signatures, num_vertices, options)?;
    let ext_balance = external_balance(signatures, &placements, num_vertices, n_external);

    match options.format {
        ReconstructDotFormat::Gammaloop => {
            Ok(render_gammaloop_dot(extracted, &placements, &ext_balance))
        }
    }
}

pub(crate) fn reconstruct_parsed_graph(
    extracted: &ExtractedSignatureExpression,
    options: &ReconstructDotOptions,
) -> Result<ParsedGraph> {
    let signatures = &extracted.signatures;
    if signatures.is_empty() {
        return Err(SignatureParseError::EmptyExpression);
    }

    let n_edges = signatures.len();
    let n_loops = signatures[0].loop_signature.len();
    let n_external = signatures[0].external_signature.len();
    for signature in signatures {
        if signature.loop_signature.len() != n_loops
            || signature.external_signature.len() != n_external
        {
            return Err(SignatureParseError::MalformedExpression(
                "inconsistent signature dimensions".to_string(),
            ));
        }
    }

    let num_vertices = resolve_num_vertices(n_edges, n_loops, options)?;

    let placements = find_edge_placements(signatures, num_vertices, options)?;
    let ext_balance = external_balance(signatures, &placements, num_vertices, n_external);
    let node_name_to_internal = (0..num_vertices)
        .map(|node| (format!("v{node}"), node))
        .collect::<BTreeMap<_, _>>();

    let internal_edges = extracted
        .signatures
        .iter()
        .zip(&placements)
        .enumerate()
        .map(
            |(edge_id, (signature, placement))| ParsedGraphInternalEdge {
                edge_id,
                tail: placement.tail,
                head: placement.head,
                label: format_momentum(signature, &extracted.loop_names, &extracted.external_names),
                mass_key: extracted
                    .masses
                    .get(edge_id)
                    .filter(|mass| !mass.is_empty())
                    .cloned(),
                signature: signature.clone(),
                had_pow: false,
            },
        )
        .collect();

    let mut external_edges = Vec::new();
    let mut next_external_id = 10_000_000usize;
    for (vertex, balance) in ext_balance.iter().enumerate() {
        for (external_id, coeff) in balance.iter().enumerate() {
            let name = &extracted.external_names[external_id];
            if *coeff > 0 {
                for _ in 0..*coeff {
                    let mut external_coefficients = vec![0; n_external];
                    external_coefficients[external_id] = -1;
                    external_edges.push(ParsedGraphExternalEdge {
                        edge_id: next_external_id,
                        source: Some(vertex),
                        destination: None,
                        label: format!("-{name}"),
                        external_coefficients,
                    });
                    next_external_id += 1;
                }
            } else if *coeff < 0 {
                for _ in 0..(-*coeff) {
                    let mut external_coefficients = vec![0; n_external];
                    external_coefficients[external_id] = 1;
                    external_edges.push(ParsedGraphExternalEdge {
                        edge_id: next_external_id,
                        source: None,
                        destination: Some(vertex),
                        label: name.clone(),
                        external_coefficients,
                    });
                    next_external_id += 1;
                }
            }
        }
    }

    Ok(ParsedGraph {
        internal_edges,
        external_edges,
        initial_state_cut_edges: Vec::new(),
        loop_names: extracted.loop_names.clone(),
        external_names: extracted.external_names.clone(),
        node_name_to_internal,
    })
}

fn resolve_num_vertices(
    n_edges: usize,
    n_loops: usize,
    options: &ReconstructDotOptions,
) -> Result<usize> {
    let expected_connected = n_edges.saturating_sub(n_loops).saturating_add(1).max(2);
    let Some(requested) = options.num_vertices else {
        return Ok(expected_connected);
    };
    let requested = requested.max(2);
    if options.require_connected && requested != expected_connected {
        return Err(SignatureParseError::MalformedExpression(format!(
            "connected explicit-edge realization with {n_edges} propagator signatures and {n_loops} loop momentum basis vector(s) requires {expected_connected} internal vertices, but --num-vertices requested {requested}; repeated/dotted propagators are currently represented as explicit repeated edges"
        )));
    }
    Ok(requested)
}

fn collect_prop_calls(
    view: AtomView<'_>,
    prop_name: &str,
    out: &mut Vec<(String, String)>,
) -> Result<()> {
    match view {
        AtomView::Fun(fun) => {
            if fun.get_symbol().get_stripped_name() == prop_name {
                let args = fun.iter().collect::<Vec<_>>();
                if args.len() != 2 {
                    return Err(SignatureParseError::InvalidPropagatorArity {
                        call: view.to_plain_string(),
                        actual: args.len(),
                    });
                }
                out.push((
                    atom_view_user_string(args[0]),
                    atom_view_user_string(args[1]),
                ));
            } else {
                for arg in fun.iter() {
                    collect_prop_calls(arg, prop_name, out)?;
                }
            }
        }
        AtomView::Add(add) => {
            for arg in add.iter() {
                collect_prop_calls(arg, prop_name, out)?;
            }
        }
        AtomView::Mul(mul) => {
            for arg in mul.iter() {
                collect_prop_calls(arg, prop_name, out)?;
            }
        }
        AtomView::Pow(pow) => {
            for arg in pow.iter() {
                collect_prop_calls(arg, prop_name, out)?;
            }
        }
        AtomView::Num(_) | AtomView::Var(_) => {}
    }
    Ok(())
}

fn atom_view_user_string(view: AtomView<'_>) -> String {
    strip_symbolica_namespaces_from_text(&view.to_plain_string())
}

fn strip_symbolica_namespaces_from_text(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    let mut token = String::new();

    let flush_token = |token: &mut String, out: &mut String| {
        if token.is_empty() {
            return;
        }
        if token.contains("::") {
            out.push_str(token.rsplit("::").next().unwrap_or(token));
        } else {
            out.push_str(token);
        }
        token.clear();
    };

    for ch in text.chars() {
        if ch == ':' || ch == '_' || ch.is_ascii_alphanumeric() {
            token.push(ch);
        } else {
            flush_token(&mut token, &mut out);
            out.push(ch);
        }
    }
    flush_token(&mut token, &mut out);
    out
}

fn build_signatures_from_qm_pairs(
    qm_pairs: Vec<(String, String)>,
    loop_prefix: &str,
    external_prefixes: &[String],
) -> Result<ExtractedSignatureExpression> {
    let mut edge_loop = Vec::<Vec<(String, i32)>>::new();
    let mut edge_external = Vec::<Vec<(String, i32)>>::new();
    let mut loop_names = BTreeSet::new();
    let mut external_names = BTreeSet::new();
    let mut masses = Vec::new();

    for (momentum_expression, mass_expression) in qm_pairs {
        let normalized = normalize_momentum_label(&momentum_expression);
        let mut loop_coeffs = Vec::new();
        let mut external_coeffs = Vec::new();
        for (sign, symbol) in split_linear_terms(&normalized)? {
            if symbol.starts_with(loop_prefix) {
                loop_names.insert(symbol.clone());
                loop_coeffs.push((symbol, sign));
            } else if external_prefixes
                .iter()
                .any(|prefix| symbol.starts_with(prefix))
            {
                external_names.insert(symbol.clone());
                external_coeffs.push((symbol, sign));
            } else {
                return Err(SignatureParseError::UnknownPrefixedSymbol {
                    symbol,
                    loop_prefix: loop_prefix.to_string(),
                    external_prefixes: external_prefixes.to_vec(),
                    expression: momentum_expression,
                });
            }
        }
        edge_loop.push(loop_coeffs);
        edge_external.push(external_coeffs);
        masses.push(mass_expression.trim().to_string());
    }

    let mut loop_names = loop_names.into_iter().collect::<Vec<_>>();
    let mut external_names = external_names.into_iter().collect::<Vec<_>>();
    loop_names.sort_by_key(|name| symbol_sort_key(name));
    external_names.sort_by_key(|name| symbol_sort_key(name));

    let mut signatures = Vec::new();
    for (loop_coeffs, external_coeffs) in edge_loop.into_iter().zip(edge_external) {
        let mut loop_signature = vec![0; loop_names.len()];
        let mut external_signature = vec![0; external_names.len()];
        for (symbol, coeff) in loop_coeffs {
            let index = loop_names
                .iter()
                .position(|name| name == &symbol)
                .expect("loop name was inserted before signature build");
            loop_signature[index] += coeff;
        }
        for (symbol, coeff) in external_coeffs {
            let index = external_names
                .iter()
                .position(|name| name == &symbol)
                .expect("external name was inserted before signature build");
            external_signature[index] += coeff;
        }
        for (symbol, coefficient) in loop_names
            .iter()
            .chain(external_names.iter())
            .zip(loop_signature.iter().chain(external_signature.iter()))
        {
            if !matches!(*coefficient, -1..=1) {
                return Err(SignatureParseError::UnsupportedCoefficient {
                    symbol: symbol.clone(),
                    coefficient: *coefficient,
                });
            }
        }
        signatures.push(MomentumSignature {
            loop_signature,
            external_signature,
        });
    }

    Ok(ExtractedSignatureExpression {
        signatures,
        loop_names,
        external_names,
        masses,
    })
}

fn find_edge_placements(
    signatures: &[MomentumSignature],
    num_vertices: usize,
    options: &ReconstructDotOptions,
) -> Result<Vec<EdgePlacement>> {
    let n_edges = signatures.len();
    let n_loops = signatures[0].loop_signature.len();
    let mut remaining_abs = vec![vec![0; n_edges + 1]; n_loops];
    for (loop_id, remaining) in remaining_abs.iter_mut().enumerate().take(n_loops) {
        for edge_id in (0..n_edges).rev() {
            remaining[edge_id] =
                remaining[edge_id + 1] + signatures[edge_id].loop_signature[loop_id].abs();
        }
    }

    let mut placements = vec![None; n_edges];
    let mut loop_balance = vec![vec![0i32; n_loops]; num_vertices];
    let n_external = signatures[0].external_signature.len();
    let mut external_balance = vec![vec![0i32; n_external]; num_vertices];
    let mut degree = vec![0usize; num_vertices];
    let mut dsu = RollbackDsu::new(num_vertices);
    let mut best_placements = None;
    let mut best_external_cost = ExternalCost::worst();

    struct Search<'a> {
        signatures: &'a [MomentumSignature],
        options: &'a ReconstructDotOptions,
        remaining_abs: &'a [Vec<i32>],
        placements: &'a mut [Option<EdgePlacement>],
        loop_balance: &'a mut [Vec<i32>],
        external_balance: &'a mut [Vec<i32>],
        degree: &'a mut [usize],
        dsu: &'a mut RollbackDsu,
        best_placements: &'a mut Option<Vec<EdgePlacement>>,
        best_external_cost: &'a mut ExternalCost,
        num_vertices: usize,
        n_loops: usize,
        n_external: usize,
    }

    impl Search<'_> {
        fn recurse(
            &mut self,
            edge_index: usize,
            components_now: usize,
            used_vertices_now: usize,
        ) -> bool {
            if self.prune(edge_index, components_now) {
                return false;
            }
            if edge_index == self.signatures.len() {
                let is_valid = self
                    .loop_balance
                    .iter()
                    .all(|row| row.iter().all(|value| *value == 0))
                    && (!self.options.require_connected || components_now == 1);
                if !is_valid {
                    return false;
                }
                if !self.options.minimize_external_legs {
                    return true;
                }

                let cost = ExternalCost::from_balance(self.external_balance);
                if cost < *self.best_external_cost {
                    *self.best_external_cost = cost;
                    *self.best_placements = Some(
                        self.placements
                            .iter()
                            .map(|placement| placement.expect("leaf search has all edges placed"))
                            .collect(),
                    );
                }
                return false;
            }

            let candidate_limit = (used_vertices_now
                + usize::from(used_vertices_now < self.num_vertices))
            .min(self.num_vertices);
            let candidates = if edge_index == 0 {
                vec![(0, 1)]
            } else {
                (0..candidate_limit)
                    .flat_map(|tail| {
                        (0..candidate_limit)
                            .filter(move |head| *head != tail)
                            .map(move |head| (tail, head))
                    })
                    .collect::<Vec<_>>()
            };

            for (tail, head) in candidates {
                let snapshot = self.dsu.snapshot();
                let mut next_components = components_now;
                if self.dsu.union(tail, head) {
                    next_components -= 1;
                }

                self.placements[edge_index] = Some(EdgePlacement { tail, head });
                for loop_id in 0..self.n_loops {
                    let coeff = self.signatures[edge_index].loop_signature[loop_id];
                    self.loop_balance[tail][loop_id] -= coeff;
                    self.loop_balance[head][loop_id] += coeff;
                }
                for external_id in 0..self.n_external {
                    let coeff = self.signatures[edge_index].external_signature[external_id];
                    self.external_balance[tail][external_id] -= coeff;
                    self.external_balance[head][external_id] += coeff;
                }
                self.degree[tail] += 1;
                self.degree[head] += 1;

                let next_used_vertices = if used_vertices_now < self.num_vertices
                    && (tail == used_vertices_now || head == used_vertices_now)
                {
                    used_vertices_now + 1
                } else {
                    used_vertices_now
                };

                if self.recurse(edge_index + 1, next_components, next_used_vertices) {
                    return true;
                }

                self.degree[tail] -= 1;
                self.degree[head] -= 1;
                for external_id in 0..self.n_external {
                    let coeff = self.signatures[edge_index].external_signature[external_id];
                    self.external_balance[tail][external_id] += coeff;
                    self.external_balance[head][external_id] -= coeff;
                }
                for loop_id in 0..self.n_loops {
                    let coeff = self.signatures[edge_index].loop_signature[loop_id];
                    self.loop_balance[tail][loop_id] += coeff;
                    self.loop_balance[head][loop_id] -= coeff;
                }
                self.placements[edge_index] = None;
                self.dsu.rollback(snapshot);
            }

            false
        }

        fn prune(&self, edge_index: usize, components_now: usize) -> bool {
            let remaining_edges = self.signatures.len() - edge_index;
            if self.options.require_connected && remaining_edges < components_now.saturating_sub(1)
            {
                return true;
            }
            for vertex in 0..self.num_vertices {
                for loop_id in 0..self.n_loops {
                    if self.loop_balance[vertex][loop_id].abs()
                        > self.remaining_abs[loop_id][edge_index]
                    {
                        return true;
                    }
                }
            }
            if let Some(max_degree) = self.options.max_degree
                && self.degree.iter().any(|degree| *degree > max_degree)
            {
                return true;
            }
            false
        }
    }

    let mut search = Search {
        signatures,
        options,
        remaining_abs: &remaining_abs,
        placements: &mut placements,
        loop_balance: &mut loop_balance,
        external_balance: &mut external_balance,
        degree: &mut degree,
        dsu: &mut dsu,
        best_placements: &mut best_placements,
        best_external_cost: &mut best_external_cost,
        num_vertices,
        n_loops,
        n_external,
    };

    let found = search.recurse(0, num_vertices, 2.min(num_vertices));
    if options.minimize_external_legs
        && let Some(best) = best_placements
    {
        return Ok(best);
    }
    if !found {
        return Err(SignatureParseError::MalformedExpression(
            "no connected DOT realization found for signatures".to_string(),
        ));
    }

    Ok(placements
        .into_iter()
        .map(|placement| placement.expect("search succeeded with all edges placed"))
        .collect())
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct ExternalCost {
    total_legs: i32,
    active_vertices: usize,
}

impl ExternalCost {
    fn worst() -> Self {
        Self {
            total_legs: i32::MAX,
            active_vertices: usize::MAX,
        }
    }

    fn from_balance(balance: &[Vec<i32>]) -> Self {
        Self {
            total_legs: balance.iter().flatten().map(|coeff| coeff.abs()).sum(),
            active_vertices: balance
                .iter()
                .filter(|row| row.iter().any(|coeff| *coeff != 0))
                .count(),
        }
    }
}

#[derive(Debug, Clone)]
struct RollbackDsu {
    parent: Vec<usize>,
    size: Vec<usize>,
    changes: Vec<(usize, usize, usize)>,
}

impl RollbackDsu {
    fn new(size: usize) -> Self {
        Self {
            parent: (0..size).collect(),
            size: vec![1; size],
            changes: Vec::new(),
        }
    }

    fn find(&self, mut value: usize) -> usize {
        while self.parent[value] != value {
            value = self.parent[value];
        }
        value
    }

    fn snapshot(&self) -> usize {
        self.changes.len()
    }

    fn rollback(&mut self, snapshot: usize) {
        while self.changes.len() > snapshot {
            let (child_root, old_parent, old_size_parent) = self.changes.pop().unwrap();
            let parent_root = self.parent[child_root];
            self.parent[child_root] = old_parent;
            self.size[parent_root] = old_size_parent;
        }
    }

    fn union(&mut self, lhs: usize, rhs: usize) -> bool {
        let mut lhs_root = self.find(lhs);
        let mut rhs_root = self.find(rhs);
        if lhs_root == rhs_root {
            return false;
        }
        if self.size[lhs_root] < self.size[rhs_root] {
            std::mem::swap(&mut lhs_root, &mut rhs_root);
        }
        self.changes
            .push((rhs_root, self.parent[rhs_root], self.size[lhs_root]));
        self.parent[rhs_root] = lhs_root;
        self.size[lhs_root] += self.size[rhs_root];
        true
    }
}

fn external_balance(
    signatures: &[MomentumSignature],
    placements: &[EdgePlacement],
    num_vertices: usize,
    n_external: usize,
) -> Vec<Vec<i32>> {
    let mut balance = vec![vec![0; n_external]; num_vertices];
    for (signature, placement) in signatures.iter().zip(placements) {
        for (external_id, coeff) in signature
            .external_signature
            .iter()
            .enumerate()
            .take(n_external)
        {
            balance[placement.tail][external_id] -= *coeff;
            balance[placement.head][external_id] += *coeff;
        }
    }
    balance
}

fn render_gammaloop_dot(
    extracted: &ExtractedSignatureExpression,
    placements: &[EdgePlacement],
    ext_balance: &[Vec<i32>],
) -> String {
    let mut out = String::new();
    writeln!(&mut out, "digraph reconstructed_3drep {{").unwrap();
    writeln!(&mut out, "  num=\"1\";").unwrap();
    writeln!(&mut out, "  overall_factor=\"1\";").unwrap();
    writeln!(&mut out, "  edge [num=\"1\"];").unwrap();
    writeln!(&mut out, "  node [num=\"1\"];").unwrap();
    writeln!(&mut out, "  ext [style=invis];").unwrap();

    for vertex in 0..ext_balance.len() {
        writeln!(&mut out, "  v{vertex};").unwrap();
    }

    let mut next_hedge_id = 0usize;
    let mut next_edge_id = 0usize;
    for (external_id, name) in extracted.external_names.iter().enumerate() {
        writeln!(
            &mut out,
            "  ext_{name} [style=invis,label=\"{}\"];",
            dot_escape(name)
        )
        .unwrap();
        let _ = external_id;
    }

    for (vertex, balance) in ext_balance.iter().enumerate() {
        for (external_id, coeff) in balance.iter().enumerate() {
            let name = &extracted.external_names[external_id];
            if *coeff > 0 {
                for _ in 0..*coeff {
                    writeln!(
                        &mut out,
                        "  v{vertex}:{} -> ext_{name} [id={next_edge_id} dir=none name=\"{}\" particle=\"scalar_0\" mass=\"0\"];",
                        next_hedge_id,
                        dot_escape(name),
                    )
                    .unwrap();
                    next_edge_id += 1;
                    next_hedge_id += 1;
                }
            } else if *coeff < 0 {
                for _ in 0..(-*coeff) {
                    writeln!(
                        &mut out,
                        "  ext_{name} -> v{vertex}:{} [id={next_edge_id} dir=none name=\"{}\" particle=\"scalar_0\" mass=\"0\"];",
                        next_hedge_id,
                        dot_escape(name),
                    )
                    .unwrap();
                    next_edge_id += 1;
                    next_hedge_id += 1;
                }
            }
        }
    }

    let loop_carriers = loop_carrier_edges(&extracted.signatures);
    for (input_edge_id, (_signature, placement)) in
        extracted.signatures.iter().zip(placements).enumerate()
    {
        let mass = &extracted.masses[input_edge_id];
        let lmb_id_attr = loop_carriers
            .get(&input_edge_id)
            .map(|loop_id| format!(" lmb_id=\"{loop_id}\""))
            .unwrap_or_default();
        writeln!(
            &mut out,
            "  v{}:{} -> v{}:{} [id={next_edge_id} dir=none{} name=\"e{next_edge_id}\" particle=\"scalar_1\" mass=\"{}\"];",
            placement.tail,
            next_hedge_id,
            placement.head,
            next_hedge_id + 1,
            lmb_id_attr,
            dot_escape(mass),
        )
        .unwrap();
        next_edge_id += 1;
        next_hedge_id += 2;
    }

    writeln!(&mut out, "}}").unwrap();
    out
}

fn loop_carrier_edges(signatures: &[MomentumSignature]) -> BTreeMap<usize, usize> {
    let mut carriers = BTreeMap::new();
    let Some(first) = signatures.first() else {
        return carriers;
    };
    for loop_id in 0..first.loop_signature.len() {
        if let Some((edge_id, _)) = signatures.iter().enumerate().find(|(_, signature)| {
            signature.external_signature.iter().all(|coeff| *coeff == 0)
                && signature
                    .loop_signature
                    .iter()
                    .enumerate()
                    .all(|(other, coeff)| {
                        if other == loop_id {
                            *coeff == 1
                        } else {
                            *coeff == 0
                        }
                    })
        }) {
            carriers.insert(edge_id, loop_id);
        }
    }
    carriers
}

fn format_momentum(
    signature: &MomentumSignature,
    loop_names: &[String],
    external_names: &[String],
) -> String {
    let mut parts = Vec::<(String, i32)>::new();
    for (name, coeff) in loop_names.iter().zip(&signature.loop_signature) {
        if *coeff != 0 {
            parts.push((name.clone(), *coeff));
        }
    }
    for (name, coeff) in external_names.iter().zip(&signature.external_signature) {
        if *coeff != 0 {
            parts.push((name.clone(), *coeff));
        }
    }
    if parts.is_empty() {
        return "0".to_string();
    }

    format_signed_terms(parts)
}

fn format_signed_terms(parts: Vec<(String, i32)>) -> String {
    if parts.is_empty() {
        return "0".to_string();
    }

    let mut out = String::new();
    for (index, (name, coeff)) in parts.iter().enumerate() {
        let abs_coeff = coeff.abs();
        if index == 0 {
            if *coeff < 0 {
                out.push('-');
            }
        } else if *coeff < 0 {
            out.push('-');
        } else {
            out.push('+');
        }
        if abs_coeff != 1 {
            write!(&mut out, "{abs_coeff}*").unwrap();
        }
        out.push_str(name);
    }
    out
}

fn dot_escape(value: &str) -> String {
    value.replace('\\', "\\\\").replace('"', "\\\"")
}

fn normalize_momentum_label(label: &str) -> String {
    let mut text = parse_lmb_rep(label);
    text = text.replace(' ', "");
    loop {
        let normalized = text
            .replace("+-", "-")
            .replace("-+", "-")
            .replace("--", "+")
            .replace("++", "+");
        if normalized == text {
            return if normalized.is_empty() {
                "0".to_string()
            } else {
                normalized
            };
        }
        text = normalized;
    }
}

fn parse_lmb_rep(label: &str) -> String {
    let mut text = trim_quotes(label).trim().to_string();
    if text.is_empty() {
        return text;
    }

    text = replace_indexed_symbol(&text, "K", "k");
    text = replace_indexed_symbol(&text, "P", "p");
    text = text.replace("*", "");
    text = text.replace("-1", "-");
    text
}

fn replace_indexed_symbol(text: &str, head: &str, replacement_prefix: &str) -> String {
    let mut out = String::with_capacity(text.len());
    let mut rest = text;

    while let Some(start) = rest.find(&format!("{head}(")) {
        out.push_str(&rest[..start]);
        let after = &rest[start + head.len() + 1..];
        let Some(end) = after.find(')') else {
            out.push_str(&rest[start..]);
            return out;
        };
        let args = &after[..end];
        let first_arg = args.split(',').next().unwrap_or("").trim();
        if let Ok(index) = first_arg.parse::<usize>() {
            out.push_str(replacement_prefix);
            out.push_str(&(index + 1).to_string());
        } else {
            out.push_str(&rest[start..start + head.len() + 1 + end + 1]);
        }
        rest = &after[end + 1..];
    }

    out.push_str(rest);
    out
}

fn split_linear_terms(expression: &str) -> Result<Vec<(i32, String)>> {
    let expression = normalize_sign_runs(expression);
    if expression.is_empty() || expression == "0" {
        return Ok(Vec::new());
    }

    let mut text = expression.clone();
    if !text.starts_with('+') && !text.starts_with('-') {
        text.insert(0, '+');
    }

    let bytes = text.as_bytes();
    let mut index = 0;
    let mut terms = Vec::new();
    while index < bytes.len() {
        let sign = match bytes[index] {
            b'+' => 1,
            b'-' => -1,
            _ => {
                return Err(SignatureParseError::MalformedExpression(
                    expression.to_string(),
                ));
            }
        };
        index += 1;
        let start = index;
        while index < bytes.len() && !matches!(bytes[index], b'+' | b'-') {
            index += 1;
        }
        let token = &text[start..index];
        if token.is_empty() {
            return Err(SignatureParseError::MalformedExpression(
                expression.to_string(),
            ));
        }
        if token
            .chars()
            .any(|ch| matches!(ch, '*' | '/' | '^' | '(' | ')' | '[' | ']'))
        {
            return Err(SignatureParseError::UnsupportedToken {
                token: token.to_string(),
                expression: expression.to_string(),
            });
        }
        let symbol = token.rsplit("::").next().unwrap_or(token);
        if !is_identifier(symbol) {
            return Err(SignatureParseError::UnsupportedToken {
                token: token.to_string(),
                expression: expression.to_string(),
            });
        }
        terms.push((sign, symbol.to_string()));
    }

    Ok(terms)
}

fn normalize_sign_runs(input: &str) -> String {
    let mut text = input.replace(' ', "");
    loop {
        let normalized = text
            .replace("+-", "-")
            .replace("-+", "-")
            .replace("--", "+")
            .replace("++", "+");
        if normalized == text {
            return normalized;
        }
        text = normalized;
    }
}

fn is_identifier(value: &str) -> bool {
    let mut chars = value.chars();
    matches!(chars.next(), Some(ch) if ch == '_' || ch.is_ascii_alphabetic())
        && chars.all(|ch| ch == '_' || ch.is_ascii_alphanumeric())
}

fn symbol_sort_key(name: &str) -> (String, usize, String) {
    let split = name
        .char_indices()
        .rev()
        .find(|(_, ch)| !ch.is_ascii_digit())
        .map(|(index, ch)| index + ch.len_utf8())
        .unwrap_or(0);
    if split < name.len()
        && split > 0
        && let Ok(value) = name[split..].parse::<usize>()
    {
        return (name[..split].to_string(), value, name.to_string());
    }
    (name.to_string(), usize::MAX, name.to_string())
}

fn trim_quotes(value: &str) -> &str {
    let trimmed = value.trim();
    if trimmed.len() >= 2 {
        let bytes = trimmed.as_bytes();
        if matches!(
            (bytes[0], bytes[trimmed.len() - 1]),
            (b'"', b'"') | (b'\'', b'\'')
        ) {
            return &trimmed[1..trimmed.len() - 1];
        }
    }
    trimmed
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn symbolica_prop_expression_extracts_signatures_and_masses() {
        let extracted = extract_signatures_and_masses_from_symbolica_expression(
            "prop(k1+p1,mA)*prop(k1+p1-q1,mB)*prop(k1-p2+q2,mC)*prop(k1,mD)",
            "k",
            &["p".to_string(), "q".to_string()],
            "prop(q_,m_)",
        )
        .unwrap();

        assert_eq!(extracted.loop_names, vec!["k1"]);
        assert_eq!(extracted.external_names, vec!["p1", "p2", "q1", "q2"]);
        assert_eq!(extracted.masses, vec!["mD", "mA", "mB", "mC"]);
        assert_eq!(extracted.signatures.len(), 4);
        assert_eq!(extracted.signatures[2].loop_signature, vec![1]);
        assert_eq!(
            extracted.signatures[2].external_signature,
            vec![1, 0, -1, 0]
        );
    }

    #[test]
    fn symbolica_prop_expression_extracts_two_loop_signatures() {
        let extracted = extract_signatures_and_masses_from_symbolica_expression(
            "prop(k1+p1,0)*prop(k1+k2+p1,0)*prop(k1+k2+p1-q1,0)*prop(k1+k2-p2+q2,0)*prop(k1+k2-p2,0)*prop(k1-p2,0)*prop(k1,m1)*prop(k2,m1)",
            "k",
            &["p".to_string(), "q".to_string()],
            "prop(q_,m_)",
        )
        .unwrap();

        assert_eq!(extracted.loop_names, vec!["k1", "k2"]);
        assert_eq!(extracted.external_names, vec!["p1", "p2", "q1", "q2"]);
        assert_eq!(
            extracted.masses,
            vec!["m1", "m1", "0", "0", "0", "0", "0", "0"]
        );
        assert_eq!(extracted.signatures.len(), 8);
    }

    #[test]
    fn reconstruct_dot_rejects_unsatisfiable_signature_set() {
        let extracted = ExtractedSignatureExpression {
            signatures: vec![
                MomentumSignature {
                    loop_signature: vec![-1, 1],
                    external_signature: vec![],
                },
                MomentumSignature {
                    loop_signature: vec![0, 1],
                    external_signature: vec![],
                },
                MomentumSignature {
                    loop_signature: vec![1, 1],
                    external_signature: vec![],
                },
                MomentumSignature {
                    loop_signature: vec![-1, 0],
                    external_signature: vec![],
                },
            ],
            loop_names: vec!["k1".to_string(), "k2".to_string()],
            external_names: vec![],
            masses: vec!["0".to_string(); 4],
        };

        let error = reconstruct_dot(&extracted, &ReconstructDotOptions::default()).unwrap_err();
        assert!(matches!(error, SignatureParseError::MalformedExpression(_)));
    }

    #[test]
    fn reconstruct_parsed_graph_keeps_mass_attributes_without_dot_roundtrip() {
        let extracted = extract_signatures_and_masses_from_symbolica_expression(
            "prop(k1+p1,mA)*prop(k1+p1-q1,mB)*prop(k1-p2+q2,mC)*prop(k1,mD)",
            "k",
            &["p".to_string(), "q".to_string()],
            "prop(q_,m_)",
        )
        .unwrap();
        let parsed =
            reconstruct_parsed_graph(&extracted, &ReconstructDotOptions::default()).unwrap();

        assert_eq!(parsed.internal_edges.len(), extracted.signatures.len());
        assert_eq!(
            parsed
                .internal_edges
                .iter()
                .map(|edge| edge.mass_key.clone().unwrap_or_else(|| "0".to_string()))
                .collect::<Vec<_>>(),
            extracted.masses
        );
        assert!(crate::validator::validate_parsed_graph(&parsed).ok);
    }

    #[test]
    fn graph_from_signatures_default_outputs_gammaloop_style_dot() {
        let (dot, _) = reconstruct_dot_from_expression(
            "prop(k1,0)*prop(k1+p1,mA)*prop(k1-p2,mB)",
            "k",
            &["p".to_string()],
            "prop(q_,m_)",
            &ReconstructDotOptions::default(),
        )
        .unwrap();

        assert!(dot.contains("num=\"1\";"));
        assert!(dot.contains("edge [num=\"1\"];"));
        assert!(dot.contains("node [num=\"1\"];"));
        assert!(dot.contains("lmb_id=\"0\""));
        assert!(!dot.contains("lmb_rep="));
        assert!(dot.contains("mass=\"mA\""));
    }
}
