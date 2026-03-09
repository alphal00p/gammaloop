//! This module represents a flat graph, i.e. a graph without subgraph.
//! To flatten the graph, existing subgraphs are transformed into their individual nodes: e.g.
//! the statement `a -> { b c }` is transformed into two statements: `a -> b` and `a -> c`.
//! Edge attributes are cloned: therefore, a graph can be flattened only if the attribute type
//! `A` implements `Clone`.
use ahash::HashSet;
use ahash::HashSetExt;
use dot_parser::ast::AList;
pub use dot_parser::ast::AttrList;
pub use dot_parser::ast::AttrStmt;
use dot_parser::ast::CompassPt;
pub use dot_parser::ast::NodeID;
pub use dot_parser::ast::NodeStmt;
use dot_parser::ast::Port;
use dot_parser::canonical::IDEq;
use dot_parser::canonical::Node;
use itertools::Either;

use std::collections::BTreeMap;
use std::ops::AddAssign;

use crate::half_edge::involution::Hedge;

/// A graph, very similar to Ast::Graph, but can not contain subgraphs.
pub struct SubGraphFreeGraph {
    /// Specifies if the `Graph` is strict or not. A "strict" graph must not
    /// contain the same edge multiple times. Notice that, for undirected edge,
    /// an edge from `A` to `B` and an edge from `B` to `A` are equals.
    pub strict: bool,
    /// Specifies if the `Graph` is directed.
    pub is_digraph: bool,
    /// The name of the `Graph`, if any.
    pub name: Option<String>,
    /// The statements that describe the graph.
    pub stmts: StmtList,
}

impl SubGraphFreeGraph {
    /// Returns all `NodeID`s that appear in the graph.
    pub fn get_node_ids(&self) -> HashSet<NodeID> {
        self.stmts.get_node_ids()
    }

    /// Returns all [EdgeStmt]s that appear in the graph.
    /// Notice that it does not recurse into nested subgraphs.
    pub fn get_edges_stmts(self) -> Vec<EdgeStmt> {
        self.stmts.get_edges_stmts()
    }

    #[allow(clippy::type_complexity)]
    pub fn nodes_and_edges(
        self,
    ) -> (
        Vec<dot_parser::ast::AttrStmt<(String, String)>>,
        Vec<IDEq>,
        BTreeMap<String, Node<(String, String)>>,
        EdgeSet,
    ) {
        let mut attrs: Vec<dot_parser::ast::AttrStmt<(String, String)>> = Vec::new();
        let mut nodes = Vec::new();
        let mut edges = Vec::new();
        let mut ideqs = Vec::new();

        for stmt in self.stmts {
            match stmt {
                Stmt::Node(node) => nodes.push(node),
                Stmt::Edge(edge) => edges.push(edge),
                Stmt::Attr(attr) => attrs.push(attr),
                Stmt::IDEq(lhs, rhs) => ideqs.push(IDEq { lhs, rhs }),
            }
        }

        let mut set: BTreeMap<String, Node<(String, String)>> = nodes
            .into_iter()
            .map(|node| (node.node.id.to_string(), node.into()))
            .collect();

        // let mut rhs = stmt.next;

        let edges: EdgeSet = (edges, &mut set).into();

        (attrs, ideqs, set, edges)
    }
}

impl From<dot_parser::ast::Graph<(String, String)>> for SubGraphFreeGraph {
    fn from(g: dot_parser::ast::Graph<(String, String)>) -> Self {
        SubGraphFreeGraph {
            strict: g.strict,
            is_digraph: g.is_digraph,
            name: g.name,
            stmts: g.stmts.into(),
        }
    }
}

/// A list of statements. This corresponds to the `stmt_list` non-terminal of the
/// grammar.
#[derive(Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct StmtList {
    /// The list of statements.
    pub stmts: Vec<Stmt>,
}

impl StmtList {
    fn get_node_ids(&self) -> HashSet<NodeID> {
        let mut hs = HashSet::new();
        for stmt in self {
            hs = hs.union(&stmt.get_node_ids()).cloned().collect();
        }
        hs
    }

    /// Returns a clone of all the EdgeStmt contained in the list.
    fn get_edges_stmts(self) -> Vec<EdgeStmt> {
        let mut v = Vec::new();
        for stmt in self {
            if let Some(edge) = stmt.get_edge() {
                v.push(edge)
            }
        }
        v
    }
}

impl IntoIterator for StmtList {
    type Item = Stmt;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.stmts.into_iter()
    }
}

impl<'a> IntoIterator for &'a StmtList {
    type Item = &'a Stmt;
    type IntoIter = std::slice::Iter<'a, Stmt>;

    fn into_iter(self) -> Self::IntoIter {
        self.stmts.iter()
    }
}

impl FromIterator<Stmt> for StmtList {
    fn from_iter<T>(iter: T) -> Self
    where
        T: IntoIterator<Item = Stmt>,
    {
        Self {
            stmts: iter.into_iter().collect(),
        }
    }
}

impl From<dot_parser::ast::StmtList<(String, String)>> for StmtList {
    fn from(stmts: dot_parser::ast::StmtList<(String, String)>) -> Self {
        let mut v = Vec::new();
        for stmt in stmts {
            match stmt {
                dot_parser::ast::Stmt::NodeStmt(n) => {
                    v.push(Stmt::Node(n));
                }
                dot_parser::ast::Stmt::EdgeStmt(e) => {
                    let edges: Vec<EdgeStmt> = EdgeStmt::from_dot_parser(e);
                    for stmt in edges {
                        v.push(Stmt::Edge(stmt));
                    }
                }
                dot_parser::ast::Stmt::AttrStmt(a) => {
                    v.push(Stmt::Attr(a));
                }
                dot_parser::ast::Stmt::IDEq(s1, s2) => {
                    v.push(Stmt::IDEq(s1, s2));
                }
                dot_parser::ast::Stmt::Subgraph(graph) => {
                    let stmts: StmtList = graph.stmts.into();
                    for stmt in stmts {
                        v.push(stmt)
                    }
                }
            }
        }
        Self { stmts: v }
    }
}

/// A statement of the graph. This corresponds to the `stmt` non-terminal of the
/// grammar.
#[derive(Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub enum Stmt {
    /// A node statement.
    Node(NodeStmt<(String, String)>),
    /// An edge statement.
    Edge(EdgeStmt),
    /// An attribute statement.
    Attr(AttrStmt<(String, String)>),
    /// An alias statement.
    IDEq(String, String),
}

impl Stmt {
    /// Returns true if `self` is a `NodeStmt` variant.
    pub fn is_node_stmt(&self) -> bool {
        matches!(self, Stmt::Node(_))
    }

    /// Returns `Some(&node)` if `&self` if a `&NodeStmt(node)`, and `None`
    /// otherwise.
    pub fn get_node_ref(&self) -> Option<&NodeStmt<(String, String)>> {
        if let Stmt::Node(node) = self {
            Some(node)
        } else {
            None
        }
    }

    /// Returns `Some(node)` if `self` if a `NodeStmt(node)`, and `None`
    /// otherwise.
    pub fn get_node(self) -> Option<NodeStmt<(String, String)>> {
        if let Stmt::Node(node) = self {
            Some(node)
        } else {
            None
        }
    }

    /// Returns true if `self` is a `EdgeStmt` variant.
    pub fn is_edge_stmt(&self) -> bool {
        matches!(self, Stmt::Edge(_))
    }

    /// Returns `Some(&edge)` if `&self` if a `&EdgeStmt(edge)`, and `None`
    /// otherwise.
    pub fn get_edge_ref(&self) -> Option<&EdgeStmt> {
        if let Stmt::Edge(edge) = self {
            Some(edge)
        } else {
            None
        }
    }

    /// Returns `Some(edge)` if `self` if a `EdgeStmt(edge)`, and `None`
    /// otherwise.
    pub fn get_edge(self) -> Option<EdgeStmt> {
        if let Stmt::Edge(edge) = self {
            Some(edge)
        } else {
            None
        }
    }

    /// Returns true if `self` is a `AttrStmt` variant.
    pub fn is_attr_stmt(&self) -> bool {
        matches!(self, Stmt::Attr(_))
    }

    /// Returns `Some(&attr)` if `&self` if a `&AttrStmt(attr)`, and `None`
    /// otherwise.
    pub fn get_attr_ref(&self) -> Option<&AttrStmt<(String, String)>> {
        if let Stmt::Attr(attr) = self {
            Some(attr)
        } else {
            None
        }
    }

    /// Returns `Some(attr)` if `self` if a `AttrStmt(attr)`, and `None`
    /// otherwise.
    pub fn get_attr(self) -> Option<AttrStmt<(String, String)>> {
        if let Stmt::Attr(attr) = self {
            Some(attr)
        } else {
            None
        }
    }

    /// Returns true if `self` is a `IDEq` variant.
    pub fn is_ideq_stmt(&self) -> bool {
        matches!(self, Stmt::IDEq(..))
    }

    /// Returns `Some((&id1, &id2))` if `&self` if a `&IDEq(id1, id2)` and `None`
    /// otherwise.
    pub fn get_ideq_ref(&self) -> Option<(&str, &str)> {
        if let Stmt::IDEq(id1, id2) = self {
            Some((id1, id2))
        } else {
            None
        }
    }

    /// Returns all `NodeID`s that appear in this statement.
    fn get_node_ids(&self) -> HashSet<NodeID> {
        match self {
            Stmt::Edge(e) => e.get_node_ids(),
            Stmt::Node(n) => HashSet::from_iter([n.node.clone()]),
            _ => HashSet::new(),
        }
    }
}

/// The description of an edge. This corresponds to the `edge_stmt` non-terminal
/// of the grammar.
#[derive(Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct EdgeStmt {
    /// The origin of the edge.
    pub from: NodeID,
    /// The destination of the edge.
    pub next: EdgeRHS,
    /// The attributes of the edge.
    pub attr: Option<AttrList<(String, String)>>,
}

impl EdgeStmt {
    fn get_node_ids(&self) -> HashSet<NodeID> {
        let mut nexts = self.next.get_node_ids();
        nexts.insert(self.from.clone());
        nexts
    }

    fn from_dot_parser(edge: dot_parser::ast::EdgeStmt<(String, String)>) -> Vec<Self> {
        let edges = edge.flatten();
        let mut v: Vec<EdgeStmt> = Vec::new();

        for edge in edges {
            // `edge` has just one destination (`edge.next.next` is `None`), as it comes from a
            // "flattened" set of edges.
            let from = edge.from;
            let to = edge.next.to;
            let attr = edge.attr;

            let (from_ids, mut extra_edges_from) = match from {
                Either::Left(node_from) => (HashSet::from_iter([node_from]), Vec::new()),
                Either::Right(subgraph) => {
                    let g = SubGraphFreeGraph {
                        strict: false,
                        is_digraph: false,
                        name: subgraph.id,
                        stmts: subgraph.stmts.into(),
                    };
                    (g.get_node_ids(), g.get_edges_stmts())
                }
            };

            let (to_ids, mut extra_edges_to) = match to {
                Either::Left(node_from) => (HashSet::from_iter([node_from]), Vec::new()),
                Either::Right(subgraph) => {
                    let g = SubGraphFreeGraph {
                        strict: false,
                        is_digraph: false,
                        name: subgraph.id,
                        stmts: subgraph.stmts.into(),
                    };
                    (g.get_node_ids(), g.get_edges_stmts())
                }
            };

            for from in from_ids {
                for to in &to_ids {
                    v.push(EdgeStmt {
                        from: from.clone(),
                        next: EdgeRHS {
                            to: to.clone(),
                            next: None,
                        },
                        attr: attr.clone(),
                    });
                }
            }
            v.append(&mut extra_edges_from);
            v.append(&mut extra_edges_to);
        }
        v
    }
}

/// The Right hand side of an edge description. This corresponds to the
/// `EdgeRHS` non-terminal of the grammar.
/// Notice that the grammar allows multiple EdgeRHS in sequence, to chain edges:
/// `A -> B -> C`.
/// Notice that, contrary to Ast::EdgeRHS, Ast::FlatGraph::EdgeRHS can not contain a subgraph.
#[derive(Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct EdgeRHS {
    /// The identifier of the destination of the edge.
    pub to: NodeID,
    /// A possible chained RHS.
    pub next: Option<Box<EdgeRHS>>,
}

impl EdgeRHS {
    fn get_node_ids(&self) -> HashSet<NodeID> {
        let mut nexts: HashSet<NodeID> = self
            .next
            .as_ref()
            .map(|n| n.get_node_ids())
            .unwrap_or_default();
        nexts.insert(self.to.clone());
        nexts
    }
}

/// A set of `Edge`s.
#[derive(Debug, Clone)]
pub struct EdgeSet {
    /// `Edge`s of the set.
    pub set: Vec<Edge>,
}

impl EdgeSet {
    fn empty() -> Self {
        Self { set: Vec::new() }
    }
}

impl AddAssign for EdgeSet {
    fn add_assign(&mut self, mut rhs: Self) {
        self.set.append(&mut rhs.set)
    }
}

/// A single edge of the graph.
// TODO: replace Strings by a reference to the relevant NodeID
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct Edge {
    /// The name of the origin of the edge.
    pub from: NodeID,
    /// The name of the destination of the edge.
    pub to: NodeID,
    /// A list of attributes that apply to this specific edge.
    pub attr: Vec<(String, String)>,
}

pub trait PortExt {
    fn hedge(&self) -> Option<Hedge>;
    fn compass(&self) -> Option<&CompassPt>;
    fn id(&self) -> Option<&str>;
}

impl PortExt for Port {
    fn hedge(&self) -> Option<Hedge> {
        match self {
            Port::ID(id, _) => id.parse().ok(),
            _ => None,
        }
    }

    fn id(&self) -> Option<&str> {
        match self {
            Port::ID(id, _) => Some(id.as_str()),
            _ => None,
        }
    }

    fn compass(&self) -> Option<&CompassPt> {
        match self {
            Port::Compass(c) => Some(c),
            Port::ID(_, c) => c.as_ref(),
        }
    }
}

impl Edge {
    pub fn source_port(&self) -> Option<&Port> {
        self.from.port.as_ref()
    }

    pub fn sink_port(&self) -> Option<&Port> {
        self.to.port.as_ref()
    }
}

impl From<(EdgeStmt, &mut BTreeMap<String, Node<(String, String)>>)> for EdgeSet {
    fn from(tuple: (EdgeStmt, &mut BTreeMap<String, Node<(String, String)>>)) -> Self {
        let (stmt, nodes) = tuple;
        let mut from = stmt.from;
        let mut rhs = stmt.next;
        let mut set = Vec::new();
        let attr = stmt
            .attr
            .map(|list| AList::from(list).elems)
            .unwrap_or(Vec::new());

        loop {
            let to = rhs.to;
            let from_id = from.id.clone();
            let to_id = to.id.clone();

            if nodes.get(&from_id).is_none() {
                nodes.insert(from_id.clone(), (&from).into());
            }

            if nodes.get(&to_id).is_none() {
                nodes.insert(to_id.clone(), (&to).into());
            }

            let edge = Edge {
                from,
                to: to.clone(),
                attr: attr.clone(),
            };

            set.push(edge);

            if rhs.next.is_none() {
                return EdgeSet { set };
            }
            from = to;
            rhs = *(rhs.next.unwrap());
        }
    }
}

impl<I> From<(I, &mut BTreeMap<String, Node<(String, String)>>)> for EdgeSet
where
    I: IntoIterator<Item = EdgeStmt>,
{
    fn from(tuple: (I, &mut BTreeMap<String, Node<(String, String)>>)) -> Self {
        let (stmts, nodes) = tuple;
        let mut set = EdgeSet::empty();
        for stmt in stmts {
            set += (stmt, &mut *nodes).into();
        }
        set
    }
}
