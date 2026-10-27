//! This module represents a flat graph, i.e. a graph without subgraph.
//! To flatten the graph, existing subgraphs are transformed into their individual nodes: e.g.
//! the statement `a -> { b c }` is transformed into two statements: `a -> b` and `a -> c`.
//! Edge attributes are cloned: therefore, a graph can be flattened only if the attribute type
//! `A` implements `Clone`.
pub use crate::ast::AttrList;
pub use crate::ast::AttrStmt;
pub use crate::ast::NodeID;
pub use crate::ast::NodeStmt;
pub use crate::ast::Port;
use either::Either;
use std::collections::HashSet;

/// A graph, very similar to Ast::Graph, but can not contain subgraphs.
pub struct Graph<A> {
    /// Specifies if the `Graph` is strict or not. A "strict" graph must not
    /// contain the same edge multiple times. Notice that, for undirected edge,
    /// an edge from `A` to `B` and an edge from `B` to `A` are equals.
    pub strict: bool,
    /// Specifies if the `Graph` is directed.
    pub is_digraph: bool,
    /// The name of the `Graph`, if any.
    pub name: Option<String>,
    /// The statements that describe the graph.
    pub stmts: StmtList<A>,
}

impl<A> Graph<A> {
    /// Filter and map attributes. The main intended usage of this function is
    /// to convert attributes as `&'a str` into an enum, e.g.
    /// to convert `["label"="whatever", "color"="foo"]` into
    /// `[Attr::Label(whatever), Attr::Color(foo)]`.
    ///
    /// To take into account non-standard attributes, the `Attr` enum has to be
    /// provided by the user.
    pub fn filter_map<B>(self, f: &dyn Fn(A) -> Option<B>) -> Graph<B> {
        let new_stmts: StmtList<B> = self.stmts.filter_map_attr(f);
        Graph {
            strict: self.strict,
            is_digraph: self.is_digraph,
            name: self.name,
            stmts: new_stmts,
        }
    }

    /// Returns all `NodeID`s that appear in the graph.
    pub fn get_node_ids(&self) -> HashSet<NodeID> {
        self.stmts.get_node_ids()
    }

    /// Returns all [EdgeStmt]s that appear in the graph.
    /// Notice that it does not recurse into nested subgraphs.
    pub fn get_edges_stmts(self) -> Vec<EdgeStmt<A>> {
        self.stmts.get_edges_stmts()
    }
}

impl<A> From<crate::ast::Graph<A>> for Graph<A>
where
    A: Clone,
{
    fn from(g: crate::ast::Graph<A>) -> Self {
        Graph {
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
pub struct StmtList<A> {
    /// The list of statements.
    pub stmts: Vec<Stmt<A>>,
}

impl<'a, A> StmtList<A> {
    fn filter_map_attr<B>(self, f: &'a dyn Fn(A) -> Option<B>) -> StmtList<B> {
        self.stmts
            .into_iter()
            .map(|stmt| stmt.filter_map_attr(f))
            .collect()
    }

    fn get_node_ids(&self) -> HashSet<NodeID> {
        let mut hs = HashSet::new();
        for stmt in self {
            hs = hs.union(&stmt.get_node_ids()).cloned().collect();
        }
        hs
    }

    /// Returns a clone of all the EdgeStmt contained in the list.
    fn get_edges_stmts(self) -> Vec<EdgeStmt<A>> {
        let mut v = Vec::new();
        for stmt in self {
            if let Some(edge) = stmt.get_edge() {
                v.push(edge)
            }
        }
        v
    }
}

impl<A> IntoIterator for StmtList<A> {
    type Item = Stmt<A>;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.stmts.into_iter()
    }
}

impl<'a, A> IntoIterator for &'a StmtList<A> {
    type Item = &'a Stmt<A>;
    type IntoIter = std::slice::Iter<'a, Stmt<A>>;

    fn into_iter(self) -> Self::IntoIter {
        self.stmts.iter()
    }
}

impl<A> FromIterator<Stmt<A>> for StmtList<A> {
    fn from_iter<T>(iter: T) -> Self
    where
        T: IntoIterator<Item = Stmt<A>>,
    {
        Self {
            stmts: iter.into_iter().collect(),
        }
    }
}

impl<A> From<crate::ast::StmtList<A>> for StmtList<A>
where
    A: Clone,
{
    fn from(stmts: crate::ast::StmtList<A>) -> Self {
        let mut v = Vec::new();
        for stmt in stmts {
            match stmt {
                crate::ast::Stmt::NodeStmt(n) => {
                    v.push(Stmt::NodeStmt(n));
                }
                crate::ast::Stmt::EdgeStmt(e) => {
                    let edges: Vec<EdgeStmt<A>> = e.into();
                    for stmt in edges {
                        v.push(Stmt::EdgeStmt(stmt));
                    }
                }
                crate::ast::Stmt::AttrStmt(a) => {
                    v.push(Stmt::AttrStmt(a));
                }
                crate::ast::Stmt::IDEq(s1, s2) => {
                    v.push(Stmt::IDEq(s1, s2));
                }
                crate::ast::Stmt::Subgraph(graph) => {
                    let stmts: StmtList<A> = graph.stmts.into();
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
pub enum Stmt<A> {
    /// A node statement.
    NodeStmt(NodeStmt<A>),
    /// An edge statement.
    EdgeStmt(EdgeStmt<A>),
    /// An attribute statement.
    AttrStmt(AttrStmt<A>),
    /// An alias statement.
    IDEq(String, String),
}

impl<'a, A> Stmt<A> {
    /// Convert a statement with attributes of type `A` into a statement with
    /// attributes of type `B`.
    fn filter_map_attr<B>(self, f: &'a dyn Fn(A) -> Option<B>) -> Stmt<B> {
        match self {
            Stmt::NodeStmt(node) => Stmt::NodeStmt(node.filter_map_attr(f)),
            Stmt::EdgeStmt(edge) => Stmt::EdgeStmt(edge.filter_map_attr(f)),
            Stmt::AttrStmt(attr) => Stmt::AttrStmt(attr.filter_map_attr(f)),
            Stmt::IDEq(a, b) => Stmt::IDEq(a, b),
        }
    }

    /// Returns true if `self` is a `NodeStmt` variant.
    pub fn is_node_stmt(&self) -> bool {
        matches!(self, Stmt::NodeStmt(_))
    }

    /// Returns `Some(&node)` if `&self` if a `&NodeStmt(node)`, and `None`
    /// otherwise.
    pub fn get_node_ref(&self) -> Option<&NodeStmt<A>> {
        if let Stmt::NodeStmt(node) = self {
            Some(node)
        } else {
            None
        }
    }

    /// Returns `Some(node)` if `self` if a `NodeStmt(node)`, and `None`
    /// otherwise.
    pub fn get_node(self) -> Option<NodeStmt<A>> {
        if let Stmt::NodeStmt(node) = self {
            Some(node)
        } else {
            None
        }
    }

    /// Returns true if `self` is a `EdgeStmt` variant.
    pub fn is_edge_stmt(&self) -> bool {
        matches!(self, Stmt::EdgeStmt(_))
    }

    /// Returns `Some(&edge)` if `&self` if a `&EdgeStmt(edge)`, and `None`
    /// otherwise.
    pub fn get_edge_ref(&self) -> Option<&EdgeStmt<A>> {
        if let Stmt::EdgeStmt(edge) = self {
            Some(edge)
        } else {
            None
        }
    }

    /// Returns `Some(edge)` if `self` if a `EdgeStmt(edge)`, and `None`
    /// otherwise.
    pub fn get_edge(self) -> Option<EdgeStmt<A>> {
        if let Stmt::EdgeStmt(edge) = self {
            Some(edge)
        } else {
            None
        }
    }

    /// Returns true if `self` is a `AttrStmt` variant.
    pub fn is_attr_stmt(&self) -> bool {
        matches!(self, Stmt::AttrStmt(_))
    }

    /// Returns `Some(&attr)` if `&self` if a `&AttrStmt(attr)`, and `None`
    /// otherwise.
    pub fn get_attr_ref(&self) -> Option<&AttrStmt<A>> {
        if let Stmt::AttrStmt(attr) = self {
            Some(attr)
        } else {
            None
        }
    }

    /// Returns `Some(attr)` if `self` if a `AttrStmt(attr)`, and `None`
    /// otherwise.
    pub fn get_attr(self) -> Option<AttrStmt<A>> {
        if let Stmt::AttrStmt(attr) = self {
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
            Stmt::EdgeStmt(e) => e.get_node_ids(),
            Stmt::NodeStmt(n) => HashSet::from_iter([n.get_node_id().clone()]),
            _ => HashSet::new(),
        }
    }
}

/// The description of an edge. This corresponds to the `edge_stmt` non-terminal
/// of the grammar.
#[derive(Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct EdgeStmt<A> {
    /// The origin of the edge.
    pub from: NodeID,
    /// The destination of the edge.
    pub next: EdgeRHS,
    /// The attributes of the edge.
    pub attr: Option<AttrList<A>>,
}

impl<A> EdgeStmt<A> {
    fn filter_map_attr<B>(self, f: &dyn Fn(A) -> Option<B>) -> EdgeStmt<B> {
        EdgeStmt {
            from: self.from,
            next: self.next,
            attr: self.attr.map(|a| a.filter_map_attr(f)),
        }
    }

    fn get_node_ids(&self) -> HashSet<NodeID> {
        let mut nexts = self.next.get_node_ids();
        nexts.insert(self.from.clone());
        nexts
    }
}

impl<A> From<crate::ast::EdgeStmt<A>> for Vec<EdgeStmt<A>>
where
    A: Clone,
{
    fn from(edge: crate::ast::EdgeStmt<A>) -> Vec<EdgeStmt<A>> {
        let edges = edge.flatten();
        let mut v: Vec<EdgeStmt<A>> = Vec::new();

        for edge in edges {
            // `edge` has just one destination (`edge.next.next` is `None`), as it comes from a
            // "flattened" set of edges.
            let from = edge.from;
            let to = edge.next.to;
            let attr = edge.attr;

            let (from_ids, mut extra_edges_from) = match from {
                Either::Left(node_from) => (HashSet::from_iter([node_from]), Vec::new()),
                Either::Right(subgraph) => {
                    let g: Graph<A> = subgraph.into_graph(false, false).into();
                    (g.get_node_ids(), g.get_edges_stmts())
                }
            };

            let (to_ids, mut extra_edges_to) = match to {
                Either::Left(node_from) => (HashSet::from_iter([node_from]), Vec::new()),
                Either::Right(subgraph) => {
                    let g: Graph<A> = subgraph.into_graph(false, false).into();
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
