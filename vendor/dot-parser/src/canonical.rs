//! This modules implement "Canonical graphs". This is motivated by the fact
//! that graphs parsed naively may be difficult to work with, from a programming
//! viewpoint: typically, `attr_list` in the grammar are defined as
//! "[ID = ID, ...][ID = ID, ...]". Therefore, we may want to flatten such
//! structures. This is done in the [`canonical`][crate::canonical] module.
//!
//! The main structure of the module is [`Graph`], which implements such a flatten
//! graph.

use std::collections::HashMap;
use std::ops::AddAssign;

#[cfg(feature = "display")]
use std::fmt::{Display, Formatter};

use crate::subgraph_free::{
    AttrStmt as AstAttrStmt, EdgeStmt, Graph as SFGraph, NodeID, NodeStmt, Port, Stmt,
};

pub use crate::ast::AList;

/// A `Graph` is a structure that can be created from a regular
/// `Graph`, but that is more friendly to work with. For instance, in a `Graph`,
/// attributes are most often given as a list of list of `Attr`, while in a
/// `Graph`, the lists are flatten.
#[derive(Debug, Clone)]
pub struct Graph<A> {
    /// Specifies if the `Graph` is strict or not. A "strict" graph must not
    /// contain the same edge multiple times. Notice that, for undirected edge,
    /// an edge from `A` to `B` and an edge from `B` to `A` are equals.
    pub strict: bool,
    /// Specifies if the `Graph` is directed.
    pub is_digraph: bool,
    /// The name of the `Graph`, if any.
    pub name: Option<String>,
    /// The global attributes of the graph.
    pub attr: Vec<AttrStmt<A>>,
    /// The nodes of the graph.
    pub nodes: NodeSet<A>,
    /// The edges of the graph.
    pub edges: EdgeSet<A>,
    /// The ID equalities declared in the graph.
    pub ideqs: Vec<IDEq>,
}

#[cfg(feature = "display")]
impl<A> Display for Graph<A>
where
    A: Display,
{
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), std::fmt::Error> {
        if self.strict {
            write!(f, "strict ")?;
        }
        if self.is_digraph {
            write!(f, "digraph ")?;
        } else {
            write!(f, "graph ")?;
        }
        if let Some(name) = &self.name {
            write!(f, "{}", name)?;
        }

        writeln!(f, "{{")?;
        for stmt in &self.attr {
            writeln!(f, "{}", stmt)?;
        }

        for ideq in &self.ideqs {
            writeln!(f, "{}", ideq)?;
        }

        write!(f, "{}", self.nodes)?;
        write!(f, "{}", self.edges)?;
        writeln!(f, "}}")?;

        Result::Ok(())
    }
}

impl<A> Graph<A> {
    /// Filter and map attributes. The main intended usage of this function is
    /// to convert attributes as `&'a str` into an enum, e.g.
    /// to convert `["label"="whatever", "color"="foo"]` into
    /// `[Attr::Label(whatever), Attr::Color(foo)]`.
    ///
    /// To take into account non-standard attributes, the `Attr` enum has to be
    /// provided by the user.
    pub fn filter_map<F, B>(self, f: F) -> Graph<B>
    where
        F: Fn(A) -> Option<B>,
    {
        let new_attr = self
            .attr
            .into_iter()
            .filter_map(|attr_stmt| attr_stmt.filter_map(&f))
            .collect();
        let new_nodes = self.nodes.map(&f);
        let new_edges = self.edges.map(&f);

        Graph {
            strict: self.strict,
            is_digraph: self.is_digraph,
            name: self.name,
            attr: new_attr,
            nodes: new_nodes,
            edges: new_edges,
            ideqs: self.ideqs,
        }
    }
}

impl<A, G> From<G> for Graph<A>
where
    G: Into<SFGraph<A>>,
    A: Clone,
{
    fn from(graph: G) -> Self {
        let graph = graph.into();
        let mut attrs: Vec<AstAttrStmt<_>> = Vec::new();
        let mut nodes = Vec::new();
        let mut edges = Vec::new();
        let mut ideqs = Vec::new();

        for stmt in graph.stmts {
            match stmt {
                Stmt::NodeStmt(node) => nodes.push(node),
                Stmt::EdgeStmt(edge) => edges.push(edge),
                Stmt::AttrStmt(attr) => attrs.push(attr),
                Stmt::IDEq(lhs, rhs) => ideqs.push(IDEq {lhs, rhs}),
            }
        }

        let mut nodes: NodeSet<A> = nodes.into();
        let edges: EdgeSet<A> = (edges, &mut nodes).into();
        let attr: Vec<AttrStmt<A>> = attrs
            .into_iter()
            .flat_map(|stmt: AstAttrStmt<_>| {
                let attr_l: Vec<AttrStmt<A>> = stmt.into();
                attr_l
            })
            .collect();

        Graph {
            strict: graph.strict,
            is_digraph: graph.is_digraph,
            name: graph.name,
            attr,
            nodes,
            edges,
            ideqs,
        }
    }
}

/// A single node of the graph.
#[derive(Debug, Clone)]
pub struct Node<A> {
    /// The identifier of the node.
    pub id: String,
    /// The port of the node.
    pub port: Option<Port>,
    /// The attributes that apply to this node.
    pub attr: AList<A>,
}

impl<A> From<NodeStmt<A>> for Node<A> {
    fn from(stmt: NodeStmt<A>) -> Self {
        Node {
            id: stmt.node.id.to_string(),
            port: stmt.node.port,
            attr: stmt.attr.map(|list| list.into()).unwrap_or(AList::empty()),
        }
    }
}

impl<A> From<&NodeID> for Node<A> {
    fn from(node: &NodeID) -> Self {
        Node {
            id: node.id.to_string(),
            port: node.port.clone(),
            attr: AList::empty(),
        }
    }
}

impl<A> Node<A> {
    fn map<F, B>(self, f: F) -> Node<B>
    where
        F: Fn(A) -> Option<B>,
    {
        Node {
            id: self.id,
            port: self.port,
            attr: self.attr.filter_map_attr(&f),
        }
    }
}

#[cfg(feature = "display")]
impl<A> Display for Node<A>
where
    A: Display,
{
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), std::fmt::Error> {
        write!(f, "\"{}\" ", self.id)?;
        if let Some(port) = &self.port {
            write!(f, "{}", port)?;
        }
        if !self.attr.is_empty() {
            write!(f, "[{}]", self.attr)?;
        }

        Ok(())
    }
}

/// A set of `Node`s.
#[derive(Debug, Clone)]
pub struct NodeSet<A> {
    /// The set of nodes in the NodeSet. They are indexed by their identifier.
    /// Note that this field being public is experimental.
    pub set: HashMap<String, Node<A>>,
}

impl<A> NodeSet<A> {
    fn insert_if_absent(&mut self, id: String, or: Node<A>) {
        // TODO: clarify what happens if id != or.id
        if self.set.get(&id).is_none() {
            self.set.insert(id, or);
        }
    }

    fn map<F, B>(self, f: F) -> NodeSet<B>
    where
        F: Fn(A) -> Option<B>,
    {
        let new_set = self
            .set
            .into_iter()
            .map(|(name, node)| (name, node.map(&f)))
            .collect();
        NodeSet { set: new_set }
    }
}

impl<A, I> From<I> for NodeSet<A>
where
    I: IntoIterator<Item = NodeStmt<A>>,
{
    fn from(nodes: I) -> Self {
        let set: HashMap<_, _> = nodes
            .into_iter()
            .map(|node| (node.node.id.to_string(), node.into()))
            .collect();
        NodeSet { set }
    }
}

#[cfg(feature = "display")]
impl<A> Display for NodeSet<A>
where
    A: Display,
{
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), std::fmt::Error> {
        for node in self.set.values() {
            writeln!(f, "{}", node)?;
        }

        Ok(())
    }
}

/// A set of `Edge`s.
#[derive(Debug, Clone)]
pub struct EdgeSet<A> {
    /// `Edge`s of the set.
    pub set: Vec<Edge<A>>,
}

impl<A> EdgeSet<A> {
    fn empty() -> Self {
        Self { set: Vec::new() }
    }

    fn map<F, B>(self, f: F) -> EdgeSet<B>
    where
        F: Fn(A) -> Option<B>,
    {
        let new_set = self.set.into_iter().map(|edge| edge.map(&f)).collect();
        EdgeSet { set: new_set }
    }
}

impl<A> AddAssign for EdgeSet<A> {
    fn add_assign(&mut self, mut rhs: Self) {
        self.set.append(&mut rhs.set)
    }
}

#[cfg(feature = "display")]
impl<A> Display for EdgeSet<A>
where
    A: Display,
{
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), std::fmt::Error> {
        for edge in &self.set {
            writeln!(f, "{}", edge)?;
        }
        Ok(())
    }
}

/// A single edge of the graph.
// TODO: replace Strings by a reference to the relevant NodeID
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct Edge<A> {
    /// The name of the origin of the edge.
    pub from: String,
    /// The name of the destination of the edge.
    pub to: String,
    /// A list of attributes that apply to this specific edge.
    pub attr: AList<A>,
}

#[cfg(feature = "display")]
impl<A> Display for Edge<A>
where
    A: Display,
{
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), std::fmt::Error> {
        write!(f, "\"{}\" -> \"{}\" [{}]", self.from, self.to, self.attr)
    }
}

impl<A> Edge<A> {
    fn map<F, B>(self, f: F) -> Edge<B>
    where
        F: Fn(A) -> Option<B>,
    {
        Edge {
            from: self.from,
            to: self.to,
            attr: self.attr.filter_map_attr(&f),
        }
    }
}

impl<A> From<(EdgeStmt<A>, &mut NodeSet<A>)> for EdgeSet<A>
where
    A: Clone,
{
    fn from(tuple: (EdgeStmt<A>, &mut NodeSet<A>)) -> Self {
        let (stmt, nodes) = tuple;
        let mut from = stmt.from;
        let mut rhs = stmt.next;
        let mut set = Vec::new();
        let attr = stmt.attr.map(|list| list.into()).unwrap_or(AList::empty());

        loop {
            let to = rhs.to;
            let from_id = from.id.clone();
            let to_id = to.id.clone();

            nodes.insert_if_absent(from.id.to_string(), (&from).into());
            nodes.insert_if_absent(to.id.to_string(), (&to).into());

            let edge = Edge {
                from: from_id.to_string(),
                to: to_id.to_string(),
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

impl<A, I> From<(I, &mut NodeSet<A>)> for EdgeSet<A>
where
    I: IntoIterator<Item = EdgeStmt<A>>,
    A: Clone,
{
    fn from(tuple: (I, &mut NodeSet<A>)) -> Self {
        let (stmts, nodes) = tuple;
        let mut set = EdgeSet::empty();
        for stmt in stmts {
            set += (stmt, &mut *nodes).into();
        }
        set
    }
}

#[derive(Debug, Copy, Clone)]
/// An `AttrStmt`, i.e. a statement that applies to either the whole graph, all
/// edges, or all nodes. Note that, in a canonical graph, `AttrStmt`s contain a
/// single statement.
pub enum AttrStmt<A> {
    /// An `AttrStmt` that applies to the whole graph.
    Graph(A),
    /// An `AttrStmt` that applies to all nodes of the graph.
    Node(A),
    /// An `AttrStmt` that applies to all edges of the graph.
    Edge(A),
}

impl<A> AttrStmt<A> {
    fn filter_map<F, B>(self, f: F) -> Option<AttrStmt<B>>
    where
        F: Fn(A) -> Option<B>,
    {
        match self {
            AttrStmt::Graph(a) => {
                f(a).map(AttrStmt::Graph)
            }
            AttrStmt::Node(a) => {
                f(a).map(AttrStmt::Node)
            }
            AttrStmt::Edge(a) => {
                f(a).map(AttrStmt::Edge)
            }
        }
    }
}

impl<A> From<AstAttrStmt<A>> for Vec<AttrStmt<A>> {
    fn from(val: AstAttrStmt<A>) -> Self {
        match val {
            AstAttrStmt::Graph(list) => {
                let alist: AList<A> = list.into();
                alist
                    .into_iter()
                    .map(|attr| AttrStmt::Graph(attr))
                    .collect()
            }
            AstAttrStmt::Node(list) => {
                let alist: AList<A> = list.into();
                alist
                    .into_iter()
                    .map(|attr| AttrStmt::Node(attr))
                    .collect()
            }
            AstAttrStmt::Edge(list) => {
                let alist: AList<A> = list.into();
                alist
                    .into_iter()
                    .map(|attr| AttrStmt::Edge(attr))
                    .collect()
            }
        }
    }
}

#[cfg(feature = "display")]
impl<A> Display for AttrStmt<A>
where
    A: Display,
{
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), std::fmt::Error> {
        match self {
            AttrStmt::Graph(attr) => write!(f, "graph [{}]", attr),
            AttrStmt::Edge(attr) => write!(f, "edge [{}]", attr),
            AttrStmt::Node(attr) => write!(f, "node [{}]", attr),
        }
    }
}

#[derive(Clone, Debug)]
/// An identifier equality, i.e. a statement that has the form `ID '=' ID`.
pub struct IDEq {
    /// The left hand side ID,
    pub lhs: String,
    /// The right hand side ID,
    pub rhs: String,
}

#[cfg(feature = "display")]
impl Display for IDEq { 
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), std::fmt::Error> {
        write!(f, "{} = {}", self.lhs, self.rhs)
    }
}
