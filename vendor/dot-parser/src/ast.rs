//! This module implements an Abstract Syntax Tree for Dot graphs. The main
//! structure is `Graph`, which corresponds to the `graph` non-terminal in the
//! grammar.

use pest::iterators::Pair;
use pest::Parser;

pub use either;
use either::Either;
use std::collections::HashSet;
use std::error::Error;
use std::fmt::{Display, Formatter};
use std::iter::FromIterator;
use std::path::Path;

#[cfg(feature = "to_tokens")]
use quote::{ToTokens, quote, TokenStreamExt};
#[cfg(feature = "to_tokens")]
use proc_macro2::*;

mod parser {
    use pest_derive::Parser;

    #[derive(Parser)]
    #[grammar = "parser/dot.pest"]
    pub struct DotParser;
}

use self::parser::DotParser;
use self::parser::Rule;

/// Type for errors that occur when parsing.
pub type PestError = pest::error::Error<Rule>;
/// Type for I/O related errors. Those may occur when reading a file for parsing.
pub type IOError = std::io::Error;

/// This enum type contains errors that can occur when using a DotParser. In principle, those
/// errors should never occur, as all parsing errors should be caught by DotParser::parse.
/// Therefore, if such error occurs, it is a bug, probably a missing feature.
#[derive(Debug, Clone)]
pub enum ParseError<'a> {
    /// This variant represents the case where we expect one of several rules, but we actually find
    /// another.
    ExpectRule {
        /// The list of accepted `Rule`s.
        expect: Vec<Rule>,
        /// The `Rule` actually found.
        found: Rule,
    },
    /// This variant represents the case where we expect a `Pair` but none is present.
    MissingPair {
        /// The parent pair, i.e. the one that terminates too early.
        parent: Pair<'a, Rule>,
        /// The missing pair should have one of the `Rule` in `expect`.
        expect: Vec<Rule>,
    },
}

impl<'a> ParseError<'a> {
    fn expect_rule(expect: Vec<Rule>, found: Rule) -> Self {
        ParseError::ExpectRule { expect, found }
    }

    fn missing_pair(parent: Pair<'a, Rule>, expect: Vec<Rule>) -> Self {
        ParseError::MissingPair { parent, expect }
    }
}

impl Display for ParseError<'_> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), std::fmt::Error> {
        match self {
            ParseError::ExpectRule { expect, found } => {
                let expected = expect
                    .iter()
                    .fold(String::new(), |acc, r| format!("{}\"{:?}\", ", acc, r));
                write!(
                    f,
                    "Expect one rule of {}but \"{:?}\" found.",
                    expected, found
                )?;
            }
            ParseError::MissingPair { parent, expect } => {
                let mut expected = expect
                    .iter()
                    .fold(String::new(), |acc, r| format!("{}\"{:?}\", ", acc, r));
                // Removing the final ", "
                expected.pop();
                expected.pop();
                write!(
                    f,
                    "The Pair:\n{}\nterminates early. Expected one of {}.",
                    parent.as_str(),
                    expected
                )?;
            }
        }
        Ok(())
    }
}

impl Error for ParseError<'_> {}

/// An error that can occur when reading from a file.
#[derive(Debug)]
pub enum GraphFromFileError<'a> {
    /// The error occured when manipulating the file (e.g. the file does not exist).
    FileError(IOError),
    /// The error occured when parsing the file (e.g. the pest parser returned an error).
    PestParseError(PestError),
    /// The error occured when traversing the `Rule`s tree returned by the pest parser.
    /// Such error occuring is likely a bug or a missing feature of the library.
    ParseError(ParseError<'a>),
}

impl Display for GraphFromFileError<'_> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), std::fmt::Error> {
        match self {
            Self::FileError(e) => write!(f, "{}", e),
            Self::PestParseError(e) => write!(f, "{}", e),
            Self::ParseError(e) => write!(f, "{}", e),
        }
    }
}

impl Error for GraphFromFileError<'_> {}

impl From<IOError> for GraphFromFileError<'_> {
    fn from(e: IOError) -> Self {
        Self::FileError(e)
    }
}

impl From<PestError> for GraphFromFileError<'_> {
    fn from(e: PestError) -> Self {
        Self::PestParseError(e)
    }
}

impl<'a> From<ParseError<'a>> for GraphFromFileError<'a> {
    fn from(e: ParseError<'a>) -> Self {
        Self::ParseError(e)
    }
}

#[derive(Debug, Eq, PartialEq, Ord, PartialOrd, Hash, Clone)]
/// This structure is an AST for the DOT graph language.  The generic `A`
/// denotes the type of attributes. By default (and when parsing), it is
/// `(&'a str, &'a str)`, i.e. two strings, one for the key and one for the
/// value of the attribute. The library provides functions to map from one type
/// to an other.
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

#[cfg(feature = "to_tokens")]
impl ToTokens for Graph<(String, String)> {
    fn to_tokens(&self, ts: &mut proc_macro2::TokenStream) {
        let name = match &self.name {
            Some(name) => quote! { std::option::Option::Some( std::string::ToString::to_string(#name) ) },
            None => quote! { std::option::Option::None },
        };
        let strict = self.strict;
        let is_digraph = self.is_digraph;
        let stmts = &self.stmts;
        let tokens = quote! {
            dot_parser::ast::Graph::<(&'static str, &'static str)> {
                strict: #strict,
                is_digraph: #is_digraph,
                name: #name,
                stmts: #stmts
            }
        };
        ts.append_all(tokens);
    }
}

impl<'a, A> TryFrom<Pair<'a, Rule>> for Graph<A>
where
    AList<A>: TryFrom<Pair<'a, Rule>, Error = ParseError<'a>>
{
    type Error = ParseError<'a>;
    fn try_from(p: Pair<'a, Rule>) -> Result<Self, ParseError<'a>> {
        let mut inner = p.clone().into_inner();
        let mut strict = false;
        let mut name = None;
        let mut pair = inner.next().ok_or(ParseError::missing_pair(
            p.clone(),
            vec![Rule::strict, Rule::digraph, Rule::graph],
        ))?;
        if let Rule::strict = pair.as_rule() {
            strict = true;
            pair = inner.next().ok_or(ParseError::missing_pair(
                p.clone(),
                vec![Rule::digraph, Rule::graph],
            ))?;
        }
        let is_digraph = match pair.as_rule() {
            Rule::digraph => true,
            Rule::graph => false,
            r => {
                return Err(ParseError::expect_rule(vec![Rule::digraph, Rule::graph], r));
            }
        };
        pair = inner.next().ok_or(ParseError::missing_pair(
            p.clone(),
            vec![Rule::ident, Rule::stmt_list],
        ))?;
        if let Rule::ident = pair.as_rule() {
            name = Some(String::from(pair.as_str()));
            pair = inner
                .next()
                .ok_or(ParseError::missing_pair(p.clone(), vec![Rule::stmt_list]))?;
        }
        let stmts = StmtList::try_from(pair)?;
        Ok(Graph {
            strict,
            is_digraph,
            name,
            stmts,
        })
    }
}

impl Graph<(String, String)> {
// Those errors are rarely return, so we can afford the size.
#[allow(clippy::result_large_err)]
    /// Parses a graph from a file, which path is given. Notice that when doing so, attributes are
    /// of type `(String, String)`, not the default `(&'a str, &'a str)`, since we are reading and taking
    /// ownership of the content of the file.
    pub fn from_file<'a, P>(p: P) -> Result<Self, GraphFromFileError<'a>>
    where
        P: AsRef<Path>,
    {
        let s = std::fs::read_to_string(p)?;
        let mut pairs = DotParser::parse(Rule::dotgraph, &s)?;
        let pair = pairs.next().expect("The toplevel `Pairs` is empty.");
        match Graph::try_from(pair) {
            Ok(g) => Ok(g),
            Err(e) => panic!("{}", e),
        }
    }
}

impl<'a> TryFrom<&'a str> for Graph<(&'a str, &'a str)> {
    type Error = PestError;
    fn try_from(s: &'a str) -> Result<Self, PestError> {
        let mut pairs = DotParser::parse(Rule::dotgraph, s)?;
        match pairs.next() {
            None => {
                panic!("The toplevel `Pairs` is empty.")
            }
            Some(pair) => Ok(Graph::try_from(pair).unwrap()),
        }
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
    pub fn get_node_ids(self) -> HashSet<NodeID> {
        let clone = self.stmts;
        clone.get_node_ids()
    }
}

/// This is a thin wrapper over vectors of `[Graph]` in order to read files with multiple dotgraph descriptions.
pub struct Graphs<A> {
    /// The set of graphs.
    pub graphs: Vec::<Graph<A>>
}

impl <'a, A> TryFrom<Pair<'a, Rule>> for Graphs<A>
where
    AList<A>: TryFrom<Pair<'a, Rule>, Error = ParseError<'a>>
{
    type Error = ParseError<'a>;
    fn try_from(p: Pair<'a, Rule>) -> Result<Self, ParseError<'a>> {
        let mut inner = p.clone().into_inner().into_iter();
        let first_graph_pair = inner.next().ok_or(ParseError::missing_pair(
            p.clone(),
            vec![Rule::dotfile],
        ))?;
        let mut graphs = vec!(Graph::try_from(first_graph_pair)?);
        for p in inner {
            graphs.push(Graph::try_from(p)?);
        }
        Ok(Graphs { graphs })
    }
}

impl<'a> TryFrom<&'a str> for Graphs<(&'a str, &'a str)> {
    type Error = PestError;
    fn try_from(s: &'a str) -> Result<Self, PestError> {
        DotParser::parse(Rule::dotfile, s).map(|mut p| match p.next() {
            None => {
                panic!("The toplevel `Pairs` is empty.")
            }
            Some(pair) => match Graphs::try_from(pair) {
                Ok(g) => g,
                Err(e) => panic!("{}", e),
            },
        })
    }
}

impl Graphs<(String, String)> {
// Those errors are rarely return, so we can afford the size.
#[allow(clippy::result_large_err)]
    /// Parses multiple graphs from a file, which path is given. Notice that when doing so, attributes are
    /// of type `(String, String)`, not the default `(&'a str, &'a str)`, since we are reading and taking
    /// ownership of the content of the file.
    pub fn from_file<'a, P>(p: P) -> Result<Self, GraphFromFileError<'a>>
    where
        P: AsRef<Path>,
    {
        let s = std::fs::read_to_string(p)?;
        let mut pairs = DotParser::parse(Rule::dotfile, &s)?;
        let pair = pairs.next().expect("The toplevel `Pairs` is empty.");
        match Graphs::try_from(pair) {
            Ok(g) => Ok(g),
            Err(e) => panic!("{}", e),
        }
    }
}

/// A list of statements. This corresponds to the `stmt_list` non-terminal of the
/// grammar.
#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct StmtList<A> {
    /// The list of statements.
    pub stmts: Vec<Stmt<A>>,
}

#[cfg(feature = "to_tokens")]
impl ToTokens for StmtList<(String, String)> {
    fn to_tokens(&self, ts: &mut proc_macro2::TokenStream) {
        let stmts = &self.stmts;
        let tokens = quote! {
            dot_parser::ast::StmtList {
                stmts: std::vec![ #( #stmts ),* ],
            }
        };
        ts.append_all(tokens);
    }
}

impl<'a, A> TryFrom<Pair<'a, Rule>> for StmtList<A>
where
    AList<A>: TryFrom<Pair<'a, Rule>, Error = ParseError<'a>>
{
    type Error = ParseError<'a>;
    fn try_from(p: Pair<'a, Rule>) -> Result<Self, ParseError<'a>> {
        let inner = p.into_inner();
        let mut stmts = Vec::new();
        for stmt in inner {
            stmts.push(Stmt::try_from(stmt)?);
        }
        Ok(StmtList { stmts })
    }
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

/// A statement of the graph. This corresponds to the `stmt` non-terminal of the
/// grammar.
#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub enum Stmt<A> {
    /// A node statement.
    NodeStmt(NodeStmt<A>),
    /// An edge statement.
    EdgeStmt(EdgeStmt<A>),
    /// An attribute statement.
    AttrStmt(AttrStmt<A>),
    /// An alias statement.
    IDEq(String, String),
    /// A subgraph.
    Subgraph(Subgraph<A>),
}

#[cfg(feature = "to_tokens")]
impl ToTokens for Stmt<(String, String)> {
    fn to_tokens(&self, ts: &mut proc_macro2::TokenStream) {
        let tokens = match &self {
            Self::NodeStmt(stmt) => {
                quote! {
                    dot_parser::ast::Stmt::NodeStmt( #stmt )
                }
            },
            Self::EdgeStmt(stmt) => {
                quote! {
                    dot_parser::ast::Stmt::EdgeStmt( #stmt )
                }
            },
            Self::AttrStmt(stmt) => {
                quote! {
                    dot_parser::ast::Stmt::AttrStmt( #stmt )
                }
            },
            Self::IDEq(s1, s2) => {
                quote!{
                    dot_parser::ast::Stmt::IDEq( std::string::ToString::to_string(#s1), std::string::ToString::to_string(#s2) )
                }
            },
            Self::Subgraph(sub) => {
                quote!{
                    dot_parser::ast::Stmt::Subgraph( #sub )
                }
            },
        };
        ts.append_all(tokens);
    }
}

impl<'a, A> TryFrom<Pair<'a, Rule>> for Stmt<A>
where
    AList<A>: TryFrom<Pair<'a, Rule>, Error = ParseError<'a>>
{
    type Error = ParseError<'a>;
    fn try_from(p: Pair<'a, Rule>) -> Result<Self, ParseError<'a>> {
        let inner = p
            .clone()
            .into_inner()
            .next()
            .ok_or(ParseError::missing_pair(
                p.clone(),
                std::vec![
                    Rule::node_stmt,
                    Rule::edge_stmt,
                    Rule::attr_stmt,
                    Rule::id_eq,
                    Rule::subgraph,
                ],
            ))?;
        match inner.as_rule() {
            Rule::node_stmt => NodeStmt::try_from(inner).map(Stmt::NodeStmt),
            Rule::edge_stmt => EdgeStmt::try_from(inner).map(Stmt::EdgeStmt),
            Rule::attr_stmt => AttrStmt::try_from(inner).map(Stmt::AttrStmt),
            Rule::id_eq => {
                let mut inners = inner.into_inner();
                let id1 = inners
                    .next()
                    .ok_or(ParseError::missing_pair(p.clone(), vec![Rule::ident]))?
                    .as_str().into();
                let id2 = inners
                    .next()
                    .ok_or(ParseError::missing_pair(p.clone(), vec![Rule::ident]))?
                    .as_str().into();
                Ok(Stmt::IDEq(id1, id2))
            }
            Rule::subgraph => Subgraph::try_from(inner).map(Stmt::Subgraph),
            other => {
                let error = ParseError::expect_rule(
                    vec![
                        Rule::node_stmt,
                        Rule::edge_stmt,
                        Rule::attr_stmt,
                        Rule::id_eq,
                        Rule::subgraph,
                    ],
                    other,
                );
                Err(error)
            }
        }
    }
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
            Stmt::Subgraph(sub) => Stmt::Subgraph(sub.filter_map_attr(f)),
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

    /// Returns `true` if `self` is a `Subgraph` variant.
    pub fn is_subgraph(&self) -> bool {
        matches!(self, Stmt::Subgraph(..))
    }

    /// Returns all `NodeID`s that appear in this statement.
    fn get_node_ids(&self) -> HashSet<NodeID> {
        match self {
            Stmt::Subgraph(g) => g.get_node_ids(),
            Stmt::EdgeStmt(e) => e.get_node_ids(),
            Stmt::NodeStmt(n) => HashSet::from_iter([n.get_node_id().clone()]),
            _ => HashSet::new(),
        }
    }
}

/// An attribute statement. This corresponds to the rule `attr_stmt`
/// non-terminal of the grammar.
#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub enum AttrStmt<A> {
    /// An `AttrStmt` on the whole graph.
    Graph(AttrList<A>),
    /// An `AttrStmt` on each nodes of the graph.
    Node(AttrList<A>),
    /// An `AttrStmt` on each edges of the graph.
    Edge(AttrList<A>),
}

#[cfg(feature = "to_tokens")]
impl ToTokens for AttrStmt<(String, String)> {
    fn to_tokens(&self, ts: &mut proc_macro2::TokenStream) {
          let new_tokens = match self {
            AttrStmt::Graph(attrs) => quote!{ dot_parser::ast::AttrStmt::Graph(#attrs) },
            AttrStmt::Node(attrs) => quote!{ dot_parser::ast::AttrStmt::Node(#attrs) },
            AttrStmt::Edge(attrs) => quote!{ dot_parser::ast::AttrStmt::Edge(#attrs) },
          };
          ts.append_all(new_tokens);
    }
}

impl<'a, A> TryFrom<Pair<'a, Rule>> for AttrStmt<A>
where
    AList<A>: TryFrom<Pair<'a, Rule>, Error = ParseError<'a>>
{
    type Error = ParseError<'a>;
    fn try_from(p: Pair<'a, Rule>) -> Result<Self, ParseError<'a>> {
        let mut inners = p.clone().into_inner();
        let kind = inners
            .next()
            .ok_or(ParseError::missing_pair(
                p.clone(),
                vec![Rule::graph, Rule::node, Rule::edge],
            ))?
            .as_rule();
        let attr_list_pair = inners
            .next()
            .ok_or(ParseError::missing_pair(p, vec![Rule::attr_list]))?;
        let attr = AttrList::try_from(attr_list_pair)?;
        match kind {
            Rule::graph => Ok(AttrStmt::Graph(attr)),
            Rule::node => Ok(AttrStmt::Node(attr)),
            Rule::edge => Ok(AttrStmt::Edge(attr)),
            r => Err(ParseError::expect_rule(
                vec![Rule::graph, Rule::node, Rule::edge],
                r,
            )),
        }
    }
}

impl<A> AttrStmt<A> {
    pub (in crate) fn filter_map_attr<B>(self, f: &dyn Fn(A) -> Option<B>) -> AttrStmt<B> {
        match self {
            AttrStmt::Graph(attr) => AttrStmt::Graph(attr.filter_map_attr(f)),
            AttrStmt::Node(attr) => AttrStmt::Node(attr.filter_map_attr(f)),
            AttrStmt::Edge(attr) => AttrStmt::Edge(attr.filter_map_attr(f)),
        }
    }
}

/// A list of `AList`s, i.e. a list of list of attributes. This (strange)
/// indirection is induced by the grammar. This structure corresponds to the
/// `attr_list` non-terminal of the grammar.
///
/// Notice methods `flatten` and `flatten_ref` to remove the indirection.
#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct AttrList<A> {
    /// The list of `AList`s.
    pub elems: Vec<AList<A>>,
}

#[cfg(feature = "to_tokens")]
impl ToTokens for AttrList<(String, String)> {
    fn to_tokens(&self, ts: &mut proc_macro2::TokenStream) {
        let elems = &self.elems;
        let tokens = quote! {
            dot_parser::ast::AttrList {
                elems: std::vec![ #( #elems ),* ]
            }
        };
        ts.append_all(tokens);
    }
}

impl<'a, A> TryFrom<Pair<'a, Rule>> for AttrList<A>
where
    AList<A>: TryFrom<Pair<'a, Rule>, Error = ParseError<'a>>
{
    type Error = ParseError<'a>;
    fn try_from(p: Pair<'a, Rule>) -> Result<Self, ParseError<'a>> {
        let mut v: Vec<AList<A>> = Vec::new();
        let mut inners = p.clone().into_inner();
        let alist_pair = inners
            .next()
            .ok_or(ParseError::missing_pair(p, vec![Rule::a_list]))?;
        let alist = AList::try_from(alist_pair)?;
        let mut tail = inners
            .next()
            .map(|p| {
                AttrList::try_from(p)
                    .map(|alist| alist.elems)
                    .unwrap_or_default()
            })
            .unwrap_or_default();
        v.push(alist);
        v.append(&mut tail);

        Ok(AttrList { elems: v })
    }
}

impl<A> AttrList<A> {
    pub (in crate) fn filter_map_attr<B>(self, f: &dyn Fn(A) -> Option<B>) -> AttrList<B> {
        AttrList {
            elems: self
                .into_iter()
                .map(|alist| alist.filter_map_attr(f))
                .collect(),
        }
    }

    /// Flatten the nested `AList`s: returns a single `AList` that contains all
    /// attributes contained in the `AttrList`.
    pub fn flatten(self) -> AList<A> {
        self.into()
    }

    /// Flatten the nested `AList`s: returns a single `AList` that contains all
    /// attributes contained in the `AttrList`.
    pub fn flatten_ref(&self) -> AList<&A> {
        self.into()
    }
}

impl<A> FromIterator<AList<A>> for AttrList<A> {
    fn from_iter<T>(iter: T) -> Self
    where
        T: IntoIterator<Item = AList<A>>,
    {
        Self {
            elems: iter.into_iter().map(|u| u.into_iter().collect()).collect(),
        }
    }
}

impl<A> IntoIterator for AttrList<A> {
    type Item = AList<A>;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.elems.into_iter()
    }
}

impl<'a, A> IntoIterator for &'a AttrList<A> {
    type Item = &'a AList<A>;
    type IntoIter = std::slice::Iter<'a, AList<A>>;

    fn into_iter(self) -> Self::IntoIter {
        self.elems.iter()
    }
}

/// A list of attributes. This corresponds to the `a_list` non-terminal of the
/// grammar.
#[derive(Debug, Eq, PartialEq, Ord, PartialOrd, Hash, Clone)]
pub struct AList<A> {
    /// The attributes in the list.
    pub elems: Vec<A>,
}

#[cfg(feature = "to_tokens")]
impl ToTokens for AList<(String, String)> {
    fn to_tokens(&self, ts: &mut proc_macro2::TokenStream) {
        let elems = &self.elems;
        let elems = elems.iter().map(|(s1, s2)| quote!{ (#s1, #s2) });
        let tokens = quote! {
            dot_parser::ast::AList {
                elems: std::vec![ #( #elems ),* ]
            }
        };
        ts.append_all(tokens);
    }
}

#[cfg(feature = "display")]
impl<A> Display for AList<A>
where
    A: Display,
{
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), std::fmt::Error> {
        for attr in &self.elems {
            write!(f, "{}; ", attr)?;
        }
        Ok(())
    }
}

impl<'a> TryFrom<Pair<'a, Rule>> for AList<(&'a str, &'a str)> {
    type Error = ParseError<'a>;
    fn try_from(p: Pair<'a, Rule>) -> Result<Self, ParseError<'a>> {
        let mut v = Vec::new();
        let mut inners = p.clone().into_inner();
        let id1 = inners
            .next()
            .ok_or(ParseError::missing_pair(p.clone(), vec![Rule::ident]))?
            .as_str();
        let id2 = inners
            .next()
            .ok_or(ParseError::missing_pair(p.clone(), vec![Rule::ident]))?
            .as_str();
        let mut tail = inners
            .next()
            .map(|p| AList::try_from(p).map(|alist| alist.elems).unwrap_or_default())
            .unwrap_or_default();
        v.push((id1, id2));
        v.append(&mut tail);

        Ok(AList {
            elems: v,
        })
    }
}

impl<'a> TryFrom<Pair<'a, Rule>> for AList<(String, String)> {
    type Error = ParseError<'a>;
    fn try_from(p: Pair<'a, Rule>) -> Result<Self, ParseError<'a>> {
        let mut v: Vec<(String, String)> = Vec::new();
        let mut inners = p.clone().into_inner();
        let id1 = inners
            .next()
            .ok_or(ParseError::missing_pair(p.clone(), vec![Rule::ident]))?
            .as_str().into();
        let id2 = inners
            .next()
            .ok_or(ParseError::missing_pair(p.clone(), vec![Rule::ident]))?
            .as_str().into();
        let mut tail = inners
            .next()
            .map(|p| AList::try_from(p).map(|alist| alist.elems).unwrap_or_default())
            .unwrap_or_default();
        v.push((id1, id2));
        v.append(&mut tail);

        Ok(AList {
            elems: v,
        })
    }
}

impl<A> AList<A> {
    pub(crate) fn filter_map_attr<B>(self, f: &dyn Fn(A) -> Option<B>) -> AList<B> {
        AList {
            elems: self.into_iter().filter_map(f).collect(),
        }
    }

    pub(crate) fn empty() -> Self {
        AList {
            elems: Vec::new(),
        }
    }

    #[cfg(feature = "display")]
    /// Returns `true` if the list of attributes is empty
    pub(crate) fn is_empty(&self) -> bool {
        self.elems.is_empty()
    }
}

impl<A> FromIterator<A> for AList<A> {
    fn from_iter<T>(iter: T) -> Self
    where
        T: IntoIterator<Item = A>,
    {
        Self {
            elems: iter.into_iter().collect(),
        }
    }
}

impl<A> IntoIterator for AList<A> {
    type Item = A;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.elems.into_iter()
    }
}

impl<'a, A> IntoIterator for &'a AList<A> {
    type Item = &'a A;
    type IntoIter = std::slice::Iter<'a, A>;

    fn into_iter(self) -> Self::IntoIter {
        self.elems.iter()
    }
}

impl<A> From<AttrList<A>> for AList<A> {
    fn from(attr: AttrList<A>) -> Self {
        attr.into_iter().flatten().collect()
    }
}

impl<'a, A> From<&'a AttrList<A>> for AList<&'a A> {
    fn from(attr: &'a AttrList<A>) -> Self {
        attr.into_iter().flatten().collect()
    }
}

/// The description of an edge. This corresponds to the `edge_stmt` non-terminal
/// of the grammar.
#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct EdgeStmt<A> {
    /// The origin of the edge.
    pub from: Either<NodeID, Subgraph<A>>,
    /// The destination of the edge.
    pub next: EdgeRHS<A>,
    /// The attributes of the edge.
    pub attr: Option<AttrList<A>>,
}

#[cfg(feature = "to_tokens")]
impl ToTokens for EdgeStmt<(String, String)> {
    fn to_tokens(&self, ts: &mut proc_macro2::TokenStream) {
        let from = match &self.from {
            Either::Left(id) => {
                quote!{ dot_parser::ast::either::Either::Left( #id ) }
            },
            Either::Right(sub) => {
                quote!{ dot_parser::ast::either::Either::Right( #sub ) }
            },
        };
        let attr = match &self.attr {
            Some(attr) => {
                quote!{ std::option::Option::Some( #attr ) }
            },
            None => {
                quote!{ std::option::Option::None }
            },
        };
        let next = &self.next;
        let tokens = quote! {
            dot_parser::ast::EdgeStmt {
                from: #from,
                next: #next,
                attr: #attr,
            }
        };
        ts.append_all(tokens);
    }
}

impl<'a, A> TryFrom<Pair<'a, Rule>> for EdgeStmt<A>
where
    AList<A>: TryFrom<Pair<'a, Rule>, Error = ParseError<'a>>
{
    type Error = ParseError<'a>;
    fn try_from(p: Pair<'a, Rule>) -> Result<Self, ParseError<'a>> {
        let mut inners = p.clone().into_inner();
        let from_pair = inners.next().ok_or(ParseError::missing_pair(
            p.clone(),
            vec![Rule::node_id, Rule::subgraph],
        ))?;
        let rule = from_pair.as_rule();
        let from = match rule {
            Rule::node_id => Either::Left(NodeID::try_from(from_pair)?),
            Rule::subgraph => Either::Right(Subgraph::try_from(from_pair)?),
            r => {
                return Err(ParseError::expect_rule(
                    vec![Rule::node_id, Rule::subgraph],
                    r,
                ));
            }
        };
        let next = inners
            .next()
            .map(EdgeRHS::try_from)
            .transpose()?
            .ok_or(ParseError::missing_pair(p, vec![Rule::edge_rhs]))?;
        let attr = inners.next().map(AttrList::try_from).transpose()?;

        Ok(EdgeStmt { from, next, attr })
    }
}

impl<'a, A> EdgeStmt<A> {
    fn filter_map_attr<B>(self, f: &'a dyn Fn(A) -> Option<B>) -> EdgeStmt<B> {
        let new_from = match self.from {
            Either::Left(node_id) => Either::Left(node_id),
            Either::Right(subgraph) => Either::Right(subgraph.filter_map_attr(f)),
        };

        let new_next = self.next.filter_map_attr(f);

        EdgeStmt {
            from: new_from,
            next: new_next,
            attr: self.attr.map(|a| a.filter_map_attr(f)),
        }
    }

    /// Flatten the EdgeStmt, i.e. removes cases where multiple EdgeRHS are in a single statement.
    pub fn flatten(self) -> Vec<EdgeStmt<A>>
    where
        A: Clone,
    {
        let mut from = self.from;
        let mut to = self.next;
        let attr = self.attr;

        let mut v = Vec::new();

        loop {
            let next_step = EdgeStmt {
                from: from.clone(),
                next: EdgeRHS {
                    to: to.to.clone(),
                    next: None,
                },
                attr: attr.clone(),
            };
            v.push(next_step);
            match to.next {
                None => return v,
                Some(rhs) => {
                    from = to.to;
                    to = *rhs;
                }
            }
        }
    }

    fn get_node_ids(&self) -> HashSet<NodeID> {
        let mut nexts = self.next.get_node_ids();
        match &self.from {
            Either::Left(node_id) => {
                nexts.insert(node_id.clone());
            }
            Either::Right(subgraph) => {
                return nexts.union(&subgraph.get_node_ids()).cloned().collect();
            }
        };
        nexts
    }
}

/// The Right hand side of an edge description. This corresponds to the
/// `EdgeRHS` non-terminal of the grammar.
/// Notice that the grammar allows multiple EdgeRHS in sequence, to chain edges:
/// `A -> B -> C`.
#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct EdgeRHS<A> {
    /// The identifier of the destination of the edge.
    pub to: Either<NodeID, Subgraph<A>>,
    /// A possible chained RHS.
    pub next: Option<Box<EdgeRHS<A>>>,
}

#[cfg(feature = "to_tokens")]
impl ToTokens for EdgeRHS<(String, String)> {
    fn to_tokens(&self, ts: &mut proc_macro2::TokenStream) {
        let to = match &self.to {
            Either::Left(id) => quote!{ dot_parser::ast::either::Either::Left( #id ) },
            Either::Right(sub) => quote!{ dot_parser::ast::either::Either::Right( #sub ) },
        };
        let next = match &self.next {
            Some(next) => quote! {
                std::option::Option::Some(
                    std::boxed::Box::new( #next )
                )
            },
            None => quote! {
                std::option::Option::None
            }
        };
        let tokens = quote! {
            dot_parser::ast::EdgeRHS {
                to: #to,
                next: #next,
            }
        };
        ts.append_all(tokens);
    }
}

impl<'a, A> TryFrom<Pair<'a, Rule>> for EdgeRHS<A>
where
    AList<A>: TryFrom<Pair<'a, Rule>, Error = ParseError<'a>>
{
    type Error = ParseError<'a>;
    fn try_from(p: Pair<'a, Rule>) -> Result<Self, ParseError<'a>> {
        let mut inners = p.clone().into_inner();
        let to_pair = inners.next().ok_or(ParseError::missing_pair(
            p,
            vec![Rule::node_id, Rule::subgraph],
        ))?;
        let to = match to_pair.as_rule() {
            Rule::node_id => Either::Left(NodeID::try_from(to_pair)?),
            Rule::subgraph => {
                // p is necessarily a subgraph.
                Either::Right(Subgraph::try_from(to_pair)?)
            }
            r => {
                return Err(ParseError::expect_rule(
                    vec![Rule::node_id, Rule::subgraph],
                    r,
                ));
            }
        };
        let next = inners.next().map(EdgeRHS::try_from).transpose()?.map(Box::new);
        Ok(EdgeRHS { to, next })
    }
}

impl<'a, A> EdgeRHS<A> {
    fn filter_map_attr<B>(self, f: &'a dyn Fn(A) -> Option<B>) -> EdgeRHS<B> {
        let to = match self.to {
            Either::Left(node_id) => Either::Left(node_id),
            Either::Right(subgraph) => Either::Right(subgraph.filter_map_attr(f)),
        };

        let next = self.next.map(|e| Box::new(e.filter_map_attr(f)));

        EdgeRHS { to, next }
    }

    fn get_node_ids(&self) -> HashSet<NodeID> {
        let mut nexts: HashSet<NodeID> = self
            .next
            .as_ref()
            .map(|n| n.get_node_ids())
            .unwrap_or_default();
        match &self.to {
            Either::Left(node_id) => {
                nexts.insert(node_id.clone());
            }
            Either::Right(subgraph) => {
                return nexts.union(&subgraph.get_node_ids()).cloned().collect();
            }
        };
        nexts
    }
}

/// This structure corresponds to the `node_stmt` non-terminal of the grammar.
/// It is basically a node identifier attached to some attributes.
#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct NodeStmt<A> {
    /// The identifier of the node.
    pub node: NodeID,
    /// The possible list of attributes.
    pub attr: Option<AttrList<A>>,
}

#[cfg(feature = "to_tokens")]
impl ToTokens for NodeStmt<(String, String)> {
    fn to_tokens(&self, ts: &mut proc_macro2::TokenStream) {
        let node = &self.node;
        let attrs = match &self.attr {
            Some(a) => quote! { std::option::Option::Some(#a) },
            None => quote! { std::option::Option::None }
        };
        let tokens = quote! {
            dot_parser::ast::NodeStmt {
                node: #node,
                attr: #attrs
            }
        };
        ts.append_all(tokens);
    }
}

impl<'a, A> TryFrom<Pair<'a, Rule>> for NodeStmt<A>
where
    AList<A>: TryFrom<Pair<'a, Rule>, Error = ParseError<'a>>
{
    type Error = ParseError<'a>;
    fn try_from(p: Pair<'a, Rule>) -> Result<Self, ParseError<'a>> {
        let mut inners = p.clone().into_inner();
        let node = inners
            .next()
            .map(NodeID::try_from)
            .transpose()?
            .ok_or(ParseError::missing_pair(p, vec![Rule::node_id]))?;
        let attr = inners.next().map(AttrList::try_from).transpose()?;
        Ok(NodeStmt { node, attr })
    }
}

impl<A> NodeStmt<A> {
    pub (in crate) fn filter_map_attr<B>(self, f: &dyn Fn(A) -> Option<B>) -> NodeStmt<B> {
        NodeStmt {
            node: self.node,
            attr: self.attr.map(|a| a.filter_map_attr(f)),
        }
    }

    /// Get the name of the `NodeStmt`, i.e. the identifier of the
    /// `NodeID` contained in the `NodeStmt`.
    pub fn name(&self) -> &str {
        &self.node.id
    }

    pub (in crate) fn get_node_id(&self) -> &NodeID {
        &self.node
    }
}

/// This structure corresponds to the `node_id` non-terminal of the grammar.
/// If contains the identifier of the node, and possibly a port.
#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct NodeID {
    /// The identifier of the node.
    pub id: String,
    /// The port of the node, if any.
    pub port: Option<Port>,
}

#[cfg(feature = "to_tokens")]
impl ToTokens for NodeID {
    fn to_tokens(&self, ts: &mut proc_macro2::TokenStream) {
        let id = &self.id;
        let port = match &self.port {
            Some(port) => quote!{ std::option::Option::Some(#port) },
            None => quote! { std::option::Option::None },
        };
        let tokens = quote!{
            dot_parser::ast::NodeID {
                id: std::string::ToString::to_string(#id),
                port: #port,
            }
        };
        ts.append_all(tokens);
    }
}

impl<'a> TryFrom<Pair<'a, Rule>> for NodeID
{
    type Error = ParseError<'a>;
    fn try_from(p: Pair<'a, Rule>) -> Result<Self, ParseError<'a>> {
        let mut inners = p.clone().into_inner();
        let id = inners
            .next()
            .ok_or(ParseError::missing_pair(p, vec![Rule::node_id]))?
            .as_str()
            .to_string();
        let port = inners.next().map(Port::try_from).transpose()?;
        Ok(NodeID { id, port })
    }
}

/// This enum corresponds to the `port` non-terminal of the grammar.
#[derive(Debug, Eq, PartialEq, Ord, PartialOrd, Hash, Clone)]
pub enum Port {
    /// The variant in which an ID is given, and possibly a compass point.
    ID(String, Option<CompassPt>),
    /// The variant in which only a compass point is given.
    Compass(CompassPt),
}

#[cfg(feature = "to_tokens")]
impl ToTokens for Port {
    fn to_tokens(&self, ts: &mut proc_macro2::TokenStream) {
        let tokens = match self {
            Port::ID(s, cpss) => {
                match cpss {
                    Some(cpss) => {
                        quote!{ dot_parser::ast::Port::ID(std::string::ToString::to_string(#s), std::option::Option::Some(#cpss)) }
                    },
                    None => {
                        quote!{ dot_parser::ast::Port::ID(std::string::ToString::to_string(#s), std::option::Option::None) }
                    }
                }
            },
            Port::Compass(cpss) => {
                quote!{ dot_parser::ast::Port::Compass(#cpss) }
            }
        };
        ts.append_all(tokens);
    }
}

impl<'a> TryFrom<Pair<'a, Rule>> for Port
{
    type Error = ParseError<'a>;
    fn try_from(p: Pair<'a, Rule>) -> Result<Self, ParseError<'a>> {
        let mut inners = p.clone().into_inner();
        let inner = inners.next().ok_or(ParseError::missing_pair(
            p,
            vec![Rule::compass_pt, Rule::ident],
        ))?;
        match inner.as_rule() {
            Rule::compass_pt => Ok(Port::Compass(CompassPt::try_from(inner)?)),
            Rule::ident => {
                let opt_comp = inners.next().map(CompassPt::try_from).transpose()?;
                Ok(Port::ID(inner.as_str().to_string(), opt_comp))
            }
            r => Err(ParseError::expect_rule(
                vec![Rule::compass_pt, Rule::ident],
                r,
            )),
        }
    }
}

#[cfg(feature = "display")]
impl Display for Port {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), std::fmt::Error> {
        match self {
            Port::ID(name, Some(cpss)) => {
                write!(f, ": {name} : {cpss}")
            }
            Port::ID(name, None) => {
                write!(f, ": {name}")
            }
            Port::Compass(cpss) => {
                write!(f, ": {}", cpss)
            }
        }
    }
}

/// A subgraph. This corresponds to the `subgraph` non-terminal of the grammar.
#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct Subgraph<A> {
    /// The (optional) identifier of the subgraph.
    pub id: Option<String>,
    /// The statements that describe the subgraph.
    pub stmts: StmtList<A>,
}

#[cfg(feature = "to_tokens")]
impl ToTokens for Subgraph<(String, String)> {
    fn to_tokens(&self, ts: &mut proc_macro2::TokenStream) {
        let id = match &self.id {
            Some(s) => quote!{ std::option::Option::Some(#s) },
            None => quote! { std::option::Option::None },
        };
        let stmts = &self.stmts;
        let tokens = quote! {
            dot_parser::ast::Subgraph {
                id: #id,
                stmts: #stmts,
            }
        };
        ts.append_all(tokens);
    }
}

impl<'a, A> TryFrom<Pair<'a, Rule>> for Subgraph<A>
where
    AList<A>: TryFrom<Pair<'a, Rule>, Error = ParseError<'a>>
{
    type Error = ParseError<'a>;
    fn try_from(p: Pair<'a, Rule>) -> Result<Self, ParseError<'a>> {
        let mut inners = p.clone().into_inner();
        let mut inner = inners.next().ok_or(ParseError::missing_pair(
            p.clone(),
            vec![Rule::ident, Rule::stmt_list],
        ))?;
        let id = if let Rule::ident = inner.as_rule() {
            let id_str = inner.as_str().to_string();
            inner = inners
                .next()
                .ok_or(ParseError::missing_pair(p.clone(), vec![Rule::stmt_list]))?;
            Some(id_str)
        } else {
            None
        };
        let stmts = StmtList::try_from(inner)?;
        Ok(Subgraph { id, stmts })
    }
}

impl<A> Subgraph<A> {
    fn filter_map_attr<B>(self, f: &dyn Fn(A) -> Option<B>) -> Subgraph<B> {
        Subgraph {
            id: self.id,
            stmts: self
                .stmts
                .into_iter()
                .map(|stmt| stmt.filter_map_attr(f))
                .collect(),
        }
    }

    /// Extract a subgraph as a standalone graph.
    pub (in crate) fn into_graph(self, strict: bool, is_digraph: bool) -> Graph<A> {
        Graph {
            strict,
            is_digraph,
            name: self.id.map(String::from),
            stmts: self.stmts,
        }
    }

    fn get_node_ids(&self) -> HashSet<NodeID> {
        self.stmts.get_node_ids()
    }
}

/// An enum that corresponds to the `compass_pt` non-terminal of the grammar.
#[derive(Debug, Eq, PartialEq, Ord, PartialOrd, Hash, Copy, Clone)]
pub enum CompassPt {
    /// A North orientation
    N,
    /// A North-East orientation
    NE,
    /// An East orientation
    E,
    /// A South-East orientation
    SE,
    /// A South orientation
    S,
    /// A South-West orientation
    SW,
    /// A West orientation
    W,
    /// A North-West orientation
    NW,
    /// A Central orientation
    C,
    /// An unspecified orientation
    Underscore,
}

#[cfg(feature = "display")]
impl Display for CompassPt {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), std::fmt::Error> {
        match self {
            CompassPt::N => write!(f, "n"),
            CompassPt::NE => write!(f, "ne"),
            CompassPt::E => write!(f, "e"),
            CompassPt::SE => write!(f, "se"),
            CompassPt::S => write!(f, "s"),
            CompassPt::SW => write!(f, "sw"),
            CompassPt::W => write!(f, "w"),
            CompassPt::NW => write!(f, "nw"),
            CompassPt::C => write!(f, "c"),
            CompassPt::Underscore => write!(f, "_"),
        }
    }
}

impl<'a> TryFrom<Pair<'a, Rule>> for CompassPt { 
    type Error = ParseError<'a>;
    fn try_from(p: Pair<'a, Rule>) -> Result<Self, ParseError<'a>> {
        match p
            .clone()
            .into_inner()
            .next()
            .ok_or(ParseError::missing_pair(
                p,
                vec![
                    Rule::n,
                    Rule::ne,
                    Rule::e,
                    Rule::se,
                    Rule::s,
                    Rule::sw,
                    Rule::w,
                    Rule::nw,
                    Rule::c,
                    Rule::underscore,
                ],
            ))?
            .as_rule()
        {
            Rule::n => Ok(CompassPt::N),
            Rule::ne => Ok(CompassPt::NE),
            Rule::e => Ok(CompassPt::E),
            Rule::se => Ok(CompassPt::SE),
            Rule::s => Ok(CompassPt::S),
            Rule::sw => Ok(CompassPt::SW),
            Rule::w => Ok(CompassPt::W),
            Rule::nw => Ok(CompassPt::NW),
            Rule::c => Ok(CompassPt::C),
            Rule::underscore => Ok(CompassPt::Underscore),
            r => Err(ParseError::expect_rule(
                vec![
                    Rule::n,
                    Rule::ne,
                    Rule::e,
                    Rule::se,
                    Rule::s,
                    Rule::sw,
                    Rule::w,
                    Rule::nw,
                    Rule::c,
                    Rule::underscore,
                ],
                r,
            )),
        }
    }
}

#[cfg(feature = "to_tokens")]
impl ToTokens for CompassPt {
      fn to_tokens(&self, tokens: &mut TokenStream) {
          let new_tokens = match self {
            CompassPt::N => quote!{ dot_parser::ast::CompassPt::N},
            CompassPt::NE => quote!{ dot_parser::ast::CompassPt::NE},
            CompassPt::E => quote!{ dot_parser::ast::CompassPt::E},
            CompassPt::SE => quote!{ dot_parser::ast::CompassPt::SE},
            CompassPt::S => quote!{ dot_parser::ast::CompassPt::S},
            CompassPt::SW => quote!{ dot_parser::ast::CompassPt::SW},
            CompassPt::W => quote!{ dot_parser::ast::CompassPt::W},
            CompassPt::NW => quote!{ dot_parser::ast::CompassPt::NW},
            CompassPt::C => quote!{ dot_parser::ast::CompassPt::C},
            CompassPt::Underscore => quote!{ dot_parser::ast::CompassPt::Underscore},
          };
          tokens.append_all(new_tokens);
      }
}
