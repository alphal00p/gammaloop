//! DOT syntax: a dependency-free reader for the
//! [DOT language](https://graphviz.org/doc/info/lang.html) that flattens each
//! graph into the statements [`DotGraph`](super::DotGraph) builds half-edge
//! graphs from.
//!
//! Lexing follows Graphviz (`lib/cgraph/scan.l`). Keywords are
//! case-insensitive, bare identifiers may contain non-ASCII characters, and a
//! number directly followed by letters (`1H`) is two identifiers. `/* */`
//! comments are skipped, and so is the rest of a line from `//` or `#` outside
//! a string. HTML strings nest and keep their angle brackets. Quoted strings
//! decode `\"` to a quote and drop a backslash-newline line continuation, as
//! Graphviz does; unlike Graphviz, they also decode `\\` to one backslash, so
//! that every string, including one that ends in a backslash, can be written
//! back. Other escapes, such as `\N` or `\l`, are kept for the consumer. String
//! concatenation with `+` is rejected.
//!
//! Flattening follows linnet's graph model. Subgraphs are inlined into the
//! enclosing graph. Edge chains and subgraph endpoints expand into single
//! edges; the statements inside an endpoint subgraph come first, then the edges
//! to and from its nodes. Attribute statements apply to the whole graph,
//! whatever their position, and `ID = ID` is a graph attribute like
//! `graph [ID = ID]`. Nodes are keyed by their decoded identifier, and
//! repeated node statements merge their attributes. `strict` is accepted and
//! dropped, because half-edge graphs keep parallel edges.

use std::{borrow::Cow, collections::BTreeMap, fmt, str::FromStr};

use ahash::{HashSet, HashSetExt};
use thiserror::Error;

/// Subgraphs nest at most this deep. Parsing recurses once per level, and the
/// Typst plugin parses user-supplied DOT on a small WebAssembly stack.
const MAX_SUBGRAPH_DEPTH: usize = 256;

/// A DOT compass point.
///
/// The variant order is the order in which linnet sorts edge endpoints (see
/// the DOT input guide); the discriminant is the archived byte value.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum CompassPt {
    N,
    NE,
    E,
    SE,
    S,
    SW,
    W,
    NW,
    C,
    /// `_`, which leaves the side to the renderer.
    Underscore,
}

impl CompassPt {
    const ALL: [Self; 10] = [
        Self::N,
        Self::NE,
        Self::E,
        Self::SE,
        Self::S,
        Self::SW,
        Self::W,
        Self::NW,
        Self::C,
        Self::Underscore,
    ];

    /// The DOT spelling, such as `ne` or `_`.
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::N => "n",
            Self::NE => "ne",
            Self::E => "e",
            Self::SE => "se",
            Self::S => "s",
            Self::SW => "sw",
            Self::W => "w",
            Self::NW => "nw",
            Self::C => "c",
            Self::Underscore => "_",
        }
    }
}

impl fmt::Display for CompassPt {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// A string that is not one of the ten DOT compass point spellings.
#[derive(Clone, Debug, PartialEq, Eq, Error)]
#[error("`{0}` is not a compass point (n, ne, e, se, s, sw, w, nw, c or _)")]
pub struct UnknownCompassPt(pub String);

impl FromStr for CompassPt {
    type Err = UnknownCompassPt;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::ALL
            .into_iter()
            .find(|compass| compass.as_str() == s)
            .ok_or_else(|| UnknownCompassPt(s.to_owned()))
    }
}

impl TryFrom<u8> for CompassPt {
    type Error = u8;

    /// Inverts `compass as u8`, returning an invalid byte unchanged.
    fn try_from(value: u8) -> Result<Self, u8> {
        Self::ALL.get(usize::from(value)).copied().ok_or(value)
    }
}

/// Invalid DOT input, located by a 1-based line and a 1-based column counted
/// in characters.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DotSyntaxError {
    /// Byte offset of the error in the input.
    pub offset: usize,
    pub line: usize,
    pub column: usize,
    /// What was expected, and what was found instead.
    pub message: String,
    /// The input line containing the error.
    pub source_line: String,
}

impl fmt::Display for DotSyntaxError {
    /// Shows the message, the source line and a caret under the column. A long
    /// line is cut to a window around the column.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        const CONTEXT: usize = 60;
        let (line, column) = (self.line, self.column);
        writeln!(f, "line {line}, column {column}: {}", self.message)?;
        let skip = column.saturating_sub(CONTEXT + 1);
        let mut shown = self.source_line.chars().skip(skip);
        let window: String = shown.by_ref().take(2 * CONTEXT).collect();
        let (lead, trail) = (
            if skip > 0 { "..." } else { "" },
            if shown.next().is_some() { "..." } else { "" },
        );
        // Tabs are kept so that the caret lines up under them.
        let pad: String = (lead.chars().map(|_| ' '))
            .chain(window.chars().take(column - 1 - skip))
            .map(|c| if c == '\t' { '\t' } else { ' ' })
            .collect();
        write!(f, "  | {lead}{window}{trail}\n  | {pad}^")
    }
}

impl std::error::Error for DotSyntaxError {}

/// An attribute `key = value` with decoded strings.
pub(crate) type Attr<'a> = (Cow<'a, str>, Cow<'a, str>);

/// The port of an edge endpoint. The derived order is part of linnet's
/// documented edge order: named ports sort before bare compass points.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) enum Port<'a> {
    /// `node:name` or `node:name:compass`.
    Named(Cow<'a, str>, Option<CompassPt>),
    /// `node:compass`.
    Compass(CompassPt),
}

impl Port<'_> {
    pub(crate) fn name(&self) -> Option<&str> {
        match self {
            Self::Named(name, _) => Some(name),
            Self::Compass(_) => None,
        }
    }

    pub(crate) fn compass(&self) -> Option<CompassPt> {
        match *self {
            Self::Named(_, compass) => compass,
            Self::Compass(compass) => Some(compass),
        }
    }
}

/// A node identifier with an optional port. Edges sort by their endpoints in
/// this derived order.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) struct Endpoint<'a> {
    pub(crate) id: Cow<'a, str>,
    pub(crate) port: Option<Port<'a>>,
}

/// A single edge after chains and subgraph endpoints have been expanded.
#[derive(Debug, PartialEq, Eq)]
pub(crate) struct Edge<'a> {
    pub(crate) from: Endpoint<'a>,
    pub(crate) to: Endpoint<'a>,
    pub(crate) attrs: Vec<Attr<'a>>,
}

/// One DOT graph flattened into the statements linnet consumes. Attribute
/// lists keep their source order, so later entries win.
#[derive(Debug, Default, PartialEq, Eq)]
pub(crate) struct FlatGraph<'a> {
    pub(crate) directed: bool,
    pub(crate) name: Option<Cow<'a, str>>,
    /// `graph [...]` attributes and `ID = ID` statements.
    pub(crate) graph_attrs: Vec<Attr<'a>>,
    /// `node [...]` attributes.
    pub(crate) node_defaults: Vec<Attr<'a>>,
    /// `edge [...]` attributes.
    pub(crate) edge_defaults: Vec<Attr<'a>>,
    /// Every node, including those named only by edges, with the attributes
    /// of all its node statements.
    pub(crate) nodes: BTreeMap<Cow<'a, str>, Vec<Attr<'a>>>,
    pub(crate) edges: Vec<Edge<'a>>,
}

impl<'a> FlatGraph<'a> {
    /// Parses input that holds exactly one graph.
    pub(crate) fn parse(src: &'a str) -> Result<Self, DotSyntaxError> {
        let mut parser = Parser::new(src)?;
        let graph = parser.graph()?;
        if parser.token != Token::End {
            return Err(parser.expected("the end of input after the graph"));
        }
        Ok(graph)
    }

    /// Parses input that holds one or more graphs.
    pub(crate) fn parse_all(src: &'a str) -> Result<Vec<Self>, DotSyntaxError> {
        let mut parser = Parser::new(src)?;
        let mut graphs = vec![parser.graph()?];
        while parser.token != Token::End {
            graphs.push(parser.graph()?);
        }
        Ok(graphs)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum IdKind {
    Bare,
    Numeral,
    Quoted,
    Html,
}

impl IdKind {
    /// The value of an identifier token spelled `text`; see the module
    /// documentation for the decoding of quoted strings.
    fn decode(self, text: &str) -> Cow<'_, str> {
        let inner = match self {
            Self::Quoted => &text[1..text.len() - 1],
            Self::Bare | Self::Numeral | Self::Html => return Cow::Borrowed(text),
        };
        if !inner.contains('\\') {
            return Cow::Borrowed(inner);
        }
        let mut decoded = String::with_capacity(inner.len());
        let mut chars = inner.chars();
        while let Some(c) = chars.next() {
            if c != '\\' {
                decoded.push(c);
                continue;
            }
            // The lexer never ends a quoted string in an unpaired backslash.
            match chars.next() {
                Some(escaped @ ('"' | '\\')) => decoded.push(escaped),
                Some('\n') => {}
                other => decoded.extend(std::iter::once('\\').chain(other)),
            }
        }
        Cow::Owned(decoded)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Token {
    Id(IdKind),
    Strict,
    Graph,
    Digraph,
    Node,
    Edge,
    Subgraph,
    /// `->`
    Arrow,
    /// `--`
    Line,
    LBrace,
    RBrace,
    LBracket,
    RBracket,
    Equals,
    Colon,
    Semicolon,
    Comma,
    End,
}

impl Token {
    fn keyword(word: &[u8]) -> Option<Self> {
        let is = |keyword: &[u8]| word.eq_ignore_ascii_case(keyword);
        Some(match word.len() {
            4 if is(b"node") => Self::Node,
            4 if is(b"edge") => Self::Edge,
            5 if is(b"graph") => Self::Graph,
            6 if is(b"strict") => Self::Strict,
            7 if is(b"digraph") => Self::Digraph,
            8 if is(b"subgraph") => Self::Subgraph,
            _ => return None,
        })
    }
}

/// The endpoints on one side of an edge operator, or of a node or subgraph
/// statement.
enum Side<'a> {
    One(Endpoint<'a>),
    Group(Vec<Endpoint<'a>>),
}

impl<'a> Side<'a> {
    fn endpoints(&self) -> &[Endpoint<'a>] {
        match self {
            Side::One(endpoint) => std::slice::from_ref(endpoint),
            Side::Group(endpoints) => endpoints,
        }
    }
}

/// A recursive-descent parser over a one-token-lookahead lexer.
struct Parser<'a> {
    src: &'a str,
    token: Token,
    /// Byte span of `token`; lexing resumes at `end`.
    start: usize,
    end: usize,
    directed: bool,
    depth: usize,
    /// Positive while `token` is, or directly follows, a bare identifier glued
    /// to a preceding number (`18s`), which reads as two identifiers.
    glued_number: u8,
}

impl<'a> Parser<'a> {
    fn new(src: &'a str) -> Result<Self, DotSyntaxError> {
        let mut parser = Parser {
            src,
            token: Token::End,
            start: 0,
            end: 0,
            directed: true,
            depth: 0,
            glued_number: 0,
        };
        parser.bump()?;
        Ok(parser)
    }

    fn error_at(&self, offset: usize, message: String) -> DotSyntaxError {
        let src = self.src;
        let line_start = src[..offset].rfind('\n').map_or(0, |i| i + 1);
        let line_end = src[offset..].find('\n').map_or(src.len(), |i| offset + i);
        DotSyntaxError {
            offset,
            line: src[..offset].bytes().filter(|&b| b == b'\n').count() + 1,
            column: src[line_start..offset].chars().count() + 1,
            message,
            source_line: src[line_start..line_end].trim_end_matches('\r').to_owned(),
        }
    }

    fn expected(&self, what: &str) -> DotSyntaxError {
        let found = match self.token {
            Token::End => "end of input".to_owned(),
            _ => {
                let text = &self.src[self.start..self.end];
                match text.char_indices().nth(24) {
                    Some((cut, _)) => format!("`{}...`", &text[..cut]),
                    None => format!("`{text}`"),
                }
            }
        };
        let mut message = format!("expected {what}, found {found}");
        if self.glued_number > 0 {
            message.push_str(
                " (a number directly followed by letters, such as `18s`, reads as two \
                 identifiers; quote it)",
            );
        }
        self.error_at(self.start, message)
    }

    /// Skips whitespace and comments, leaving `end` at the next token.
    fn skip_trivia(&mut self) -> Result<(), DotSyntaxError> {
        let bytes = self.src.as_bytes();
        loop {
            while bytes.get(self.end).is_some_and(u8::is_ascii_whitespace) {
                self.end += 1;
            }
            let rest = &bytes[self.end..];
            if rest.starts_with(b"//") || rest.first() == Some(&b'#') {
                self.end += rest.iter().position(|&b| b == b'\n').unwrap_or(rest.len());
            } else if rest.starts_with(b"/*") {
                let Some(length) = self.src[self.end + 2..].find("*/") else {
                    return Err(self.error_at(self.end, "unterminated `/*` comment".to_owned()));
                };
                self.end += length + 4;
            } else {
                return Ok(());
            }
        }
    }

    /// Advances to the next token.
    fn bump(&mut self) -> Result<(), DotSyntaxError> {
        let previous_end = self.end;
        self.skip_trivia()?;
        let bytes = self.src.as_bytes();
        let start = self.end;
        let single = |token| (token, start + 1);
        let (token, end) = match bytes.get(start).copied() {
            None => (Token::End, start),
            Some(b'{') => single(Token::LBrace),
            Some(b'}') => single(Token::RBrace),
            Some(b'[') => single(Token::LBracket),
            Some(b']') => single(Token::RBracket),
            Some(b'=') => single(Token::Equals),
            Some(b':') => single(Token::Colon),
            Some(b';') => single(Token::Semicolon),
            Some(b',') => single(Token::Comma),
            Some(b'-') if bytes.get(start + 1) == Some(&b'>') => (Token::Arrow, start + 2),
            Some(b'-') if bytes.get(start + 1) == Some(&b'-') => (Token::Line, start + 2),
            Some(b'"') => {
                let mut i = start + 1;
                loop {
                    match bytes.get(i) {
                        None => {
                            return Err(
                                self.error_at(start, "unterminated quoted string".to_owned())
                            )
                        }
                        Some(b'"') => break,
                        Some(b'\\') => i += 2,
                        Some(_) => i += 1,
                    }
                }
                (Token::Id(IdKind::Quoted), i + 1)
            }
            Some(b'<') => {
                let (mut depth, mut i) = (0usize, start);
                loop {
                    match bytes.get(i) {
                        None => {
                            return Err(self.error_at(start, "unterminated HTML string".to_owned()))
                        }
                        Some(b'<') => depth += 1,
                        Some(b'>') => {
                            depth -= 1;
                            if depth == 0 {
                                break;
                            }
                        }
                        Some(_) => {}
                    }
                    i += 1;
                }
                (Token::Id(IdKind::Html), i + 1)
            }
            Some(c @ (b'-' | b'.' | b'0'..=b'9')) => {
                let digits = |from: usize| {
                    from + bytes[from..]
                        .iter()
                        .take_while(|b| b.is_ascii_digit())
                        .count()
                };
                let integer_start = start + usize::from(c == b'-');
                let mut i = digits(integer_start);
                let mut has_digits = i > integer_start;
                if bytes.get(i) == Some(&b'.') {
                    let fraction_end = digits(i + 1);
                    has_digits |= fraction_end > i + 1;
                    i = fraction_end;
                }
                if !has_digits {
                    return Err(self.error_at(start, format!("unexpected `{}`", char::from(c))));
                }
                (Token::Id(IdKind::Numeral), i)
            }
            Some(c) if c.is_ascii_alphabetic() || c == b'_' || c >= 0x80 => {
                let end = start
                    + bytes[start..]
                        .iter()
                        .take_while(|&&b| b.is_ascii_alphanumeric() || b == b'_' || b >= 0x80)
                        .count();
                let token = Token::keyword(&bytes[start..end]).unwrap_or(Token::Id(IdKind::Bare));
                (token, end)
            }
            Some(b'+') => {
                return Err(self.error_at(
                    start,
                    "string concatenation with `+` is not supported".to_owned(),
                ))
            }
            Some(_) => {
                let c = self.src[start..].chars().next().unwrap_or_default();
                return Err(self.error_at(start, format!("unexpected character `{c}`")));
            }
        };
        let glued = token == Token::Id(IdKind::Bare)
            && self.token == Token::Id(IdKind::Numeral)
            && previous_end == start;
        self.glued_number = if glued {
            2
        } else {
            self.glued_number.saturating_sub(1)
        };
        (self.token, self.start, self.end) = (token, start, end);
        Ok(())
    }

    fn eat(&mut self, token: Token) -> Result<bool, DotSyntaxError> {
        let found = self.token == token;
        if found {
            self.bump()?;
        }
        Ok(found)
    }

    fn take_id(&mut self) -> Result<Option<Cow<'a, str>>, DotSyntaxError> {
        let Token::Id(kind) = self.token else {
            return Ok(None);
        };
        let id = kind.decode(&self.src[self.start..self.end]);
        self.bump()?;
        Ok(Some(id))
    }

    fn expect_id(&mut self, what: &str) -> Result<Cow<'a, str>, DotSyntaxError> {
        match self.take_id()? {
            Some(id) => Ok(id),
            None => Err(self.expected(what)),
        }
    }

    /// `strict? (graph | digraph) ID? { stmt_list }`
    fn graph(&mut self) -> Result<FlatGraph<'a>, DotSyntaxError> {
        self.eat(Token::Strict)?;
        self.directed = match self.token {
            Token::Digraph => true,
            Token::Graph => false,
            _ => return Err(self.expected("`digraph` or `graph`")),
        };
        self.bump()?;
        let mut graph = FlatGraph {
            directed: self.directed,
            name: self.take_id()?,
            ..FlatGraph::default()
        };
        if !self.eat(Token::LBrace)? {
            return Err(self.expected("`{`"));
        }
        self.stmt_list(&mut graph, None)?;
        let FlatGraph { nodes, edges, .. } = &mut graph;
        for edge in edges.iter() {
            nodes.entry(edge.from.id.clone()).or_default();
            nodes.entry(edge.to.id.clone()).or_default();
        }
        Ok(graph)
    }

    /// Statements up to and including the closing `}`. Inside a subgraph,
    /// `members` collects the endpoints of its statements, repeats included,
    /// in case the subgraph is an edge endpoint.
    fn stmt_list(
        &mut self,
        graph: &mut FlatGraph<'a>,
        mut members: Option<&mut Vec<Endpoint<'a>>>,
    ) -> Result<(), DotSyntaxError> {
        loop {
            match self.token {
                Token::RBrace => return self.bump(),
                Token::End => return Err(self.expected("a statement or `}`")),
                _ => self.stmt(graph, members.as_deref_mut())?,
            }
            self.eat(Token::Semicolon)?;
        }
    }

    fn stmt(
        &mut self,
        graph: &mut FlatGraph<'a>,
        members: Option<&mut Vec<Endpoint<'a>>>,
    ) -> Result<(), DotSyntaxError> {
        let defaults = match self.token {
            Token::Graph => Some(&mut graph.graph_attrs),
            Token::Node => Some(&mut graph.node_defaults),
            Token::Edge => Some(&mut graph.edge_defaults),
            _ => None,
        };
        if let Some(defaults) = defaults {
            self.bump()?;
            if self.token != Token::LBracket {
                return Err(self.expected("`[` after an attribute statement keyword"));
            }
            defaults.extend(self.attr_lists()?);
            return Ok(());
        }
        let first = match self.token {
            Token::Subgraph | Token::LBrace => Side::Group(self.subgraph(graph)?),
            Token::Id(_) => {
                let id = self.expect_id("a node")?;
                if self.eat(Token::Equals)? {
                    let value = self.expect_id("a value after `=`")?;
                    graph.graph_attrs.push((id, value));
                    return Ok(());
                }
                Side::One(Endpoint {
                    id,
                    port: self.port()?,
                })
            }
            _ => return Err(self.expected("a statement or `}`")),
        };
        if matches!(self.token, Token::Arrow | Token::Line) {
            return self.edge_stmt(graph, first, members);
        }
        if let Side::One(node) = &first {
            let attrs = self.attr_lists()?;
            graph
                .nodes
                .entry(node.id.clone())
                .or_default()
                .extend(attrs);
        }
        if let Some(members) = members {
            members.extend_from_slice(first.endpoints());
        }
        Ok(())
    }

    /// The rest of an edge statement after its first endpoint group.
    fn edge_stmt(
        &mut self,
        graph: &mut FlatGraph<'a>,
        first: Side<'a>,
        members: Option<&mut Vec<Endpoint<'a>>>,
    ) -> Result<(), DotSyntaxError> {
        let mut sides = vec![first];
        while matches!(self.token, Token::Arrow | Token::Line) {
            if (self.token == Token::Arrow) != self.directed {
                let (found, kind, operator) = match self.directed {
                    true => ("--", "digraph", "->"),
                    false => ("->", "graph", "--"),
                };
                return Err(self.error_at(
                    self.start,
                    format!("`{found}` is not an edge operator in a {kind}; use `{operator}`"),
                ));
            }
            self.bump()?;
            sides.push(match self.token {
                Token::Subgraph | Token::LBrace => Side::Group(self.subgraph(graph)?),
                _ => {
                    let id = self.expect_id("an edge endpoint (a node or a subgraph)")?;
                    Side::One(Endpoint {
                        id,
                        port: self.port()?,
                    })
                }
            });
        }
        let mut attrs = self.attr_lists()?;
        let mut pairs = sides
            .windows(2)
            .flat_map(|pair| {
                let targets = pair[1].endpoints();
                (pair[0].endpoints().iter())
                    .flat_map(move |from| targets.iter().map(move |to| (from, to)))
            })
            .peekable();
        while let Some((from, to)) = pairs.next() {
            graph.edges.push(Edge {
                from: from.clone(),
                to: to.clone(),
                attrs: match pairs.peek() {
                    Some(_) => attrs.clone(),
                    None => std::mem::take(&mut attrs),
                },
            });
        }
        if let Some(members) = members {
            members.extend(sides.iter().flat_map(Side::endpoints).cloned());
        }
        Ok(())
    }

    /// `subgraph ID? { ... }` or `{ ... }`. Its statements are flattened into
    /// the enclosing graph; returns the endpoints it contributes to an edge.
    fn subgraph(&mut self, graph: &mut FlatGraph<'a>) -> Result<Vec<Endpoint<'a>>, DotSyntaxError> {
        if self.eat(Token::Subgraph)? {
            self.take_id()?;
        }
        let open = self.start;
        if !self.eat(Token::LBrace)? {
            return Err(self.expected("`{`"));
        }
        if self.depth == MAX_SUBGRAPH_DEPTH {
            return Err(self.error_at(
                open,
                format!("subgraphs are nested more than {MAX_SUBGRAPH_DEPTH} levels deep"),
            ));
        }
        self.depth += 1;
        let mut members = Vec::new();
        self.stmt_list(graph, Some(&mut members))?;
        self.depth -= 1;
        // Each member once, in order of first appearance. Removing repeats
        // once here rather than on every insertion keeps large subgraphs
        // linear.
        let mut seen = HashSet::with_capacity(members.len());
        members.retain(|member| seen.insert(member.clone()));
        Ok(members)
    }

    /// `(: ID)? (: compass)?`. A single component is a compass point when it is
    /// spelled like one, whatever its quoting, as in Graphviz. HTML strings
    /// keep their angle brackets, so they are never compass points.
    fn port(&mut self) -> Result<Option<Port<'a>>, DotSyntaxError> {
        if !self.eat(Token::Colon)? {
            return Ok(None);
        }
        let name = self.expect_id("a port name or compass point after `:`")?;
        if self.token == Token::Colon {
            self.bump()?;
            let at = self.start;
            let compass = (self.expect_id("a compass point after `:`")?.parse())
                .map_err(|error: UnknownCompassPt| self.error_at(at, error.to_string()))?;
            return Ok(Some(Port::Named(name, Some(compass))));
        }
        Ok(Some(match name.parse() {
            Ok(compass) => Port::Compass(compass),
            Err(_) => Port::Named(name, None),
        }))
    }

    /// Zero or more `[ key = value, ... ]` lists, concatenated.
    fn attr_lists(&mut self) -> Result<Vec<Attr<'a>>, DotSyntaxError> {
        let mut attrs = Vec::new();
        while self.eat(Token::LBracket)? {
            while !self.eat(Token::RBracket)? {
                let key = self.expect_id("an attribute name or `]`")?;
                if !self.eat(Token::Equals)? {
                    return Err(self.expected("`=` after an attribute name"));
                }
                let value = self.expect_id("an attribute value")?;
                attrs.push((key, value));
                if !self.eat(Token::Semicolon)? {
                    self.eat(Token::Comma)?;
                }
            }
        }
        Ok(attrs)
    }
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use super::*;

    fn one(src: &str) -> FlatGraph<'_> {
        FlatGraph::parse(src).unwrap_or_else(|error| panic!("{error}"))
    }

    fn error(src: &str) -> DotSyntaxError {
        FlatGraph::parse_all(src).expect_err(src)
    }

    fn attrs<'a>(pairs: &[(&'a str, &'a str)]) -> Vec<Attr<'a>> {
        (pairs.iter())
            .map(|&(key, value)| (key.into(), value.into()))
            .collect()
    }

    fn endpoints<'g>(graph: &'g FlatGraph<'_>) -> Vec<(&'g str, &'g str)> {
        (graph.edges.iter())
            .map(|edge| (edge.from.id.as_ref(), edge.to.id.as_ref()))
            .collect()
    }

    #[test]
    fn compass_points_round_trip_through_text_and_bytes() {
        for compass in CompassPt::ALL {
            assert_eq!(compass.to_string().parse(), Ok(compass));
            assert_eq!(CompassPt::try_from(compass as u8), Ok(compass));
        }
        assert!(CompassPt::ALL.is_sorted());
        assert_eq!(CompassPt::try_from(10), Err(10));
        assert_eq!(
            "N".parse::<CompassPt>().unwrap_err().to_string(),
            "`N` is not a compass point (n, ne, e, se, s, sw, w, nw, c or _)"
        );
    }

    /// dot-parser read `a:ne` as compass `n` followed by a node `e`.
    #[test]
    fn ports_are_compass_points_only_when_spelled_like_one() {
        let graph = one(
            r#"digraph { a:ne -> b:sw [id=0]; a:source -> b:0:s; c:"n" -> d:p:"_"; e:<n> -> f }"#,
        );
        let ports: Vec<_> = (graph.edges.iter())
            .map(|edge| (edge.from.port.clone(), edge.to.port.clone()))
            .collect();
        assert_eq!(
            ports,
            [
                (
                    Some(Port::Compass(CompassPt::NE)),
                    Some(Port::Compass(CompassPt::SW))
                ),
                (
                    Some(Port::Named("source".into(), None)),
                    Some(Port::Named("0".into(), Some(CompassPt::S)))
                ),
                (
                    Some(Port::Compass(CompassPt::N)),
                    Some(Port::Named("p".into(), Some(CompassPt::Underscore)))
                ),
                (Some(Port::Named("<n>".into(), None)), None),
            ]
        );
        assert_eq!(graph.edges[0].attrs, attrs(&[("id", "0")]));
        assert_eq!(graph.nodes.len(), 6);
        assert_eq!(
            error("digraph { a:p:<n> -> b }").message,
            "`<n>` is not a compass point (n, ne, e, se, s, sw, w, nw, c or _)"
        );
        assert!(Port::Named("z".into(), None) < Port::Compass(CompassPt::N));
    }

    /// `stringify!` puts spaces between all tokens of the `dot!` macro.
    #[test]
    fn token_spaced_input_parses_like_compact_input() {
        assert_eq!(
            one("digraph { 0 : 1 : s -> 1 : 0 [ id = 0 ] ; }"),
            one("digraph{0:1:s->1:0[id=0];}")
        );
    }

    #[test]
    fn quoted_strings_decode_quotes_backslashes_and_continuations() {
        let graph = one(concat!(
            "digraph { a [q=\"say \\\"hi\\\"\", end=\"x\\\\\", two=\"a\\\\\\\\b\", ",
            "keep=\"\\N\\l\\n\", cont=\"line\\\ncontinued\", html=<<b>x</b>>, multi=\"x\n  y\"] }"
        ));
        assert_eq!(
            graph.nodes["a"],
            attrs(&[
                ("q", "say \"hi\""),
                ("end", "x\\"),
                ("two", "a\\\\b"),
                ("keep", "\\N\\l\\n"),
                ("cont", "linecontinued"),
                ("html", "<<b>x</b>>"),
                ("multi", "x\n  y"),
            ])
        );
        assert!(matches!(graph.nodes["a"][0].1, Cow::Owned(_)));
        assert!(matches!(graph.nodes["a"][6].1, Cow::Borrowed(_)));
    }

    #[test]
    fn keywords_are_case_insensitive() {
        let graph = one("DiGraph G { Node [x=1]; EDGE [y=2]; GRAPH [z=3]; a -> b }");
        assert!(graph.directed);
        assert_eq!(graph.name.as_deref(), Some("G"));
        assert_eq!(graph.node_defaults, attrs(&[("x", "1")]));
        assert_eq!(graph.edge_defaults, attrs(&[("y", "2")]));
        assert_eq!(graph.graph_attrs, attrs(&[("z", "3")]));
        assert!(error("digraph { Node -> b }")
            .message
            .starts_with("expected `[`"));
    }

    /// Graphviz skips the rest of a line from `#` wherever it appears, not
    /// only at the start of a line.
    #[test]
    fn comments_and_numerals() {
        let graph = one(concat!(
            "# a preprocessor line\n",
            "strict digraph { // a comment\n",
            "/* a block\ncomment */ a -> -1.5; b -> .5 # -> x\n; c -> 1.; 1H; é -> \"#\" }",
        ));
        assert_eq!(
            endpoints(&graph),
            [("a", "-1.5"), ("b", ".5"), ("c", "1."), ("é", "#")]
        );
        assert!(graph.nodes.contains_key("1") && graph.nodes.contains_key("H"));
        assert!(!graph.nodes.contains_key("x"));
    }

    /// Graphviz merges repeated node statements and identifies `a` with `"a"`.
    #[test]
    fn repeated_node_statements_merge_by_decoded_identifier() {
        let graph = one(r#"digraph { a [x=1, y=1]; "a" [y=2]; "a" -> a }"#);
        assert_eq!(graph.nodes.len(), 1);
        assert_eq!(
            graph.nodes["a"],
            attrs(&[("x", "1"), ("y", "1"), ("y", "2")])
        );
        assert_eq!(endpoints(&graph), [("a", "a")]);
    }

    #[test]
    fn subgraphs_flatten_and_edges_expand() {
        let graph = one(concat!(
            "digraph { a -> {b c} -> d [w=1]; subgraph s { node [k=v]; e [z=1]; x = y } ",
            "{f -> g [inner=1]; f; {g; f -> f2}} -> h; i -> subgraph { j } }"
        ));
        assert_eq!(
            endpoints(&graph),
            [
                ("a", "b"),
                ("a", "c"),
                ("b", "d"),
                ("c", "d"),
                ("f", "g"),
                ("f", "f2"),
                ("f", "h"),
                ("g", "h"),
                ("f2", "h"),
                ("i", "j"),
            ]
        );
        assert!(graph.edges[..4]
            .iter()
            .all(|e| e.attrs == attrs(&[("w", "1")])));
        assert_eq!(graph.edges[4].attrs, attrs(&[("inner", "1")]));
        assert_eq!(graph.node_defaults, attrs(&[("k", "v")]));
        assert_eq!(graph.graph_attrs, attrs(&[("x", "y")]));
        assert_eq!(graph.nodes["e"], attrs(&[("z", "1")]));
    }

    /// A subgraph's members are collected with repeats and deduplicated once
    /// it closes; checking each insertion made this quadratic in its size.
    #[test]
    fn large_subgraphs_parse_in_linear_time() {
        let edges: String = (0..20_000).map(|i| format!("n{i} -> m{i}; ")).collect();
        let fastest = |src: String| {
            let runs = (0..3).map(|_| {
                let start = std::time::Instant::now();
                assert_eq!(FlatGraph::parse(&src).unwrap().edges.len(), 20_000);
                start.elapsed()
            });
            runs.min().unwrap()
        };
        let flat = fastest(format!("digraph {{ {edges} }}"));
        let nested = fastest(format!("digraph {{ subgraph s {{ {edges} }} }}"));
        assert!(
            nested < 5 * flat,
            "{nested:?} against {flat:?} without the subgraph"
        );
    }

    #[test]
    fn attribute_lists_concatenate() {
        let graph =
            one("digraph { node [k=1][l=2]; a [x=1] [y=2]; a -> b [p=1][q=2]; graph [g=1][h=2] }");
        assert_eq!(graph.node_defaults, attrs(&[("k", "1"), ("l", "2")]));
        assert_eq!(graph.nodes["a"], attrs(&[("x", "1"), ("y", "2")]));
        assert_eq!(graph.edges[0].attrs, attrs(&[("p", "1"), ("q", "2")]));
        assert_eq!(graph.graph_attrs, attrs(&[("g", "1"), ("h", "2")]));
    }

    /// `ID = ID` and `graph [ID = ID]` set the same graph attributes, so they
    /// share one list in source order.
    #[test]
    fn graph_attributes_keep_source_order() {
        let graph = one("digraph { k = 1; graph [k=2, j=1]; j = 2 }");
        assert_eq!(
            graph.graph_attrs,
            attrs(&[("k", "1"), ("k", "2"), ("j", "1"), ("j", "2")])
        );
    }

    #[test]
    fn edge_operators_follow_the_graph_kind() {
        assert!(!one("graph { a -- b -- c }").directed);
        let e = error("digraph {\n  a -- b }");
        assert_eq!((e.line, e.column), (2, 5));
        assert_eq!(
            e.message,
            "`--` is not an edge operator in a digraph; use `->`"
        );
        assert!(error("graph { a -> b }").message.contains("use `--`"));
    }

    #[test]
    fn single_graph_input_rejects_trailing_input() {
        let e = FlatGraph::parse("digraph a {}\ndigraph b {}").unwrap_err();
        assert_eq!((e.line, e.column, e.offset), (2, 1, 13));
        assert_eq!(
            e.message,
            "expected the end of input after the graph, found `digraph`"
        );
        assert!(FlatGraph::parse("digraph {} }").is_err());
        let graphs = FlatGraph::parse_all("digraph a {}\n graph b {}").unwrap();
        assert_eq!(graphs.len(), 2);
        assert_eq!(
            error(" ").message,
            "expected `digraph` or `graph`, found end of input"
        );
    }

    #[test]
    fn rejected_forms_have_dedicated_messages() {
        for (src, message) in [
            (
                r#"digraph { a [l="x" + "y"] }"#,
                "string concatenation with `+` is not supported",
            ),
            ("digraph { a [l=\"x] }", "unterminated quoted string"),
            ("digraph { a /* b }", "unterminated `/*` comment"),
            ("digraph { a [l=<<b>x] }", "unterminated HTML string"),
            ("digraph { a:p:x -> b }", "`x` is not a compass point"),
            ("digraph { a -> node }", "expected an edge endpoint"),
            ("digraph { a -> b [id=18s] }", "reads as two identifiers"),
            (
                "digraph { a -> }",
                "expected an edge endpoint (a node or a subgraph), found `}`",
            ),
            ("digraph { a; @ }", "unexpected character `@`"),
            (
                "digraph { a [x] }",
                "expected `=` after an attribute name, found `]`",
            ),
            (
                r#"digraph { a [x "abcdefghijklmnopqrstuvwxyz"] }"#,
                r#"found `"abcdefghijklmnopqrstuvw...`"#,
            ),
        ] {
            let e = error(src);
            assert!(e.message.contains(message), "{src}: {e}");
        }
    }

    #[test]
    fn only_glued_numbers_get_the_two_identifiers_hint() {
        let hint = |src| error(src).message.contains("reads as two identifiers");
        assert!(hint("digraph { a [w=1label] }"));
        assert!(!hint("digraph { a [w=1 label] }"));
        assert!(!hint("digraph { a [w=1\nlabel] }"));
        assert!(!hint("digraph { a [w=1/* */label] }"));
    }

    fn nested(depth: usize) -> String {
        format!(
            "digraph {{ {}a{} }}",
            "{ ".repeat(depth),
            " }".repeat(depth)
        )
    }

    #[test]
    fn subgraph_nesting_is_limited() {
        assert!(FlatGraph::parse(&nested(MAX_SUBGRAPH_DEPTH)).is_ok());
        let e = FlatGraph::parse(&nested(MAX_SUBGRAPH_DEPTH + 1)).unwrap_err();
        assert_eq!(e.message, "subgraphs are nested more than 256 levels deep");
        assert_eq!(e.offset, "digraph { ".len() + 2 * MAX_SUBGRAPH_DEPTH);
        assert_eq!(
            FlatGraph::parse(&nested(100_000)).unwrap_err().offset,
            e.offset
        );
    }

    #[test]
    fn errors_show_the_line_with_a_caret() {
        let e = error("digraph {\n\tx [a=1];\n\ta -> ;\n}");
        assert_eq!(
            e.to_string(),
            "line 3, column 7: expected an edge endpoint (a node or a subgraph), found `;`\n  \
             | \ta -> ;\n  | \t     ^"
        );
        let long = format!("digraph {{ {0} a -> ; {0} }}", "x;".repeat(100));
        let shown = error(&long).to_string();
        let lines: Vec<&str> = shown.lines().collect();
        assert!(lines[1].starts_with("  | ...") && lines[1].ends_with("..."));
        assert_eq!(lines[2].find('^'), lines[1].find("-> ;").map(|i| i + 3));
    }

    const FRAGMENTS: &[&str] = &[
        "digraph",
        "graph",
        "strict",
        "subgraph",
        "node",
        "edge",
        "{",
        "}",
        "[",
        "]",
        "=",
        ":",
        ";",
        ",",
        "->",
        "--",
        "a",
        "b1",
        "\"q\\\"x\"",
        "\"",
        "\\",
        "<",
        ">",
        "<<b>>",
        "#",
        "\n",
        " ",
        "\t",
        "//",
        "/*",
        "*/",
        "-",
        ".",
        "1",
        "-2.5",
        "é",
        "+",
        "@",
        "n",
        "sw",
        "\"\\\n\"",
    ];

    /// Error positions agree with the input around them.
    fn check_error_location(src: &str, e: &DotSyntaxError) {
        assert!(
            e.offset <= src.len() && src.is_char_boundary(e.offset),
            "{e:?}"
        );
        let before = &src[..e.offset];
        assert_eq!(e.line, before.matches('\n').count() + 1);
        let line_start = before.rfind('\n').map_or(0, |i| i + 1);
        assert_eq!(e.column, before[line_start..].chars().count() + 1);
        let line = src[line_start..].lines().next().unwrap_or_default();
        assert_eq!(e.source_line, line.trim_end_matches('\r'));
        let _ = e.to_string();
    }

    fn graph_text() -> impl Strategy<Value = (Vec<(usize, usize)>, String)> {
        (prop::collection::vec((0..6usize, 0..6usize, any::<bool>()), 1..12)).prop_map(|edges| {
            let statements: Vec<String> = (edges.iter())
                .map(|&(from, to, quoted)| match quoted {
                    true => format!("\"n{from}\" -> n{to} [w = \"{from}\\\"{to}\"]"),
                    false => format!("n{from}:p -> n{to}:s"),
                })
                .collect();
            let pairs = edges.iter().map(|&(from, to, _)| (from, to)).collect();
            (pairs, format!("digraph g {{ {} }}", statements.join(" ; ")))
        })
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(512))]

        #[test]
        fn arbitrary_input_never_panics(
            fragments in prop::collection::vec(prop::sample::select(FRAGMENTS), 0..48)
        ) {
            let src = fragments.concat();
            match FlatGraph::parse(&src) {
                Ok(graph) => prop_assert_eq!(FlatGraph::parse_all(&src).unwrap(), vec![graph]),
                Err(e) => check_error_location(&src, &e),
            }
            if let Err(e) = FlatGraph::parse_all(&src) {
                check_error_location(&src, &e);
            }
        }

        #[test]
        fn generated_graphs_keep_their_edges((pairs, src) in graph_text()) {
            let graph = FlatGraph::parse(&src).unwrap();
            let expected: Vec<(String, String)> = (pairs.iter())
                .map(|(from, to)| (format!("n{from}"), format!("n{to}")))
                .collect();
            let found: Vec<(String, String)> = (endpoints(&graph).into_iter())
                .map(|(from, to)| (from.to_owned(), to.to_owned()))
                .collect();
            prop_assert_eq!(found, expected);
            for (edge, (from, to)) in graph.edges.iter().zip(&pairs) {
                if let [(_, value)] = edge.attrs.as_slice() {
                    prop_assert_eq!(value.as_ref(), format!("{from}\"{to}"));
                }
            }
        }

        #[test]
        fn an_unexpected_character_is_located(
            (_, src) in graph_text(),
            at in any::<prop::sample::Index>(),
        ) {
            let spaces: Vec<usize> = (src.match_indices(' '))
                .map(|(i, _)| i)
                .filter(|&i| src[..i].matches('"').count() % 2 == 0 && !src[..i].ends_with('\\'))
                .collect();
            let offset = spaces[at.index(spaces.len())];
            let src = format!("{}\n@{}", &src[..offset], &src[offset..]);
            let e = FlatGraph::parse(&src).unwrap_err();
            prop_assert_eq!(e.message.as_str(), "unexpected character `@`");
            prop_assert_eq!((e.offset, e.line, e.column), (offset + 1, 2, 1));
            check_error_location(&src, &e);
        }
    }
}
