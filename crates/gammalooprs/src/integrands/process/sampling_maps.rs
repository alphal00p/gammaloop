//! Symbolica-validated descriptions of graph-aware sampling maps.
//!
//! This module intentionally contains no graph or numerical-sampling code yet.  It
//! provides the stable, topology-independent language that later map kernels can
//! consume.  Parsing is structural: a complete Symbolica expression is checked
//! recursively, so a valid nested call cannot hide an invalid root or sibling.

use color_eyre::eyre::{Result, eyre};
use serde::{Deserialize, Serialize};
use symbolica::{
    atom::{Atom, AtomView},
    symbol, try_parse,
};

/// The geometric role of a resolved sampling map.
#[derive(Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SamplingMapDefinition {
    /// A loop-momentum basis channel, identified by its ordered edge list.
    Lmb(Vec<usize>),
    /// An energy surface in the supplied edge subspace.
    Surface(Vec<usize>),
    /// A physical Cutkosky cut.
    Cut(Vec<usize>),
    /// A soft channel around one routed edge.
    Soft(usize),
    /// A qualified collinear channel between two routed edges.
    Collinear(usize, usize),
    /// Coordinates left for a complement map.
    Complement(Vec<usize>),
    /// Independent maps on a direct-sum block.
    Product(Vec<Self>),
    /// Joint coordinates for compatible constraints.
    Intersect(Vec<Self>),
    /// Ordered conditional maps.  Later maps may depend on earlier data.
    Then(Vec<Self>),
    /// Prepare physical cut phase-space data before mapping the remaining maps.
    PhaseSpace(Box<Self>),
    /// Resolve a target on the left side of a prepared cut.
    Left(Box<Self>),
    /// Resolve a target on the right side of a prepared cut.
    Right(Box<Self>),
}

impl SamplingMapDefinition {
    /// Parse and structurally validate one complete Symbolica map expression.
    pub fn parse(source: &str) -> Result<Self> {
        let atom = try_parse!(source.trim())
            .map_err(|error| eyre!("failed to parse sampling map `{source}`: {error}"))?;
        Self::from_atom(atom.as_view())
    }

    /// Build a map from an already parsed Symbolica atom.
    pub fn from_atom(atom: AtomView<'_>) -> Result<Self> {
        let AtomView::Fun(function) = atom else {
            return Err(eyre!(
                "sampling map must be a constructor call, got `{atom}`"
            ));
        };
        let symbol = function.get_symbol();
        let name = symbol.get_name();
        let arguments = function.iter().collect::<Vec<_>>();
        let short_name = name.rsplit("::").next().unwrap_or(name);
        // `try_parse!` resolves unqualified constructors in the current
        // crate namespace (`gammalooprs::surface`), while canonical atoms
        // emitted by `to_atom` use the dedicated sampling-map namespace.
        // Accept both forms but reject unrelated user namespaces.
        if name != short_name
            && name != format!("gammalooprs::{short_name}")
            && name != format!("gammalooprs::sampling_map::{short_name}")
        {
            return Err(eyre!(
                "sampling-map constructor `{name}` is outside the reserved sampling-map namespace"
            ));
        }

        match short_name {
            "lmb" => Self::ordered_edges(arguments, "lmb").map(Self::Lmb),
            "surface" => Self::edges(arguments, "surface").map(Self::Surface),
            "cut" => Self::edges(arguments, "cut").map(Self::Cut),
            "complement" => Self::edges(arguments, "complement").map(Self::Complement),
            "soft" => Self::one_edge(arguments, "soft").map(Self::Soft),
            "collinear" => {
                Self::two_edges(arguments, "collinear").map(|(a, b)| Self::Collinear(a, b))
            }
            "product" => Self::maps(arguments, "product").map(Self::Product),
            "intersect" => Self::maps(arguments, "intersect").map(Self::Intersect),
            "then" => Self::maps(arguments, "then").map(Self::Then),
            "phase_space" => {
                Self::one_map(arguments, "phase_space").map(|map| Self::PhaseSpace(Box::new(map)))
            }
            "left" => Self::one_map(arguments, "left").map(|map| Self::Left(Box::new(map))),
            "right" => Self::one_map(arguments, "right").map(|map| Self::Right(Box::new(map))),
            _ => Err(eyre!(
                "unknown sampling-map constructor `{name}`; expected lmb, surface, cut, soft, collinear, complement, product, intersect, then, phase_space, left or right"
            )),
        }
    }

    /// Return a canonical Symbolica expression for this definition.
    pub fn to_atom(&self) -> Atom {
        match self {
            Self::Lmb(edges) => call(
                "lmb",
                edges.iter().copied().map(|edge| Atom::num(edge as i64)),
            ),
            Self::Surface(edges) => call(
                "surface",
                edges.iter().copied().map(|edge| Atom::num(edge as i64)),
            ),
            Self::Cut(edges) => call(
                "cut",
                edges.iter().copied().map(|edge| Atom::num(edge as i64)),
            ),
            Self::Soft(edge) => call("soft", [Atom::num(*edge as i64)]),
            Self::Collinear(a, b) => {
                call("collinear", [Atom::num(*a as i64), Atom::num(*b as i64)])
            }
            Self::Complement(edges) => call(
                "complement",
                edges.iter().copied().map(|edge| Atom::num(edge as i64)),
            ),
            Self::Product(maps) => call("product", maps.iter().map(Self::to_atom)),
            Self::Intersect(maps) => call("intersect", maps.iter().map(Self::to_atom)),
            Self::Then(maps) => call("then", maps.iter().map(Self::to_atom)),
            Self::PhaseSpace(map) => call("phase_space", [map.to_atom()]),
            Self::Left(map) => call("left", [map.to_atom()]),
            Self::Right(map) => call("right", [map.to_atom()]),
        }
    }

    fn edges(arguments: Vec<AtomView<'_>>, constructor: &str) -> Result<Vec<usize>> {
        if arguments.is_empty() {
            return Err(eyre!("{constructor} expects at least one edge"));
        }
        let mut edges = arguments
            .into_iter()
            .map(|argument| parse_edge(argument, constructor))
            .collect::<Result<Vec<_>>>()?;
        edges.sort_unstable();
        if edges.windows(2).any(|pair| pair[0] == pair[1]) {
            return Err(eyre!(
                "{constructor} contains a duplicate edge in {edges:?}"
            ));
        }
        Ok(edges)
    }

    fn ordered_edges(arguments: Vec<AtomView<'_>>, constructor: &str) -> Result<Vec<usize>> {
        if arguments.is_empty() {
            return Err(eyre!("{constructor} expects at least one edge"));
        }
        let edges = arguments
            .into_iter()
            .map(|argument| parse_edge(argument, constructor))
            .collect::<Result<Vec<_>>>()?;
        let mut sorted_edges = edges.clone();
        sorted_edges.sort_unstable();
        if sorted_edges.windows(2).any(|pair| pair[0] == pair[1]) {
            return Err(eyre!(
                "{constructor} contains a duplicate edge in {edges:?}"
            ));
        }
        Ok(edges)
    }

    fn one_edge(arguments: Vec<AtomView<'_>>, constructor: &str) -> Result<usize> {
        match arguments.as_slice() {
            [edge] => parse_edge(*edge, constructor),
            _ => Err(eyre!(
                "{constructor} expects exactly one edge, got {} arguments",
                arguments.len()
            )),
        }
    }

    fn two_edges(arguments: Vec<AtomView<'_>>, constructor: &str) -> Result<(usize, usize)> {
        match arguments.as_slice() {
            [a, b] => {
                let a = parse_edge(*a, constructor)?;
                let b = parse_edge(*b, constructor)?;
                if a == b {
                    return Err(eyre!("{constructor} requires two distinct edges, got {a}"));
                }
                Ok((a, b))
            }
            _ => Err(eyre!(
                "{constructor} expects exactly two edges, got {} arguments",
                arguments.len()
            )),
        }
    }

    fn maps(arguments: Vec<AtomView<'_>>, constructor: &str) -> Result<Vec<Self>> {
        if arguments.is_empty() {
            return Err(eyre!("{constructor} expects at least one map"));
        }
        arguments.into_iter().map(Self::from_atom).collect()
    }

    fn one_map(arguments: Vec<AtomView<'_>>, constructor: &str) -> Result<Self> {
        match arguments.as_slice() {
            [map] => Self::from_atom(*map),
            _ => Err(eyre!(
                "{constructor} expects exactly one map, got {} arguments",
                arguments.len()
            )),
        }
    }
}

/// How a numerical kernel establishes the domain of its push-forward density.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SamplingSupport {
    /// The map covers the full raw domain, including absent-surface fallback data.
    Full,
    /// The map is valid only after an earlier conditional map has prepared data.
    Conditional,
    /// The map has explicitly enumerated regular branches.
    Branched,
}

/// Contract for the Jacobian supplied by a map kernel.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SamplingJacobian {
    /// Exact determinant of the selected forward map.
    ExactForward,
    /// Exact determinant obtained by implicit differentiation of a solved root.
    ExactImplicit,
    /// A positive partition score only; never the selected channel's integration Jacobian.
    ProxyOnly,
}

/// Static contract carried by a future compiled sampling-map kernel.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct SamplingMapContract {
    pub support: SamplingSupport,
    pub jacobian: SamplingJacobian,
}

fn parse_edge(atom: AtomView<'_>, constructor: &str) -> Result<usize> {
    let value = i64::try_from(atom)
        .map_err(|_| eyre!("{constructor} expects non-negative integer edge IDs, got `{atom}`"))?;
    usize::try_from(value)
        .map_err(|_| eyre!("{constructor} expects non-negative integer edge IDs, got {value}"))
}

fn call<I>(name: &str, arguments: I) -> Atom
where
    I: IntoIterator<Item = Atom>,
{
    symbol!(&format!("gammalooprs::sampling_map::{name}")).call_args(arguments)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_nested_definition_and_round_trips() {
        let source =
            "then(phase_space(cut(10, 2)), product(left(surface(12, 4)), right(complement(7, 8))))";
        let parsed = SamplingMapDefinition::parse(source).expect("valid map");
        assert_eq!(
            parsed,
            SamplingMapDefinition::Then(vec![
                SamplingMapDefinition::PhaseSpace(Box::new(SamplingMapDefinition::Cut(vec![
                    2, 10
                ]))),
                SamplingMapDefinition::Product(vec![
                    SamplingMapDefinition::Left(Box::new(SamplingMapDefinition::Surface(vec![
                        4, 12
                    ]))),
                    SamplingMapDefinition::Right(Box::new(SamplingMapDefinition::Complement(
                        vec![7, 8]
                    ))),
                ]),
            ])
        );
        let reparsed =
            SamplingMapDefinition::from_atom(parsed.to_atom().as_view()).expect("round trip");
        assert_eq!(parsed, reparsed);
    }

    #[test]
    fn rejects_unknown_root_and_invalid_nested_map() {
        for source in [
            "foo(surface(1))",
            "product(surface(1), foo(2))",
            "surface(1) + soft(2)",
            "user_defined::surface(1)",
        ] {
            assert!(
                SamplingMapDefinition::parse(source).is_err(),
                "accepted {source}"
            );
        }
    }

    #[test]
    fn rejects_bad_arities_and_duplicate_edges() {
        for source in [
            "surface()",
            "soft(1, 2)",
            "collinear(3, 3)",
            "then()",
            "phase_space(cut(1), cut(2))",
            "surface(2, 2)",
            "surface(-1)",
        ] {
            assert!(
                SamplingMapDefinition::parse(source).is_err(),
                "accepted {source}"
            );
        }
        assert!(SamplingMapDefinition::parse("lmb(1, 3, 1)").is_err());
    }

    #[test]
    fn preserves_parent_lmb_order_but_canonicalizes_edge_sets() {
        assert_eq!(
            SamplingMapDefinition::parse("lmb(4, 2, 7)").expect("valid lmb"),
            SamplingMapDefinition::Lmb(vec![4, 2, 7])
        );
        assert_eq!(
            SamplingMapDefinition::parse("surface(7, 2, 4)").expect("valid surface"),
            SamplingMapDefinition::Surface(vec![2, 4, 7])
        );
    }
}
