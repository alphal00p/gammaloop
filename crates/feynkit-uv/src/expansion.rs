//! The forest formula written term by term.
//!
//! A forest term is `± m(Γ|γ) · Op[m(S|G) · inner] ⋯`: the cofactor `m(Γ|γ)` is
//! the diagram with its outermost elements `γ` treated, and each factor
//! approximates an outermost element `S` after the elements `G` nested in it
//! were treated by the factors in `inner`. An element with nothing nested in
//! it is `Op[m(S|∅)]`. The operator `Op` is any symbol, `K` by default.
//!
//! The spinney of the whole graph is named `Γ` and the others `a`, `b`, … in
//! [`Wood::nodes_by_edges`] order, so a factor reads `K_{a|bc}`: spinney `a`
//! given `b` and `c`. [`Expansion::rename`] sets other names. Terms are numbered
//! by their stored quantities `k₀, k₁, …` in topological order: by depth, then
//! by their elements, outermost and largest first. The empty forest is
//! `K[∅] = k₀`.

use std::cmp::Reverse;
use std::collections::BTreeMap;
use std::slice;

use feynkit_graph::{FeynmanDiagram, LoopMomentumBasis};
use linnet::half_edge::NodeIndex;
use linnet::half_edge::subgraph::{Inclusion, SuBitGraph};
use thiserror::Error;

use crate::{Forests, Wood};

/// One factor `Op[m(S|G) · inner]` of a forest term.
///
/// See [`Expansion`] for an example.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TermFactor {
    /// Wood node of the approximated spinney `S`.
    pub spinney: NodeIndex,
    /// Wood nodes of the treated spinneys `G`: the elements nested directly in
    /// `S`, which `inner` approximates.
    pub given: Vec<NodeIndex>,
    /// The factors of the nested term, one per element of `given`.
    pub inner: Vec<TermFactor>,
}

/// One term of the forest formula.
///
/// See [`Expansion`] for an example.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ForestTerm {
    /// Node of the forest in its [`Forests`].
    pub forest: NodeIndex,
    /// Index of the stored quantity `k`.
    pub k: usize,
    /// Sign `(-1)^depth`.
    pub sign: i8,
    /// Wood nodes of the outermost elements, treated in the cofactor `m(Γ|γ)`.
    pub outermost: Vec<NodeIndex>,
    /// The factors, one per outermost element; none for the empty forest.
    pub factors: Vec<TermFactor>,
}

/// The forest formula of a diagram, term by term in `k` order.
///
/// ```
/// use feynkit_graph::FeynmanDiagram;
/// use feynkit_model::Model;
/// use feynkit_uv::{RenormalizationPrescriptionSettings, UvGraph, Wood};
///
/// // The vacuum sunrise of a cubic theory: three overlapping bubbles a, b, c
/// // inside the quadratically divergent whole graph Γ.
/// let sunrise = FeynmanDiagram::from_dot(
///     Model::phi3(),
///     "digraph sunrise { edge [particle=phi]; a -> b; a -> b; a -> b; }",
/// )?;
/// let spinneys = sunrise.classify_spinneys(
///     &sunrise.hedge_graph().full_filter(),
///     sunrise.model(),
///     &RenormalizationPrescriptionSettings::default(),
///     sunrise.loop_momentum_basis(),
/// )?;
/// let forests = Wood::from_spinneys(spinneys, &sunrise, sunrise.loop_momentum_basis())?.unfold();
/// let expansion = forests.expansion(&sunrise);
/// let table: Vec<_> = expansion.k_table("K");
/// assert_eq!(
///     table,
///     [
///         "k0 = K[∅]",
///         "k1 = K_{Γ}",
///         "k2 = K_{a}",
///         "k3 = K_{b}",
///         "k4 = K_{c}",
///         "k5 = K_{Γ|a}[k2]",
///         "k6 = K_{Γ|b}[k3]",
///         "k7 = K_{Γ|c}[k4]",
///     ]
/// );
/// let signs: Vec<_> = expansion.terms().iter().map(|term| term.sign).collect();
/// assert_eq!(signs, [1, -1, -1, -1, -1, 1, 1, 1]);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Clone, Debug)]
pub struct Expansion {
    names: BTreeMap<NodeIndex, String>,
    /// Wood node of the whole graph, when it is a spinney.
    whole: Option<NodeIndex>,
    terms: Vec<ForestTerm>,
}

impl Forests<LoopMomentumBasis> {
    /// The forest formula of these forests of `diagram`, term by term.
    ///
    /// `diagram` must be the diagram the forests were built for; it is not
    /// checked. See [`Expansion`] for an example.
    pub fn expansion(&self, diagram: &FeynmanDiagram) -> Expansion {
        let wood = self.wood();
        let order: Vec<NodeIndex> = wood
            .nodes_by_edges(diagram)
            .into_iter()
            .filter(|node| !wood.spinney(*node).is_empty())
            .collect();
        let rank: BTreeMap<NodeIndex, usize> = order
            .iter()
            .enumerate()
            .map(|(index, node)| (*node, index))
            .collect();
        let internal = diagram.internal_subgraph();
        let whole = order
            .iter()
            .copied()
            .find(|node| wood.spinney(*node).filter() == &internal);
        let mut letters = 0;
        let names: BTreeMap<NodeIndex, String> = order
            .iter()
            .map(|node| {
                let name = if Some(*node) == whole {
                    WHOLE_GRAPH.to_owned()
                } else {
                    letters += 1;
                    spinney_name(letters - 1)
                };
                (*node, name)
            })
            .collect();
        let node_of: BTreeMap<&SuBitGraph, NodeIndex> = wood
            .spinneys()
            .map(|(node, spinney)| (spinney.filter(), node))
            .collect();
        let mut terms: Vec<(_, ForestTerm)> = self
            .operations()
            .map(|(forest, operation)| {
                let elements: Vec<NodeIndex> = self
                    .forest_spinneys(forest)
                    .into_iter()
                    .map(|spinney| node_of[spinney.filter()])
                    .collect();
                let mut key: Vec<_> = elements
                    .iter()
                    .map(|node| {
                        let spinney = wood.spinney(*node);
                        (Reverse(spinney.loop_count()), spinney.edges(diagram))
                    })
                    .collect();
                key.sort();
                let outermost = maximal(wood, &rank, &elements, None);
                let factors = outermost
                    .iter()
                    .map(|node| factor(wood, &rank, &elements, *node))
                    .collect();
                (
                    (operation.depth(), key),
                    ForestTerm {
                        forest,
                        k: 0,
                        sign: operation.sign(),
                        outermost,
                        factors,
                    },
                )
            })
            .collect();
        terms.sort_by(|(a, _), (b, _)| a.cmp(b));
        let terms = terms
            .into_iter()
            .enumerate()
            .map(|(k, (_, mut term))| {
                term.k = k;
                term
            })
            .collect();
        Expansion {
            names,
            whole,
            terms,
        }
    }
}

/// The elements of `elements` strictly inside `within` (all when `None`) that
/// no other such element contains, in naming order.
fn maximal(
    wood: &Wood<LoopMomentumBasis>,
    rank: &BTreeMap<NodeIndex, usize>,
    elements: &[NodeIndex],
    within: Option<NodeIndex>,
) -> Vec<NodeIndex> {
    let inside = |node: NodeIndex, outer: NodeIndex| {
        node != outer
            && wood
                .spinney(outer)
                .subgraph()
                .includes(wood.spinney(node).subgraph())
    };
    let candidates: Vec<NodeIndex> = elements
        .iter()
        .copied()
        .filter(|node| within.is_none_or(|outer| inside(*node, outer)))
        .collect();
    let mut maximal: Vec<NodeIndex> = candidates
        .iter()
        .copied()
        .filter(|node| !candidates.iter().any(|other| inside(*node, *other)))
        .collect();
    maximal.sort_by_key(|node| rank.get(node).copied());
    maximal
}

fn factor(
    wood: &Wood<LoopMomentumBasis>,
    rank: &BTreeMap<NodeIndex, usize>,
    elements: &[NodeIndex],
    node: NodeIndex,
) -> TermFactor {
    let given = maximal(wood, rank, elements, Some(node));
    let inner = given
        .iter()
        .map(|child| factor(wood, rank, elements, *child))
        .collect();
    TermFactor {
        spinney: node,
        given,
        inner,
    }
}

/// Name of the spinney of the whole graph.
pub const WHOLE_GRAPH: &str = "Γ";

/// A rejected spinney name passed to [`Expansion::rename`].
#[derive(Clone, Debug, Error, PartialEq, Eq)]
pub enum RenameError {
    /// The wood node is no nonempty spinney of the expansion.
    #[error("wood node {0:?} is not a named spinney of this expansion")]
    UnknownSpinney(NodeIndex),
    /// The name is empty, holds whitespace, control characters or characters
    /// the notation uses (`|`, braces, brackets, `,`, `·`, `∅`), or reads as a
    /// stored quantity `k0`, `k1`, ….
    #[error(
        "spinney name {0:?} must be nonempty, without whitespace, control characters, \
         '|', braces, brackets, ',', '·' or '∅', and not of the form k<digits>"
    )]
    InvalidName(String),
    /// `Γ` names only the whole graph.
    #[error("spinney name \"Γ\" is reserved for the whole graph")]
    ReservedName,
    /// The same wood node is renamed twice in one call.
    #[error("wood node {0:?} is given more than one new name")]
    RenamedTwice(NodeIndex),
    /// Two spinneys would share the name.
    #[error("spinney name {0:?} is given to more than one spinney")]
    Duplicate(String),
}

/// `a`, …, `z`, then `a1`, …: names that stay short for every wood.
fn spinney_name(index: usize) -> String {
    let letter = char::from(b'a' + (index % 26) as u8);
    match index / 26 {
        0 => letter.to_string(),
        round => format!("{letter}{round}"),
    }
}

impl Expansion {
    /// The terms, in `k` order. See [`Expansion`] for an example.
    pub fn terms(&self) -> &[ForestTerm] {
        &self.terms
    }

    /// The term of the forest at `forest`, if it is one of these forests.
    ///
    /// See [`Expansion`] for an example of the terms.
    pub fn term(&self, forest: NodeIndex) -> Option<&ForestTerm> {
        self.terms.iter().find(|term| term.forest == forest)
    }

    /// The name of the spinney at wood node `spinney`, such as `a`.
    ///
    /// See [`Expansion`] for the naming order.
    pub fn name(&self, spinney: NodeIndex) -> Option<&str> {
        self.names.get(&spinney).map(String::as_str)
    }

    /// Gives the spinneys at the wood nodes in `names` their new names; the
    /// others keep theirs. Nothing changes on error.
    ///
    /// ```
    /// # use feynkit_graph::FeynmanDiagram;
    /// # use feynkit_model::Model;
    /// # use feynkit_uv::{RenameError, RenormalizationPrescriptionSettings, UvGraph, Wood};
    /// # let sunrise = FeynmanDiagram::from_dot(
    /// #     Model::phi3(),
    /// #     "digraph sunrise { edge [particle=phi]; a -> b; a -> b; a -> b; }",
    /// # )?;
    /// # let spinneys = sunrise.classify_spinneys(
    /// #     &sunrise.hedge_graph().full_filter(),
    /// #     sunrise.model(),
    /// #     &RenormalizationPrescriptionSettings::default(),
    /// #     sunrise.loop_momentum_basis(),
    /// # )?;
    /// # let forests = Wood::from_spinneys(spinneys, &sunrise, sunrise.loop_momentum_basis())?.unfold();
    /// let mut expansion = forests.expansion(&sunrise);
    /// let a = expansion.terms()[2].factors[0].spinney;
    /// expansion.rename([(a, "x".to_owned())])?;
    /// assert_eq!(expansion.k_table("K")[5], "k5 = K_{Γ|x}[k2]");
    /// let b = expansion.terms()[3].factors[0].spinney;
    /// assert_eq!(
    ///     expansion.rename([(b, "x".to_owned())]),
    ///     Err(RenameError::Duplicate("x".to_owned()))
    /// );
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn rename(
        &mut self,
        names: impl IntoIterator<Item = (NodeIndex, String)>,
    ) -> Result<(), RenameError> {
        let mut renamed = self.names.clone();
        let mut touched = std::collections::BTreeSet::new();
        for (node, name) in names {
            let slot = renamed
                .get_mut(&node)
                .ok_or(RenameError::UnknownSpinney(node))?;
            if !touched.insert(node) {
                return Err(RenameError::RenamedTwice(node));
            }
            let stored_quantity = name.strip_prefix('k').is_some_and(|digits| {
                !digits.is_empty() && digits.chars().all(|c| c.is_ascii_digit())
            });
            if name.is_empty()
                || stored_quantity
                || name
                    .chars()
                    .any(|c| c.is_whitespace() || c.is_control() || "|{}[],·∅".contains(c))
            {
                return Err(RenameError::InvalidName(name));
            }
            if name == WHOLE_GRAPH && Some(node) != self.whole {
                return Err(RenameError::ReservedName);
            }
            *slot = name;
        }
        let mut seen = std::collections::BTreeSet::new();
        for name in renamed.values() {
            if !seen.insert(name) {
                return Err(RenameError::Duplicate(name.clone()));
            }
        }
        self.names = renamed;
        Ok(())
    }

    /// The factor `Op_{S|G}` in text, such as `K_{Γ|a}`, or `K_{a}` when nothing
    /// is nested in `S`.
    ///
    /// See [`Expansion`] for an example.
    pub fn factor_notation(&self, factor: &TermFactor, op: &str) -> String {
        let given: Vec<&str> = factor
            .given
            .iter()
            .filter_map(|node| self.name(*node))
            .collect();
        // Single-character names read run together; longer ones need commas.
        let given = if given.iter().all(|name| name.chars().count() == 1) {
            given.concat()
        } else {
            given.join(",")
        };
        let name = self.name(factor.spinney).unwrap_or("?");
        if given.is_empty() {
            format!("{op}_{{{name}}}")
        } else {
            format!("{op}_{{{name}|{given}}}")
        }
    }

    /// A term in text: its factors side by side, such as `K_{Γ|a}`, or
    /// `Op[∅]` for the empty forest.
    ///
    /// See [`Expansion`] for an example.
    pub fn notation(&self, term: &ForestTerm, op: &str) -> String {
        if term.factors.is_empty() {
            return format!("{op}[∅]");
        }
        term.factors
            .iter()
            .map(|factor| self.factor_notation(factor, op))
            .collect()
    }

    /// The stored quantities, one line per term: `k5 = K_{Γ|a}[k2]`, where the
    /// bracket names the nested term, and `k7 = k1·k2` for a product of
    /// disjoint factors.
    ///
    /// See [`Expansion`] for an example.
    pub fn k_table(&self, op: &str) -> Vec<String> {
        // Each factor stands for the forest made of its own subtree.
        let k_of: BTreeMap<Vec<NodeIndex>, usize> = self
            .terms
            .iter()
            .map(|term| (subtree(&term.factors), term.k))
            .collect();
        let k_name = |factors: &[TermFactor]| {
            k_of.get(&subtree(factors))
                .map_or_else(|| "?".to_owned(), |k| format!("k{k}"))
        };
        self.terms
            .iter()
            .map(|term| {
                let value = match term.factors.as_slice() {
                    [] => format!("{op}[∅]"),
                    [factor] if factor.inner.is_empty() => self.factor_notation(factor, op),
                    [factor] => format!(
                        "{}[{}]",
                        self.factor_notation(factor, op),
                        k_name(&factor.inner)
                    ),
                    factors => factors
                        .iter()
                        .map(|factor| k_name(slice::from_ref(factor)))
                        .collect::<Vec<_>>()
                        .join("·"),
                };
                format!("k{} = {value}", term.k)
            })
            .collect()
    }
}

/// Every spinney of `factors` and their nested factors, sorted: the elements of
/// the forest they form.
fn subtree(factors: &[TermFactor]) -> Vec<NodeIndex> {
    fn collect(factors: &[TermFactor], into: &mut Vec<NodeIndex>) {
        for factor in factors {
            into.push(factor.spinney);
            collect(&factor.inner, into);
        }
    }
    let mut nodes = Vec::new();
    collect(factors, &mut nodes);
    nodes.sort();
    nodes
}
