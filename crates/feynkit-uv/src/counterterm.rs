//! Counterterms as combinatorial objects: a graph with a forest of spinneys
//! contracted, and the scheme that fixes their finite parts.
//!
//! Counterterm expressions and integrands are not built here yet; this type
//! only identifies which spinneys are contracted and in which scheme.

use std::fmt::{self, Display, Formatter};
use std::str::FromStr;

use linnet::half_edge::NodeIndex;
use linnet::half_edge::subgraph::{InternalSubGraph, SubSetLike};
use thiserror::Error;

use crate::{ApproximationType, Forests, Spinney, SpinneyRelation, UvGraph, UvGraphError};

/// The renormalization scheme of a counterterm.
///
/// ```
/// use feynkit_uv::CountertermScheme;
///
/// let scheme: CountertermScheme = "ms-BAR".parse()?;
/// assert_eq!(scheme, CountertermScheme::MsBar);
/// assert_eq!(CountertermScheme::OnShell.to_string(), "on-shell");
/// assert!("minimal".parse::<CountertermScheme>().is_err());
/// # Ok::<(), feynkit_uv::UnknownCountertermScheme>(())
/// ```
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CountertermScheme {
    /// No scheme chosen yet.
    #[default]
    Unspecified,
    /// Modified minimal subtraction.
    MsBar,
    /// On-shell subtraction.
    OnShell,
}

/// A scheme name other than `MS-bar`, `on-shell` or `unspecified`.
///
/// See [`CountertermScheme`] for an example.
#[derive(Clone, Copy, Debug, Error, PartialEq, Eq)]
#[error("counterterm scheme must be 'MS-bar', 'on-shell', or 'unspecified'")]
pub struct UnknownCountertermScheme;

impl CountertermScheme {
    /// The scheme the approximation of `spinneys` settles: on shell when every
    /// one is subtracted on shell ([`ApproximationType::OS`]), and unspecified
    /// otherwise. `MUV` does not imply MS-bar.
    ///
    /// ```
    /// use feynkit_graph::FeynmanDiagram;
    /// use feynkit_model::Model;
    /// use feynkit_uv::{CountertermScheme, RenormalizationPrescriptionSettings, UvGraph};
    ///
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
    /// // The default prescription subtracts with MUV everywhere.
    /// assert_eq!(CountertermScheme::implied_by(&spinneys), CountertermScheme::Unspecified);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn implied_by<'a, B: 'a>(spinneys: impl IntoIterator<Item = &'a Spinney<B>>) -> Self {
        let mut approximations = spinneys
            .into_iter()
            .filter(|spinney| !spinney.is_empty())
            .map(Spinney::renormalization_scheme)
            .peekable();
        if approximations.peek().is_some()
            && approximations.all(|approximation| approximation == ApproximationType::OS)
        {
            Self::OnShell
        } else {
            Self::Unspecified
        }
    }
}

impl Display for CountertermScheme {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::Unspecified => "unspecified",
            Self::MsBar => "MS-bar",
            Self::OnShell => "on-shell",
        })
    }
}

impl FromStr for CountertermScheme {
    type Err = UnknownCountertermScheme;

    /// Case-insensitive `MS-bar`, `on-shell` or `unspecified`.
    fn from_str(name: &str) -> Result<Self, Self::Err> {
        [Self::Unspecified, Self::MsBar, Self::OnShell]
            .into_iter()
            .find(|scheme| scheme.to_string().eq_ignore_ascii_case(name))
            .ok_or(UnknownCountertermScheme)
    }
}

/// A selection of spinneys that cannot be contracted together; `E` is the
/// graph's own error type.
#[derive(Debug, Error)]
pub enum CountertermError<E> {
    /// No spinney was selected, or only the empty one.
    #[error("a counterterm contracts at least one nonempty spinney")]
    Empty,
    /// Two selected spinneys overlap or touch, so they belong to no common forest.
    #[error("spinneys {first} and {second} are {relation} and cannot be contracted together")]
    Incompatible {
        first: String,
        second: String,
        relation: SpinneyRelation,
    },
    /// A spinney has a convergent connected factor, so classification would
    /// not retain it.
    #[error("spinney {spinney} has a convergent connected component")]
    ConvergentFactor { spinney: String },
    /// A spinney is selected twice, possibly as a connected component of a
    /// selected disconnected spinney, so the selection is no forest.
    #[error("spinney {spinney} is selected more than once")]
    Duplicate { spinney: String },
    /// A spinney does not belong to the destination graph.
    #[error(transparent)]
    Foreign(#[from] UvGraphError<E>),
}

/// A graph with a forest of spinneys contracted, in a renormalization scheme.
///
/// The contracted spinneys are connected, pairwise nested or disjoint, and
/// outermost first: by decreasing loop count, then by subgraph. A disconnected
/// spinney contracts as its connected components, as its product forest does.
///
/// ```
/// use feynkit_graph::{DiagramError, FeynmanDiagram};
/// use feynkit_model::Model;
/// use feynkit_uv::{Counterterm, CountertermScheme, RenormalizationPrescriptionSettings, UvGraph, Wood};
///
/// let eye = FeynmanDiagram::from_dot(
///     Model::phi4(),
///     "digraph eye { edge [particle=phi]; ext [style=invis]; ext -> a; ext -> a; \
///          a -> b; a -> c; b -> c; b -> c; b -> ext; c -> ext; }",
/// )?;
/// let spinneys = eye.classify_spinneys(
///     &eye.hedge_graph().full_filter(),
///     eye.model(),
///     &RenormalizationPrescriptionSettings::default(),
///     eye.loop_momentum_basis(),
/// )?;
/// let forests = Wood::from_spinneys(spinneys, &eye, eye.loop_momentum_basis())?.unfold();
/// let (node, _) = forests
///     .operations()
///     .find(|(_, forest)| forest.depth() == 2)
///     .expect("the bubble nests in the full graph");
/// let counterterm = Counterterm::of_forest::<DiagramError>(&forests, node, CountertermScheme::MsBar)?;
/// let loops: Vec<_> = counterterm.contracted().iter().map(|s| s.loop_count()).collect();
/// assert_eq!(loops, [2, 1]);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Clone, Debug)]
pub struct Counterterm<B> {
    contracted: Vec<Spinney<B>>,
    scheme: CountertermScheme,
}

impl<B: Clone> Counterterm<B> {
    /// Contract `spinneys` of `graph` in `scheme`.
    ///
    /// Each spinney is checked against `graph` and `parent_basis` as
    /// [`crate::Wood::from_spinneys`] checks it. A disconnected spinney is
    /// replaced by its connected components, so the contraction matches
    /// [`Self::of_forest`] of its product forest.
    ///
    /// # Errors
    ///
    /// Fails when no nonempty spinney is given, when a spinney does not belong to
    /// `graph` or has a convergent connected factor, when two of them overlap or
    /// touch, or when a spinney is selected twice, also as a component of a
    /// disconnected one.
    ///
    /// See [`Counterterm`] for an example of a forest counterterm.
    pub fn new<G: UvGraph<Basis = B> + ?Sized>(
        spinneys: impl IntoIterator<Item = Spinney<B>>,
        scheme: CountertermScheme,
        graph: &G,
        parent_basis: &B,
    ) -> Result<Self, CountertermError<G::Error>>
    where
        B: PartialEq,
    {
        let mut factors: Vec<Spinney<B>> = Vec::new();
        for spinney in spinneys {
            spinney.ensure_belongs_to(graph, parent_basis)?;
            if spinney.is_empty() {
                continue;
            }
            if !spinney.has_divergent_factors(graph)? {
                return Err(CountertermError::ConvergentFactor {
                    spinney: spinney.label(graph),
                });
            }
            if spinney.n_components() == 1 {
                factors.push(spinney);
                continue;
            }
            for component in spinney.components() {
                let subgraph = InternalSubGraph::cleaned_filter_optimist(
                    component.clone(),
                    graph.hedge_graph(),
                );
                let dod = graph
                    .degree_of_divergence(&subgraph.filter)
                    .map_err(|source| UvGraphError::DegreeOfDivergence {
                        subgraph: subgraph.string_label(),
                        source,
                    })?;
                let factor = Spinney::with_scheme(
                    subgraph,
                    graph,
                    parent_basis,
                    spinney.renormalization_scheme(),
                    dod,
                )?
                .ok_or_else(|| CountertermError::ConvergentFactor {
                    spinney: spinney.label(graph),
                })?;
                factors.push(factor);
            }
        }
        let mut contracted: Vec<Spinney<B>> = Vec::new();
        for spinney in factors {
            if contracted
                .iter()
                .any(|other| other.filter() == spinney.filter())
            {
                return Err(CountertermError::Duplicate {
                    spinney: spinney.label(graph),
                });
            }
            for other in &contracted {
                let relation = spinney.relation(other, graph);
                if !relation.is_compatible() {
                    return Err(CountertermError::Incompatible {
                        first: other.label(graph),
                        second: spinney.label(graph),
                        relation,
                    });
                }
            }
            contracted.push(spinney);
        }
        Self::outermost_first(contracted, scheme)
    }

    /// Contract the elements of the forest at `node` of `forests`, which belong
    /// to the forests' graph and are pairwise compatible by construction. A
    /// product forest contracts each factor, as [`Self::new`] does for the same
    /// spinneys.
    ///
    /// # Errors
    ///
    /// Fails with [`CountertermError::Empty`] for the empty forest.
    ///
    /// # Panics
    ///
    /// Panics if `node` is not a forest of `forests`.
    pub fn of_forest<E>(
        forests: &Forests<B>,
        node: NodeIndex,
        scheme: CountertermScheme,
    ) -> Result<Self, CountertermError<E>> {
        Self::outermost_first(
            forests.forest_spinneys(node).into_iter().cloned().collect(),
            scheme,
        )
    }

    fn outermost_first<E>(
        mut contracted: Vec<Spinney<B>>,
        scheme: CountertermScheme,
    ) -> Result<Self, CountertermError<E>> {
        if contracted.is_empty() {
            return Err(CountertermError::Empty);
        }
        contracted.sort_by(|a, b| {
            b.loop_count()
                .cmp(&a.loop_count())
                .then_with(|| a.filter().cmp(b.filter()))
        });
        Ok(Self { contracted, scheme })
    }
}

impl<B> Counterterm<B> {
    /// The contracted spinneys, outermost first.
    ///
    /// See [`Counterterm`] for an example.
    pub fn contracted(&self) -> &[Spinney<B>] {
        &self.contracted
    }

    /// The renormalization scheme.
    ///
    /// See [`Counterterm`] for an example.
    pub fn scheme(&self) -> CountertermScheme {
        self.scheme
    }

    /// The same contraction in another scheme.
    ///
    /// See [`Counterterm`] for an example.
    #[must_use]
    pub fn with_scheme(mut self, scheme: CountertermScheme) -> Self {
        self.scheme = scheme;
        self
    }
}
