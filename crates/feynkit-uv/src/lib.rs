//! Ultraviolet (UV) subtraction combinatorics for Feynman graphs.
//!
//! The UV-divergent subgraphs of a graph (spinneys) are enumerated as unions of
//! cycles and classified by power counting and a renormalization prescription.
//! Ordered by inclusion they form a [`Wood`]. Unfolding the commuting
//! subtraction operations of the wood yields its [`Forests`]: one node per
//! forest of the forest formula, in Foata normal form.
//!
//! The combinatorics only needs the queries of [`UvGraph`]. FeynKit's
//! [`feynkit_graph::FeynmanDiagram`] implements it; applications implement it
//! for their own graphs and keep any per-forest computed values themselves,
//! keyed by forest node or [`OperationNode`].
//!
//! ```
//! use feynkit_graph::FeynmanDiagram;
//! use feynkit_model::Model;
//! use feynkit_uv::{RenormalizationPrescriptionSettings, UvGraph, Wood};
//!
//! let diagram = FeynmanDiagram::from_dot(
//!     Model::phi4(),
//!     "digraph sunrise { edge [particle=phi]; ext [style=invis]; \
//!          ext -> a; a -> b; a -> b; a -> b; b -> ext; }",
//! )?;
//! let spinneys = diagram.classify_spinneys(
//!     &diagram.hedge_graph().full_filter(),
//!     diagram.model(),
//!     &RenormalizationPrescriptionSettings::default(),
//!     diagram.loop_momentum_basis(),
//! )?;
//! let forests = Wood::from_spinneys(spinneys, &diagram, diagram.loop_momentum_basis())?.unfold();
//! for (_, forest) in forests.operations() {
//!     println!("{forest}: depth {}, sign {}", forest.depth(), forest.sign());
//! }
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

#![forbid(unsafe_code)]

mod counterterm;
mod diagram;
#[cfg(feature = "svg")]
mod display;
mod expansion;
mod forests;
mod graph;
mod scheme;
mod spinney;
mod wood;

pub use counterterm::{Counterterm, CountertermError, CountertermScheme, UnknownCountertermScheme};
#[cfg(feature = "svg")]
pub use display::{FULL_FILL, ROOT_FILL, STAMP_SIZE, StampOptions, UvDisplayError, UvRegion};
pub use expansion::{Expansion, ForestTerm, RenameError, TermFactor, WHOLE_GRAPH};
pub use forests::{
    ForestEdge, ForestError, Forests, Operation, OperationNode, OperationSpinneyError,
    UnionReplayState,
};
pub use graph::{ForeignReason, UvGraph, UvGraphError};
pub use scheme::{
    ApproximationType, CTIdentifier, CTRenormalizationRule, RenormalizationPrescriptionSettings,
};
pub use spinney::{Spinney, SpinneyRelation};
pub use wood::{ForestNode, OwnedForestNode, Wood, WoodInclusion};
