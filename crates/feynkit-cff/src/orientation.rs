use linnet::half_edge::involution::{EdgeIndex, EdgeVec, Orientation};
use serde::{Deserialize, Serialize};
use symbolica::atom::{Atom, AtomCore};

pub type EdgeOrientation = Orientation;
pub type GraphOrientation = EdgeVec<Orientation>;

/// Structural operations on a complete Linnet graph orientation.
pub trait GraphOrientationExt {
    /// Compare orientations while treating undirected reference edges as wildcards.
    fn is_compatible_with(&self, reference: &Self) -> bool;
}

impl GraphOrientationExt for GraphOrientation {
    fn is_compatible_with(&self, reference: &Self) -> bool {
        orientation_entries_are_compatible(
            self.iter().map(|(edge, orientation)| (edge, *orientation)),
            reference
                .iter()
                .map(|(edge, orientation)| (edge, *orientation)),
        )
    }
}

fn orientation_entries_are_compatible(
    left: impl ExactSizeIterator<Item = (EdgeIndex, Orientation)>,
    reference: impl ExactSizeIterator<Item = (EdgeIndex, Orientation)>,
) -> bool {
    left.len() == reference.len()
        && left
            .zip(reference)
            .all(|((left_id, left), (right_id, right))| {
                left_id == right_id && (right == Orientation::Undirected || left == right)
            })
}

/// Stable index of an orientation within a generated expression.
#[derive(
    Clone,
    Copy,
    Debug,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Hash,
    Serialize,
    Deserialize,
    bincode::Encode,
    bincode::Decode,
)]
pub struct OrientationId(pub usize);

impl OrientationId {
    pub const fn index(self) -> usize {
        self.0
    }
}

impl From<usize> for OrientationId {
    fn from(value: usize) -> Self {
        Self(value)
    }
}

impl From<OrientationId> for usize {
    fn from(value: OrientationId) -> Self {
        value.0
    }
}

/// Exact weight `2^-halvings` multiplying one orientation's denominator tree.
///
/// Each of the two closures of a self-loop carries half of the tadpole
/// contour, as in the shared recursion of `three-dimensional-reps`, so an
/// orientation directing `k` self-loops has weight `2^-k`. Human-readable
/// expressions without a weight field read as the unit weight. Compact binary
/// encodings always carry the weight, so their layout differs from encodings
/// made before this weight existed.
#[derive(
    Clone,
    Copy,
    Debug,
    Default,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Hash,
    Serialize,
    Deserialize,
    bincode::Encode,
    bincode::Decode,
)]
pub struct OrientationWeight {
    halvings: u32,
}

impl OrientationWeight {
    pub const UNIT: Self = Self { halvings: 0 };

    pub const fn from_halvings(halvings: u32) -> Self {
        Self { halvings }
    }

    /// Number of factors of one half in this weight.
    pub const fn halvings(self) -> u32 {
        self.halvings
    }

    pub const fn is_unit(&self) -> bool {
        self.halvings == 0
    }

    /// The exact rational weight as a Symbolica number.
    pub fn to_atom(self) -> Atom {
        Atom::num(2).pow(-i64::from(self.halvings))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn undirected_reference_edges_are_wildcards() {
        let orientation =
            GraphOrientation::from_iter([Orientation::Default, Orientation::Reversed]);
        let reference =
            GraphOrientation::from_iter([Orientation::Default, Orientation::Undirected]);

        assert!(orientation.is_compatible_with(&reference));
        assert!(!reference.is_compatible_with(&orientation));
    }

    #[test]
    fn unequal_edge_coverage_is_incompatible_in_both_directions() {
        let short = GraphOrientation::from_iter([Orientation::Default]);
        let long = GraphOrientation::from_iter([Orientation::Default, Orientation::Undirected]);

        assert!(!short.is_compatible_with(&long));
        assert!(!long.is_compatible_with(&short));
    }

    #[test]
    fn different_edge_ids_are_incompatible() {
        let left = [(EdgeIndex(0), Orientation::Default)];
        let reference = [(EdgeIndex(1), Orientation::Default)];

        assert!(!orientation_entries_are_compatible(
            left.into_iter(),
            reference.into_iter()
        ));
    }
}
