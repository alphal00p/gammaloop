use color_eyre::Result;
use eyre::eyre;
use itertools::{EitherOrBoth, Itertools};
use linnet::half_edge::involution::{EdgeIndex, EdgeVec, EdgeVecIter, Orientation};
use symbolica::prelude::*;

use crate::utils::GS;

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct RepresentativeScore {
    undirected_count: usize,
    default_count: usize,
    leading_defaults: Vec<bool>,
}

pub trait GraphOrientation: Sized {
    fn orientation(&self) -> &EdgeVec<Orientation>;

    fn orientation_thetas(&self) -> Atom {
        let mut thetas = Atom::num(1);

        for (e, h) in self.orientation() {
            match h {
                Orientation::Default => {
                    thetas *= GS.sign_theta(GS.sign(e));
                }
                Orientation::Reversed => {
                    thetas *= GS.sign_theta(-GS.sign(e));
                }
                _ => {}
            }
        }
        thetas
    }

    fn orientation_delta(&self) -> Atom
    where
        Self: Sized,
    {
        GS.orientation_delta(self)
    }

    fn select<'a>(&self, atom: impl Into<AtomOrView<'a>>) -> Atom {
        let theta_reps = vec![
            Replacement::new(GS.sign_theta(1).to_pattern(), Atom::num(1)),
            Replacement::new(GS.sign_theta(-1).to_pattern(), Atom::Zero),
        ];

        let mut reps = Vec::new();

        for (e, h) in self.orientation() {
            match h {
                Orientation::Default => {
                    reps.push(Replacement::new(GS.sign(e).to_pattern(), Atom::num(1)));
                }
                Orientation::Reversed => {
                    reps.push(Replacement::new(GS.sign(e).to_pattern(), Atom::num(-1)));
                }
                _ => {}
            }
        }

        let orientation = self.orientation();
        atom.into()
            .replace_multiple(&reps)
            .replace_multiple(&theta_reps)
            .replace_map(|term, _ctx, out| {
                if let AtomView::Fun(f) = term
                    && f.get_symbol() == GS.orientation_delta
                {
                    if f.iter()
                        .zip_longest(orientation)
                        .all(|either| match either {
                            EitherOrBoth::Both(a, (_, o)) => {
                                if let Ok(a) = i64::try_from(a) {
                                    match o {
                                        Orientation::Default => a >= 0,
                                        Orientation::Reversed => a <= 0,
                                        Orientation::Undirected => true,
                                    }
                                } else {
                                    false
                                }
                            }
                            EitherOrBoth::Left(_) | EitherOrBoth::Right(_) => false,
                        })
                    {
                        **out = Atom::num(1);
                    } else {
                        **out = Atom::Zero;
                    }
                }
            })
    }

    fn iterate<'a>(&'a self) -> EdgeVecIter<'a, Orientation> {
        self.orientation().iter()
    }

    fn get(&self, index: EdgeIndex) -> Orientation {
        self.orientation()[index]
    }

    /// Returns `true` if the orientation matches the reference, ignoring `Undirected` edges of the reference.
    fn is_compatible_with(&self, reference: &Self) -> bool {
        self.iterate()
            .zip_eq(reference.orientation())
            .all(|((_, l), (_, r))| matches!(r, Orientation::Undirected) || l == r)
    }

    fn score(&self, internal_edges: &[EdgeIndex]) -> RepresentativeScore {
        let mut undirected_count = 0;
        let mut default_count = 0;
        let mut leading_defaults = Vec::with_capacity(internal_edges.len());

        for edge in internal_edges {
            let is_default = match self.get(*edge) {
                Orientation::Undirected => {
                    undirected_count += 1;
                    false
                }
                Orientation::Default => {
                    default_count += 1;
                    true
                }
                Orientation::Reversed => false,
            };
            leading_defaults.push(is_default);
        }

        RepresentativeScore {
            undirected_count,
            default_count,
            leading_defaults,
        }
    }

    /// Select the deterministic full-graph representative used to host the integrated UV term.
    ///
    /// The reduced orientation already fixes all edges that remain explicit after contracting the
    /// integrated subgraph. Any contracted edge appears as `Undirected` and is therefore ignored when
    /// matching compatible full orientations.
    ///
    /// Among the compatible candidates, the representative is chosen by maximizing:
    /// 1. the number of `Undirected` internal edges
    /// 2. then the number of `Default` internal edges
    /// 3. then the lexicographically maximal `is_default` pattern in ascending internal-edge order
    ///
    /// The first criterion is currently future-proof only: the canonical full-graph acyclic basis is
    /// built from `Default` / `Reversed` assignments on paired internal edges, so `Undirected`
    /// internal entries are not expected there unless the global orientation-generation step changes.
    fn select_representative_orientation<'b>(
        self,
        valid_global_orientations: &'b [Self],
        internal_edges: &[EdgeIndex],
    ) -> Result<&'b Self> {
        valid_global_orientations
            .iter()
            .filter(|candidate| candidate.is_compatible_with(&self))
            .max_by_key(|candidate| candidate.score(internal_edges))
            .ok_or_else(|| {
                eyre!(
                    "no valid global orientation matches reduced orientation {}",
                    self.orientation_delta()
                )
            })
    }

    /// Build only the selector factors for the internal edges of the integrated subgraph.
    ///
    /// This intentionally does **not** use `orientation_thetas()` on the whole representative
    /// orientation. The reduced-graph CFF term already carries selectors for the edges that remain
    /// explicit after contraction, so re-applying them here would duplicate selectors on external
    /// edges.
    fn internal_orientation_selector(&self, internal_edges: &[EdgeIndex]) -> Atom {
        let mut selector = Atom::num(1);

        for edge in internal_edges {
            match self.get(*edge) {
                Orientation::Default => selector *= GS.sign_theta(GS.sign(*edge)),
                Orientation::Reversed => selector *= GS.sign_theta(-GS.sign(*edge)),
                Orientation::Undirected => {}
            }
        }

        selector
    }
}

#[cfg(test)]
mod tests {

    use crate::{cff::orientations::GraphOrientation, utils::GS};
    use linnet::half_edge::involution::{EdgeIndex, EdgeVec, Orientation};

    fn orientation(value: i8) -> Orientation {
        match value {
            1 => Orientation::Default,
            -1 => Orientation::Reversed,
            0 => Orientation::Undirected,
            _ => panic!("invalid orientation encoding"),
        }
    }

    fn edgevec(values: impl IntoIterator<Item = i8>) -> EdgeVec<Orientation> {
        EdgeVec::from_iter(values.into_iter().map(orientation))
    }

    fn edges(values: impl IntoIterator<Item = usize>) -> Vec<EdgeIndex> {
        values.into_iter().map(EdgeIndex).collect()
    }

    #[test]
    fn picks_first_valid_internal_representatives_per_external_class() {
        let valid = vec![
            edgevec([1, -1, 1, -1, 1]),
            edgevec([1, 1, -1, -1, 1]),
            edgevec([1, 1, -1, -1, -1]),
            edgevec([-1, 1, 1, 1, -1]),
        ];
        let internal = edges([2, 3, 4]);

        let reduced_pm = edgevec([1, -1, 0, 0, 0]);
        let reduced_pp = edgevec([1, 1, 0, 0, 0]);
        let reduced_mp = edgevec([-1, 1, 0, 0, 0]);

        assert_eq!(
            reduced_pm
                .select_representative_orientation(&valid, &internal)
                .unwrap(),
            &valid[0]
        );
        assert_eq!(
            reduced_pp
                .select_representative_orientation(&valid, &internal)
                .unwrap(),
            &valid[1]
        );
        assert_eq!(
            reduced_mp
                .select_representative_orientation(&valid, &internal)
                .unwrap(),
            &valid[3]
        );
    }

    #[test]
    fn prefers_more_undirected_internal_edges_first() {
        let valid = vec![edgevec([1, 0, 1]), edgevec([1, 1, 1]), edgevec([1, -1, 1])];
        let reduced = edgevec([1, 0, 1]);
        let internal = edges([1]);

        assert_eq!(
            reduced
                .select_representative_orientation(&valid, &internal)
                .unwrap(),
            &valid[0]
        );
    }

    #[test]
    fn prefers_more_default_internal_edges_when_undirected_counts_tie() {
        let valid = vec![
            edgevec([1, 1, -1]),
            edgevec([1, -1, -1]),
            edgevec([1, -1, 1]),
        ];
        let reduced = edgevec([1, 0, 0]);
        let internal = edges([1, 2]);

        assert_eq!(
            reduced
                .select_representative_orientation(&valid, &internal)
                .unwrap(),
            &valid[0]
        );
    }

    #[test]
    fn breaks_default_ties_by_lexicographically_maximal_leading_defaults() {
        let valid = vec![
            edgevec([1, 1, -1, 1]),
            edgevec([1, -1, 1, 1]),
            edgevec([1, 1, 1, -1]),
        ];
        let reduced = edgevec([1, 0, 0, 0]);
        let internal = edges([1, 2, 3]);

        assert_eq!(
            reduced
                .select_representative_orientation(&valid, &internal)
                .unwrap(),
            &valid[2]
        );
    }

    #[test]
    fn errors_for_missing_external_orientation_class() {
        let valid = vec![edgevec([1, 1, 1]), edgevec([1, -1, 1])];
        let reduced = edgevec([-1, 0, 1]);
        let internal = edges([1]);

        assert!(
            reduced
                .select_representative_orientation(&valid, &internal)
                .is_err()
        );
    }

    #[test]
    fn selector_only_depends_on_internal_edges() {
        let representative = edgevec([1, -1, 1, -1]);
        let selector = representative.internal_orientation_selector(&edges([1, 3]));

        let expected =
            GS.sign_theta(-GS.sign(EdgeIndex(1))) * GS.sign_theta(-GS.sign(EdgeIndex(3)));
        assert_eq!(selector, expected);
    }
}
