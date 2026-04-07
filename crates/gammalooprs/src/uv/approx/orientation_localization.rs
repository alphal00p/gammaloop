use crate::{cff::expression::GraphOrientation, utils::GS};
use color_eyre::Result;
use eyre::eyre;
use linnet::half_edge::involution::{EdgeIndex, EdgeVec, Orientation};
use symbolica::atom::Atom;

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct RepresentativeScore {
    undirected_count: usize,
    default_count: usize,
    leading_defaults: Vec<bool>,
}

fn representative_score(
    orientation: &EdgeVec<Orientation>,
    internal_edges: &[EdgeIndex],
) -> RepresentativeScore {
    let mut undirected_count = 0;
    let mut default_count = 0;
    let mut leading_defaults = Vec::with_capacity(internal_edges.len());

    for edge in internal_edges {
        let is_default = match orientation[*edge] {
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

fn orientations_match_outside_integrated_subgraph(
    reduced_orientation: &EdgeVec<Orientation>,
    global_orientation: &EdgeVec<Orientation>,
) -> bool {
    reduced_orientation.iter().all(|(edge, reduced)| {
        matches!(reduced, Orientation::Undirected) || *reduced == global_orientation[edge]
    })
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
pub(crate) fn select_representative_orientation<'a>(
    reduced_orientation: &EdgeVec<Orientation>,
    valid_global_orientations: &'a [EdgeVec<Orientation>],
    internal_edges: &[EdgeIndex],
) -> Result<&'a EdgeVec<Orientation>> {
    valid_global_orientations
        .iter()
        .filter(|candidate| {
            orientations_match_outside_integrated_subgraph(reduced_orientation, candidate)
        })
        .max_by_key(|candidate| representative_score(candidate, internal_edges))
        .ok_or_else(|| {
            eyre!(
                "no valid global orientation matches reduced orientation {}",
                GS.orientation_delta(reduced_orientation)
            )
        })
}

/// Build only the selector factors for the internal edges of the integrated subgraph.
///
/// This intentionally does **not** use `orientation_thetas()` on the whole representative
/// orientation. The reduced-graph CFF term already carries selectors for the edges that remain
/// explicit after contraction, so re-applying them here would duplicate selectors on external
/// edges.
pub(crate) fn internal_orientation_selector(
    representative: &EdgeVec<Orientation>,
    internal_edges: &[EdgeIndex],
) -> Atom {
    let mut selector = Atom::num(1);

    for edge in internal_edges {
        match representative[*edge] {
            Orientation::Default => selector *= GS.sign_theta(GS.sign(*edge)),
            Orientation::Reversed => selector *= GS.sign_theta(-GS.sign(*edge)),
            Orientation::Undirected => {}
        }
    }

    selector
}

/// Localize a reduced-graph orientation term onto one representative full orientation.
///
/// The reduced term already contains the external-edge selector structure through
/// `reduced_orientation.orientation_thetas()`. We only add the selector for the integrated
/// subgraph's internal edges chosen by `select_representative_orientation`.
pub(crate) fn localize_reduced_orientation_term(
    reduced_expression: &Atom,
    reduced_orientation: &EdgeVec<Orientation>,
    valid_global_orientations: &[EdgeVec<Orientation>],
    internal_edges: &[EdgeIndex],
) -> Result<Atom> {
    let representative = select_representative_orientation(
        reduced_orientation,
        valid_global_orientations,
        internal_edges,
    )?;

    Ok(reduced_expression.clone()
        * reduced_orientation.orientation_thetas()
        * internal_orientation_selector(representative, internal_edges))
}

#[cfg(test)]
mod tests {
    use super::{
        internal_orientation_selector, localize_reduced_orientation_term,
        select_representative_orientation,
    };
    use crate::utils::GS;
    use linnet::half_edge::involution::{EdgeIndex, EdgeVec, Orientation};
    use symbolica::function;

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
            select_representative_orientation(&reduced_pm, &valid, &internal).unwrap(),
            &valid[0]
        );
        assert_eq!(
            select_representative_orientation(&reduced_pp, &valid, &internal).unwrap(),
            &valid[1]
        );
        assert_eq!(
            select_representative_orientation(&reduced_mp, &valid, &internal).unwrap(),
            &valid[3]
        );
    }

    #[test]
    fn prefers_more_undirected_internal_edges_first() {
        let valid = vec![edgevec([1, 0, 1]), edgevec([1, 1, 1]), edgevec([1, -1, 1])];
        let reduced = edgevec([1, 0, 1]);
        let internal = edges([1]);

        assert_eq!(
            select_representative_orientation(&reduced, &valid, &internal).unwrap(),
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
            select_representative_orientation(&reduced, &valid, &internal).unwrap(),
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
            select_representative_orientation(&reduced, &valid, &internal).unwrap(),
            &valid[2]
        );
    }

    #[test]
    fn errors_for_missing_external_orientation_class() {
        let valid = vec![edgevec([1, 1, 1]), edgevec([1, -1, 1])];
        let reduced = edgevec([-1, 0, 1]);
        let internal = edges([1]);

        assert!(select_representative_orientation(&reduced, &valid, &internal).is_err());
    }

    #[test]
    fn selector_only_depends_on_internal_edges() {
        let representative = edgevec([1, -1, 1, -1]);
        let selector = internal_orientation_selector(&representative, &edges([1, 3]));

        let expected =
            GS.sign_theta(-GS.sign(EdgeIndex(1))) * GS.sign_theta(-GS.sign(EdgeIndex(3)));
        assert_eq!(selector, expected);
    }

    #[test]
    fn localization_keeps_external_selectors_and_adds_internal_ones() {
        let reduced_expression = function!(GS.ose, 0);
        let reduced_orientation = edgevec([1, 0, -1]);
        let valid = vec![edgevec([1, 1, -1]), edgevec([1, -1, -1])];
        let localized = localize_reduced_orientation_term(
            &reduced_expression,
            &reduced_orientation,
            &valid,
            &edges([1]),
        )
        .unwrap();

        let expected = reduced_expression
            * GS.sign_theta(GS.sign(EdgeIndex(0)))
            * GS.sign_theta(-GS.sign(EdgeIndex(2)))
            * GS.sign_theta(GS.sign(EdgeIndex(1)));
        assert_eq!(localized, expected);
    }
}
