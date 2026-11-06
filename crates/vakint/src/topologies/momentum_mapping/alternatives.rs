//! Native graph automorphism proposals, with exact momentum-map replay.
mod edge_orbits;
#[cfg(test)]
mod tests;

use super::*;
use std::collections::BTreeSet;

/// Optional work-count diagnostics for scope-aware native witness selection.
/// Ordinary requests use the shared RAM/time budget instead of these counts.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MomentumMatchLimits {
    pub diagnostics_enabled: bool,
    pub max_graph_automorphisms: usize,
    pub max_orientation_candidates: usize,
    pub max_verified_witnesses: usize,
}
impl Default for MomentumMatchLimits {
    fn default() -> Self {
        Self {
            diagnostics_enabled: false,
            max_graph_automorphisms: 100_000,
            max_orientation_candidates: 1_000_000,
            max_verified_witnesses: 100_000,
        }
    }
}
impl MomentumMatchLimits {
    pub fn with_diagnostics(mut self) -> Self {
        self.diagnostics_enabled = true;
        self
    }
}

fn charge(
    count: &mut usize,
    limit: Option<usize>,
    resource: &'static str,
) -> Result<(), VakintError> {
    crate::execution::checkpoint()?;
    let requested = count.saturating_add(1);
    if let Some(limit) = limit
        && requested > limit
    {
        return Err(VakintError::MomentumMatchResourceLimit {
            resource,
            requested,
            limit,
        });
    }
    *count = requested;
    Ok(())
}

/// Select the first verified witness, in enumeration order, whose acceptance
/// key is minimal. Without any accepted key, `None` keeps the ordinary first
/// witness. Errors from acceptance are never confused with ordinary rejection.
pub(crate) fn first_minimal_acceptance<K: Ord>(
    witnesses: &[ReplacementRules],
    mut accept: impl FnMut(&ReplacementRules) -> Result<Option<K>, VakintError>,
) -> Result<Option<usize>, VakintError> {
    let mut best: Option<(K, usize)> = None;
    for (index, candidate) in witnesses.iter().enumerate() {
        let accepted = accept(candidate)?;
        if index < 8 {
            log::debug!(target: "vakint::five_loop_match", "phase=alternative_admission ordinal={} admitted={}", index + 1, accepted.is_some());
        }
        if let Some(key) = accepted
            && best.as_ref().is_none_or(|(old, _)| key < *old)
        {
            best = Some((key, index));
        }
    }
    Ok(best.map(|(_, index)| index))
}

impl MomentumInput {
    /// Power colors of the propagators (each an integer power in 1..=254), or
    /// `None` when some power is symbolic, non-positive or too large.
    fn power_colors(&self) -> Option<Vec<u8>> {
        self.propagators
            .iter()
            .map(|propagator| {
                get_integer_from_atom(propagator.power.as_view())
                    .filter(|power| (1..=254).contains(power))
                    .map(|power| power as u8)
            })
            .collect()
    }

    /// Canonical form of the power-colored edge-marked incidence graph: two
    /// inputs with positive integer powers can only be equal under a momentum
    /// map that sends each line to a line of equal power if these agree.
    pub fn power_colored_canonical_graph(
        &self,
    ) -> Result<Option<symbolica::graph::Graph<u8, Empty>>, VakintError> {
        let Some(colors) = self.power_colors() else {
            return Ok(None);
        };
        Ok(Some(
            edge_orbits::EdgeGraph::new_colored(&self.incidence, &colors)?
                .canonical
                .graph,
        ))
    }

    /// Every verified unimodular loop-momentum map from this input onto
    /// `target` that sends each propagator onto a propagator of equal power
    /// (all power-colored graph isomorphisms, each with every orientation of
    /// the loop basis, verified denominator by denominator with the Jacobian
    /// by FeynKit), as the images of this input's loop momenta in the target's
    /// loop momenta, in deterministic order. `None` if a power is not a
    /// positive integer. Used by the tests; evaluation stops at the first
    /// accepted map ([`Self::find_power_preserving_map`]).
    #[cfg(test)]
    pub fn power_preserving_maps(
        &self,
        target: &MomentumInput,
    ) -> Result<Option<Vec<Vec<Atom>>>, VakintError> {
        let mut maps = Vec::new();
        Ok(self
            .find_power_preserving_map(target, |images| {
                maps.push(images.to_vec());
                Ok(false)
            })?
            .map(|_| maps))
    }

    /// Every verified power-preserving unimodular map onto `target` (as
    /// described for the test helper `power_preserving_maps`, in the same
    /// order), offered to `accept` one at a time until it returns `true` (then
    /// `Some(true)`); the enumeration of the remaining isomorphisms is skipped.
    pub fn find_power_preserving_map(
        &self,
        target: &MomentumInput,
        mut accept: impl FnMut(&[Atom]) -> Result<bool, VakintError>,
    ) -> Result<Option<bool>, VakintError> {
        let (Some(source_colors), Some(target_colors)) =
            (self.power_colors(), target.power_colors())
        else {
            return Ok(None);
        };
        let loops = self.family.loop_momenta();
        if loops.len() != target.family.loop_momenta().len()
            || self.propagators.len() != target.propagators.len()
        {
            return Ok(Some(false));
        }
        let source_graph = edge_orbits::EdgeGraph::new_colored(&self.incidence, &source_colors)?;
        let target_graph = edge_orbits::EdgeGraph::new_colored(&target.incidence, &target_colors)?;
        if source_graph.canonical.graph != target_graph.canonical.graph {
            return Ok(Some(false));
        }
        let variables = (1..=loops.len())
            .map(|axis| Atom::var(vk_symbol!(format!("momentum_match_source_{axis}"))))
            .collect::<Vec<_>>();
        let substitute = loops
            .iter()
            .zip(&variables)
            .map(|(source, variable)| Replacement::new(source.clone(), variable.clone()))
            .collect::<Vec<_>>();
        let sources = self
            .basis
            .iter()
            .map(|&slot| {
                self.propagators[slot]
                    .momentum
                    .replace_multiple(&substitute)
            })
            .collect::<Vec<_>>();
        let orientations = 1usize
            .checked_shl(loops.len().try_into().map_err(invalid)?)
            .ok_or_else(|| invalid("momentum orientation count overflow"))?;
        let mut seen = BTreeSet::new();
        let mut work = 0;
        let accepted = edge_orbits::visit_until(
            &source_graph,
            &target_graph,
            None,
            &mut work,
            |permutation| {
                let images = self
                    .basis
                    .iter()
                    .map(|&slot| &target.propagators[permutation[slot]].momentum)
                    .collect::<Vec<_>>();
                for signs in 0..orientations {
                    crate::execution::checkpoint()?;
                    let system = sources
                        .iter()
                        .zip(&images)
                        .enumerate()
                        .map(|(axis, (source, image))| {
                            source
                                - Atom::num(if signs & (1 << axis) == 0 { 1 } else { -1 })
                                    * (*image)
                        })
                        .collect::<Vec<_>>();
                    let coordinates = Topology::solve_loop_coordinates(&system, &variables)?;
                    if !seen.insert(coordinates.clone()) {
                        continue;
                    }
                    let Some(mapping) = self
                        .family
                        .mapping_to(&target.family, &coordinates)
                        .map_err(invalid)?
                    else {
                        continue;
                    };
                    if mapping
                        .denominator_map()
                        .iter()
                        .enumerate()
                        .any(|(source, target)| source_colors[source] != target_colors[*target])
                    {
                        continue;
                    }
                    let images = mapping
                        .momentum_rules()
                        .iter()
                        .map(|(_, image)| image.clone())
                        .collect::<Vec<_>>();
                    if accept(&images)? {
                        return Ok(true);
                    }
                }
                Ok(false)
            },
        )?;
        Ok(Some(accepted))
    }

    /// Integer coefficient vectors of the propagator momenta over the loop
    /// momenta, or `None` when a coefficient is not an integer.
    fn line_vectors(&self) -> Option<Vec<Vec<i128>>> {
        let loops = self.family.loop_momenta();
        self.propagators
            .iter()
            .map(|propagator| integer_line_vector(&propagator.momentum, loops))
            .collect()
    }

    /// Power-colored profile of the cycle matroid of the lines (see
    /// [`MatroidProfile`]), or `None` when a power is not a positive integer
    /// or a momentum coefficient is not an integer.
    pub fn power_colored_matroid_profile(&self) -> Option<MatroidProfile> {
        let colors = self.power_colors()?;
        let vectors = self.line_vectors()?;
        let (bases, per_line) = matroid_basis_counts(&vectors, self.family.loop_momenta().len())?;
        let mut lines = colors.into_iter().zip(per_line).collect::<Vec<_>>();
        lines.sort_unstable();
        Some(MatroidProfile { bases, lines })
    }

    /// The verified power-preserving unimodular loop-momentum maps from this
    /// input onto `target` at the level of the line matroid, offered to
    /// `accept` one at a time until it returns `true` (then `Some(true)`).
    ///
    /// Unlike [`Self::find_power_preserving_map`], the two incidence graphs
    /// need not be isomorphic: a vacuum integral depends on its line momenta
    /// only, so 2-isomorphic graphs (related by Whitney twists at a 2-vertex
    /// cut, e.g. two embeddings of one contracted master) carry the same
    /// integral. The search assigns the lines of a loop-momentum basis of
    /// this input to lines of equal power of the target (with a sign; the
    /// global sign is fixed because scalar products are even), derives every
    /// other line's image from its exact basis coordinates and prunes as soon
    /// as an image is not a remaining target line of equal power. Each
    /// complete assignment gives one integer loop map, which FeynKit verifies
    /// denominator by denominator with the Jacobian; the powers are then
    /// re-checked on the verified denominator map. `None` if a power is not a
    /// positive integer, a coefficient is not an integer, or the basis of this
    /// input is not unimodular.
    pub fn find_power_preserving_matroid_map(
        &self,
        target: &MomentumInput,
        mut accept: impl FnMut(&[Atom]) -> Result<bool, VakintError>,
    ) -> Result<Option<bool>, VakintError> {
        let (Some(source_colors), Some(target_colors)) =
            (self.power_colors(), target.power_colors())
        else {
            return Ok(None);
        };
        let (Some(source), Some(destination)) = (self.line_vectors(), target.line_vectors()) else {
            return Ok(None);
        };
        let loops = self.family.loop_momenta().len();
        if loops != target.family.loop_momenta().len() || source.len() != destination.len() {
            return Ok(Some(false));
        }
        let basis_rows = self
            .basis
            .iter()
            .map(|&line| source[line].clone())
            .collect::<Vec<_>>();
        let determinant = integer_determinant(basis_rows.clone());
        if determinant.abs() != 1 {
            return Ok(None);
        }
        // Inverse of the (unimodular) basis matrix: adjugate times det.
        let inverse = adjugate(&basis_rows)
            .into_iter()
            .map(|row| row.into_iter().map(|value| value * determinant).collect())
            .collect::<Vec<Vec<i128>>>();
        let coordinates = source
            .iter()
            .map(|vector| row_times(vector, &inverse))
            .collect::<Vec<_>>();
        // A non-basis line is checked as soon as the last basis line of its
        // coordinate support has an image.
        let mut checks = vec![Vec::new(); loops];
        for (line, coordinate) in coordinates.iter().enumerate() {
            if self.basis.contains(&line) {
                continue;
            }
            let Some(last) = (0..loops).rev().find(|&axis| coordinate[axis] != 0) else {
                return Ok(None);
            };
            checks[last].push(line);
        }
        let mut classes = BTreeMap::<(Vec<i128>, u8), usize>::new();
        for (vector, color) in destination.iter().zip(&target_colors) {
            *classes.entry((normalized(vector), *color)).or_default() += 1;
        }
        let mut search = MatroidSearch {
            source: self,
            target,
            source_colors: &source_colors,
            target_colors: &target_colors,
            basis: &self.basis,
            inverse: &inverse,
            coordinates: &coordinates,
            checks: &checks,
            classes: classes.into_iter().collect(),
            images: Vec::with_capacity(loops),
            seen: BTreeSet::new(),
            accept: &mut accept,
        };
        Ok(Some(search.assign()?))
    }

    #[cfg(test)]
    pub fn match_topologies_with_acceptance<'a, K: Ord>(
        &self,
        topologies: impl IntoIterator<Item = &'a Topology>,
        limits: MomentumMatchLimits,
        accept: impl FnMut(&ReplacementRules) -> Result<Option<K>, VakintError>,
    ) -> Result<Option<ReplacementRules>, VakintError> {
        self.match_topologies_with_selection(topologies, limits, |witnesses| {
            first_minimal_acceptance(witnesses, accept)
        })
    }

    /// Enumerate every verified witness under the shared budget, then let the
    /// caller select one by index. `None` from the selector, or an index out
    /// of range, keeps the ordinary first witness, so a complete scope miss
    /// has a deterministic fallback for partial reduction. Without any
    /// verified witness the selector is not called.
    pub fn match_topologies_with_selection<'a>(
        &self,
        topologies: impl IntoIterator<Item = &'a Topology>,
        limits: MomentumMatchLimits,
        select: impl FnOnce(&[ReplacementRules]) -> Result<Option<usize>, VakintError>,
    ) -> Result<Option<ReplacementRules>, VakintError> {
        let mut witnesses = self.verified_witnesses(topologies, limits)?;
        if witnesses.is_empty() {
            log::debug!(target:"vakint::five_loop_match", "phase=admitted_witness_search verified_witnesses=0 admitted=false");
            return Ok(None);
        }
        let selected = select(&witnesses)?;
        log::debug!(target:"vakint::five_loop_match", "phase=admitted_witness_search verified_witnesses={} admitted={} selected={selected:?}", witnesses.len(), selected.is_some());
        let index = selected
            .filter(|index| *index < witnesses.len())
            .unwrap_or(0);
        Ok(Some(witnesses.swap_remove(index)))
    }

    /// Every verified exact momentum witness, in the deterministic order of
    /// the native proposals: the ordinary incidence witness of each topology
    /// first, then its edge-automorphism images. Coordinates repeated by
    /// different proposals are verified once.
    pub fn verified_witnesses<'a>(
        &self,
        topologies: impl IntoIterator<Item = &'a Topology>,
        limits: MomentumMatchLimits,
    ) -> Result<Vec<ReplacementRules>, VakintError> {
        let topologies = topologies.into_iter().collect::<Vec<_>>();
        let mut witnesses = Vec::new();
        let source_graph = edge_orbits::EdgeGraph::new(&self.incidence)?;
        let loops = self.family.loop_momenta();
        let variables = (1..=loops.len())
            .map(|axis| Atom::var(vk_symbol!(format!("momentum_match_source_{axis}"))))
            .collect::<Vec<_>>();
        let substitute = loops
            .iter()
            .zip(&variables)
            .map(|(source, variable)| Replacement::new(source.clone(), variable.clone()))
            .collect::<Vec<_>>();
        let sources = self
            .basis
            .iter()
            .map(|&slot| {
                self.propagators[slot]
                    .momentum
                    .replace_multiple(&substitute)
            })
            .collect::<Vec<_>>();
        let orientations = 1usize
            .checked_shl(loops.len().try_into().map_err(invalid)?)
            .ok_or_else(|| invalid("momentum orientation count overflow"))?;
        let mut tried = 0;
        let mut verified = 0;
        let mut graph_maps = 0;
        for topology in topologies {
            let integral = topology.get_integral();
            if integral.n_loops != loops.len()
                || integral.graph.edges.len() != self.propagators.len()
            {
                continue;
            }
            let mut unit_graph = integral.graph.clone();
            for edge in unit_graph.edges.values_mut() {
                edge.mass = Atom::one();
            }
            let target_graph = edge_orbits::EdgeGraph::new(&unit_graph.to_symbolica_graph(false))?;
            if source_graph.canonical.graph != target_graph.canonical.graph {
                continue;
            }
            let expression = integral
                .canonical_expression
                .as_ref()
                .ok_or_else(|| invalid("missing canonical topology"))?;
            let physical = (1..=integral.n_props)
                .filter_map(|id| {
                    get_prop_with_id(expression.as_view(), id)
                        .map(|propagator| (id, propagator[&vk_symbol!("q_")].clone()))
                })
                .collect::<Vec<_>>();
            if physical.len() != self.propagators.len() {
                return Err(invalid("canonical graph and propagators disagree"));
            }
            let kinematics = self.family.kinematics();
            let denominators = physical
                .iter()
                .map(|(_, momentum)| {
                    kinematics
                        .scalar_product(momentum, momentum)
                        .map(|s| s - Atom::one())
                })
                .collect::<Result<Vec<_>, _>>()
                .map_err(invalid)?;
            let target = IntegralFamily::new(loops.to_vec(), Vec::new(), denominators, kinematics)
                .map_err(invalid)?;
            let mut seen_coordinates = BTreeSet::new();
            let mut seen_permutations = BTreeSet::new();
            let mut evaluate = |permutation: &[usize]| {
                if !seen_permutations.insert(permutation.to_vec()) {
                    return Ok(());
                }
                let images = self
                    .basis
                    .iter()
                    .map(|&slot| &physical[permutation[slot]].1)
                    .collect::<Vec<_>>();
                for signs in 0..orientations {
                    charge(
                        &mut tried,
                        limits
                            .diagnostics_enabled
                            .then_some(limits.max_orientation_candidates),
                        "momentum orientation candidates",
                    )?;
                    let system = sources
                        .iter()
                        .zip(&images)
                        .enumerate()
                        .map(|(axis, (source, image))| {
                            source
                                - Atom::num(if signs & (1 << axis) == 0 { 1 } else { -1 })
                                    * (*image)
                        })
                        .collect::<Vec<_>>();
                    let coordinates = Topology::solve_loop_coordinates(&system, &variables)?;
                    if !seen_coordinates.insert(coordinates.clone()) {
                        continue;
                    }
                    let Some(mapping) = self
                        .family
                        .mapping_to(&target, &coordinates)
                        .map_err(invalid)?
                    else {
                        continue;
                    };
                    charge(
                        &mut verified,
                        limits
                            .diagnostics_enabled
                            .then_some(limits.max_verified_witnesses),
                        "verified momentum witnesses",
                    )?;
                    if verified <= 8 {
                        log::debug!(target: "vakint::five_loop_match",
                            "phase=verified_alternative ordinal={verified} topology={} source_rows={} target_rows={} denominator_map={:?} coordinates={:?}",
                            integral.name,self.propagators.len(),physical.len(),mapping.denominator_map(),coordinates);
                    }
                    witnesses.push(self.replacement(topology, &physical, &mapping));
                }
                Ok(())
            };
            // Try precisely the ordinary incidence witness first, under the
            // same budget, so a complete scope miss has a deterministic fallback.
            if let Some(first) = self.incidence_permutation(&integral.graph) {
                charge(
                    &mut graph_maps,
                    limits
                        .diagnostics_enabled
                        .then_some(limits.max_graph_automorphisms),
                    "native edge automorphism witnesses",
                )?;
                evaluate(&first)?;
            }
            edge_orbits::visit(
                &source_graph,
                &target_graph,
                limits
                    .diagnostics_enabled
                    .then_some(limits.max_graph_automorphisms),
                &mut graph_maps,
                &mut evaluate,
            )?;
        }
        log::debug!(target:"vakint::five_loop_match", "phase=verified_witness_enumeration graph_maps={graph_maps} orientation_candidates={tried} verified_witnesses={verified}");
        Ok(witnesses)
    }
}

/// Power-colored profile of the cycle matroid of a list of lines: the number
/// of loop-momentum bases (independent `L`-subsets of line momenta) and, sorted,
/// each line's (power, number of bases containing it). A momentum map sends
/// lines to lines of equal power by an invertible linear map, so inputs with
/// different profiles have no such map.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MatroidProfile {
    bases: usize,
    lines: Vec<(u8, usize)>,
}

/// Integer coefficients of a homogeneous linear `momentum` over `loops`, or
/// `None` when it is not such a combination.
pub(crate) fn integer_line_vector(momentum: &Atom, loops: &[Atom]) -> Option<Vec<i128>> {
    let variables = loops.iter().collect::<Vec<_>>();
    let mut row = vec![0i128; loops.len()];
    for (monomial, coefficient) in momentum.expand().coefficient_list::<u16>(&variables) {
        let axis = loops.iter().position(|variable| variable == &monomial)?;
        row[axis] += i128::from(get_integer_from_atom(coefficient.as_view())?);
    }
    Some(row)
}

/// Exact determinant of a square integer matrix (fraction-free Bareiss
/// elimination).
fn integer_determinant(mut matrix: Vec<Vec<i128>>) -> i128 {
    let n = matrix.len();
    let mut sign = 1;
    let mut previous = 1;
    for pivot in 0..n {
        let Some(row) = (pivot..n).find(|&row| matrix[row][pivot] != 0) else {
            return 0;
        };
        if row != pivot {
            matrix.swap(pivot, row);
            sign = -sign;
        }
        for row in pivot + 1..n {
            for column in pivot + 1..n {
                matrix[row][column] = (matrix[row][column] * matrix[pivot][pivot]
                    - matrix[row][pivot] * matrix[pivot][column])
                    / previous;
            }
            matrix[row][pivot] = 0;
        }
        previous = matrix[pivot][pivot];
    }
    sign * previous
}

/// Adjugate (transposed cofactor matrix) of a square integer matrix.
fn adjugate(matrix: &[Vec<i128>]) -> Vec<Vec<i128>> {
    let n = matrix.len();
    let cofactor = |row: usize, column: usize| {
        let minor = matrix
            .iter()
            .enumerate()
            .filter(|(index, _)| *index != row)
            .map(|(_, values)| {
                values
                    .iter()
                    .enumerate()
                    .filter(|(index, _)| *index != column)
                    .map(|(_, value)| *value)
                    .collect()
            })
            .collect();
        let sign = if (row + column).is_multiple_of(2) {
            1
        } else {
            -1
        };
        sign * integer_determinant(minor)
    };
    // Entry (column, row) is the (row, column) cofactor.
    (0..n)
        .map(|column| (0..n).map(|row| cofactor(row, column)).collect())
        .collect()
}

/// Row vector times matrix.
fn row_times(vector: &[i128], matrix: &[Vec<i128>]) -> Vec<i128> {
    (0..matrix.first().map_or(0, Vec::len))
        .map(|column| {
            vector
                .iter()
                .zip(matrix)
                .map(|(value, row)| value * row[column])
                .sum()
        })
        .collect()
}

/// The sign representative of a line momentum: first nonzero entry positive.
fn normalized(vector: &[i128]) -> Vec<i128> {
    match vector.iter().find(|value| **value != 0) {
        Some(first) if *first < 0 => vector.iter().map(|value| -value).collect(),
        _ => vector.to_vec(),
    }
}

/// The number of linearly independent `rank`-subsets of `vectors` (the bases
/// of the represented matroid) and, per vector, the number of them containing
/// it; `None` if the vectors do not have length `rank`.
pub(crate) fn matroid_basis_counts(
    vectors: &[Vec<i128>],
    rank: usize,
) -> Option<(usize, Vec<usize>)> {
    if rank == 0 || vectors.len() < rank || vectors.iter().any(|vector| vector.len() != rank) {
        return None;
    }
    let mut total = 0usize;
    let mut per_line = vec![0usize; vectors.len()];
    let mut subset = (0..rank).collect::<Vec<_>>();
    loop {
        if integer_determinant(subset.iter().map(|&line| vectors[line].clone()).collect()) != 0 {
            total += 1;
            for &line in &subset {
                per_line[line] += 1;
            }
        }
        // Next combination in lexicographic order.
        let mut index = rank;
        while index > 0 && subset[index - 1] == vectors.len() - rank + index - 1 {
            index -= 1;
        }
        if index == 0 {
            break;
        }
        subset[index - 1] += 1;
        for next in index..rank {
            subset[next] = subset[next - 1] + 1;
        }
    }
    Some((total, per_line))
}

/// Backtracking state of [`MomentumInput::find_power_preserving_matroid_map`].
struct MatroidSearch<'a> {
    source: &'a MomentumInput,
    target: &'a MomentumInput,
    source_colors: &'a [u8],
    target_colors: &'a [u8],
    /// The source loop-momentum basis (line indices) in assignment order.
    basis: &'a [usize],
    /// Inverse of the source basis matrix.
    inverse: &'a [Vec<i128>],
    /// Coordinates of every source line in the source basis.
    coordinates: &'a [Vec<i128>],
    /// Non-basis lines determined once basis position `t` has an image.
    checks: &'a [Vec<usize>],
    /// Remaining target lines by (sign representative, power).
    classes: Vec<((Vec<i128>, u8), usize)>,
    /// Images (target coefficient vectors) of the assigned basis lines.
    images: Vec<Vec<i128>>,
    seen: BTreeSet<Vec<Vec<i128>>>,
    accept: &'a mut dyn FnMut(&[Atom]) -> Result<bool, VakintError>,
}

impl MatroidSearch<'_> {
    fn take(&mut self, vector: &[i128], color: u8) -> Option<usize> {
        let key = (normalized(vector), color);
        let index = self
            .classes
            .binary_search_by(|(class, _)| class.cmp(&key))
            .ok()?;
        if self.classes[index].1 == 0 {
            return None;
        }
        self.classes[index].1 -= 1;
        Some(index)
    }

    fn assign(&mut self) -> Result<bool, VakintError> {
        crate::execution::checkpoint()?;
        let position = self.images.len();
        if position == self.basis.len() {
            return self.complete();
        }
        let color = self.source_colors[self.basis[position]];
        for class in 0..self.classes.len() {
            let ((ref vector, class_color), remaining) = self.classes[class];
            if class_color != color || remaining == 0 {
                continue;
            }
            let vector = vector.clone();
            // The global sign is immaterial (scalar products are even).
            let signs: &[i128] = if position == 0 { &[1] } else { &[1, -1] };
            for &sign in signs {
                self.classes[class].1 -= 1;
                self.images
                    .push(vector.iter().map(|value| sign * value).collect());
                let mut taken = Vec::new();
                let mut consistent = true;
                for &line in &self.checks[position] {
                    let image = row_times(&self.coordinates[line][..=position], &self.images);
                    match self.take(&image, self.source_colors[line]) {
                        Some(index) => taken.push(index),
                        None => {
                            consistent = false;
                            break;
                        }
                    }
                }
                let accepted = consistent && self.assign()?;
                for index in taken {
                    self.classes[index].1 += 1;
                }
                self.images.pop();
                self.classes[class].1 += 1;
                if accepted {
                    return Ok(true);
                }
            }
        }
        Ok(false)
    }

    /// One complete assignment: the loop map `T = B^-1 W`, FeynKit
    /// verification, power check on the verified denominator map, `accept`.
    fn complete(&mut self) -> Result<bool, VakintError> {
        if integer_determinant(self.images.clone()).abs() != 1 {
            return Ok(false);
        }
        let map = self
            .inverse
            .iter()
            .map(|row| row_times(row, &self.images))
            .collect::<Vec<_>>();
        if !self.seen.insert(map.clone()) {
            return Ok(false);
        }
        let loops = self.target.family.loop_momenta();
        let images = map
            .iter()
            .map(|row| {
                row.iter()
                    .zip(loops)
                    .map(|(coefficient, momentum)| Atom::num(*coefficient as i64) * momentum)
                    .sum::<Atom>()
            })
            .collect::<Vec<_>>();
        let Some(mapping) = self
            .source
            .family
            .mapping_to(&self.target.family, &images)
            .map_err(invalid)?
        else {
            return Ok(false);
        };
        if mapping
            .denominator_map()
            .iter()
            .enumerate()
            .any(|(source, target)| self.source_colors[source] != self.target_colors[*target])
        {
            return Ok(false);
        }
        let images = mapping
            .momentum_rules()
            .iter()
            .map(|(_, image)| image.clone())
            .collect::<Vec<_>>();
        (self.accept)(&images)
    }
}
