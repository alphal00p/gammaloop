//! Exact permutation-group operations for signed canonical symmetry sites.
//!
//! A group element always carries the complete canonical-vertex permutation.
//! The signed site projection is additional data, not a replacement for that
//! permutation: scope stabilizers must first be selected using graph vertices.

#[cfg(test)]
use std::collections::BTreeMap;
use std::collections::{BTreeSet, HashSet, VecDeque};

use linnet::permutation::Permutation;
use thiserror::Error;

/// One antisymmetric site's fixed canonical port frame.
///
/// `layout_path` distinguishes multiple symmetry groups owned by the same
/// vertex. Exchangeable Young columns instead have distinct owner vertices
/// and deliberately share a height-based path. `members` are ordered in the
/// canonical frame, so transporting a whole column is unsigned while reversing
/// its internal frame remains odd. Strict private-carrier numeric scalar sign
/// payloads use a reserved path and the singleton owner frame, whose transport
/// has no intrinsic parity.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct SignSiteFrame {
    pub(crate) owner: usize,
    pub(crate) layout_path: Vec<usize>,
    pub(crate) members: Vec<usize>,
}

impl SignSiteFrame {
    fn validate(&self, site: usize, vertex_count: usize) -> Result<(), SignedGroupError> {
        if self.members.is_empty() {
            return Err(SignedGroupError::EmptySiteFrame { site });
        }
        if self.owner >= vertex_count {
            return Err(SignedGroupError::VertexOutOfRange {
                vertex: self.owner,
                vertex_count,
            });
        }
        let mut seen = HashSet::new();
        for &member in &self.members {
            if member >= vertex_count {
                return Err(SignedGroupError::VertexOutOfRange {
                    vertex: member,
                    vertex_count,
                });
            }
            if !seen.insert(member) {
                return Err(SignedGroupError::RepeatedFrameMember { site, member });
            }
        }
        Ok(())
    }
}

/// The image of one sign site under a group element.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) struct SignedSiteImage {
    pub(crate) target: usize,
    pub(crate) odd: bool,
}

/// A paired complete-vertex and affine signed-site action.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) struct SignedAction {
    vertices: Permutation,
    sites: Vec<SignedSiteImage>,
}

#[derive(Clone, Debug, Error, PartialEq, Eq)]
pub(crate) enum SignedGroupError {
    #[error("invalid graph automorphism cycles: {0}")]
    InvalidCycles(String),
    #[error("vertex {vertex} is outside the canonical graph domain 0..{vertex_count}")]
    VertexOutOfRange { vertex: usize, vertex_count: usize },
    #[error("site {site} has an empty canonical frame")]
    EmptySiteFrame { site: usize },
    #[error("site {site} repeats canonical frame member {member}")]
    RepeatedFrameMember { site: usize, member: usize },
    #[error("site {site} does not have exactly one image under the vertex map")]
    InvalidSiteImage { site: usize },
    #[error("the site targets do not form a permutation")]
    InvalidSitePermutation,
    #[error("the canonical vertex map is not a permutation")]
    InvalidVertexPermutation,
    #[error("signed actions belong to different vertex or site domains")]
    IncompatibleDomains,
    #[error("fixed site {site} is outside the site domain 0..{site_count}")]
    SiteOutOfRange { site: usize, site_count: usize },
    #[error("signed decoration has {actual} sites, expected {expected}")]
    DecorationSize { actual: usize, expected: usize },
    #[error("active-site mask has {actual} sites, expected {expected}")]
    SiteMaskSize { actual: usize, expected: usize },
    #[error("active-site mask is not invariant under site transport {source_site} -> {target}")]
    NonInvariantSiteMask { source_site: usize, target: usize },
    #[error(
        "signed decoration orbit contains at least {observed_at_least} states across {sites} active sites, exceeding limit {limit}"
    )]
    DecorationOrbitLimit {
        observed_at_least: usize,
        limit: usize,
        sites: usize,
    },
}

impl SignedAction {
    pub(crate) fn identity(vertex_count: usize, site_count: usize) -> Self {
        Self {
            vertices: Permutation::id(vertex_count),
            sites: (0..site_count)
                .map(|target| SignedSiteImage { target, odd: false })
                .collect(),
        }
    }

    pub(crate) fn from_parts(
        vertices: Permutation,
        sites: Vec<SignedSiteImage>,
    ) -> Result<Self, SignedGroupError> {
        let mut targets = vec![false; sites.len()];
        for image in &sites {
            let Some(seen) = targets.get_mut(image.target) else {
                return Err(SignedGroupError::InvalidSitePermutation);
            };
            if std::mem::replace(seen, true) {
                return Err(SignedGroupError::InvalidSitePermutation);
            }
        }
        Ok(Self { vertices, sites })
    }

    /// Return `self ∘ earlier`, applying `earlier` first.
    pub(crate) fn compose(&self, earlier: &Self) -> Result<Self, SignedGroupError> {
        if self.vertex_count() != earlier.vertex_count()
            || self.site_count() != earlier.site_count()
        {
            return Err(SignedGroupError::IncompatibleDomains);
        }

        let sites = earlier
            .sites
            .iter()
            .map(|first| {
                let second = self.sites[first.target];
                SignedSiteImage {
                    target: second.target,
                    odd: first.odd ^ second.odd,
                }
            })
            .collect();
        Self::from_parts(self.vertices.compose(&earlier.vertices), sites)
    }

    pub(crate) fn inverse(&self) -> Self {
        let mut sites = vec![
            SignedSiteImage {
                target: 0,
                odd: false
            };
            self.site_count()
        ];
        for (source, image) in self.sites.iter().enumerate() {
            sites[image.target] = SignedSiteImage {
                target: source,
                odd: image.odd,
            };
        }
        Self {
            vertices: self.vertices.inverse(),
            sites,
        }
    }

    pub(crate) fn vertex_count(&self) -> usize {
        self.vertices.map().len()
    }

    pub(crate) fn vertex_map(&self) -> &[usize] {
        self.vertices.map()
    }

    pub(crate) fn site_count(&self) -> usize {
        self.sites.len()
    }

    #[cfg(test)]
    pub(crate) fn site_images(&self) -> &[SignedSiteImage] {
        &self.sites
    }

    #[cfg(test)]
    pub(crate) fn vertex_image(&self, vertex: usize) -> Result<usize, SignedGroupError> {
        self.vertices
            .map()
            .get(vertex)
            .copied()
            .ok_or(SignedGroupError::VertexOutOfRange {
                vertex,
                vertex_count: self.vertex_count(),
            })
    }

    #[cfg(test)]
    pub(crate) fn site_image(&self, site: usize) -> Result<SignedSiteImage, SignedGroupError> {
        self.sites
            .get(site)
            .copied()
            .ok_or(SignedGroupError::SiteOutOfRange {
                site,
                site_count: self.site_count(),
            })
    }

    /// Apply the affine site action: the phase at a source site is transported
    /// to its target and XORed with the frame-relative parity.
    pub(crate) fn apply_site_phases(&self, phases: &[bool]) -> Result<Vec<bool>, SignedGroupError> {
        if phases.len() != self.site_count() {
            return Err(SignedGroupError::DecorationSize {
                actual: phases.len(),
                expected: self.site_count(),
            });
        }
        let mut result = vec![false; phases.len()];
        for (source, image) in self.sites.iter().enumerate() {
            result[image.target] = phases[source] ^ image.odd;
        }
        Ok(result)
    }

    pub(crate) fn is_identity(&self) -> bool {
        self.vertices.is_identity()
            && self
                .sites
                .iter()
                .enumerate()
                .all(|(site, image)| image.target == site && !image.odd)
    }

    fn cmp_key(&self, other: &Self) -> std::cmp::Ordering {
        self.vertices
            .map()
            .cmp(other.vertices.map())
            .then_with(|| self.sites.cmp(&other.sites))
    }
}

/// Site transport induced by an input-to-canonical map or canonical
/// automorphism. Frame parity is always measured relative to the target's
/// canonical frame.
pub(crate) fn transport_site_frames(
    vertex_map: &[usize],
    source_frames: &[SignSiteFrame],
    target_frames: &[SignSiteFrame],
) -> Result<Vec<SignedSiteImage>, SignedGroupError> {
    let vertex_count = vertex_map.len();
    if let Some(&vertex) = vertex_map.iter().find(|&&vertex| vertex >= vertex_count) {
        return Err(SignedGroupError::VertexOutOfRange {
            vertex,
            vertex_count,
        });
    }
    let mut vertex_targets = vec![false; vertex_count];
    if vertex_map
        .iter()
        .any(|&target| std::mem::replace(&mut vertex_targets[target], true))
    {
        return Err(SignedGroupError::InvalidVertexPermutation);
    }
    for (site, frame) in source_frames.iter().enumerate() {
        frame.validate(site, vertex_count)?;
    }
    for (site, frame) in target_frames.iter().enumerate() {
        frame.validate(site, vertex_count)?;
    }

    let images = source_frames
        .iter()
        .enumerate()
        .map(|(site, source)| {
            let mapped_owner = vertex_map.get(source.owner).copied().ok_or(
                SignedGroupError::VertexOutOfRange {
                    vertex: source.owner,
                    vertex_count,
                },
            )?;
            let mapped_members = source
                .members
                .iter()
                .map(|&member| {
                    vertex_map
                        .get(member)
                        .copied()
                        .ok_or(SignedGroupError::VertexOutOfRange {
                            vertex: member,
                            vertex_count,
                        })
                })
                .collect::<Result<Vec<_>, _>>()?;

            let candidates = target_frames
                .iter()
                .enumerate()
                .filter(|(_, target)| {
                    target.owner == mapped_owner
                        && target.layout_path == source.layout_path
                        && target.members.len() == mapped_members.len()
                        && mapped_members
                            .iter()
                            .all(|member| target.members.contains(member))
                })
                .collect::<Vec<_>>();
            let [(target_site, target)] = candidates.as_slice() else {
                return Err(SignedGroupError::InvalidSiteImage { site });
            };
            let member_map = mapped_members
                .iter()
                .map(|member| {
                    target
                        .members
                        .iter()
                        .position(|target| target == member)
                        .unwrap()
                })
                .collect();
            Ok(SignedSiteImage {
                target: *target_site,
                odd: Permutation::from_map(member_map).sign() < 0,
            })
        })
        .collect::<Result<Vec<_>, _>>()?;
    if images.len() != target_frames.len() {
        return Err(SignedGroupError::InvalidSitePermutation);
    }
    let mut targets = vec![false; images.len()];
    for image in &images {
        let Some(seen) = targets.get_mut(image.target) else {
            return Err(SignedGroupError::InvalidSitePermutation);
        };
        if std::mem::replace(seen, true) {
            return Err(SignedGroupError::InvalidSitePermutation);
        }
    }
    Ok(images)
}

/// A finite group generated by paired signed actions.
#[derive(Clone, Debug)]
pub(crate) struct SignedGroup {
    vertex_count: usize,
    site_count: usize,
    generators: Vec<SignedAction>,
}

impl SignedGroup {
    pub(crate) fn new(
        vertex_count: usize,
        site_count: usize,
        generators: Vec<SignedAction>,
    ) -> Result<Self, SignedGroupError> {
        if generators.iter().any(|generator| {
            generator.vertex_count() != vertex_count || generator.site_count() != site_count
        }) {
            return Err(SignedGroupError::IncompatibleDomains);
        }
        Ok(Self {
            vertex_count,
            site_count,
            generators: Self::normalized_generators(generators),
        })
    }

    /// Decode Graphica's canonical-numbered generators. Every outer entry is
    /// one generator; its inner disjoint cycles are materialized together.
    pub(crate) fn from_graphica(
        vertex_count: usize,
        site_frames: &[SignSiteFrame],
        orbit_generators: &[Vec<Vec<usize>>],
    ) -> Result<Self, SignedGroupError> {
        let generators = orbit_generators
            .iter()
            .map(|cycles| {
                let cycle_permutation = Permutation::from_disjoint_cycles(cycles)
                    .map_err(SignedGroupError::InvalidCycles)?;
                if cycle_permutation.map().len() > vertex_count {
                    return Err(SignedGroupError::VertexOutOfRange {
                        vertex: cycle_permutation.map().len() - 1,
                        vertex_count,
                    });
                }
                let mut vertex_map = cycle_permutation.map().to_vec();
                vertex_map.extend(vertex_map.len()..vertex_count);
                let sites = transport_site_frames(&vertex_map, site_frames, site_frames)?;
                SignedAction::from_parts(Permutation::from_map(vertex_map), sites)
            })
            .collect::<Result<Vec<_>, _>>()?;
        Self::new(vertex_count, site_frames.len(), generators)
    }

    pub(crate) fn generators(&self) -> &[SignedAction] {
        &self.generators
    }

    pub(crate) fn identity(&self) -> SignedAction {
        SignedAction::identity(self.vertex_count, self.site_count)
    }

    /// Choose the lexicographically least decoration in the exact affine
    /// orbit. The uncapped test helper is bounded by the boolean decoration
    /// domain rather than by the usually much larger graph automorphism group.
    #[cfg(test)]
    pub(crate) fn canonical_decoration(
        &self,
        phases: &[bool],
    ) -> Result<Vec<bool>, SignedGroupError> {
        Ok(self
            .decoration_orbit(phases, &vec![true; self.site_count], None)?
            .into_iter()
            .next()
            .unwrap())
    }

    /// Enumerate the exact observable affine decoration orbit in lexicographic
    /// order, stopping before an optional state-count limit would be exceeded.
    pub(crate) fn decoration_orbit(
        &self,
        phases: &[bool],
        active_sites: &[bool],
        limit: Option<usize>,
    ) -> Result<Vec<Vec<bool>>, SignedGroupError> {
        if phases.len() != self.site_count {
            return Err(SignedGroupError::DecorationSize {
                actual: phases.len(),
                expected: self.site_count,
            });
        }
        if active_sites.len() != self.site_count {
            return Err(SignedGroupError::SiteMaskSize {
                actual: active_sites.len(),
                expected: self.site_count,
            });
        }
        let sites = active_sites.iter().filter(|active| **active).count();
        if limit == Some(0) {
            return Err(SignedGroupError::DecorationOrbitLimit {
                observed_at_least: 1,
                limit: 0,
                sites,
            });
        }
        if self.site_count == 0 {
            return Ok(vec![Vec::new()]);
        }

        let generators = self.symmetric_generators();
        for generator in &generators {
            for (source, image) in generator.sites.iter().enumerate() {
                if active_sites[source] != active_sites[image.target] {
                    return Err(SignedGroupError::NonInvariantSiteMask {
                        source_site: source,
                        target: image.target,
                    });
                }
            }
        }

        let erase_inactive = |mut decoration: Vec<bool>| {
            for (phase, active) in decoration.iter_mut().zip(active_sites) {
                if !active {
                    *phase = false;
                }
            }
            decoration
        };
        let initial = erase_inactive(phases.to_vec());
        let mut seen = BTreeSet::from([initial.clone()]);
        let mut queue = VecDeque::from([initial]);
        while let Some(decoration) = queue.pop_front() {
            for generator in &generators {
                let image = erase_inactive(generator.apply_site_phases(&decoration)?);
                if !seen.contains(&image) {
                    if let Some(limit) = limit
                        && seen.len() == limit
                    {
                        return Err(SignedGroupError::DecorationOrbitLimit {
                            observed_at_least: limit.saturating_add(1),
                            limit,
                            sites,
                        });
                    }
                    seen.insert(image.clone());
                    queue.push_back(image);
                }
            }
        }
        Ok(seen.into_iter().collect())
    }

    #[cfg(test)]
    pub(crate) fn vertex_orbit_transversal(
        &self,
        vertex: usize,
    ) -> Result<OrbitTransversal, SignedGroupError> {
        self.check_vertex(vertex)?;
        self.orbit_transversal(self.vertex_count, vertex, |action, point| {
            action.vertices.map()[point]
        })
    }

    #[cfg(test)]
    pub(crate) fn site_orbit_transversal(
        &self,
        site: usize,
    ) -> Result<OrbitTransversal, SignedGroupError> {
        self.check_site(site)?;
        self.orbit_transversal(self.site_count, site, |action, point| {
            action.sites[point].target
        })
    }

    /// Exact subgroup fixing every listed canonical graph vertex pointwise.
    pub(crate) fn pointwise_vertex_stabilizer(
        &self,
        vertices: &[usize],
    ) -> Result<Self, SignedGroupError> {
        let mut stabilizer = self.clone();
        for &vertex in vertices {
            self.check_vertex(vertex)?;
            stabilizer = stabilizer.stabilizer(self.vertex_count, vertex, |action, point| {
                action.vertices.map()[point]
            })?;
        }
        Ok(stabilizer)
    }

    /// Exact subgroup fixing every listed sign-site target pointwise. This is
    /// useful only after graph-boundary stabilizers have been selected.
    pub(crate) fn pointwise_site_stabilizer(
        &self,
        sites: &[usize],
    ) -> Result<Self, SignedGroupError> {
        let mut stabilizer = self.clone();
        for &site in sites {
            self.check_site(site)?;
            stabilizer = stabilizer.stabilizer(self.site_count, site, |action, point| {
                action.sites[point].target
            })?;
        }
        Ok(stabilizer)
    }

    /// Exact subgroup preserving a set of sign sites as a set.
    ///
    /// This is a Schreier stabilizer on the induced subset action. In
    /// particular, it retains words whose individual input generators move a
    /// site outside the set before a later generator moves it back.
    #[cfg(test)]
    pub(crate) fn setwise_site_stabilizer(
        &self,
        sites: &BTreeSet<usize>,
    ) -> Result<Self, SignedGroupError> {
        for &site in sites {
            self.check_site(site)?;
        }
        let base = sites.iter().copied().collect::<Vec<_>>();
        let symmetric_generators = self.symmetric_generators();
        let mut representatives = BTreeMap::from([(base.clone(), self.identity())]);
        let mut queue = VecDeque::from([base.clone()]);
        while let Some(subset) = queue.pop_front() {
            let representative = representatives[&subset].clone();
            for generator in &symmetric_generators {
                let target = subset
                    .iter()
                    .map(|site| generator.sites[*site].target)
                    .collect::<BTreeSet<_>>()
                    .into_iter()
                    .collect::<Vec<_>>();
                if !representatives.contains_key(&target) {
                    representatives.insert(target.clone(), generator.compose(&representative)?);
                    queue.push_back(target);
                }
            }
        }

        let mut generators = Vec::new();
        for representative in representatives.values() {
            for generator in &symmetric_generators {
                let candidate = generator.compose(representative)?;
                let target = base
                    .iter()
                    .map(|site| candidate.sites[*site].target)
                    .collect::<BTreeSet<_>>()
                    .into_iter()
                    .collect::<Vec<_>>();
                let target_inverse = representatives[&target].inverse();
                let schreier = target_inverse.compose(&candidate)?;
                if !schreier.is_identity() {
                    debug_assert_eq!(
                        base.iter()
                            .map(|site| schreier.sites[*site].target)
                            .collect::<BTreeSet<_>>()
                            .into_iter()
                            .collect::<Vec<_>>(),
                        base
                    );
                    generators.push(schreier);
                }
            }
        }
        Self::new(self.vertex_count, self.site_count, generators)
    }

    /// Whether the exact subgroup fixing `site` contains an odd self-action.
    /// The test is made on its Schreier generators, on which parity is a
    /// one-dimensional character because every generator fixes the site.
    pub(crate) fn has_odd_site_stabilizer(&self, site: usize) -> Result<bool, SignedGroupError> {
        let stabilizer = self.pointwise_site_stabilizer(&[site])?;
        Ok(stabilizer
            .generators
            .iter()
            .any(|generator| generator.sites[site].odd))
    }

    fn stabilizer(
        &self,
        degree: usize,
        base: usize,
        image: impl Copy + Fn(&SignedAction, usize) -> usize,
    ) -> Result<Self, SignedGroupError> {
        let orbit = self.orbit_transversal(degree, base, image)?;
        let symmetric_generators = self.symmetric_generators();
        let mut generators = Vec::new();
        for representative in orbit.representatives.iter().flatten() {
            for generator in &symmetric_generators {
                let candidate = generator.compose(representative)?;
                let target = image(&candidate, base);
                let target_inverse = orbit.representatives[target]
                    .as_ref()
                    .expect("an orbit image has a transversal")
                    .inverse();
                let schreier = target_inverse.compose(&candidate)?;
                if !schreier.is_identity() {
                    debug_assert_eq!(image(&schreier, base), base);
                    generators.push(schreier);
                }
            }
        }
        Self::new(self.vertex_count, self.site_count, generators)
    }

    fn orbit_transversal(
        &self,
        degree: usize,
        base: usize,
        image: impl Copy + Fn(&SignedAction, usize) -> usize,
    ) -> Result<OrbitTransversal, SignedGroupError> {
        debug_assert!(base < degree);

        let generators = self.symmetric_generators();
        let mut representatives = vec![None; degree];
        representatives[base] = Some(self.identity());
        let mut queue = VecDeque::from([base]);
        while let Some(point) = queue.pop_front() {
            for generator in &generators {
                let target = image(generator, point);
                if representatives[target].is_none() {
                    let representative = generator.compose(
                        representatives[point]
                            .as_ref()
                            .expect("queued orbit points have transversals"),
                    )?;
                    representatives[target] = Some(representative);
                    queue.push_back(target);
                }
            }
        }
        Ok(OrbitTransversal {
            #[cfg(test)]
            base,
            representatives,
        })
    }

    fn check_vertex(&self, vertex: usize) -> Result<(), SignedGroupError> {
        if vertex < self.vertex_count {
            Ok(())
        } else {
            Err(SignedGroupError::VertexOutOfRange {
                vertex,
                vertex_count: self.vertex_count,
            })
        }
    }

    fn check_site(&self, site: usize) -> Result<(), SignedGroupError> {
        if site < self.site_count {
            Ok(())
        } else {
            Err(SignedGroupError::SiteOutOfRange {
                site,
                site_count: self.site_count,
            })
        }
    }

    pub(crate) fn symmetric_generators(&self) -> Vec<SignedAction> {
        Self::normalized_generators(
            self.generators
                .iter()
                .flat_map(|generator| [generator.clone(), generator.inverse()])
                .collect(),
        )
    }

    fn normalized_generators(generators: Vec<SignedAction>) -> Vec<SignedAction> {
        let mut seen = HashSet::new();
        let mut generators = generators
            .into_iter()
            .filter(|generator| !generator.is_identity() && seen.insert(generator.clone()))
            .collect::<Vec<_>>();
        generators.sort_by(SignedAction::cmp_key);
        generators
    }
}

/// A deterministic transversal: the representative at point `p` maps `base`
/// to `p` while retaining both paired actions.
#[derive(Clone, Debug)]
pub(crate) struct OrbitTransversal {
    #[cfg(test)]
    base: usize,
    representatives: Vec<Option<SignedAction>>,
}

impl OrbitTransversal {
    #[cfg(test)]
    pub(crate) fn base(&self) -> usize {
        self.base
    }

    #[cfg(test)]
    pub(crate) fn orbit(&self) -> impl Iterator<Item = usize> + '_ {
        self.representatives
            .iter()
            .enumerate()
            .filter_map(|(point, representative)| representative.as_ref().map(|_| point))
    }

    #[cfg(test)]
    pub(crate) fn representative(&self, point: usize) -> Option<&SignedAction> {
        self.representatives.get(point)?.as_ref()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use symbolica::graph::Graph;

    fn permutation(map: &[usize]) -> Permutation {
        Permutation::from_map(map.to_vec())
    }

    fn action(vertices: &[usize], sites: &[(usize, bool)]) -> SignedAction {
        SignedAction::from_parts(
            permutation(vertices),
            sites
                .iter()
                .map(|&(target, odd)| SignedSiteImage { target, odd })
                .collect(),
        )
        .unwrap()
    }

    fn generated_elements(group: &SignedGroup) -> HashSet<SignedAction> {
        let generators = group.symmetric_generators();
        let identity = group.identity();
        let mut elements = HashSet::from([identity.clone()]);
        let mut queue = VecDeque::from([identity]);
        while let Some(element) = queue.pop_front() {
            for generator in &generators {
                let next = generator.compose(&element).unwrap();
                if elements.insert(next.clone()) {
                    queue.push_back(next);
                }
            }
        }
        elements
    }

    #[test]
    fn composition_and_inverse_preserve_the_paired_law() {
        let first = action(&[1, 2, 0], &[(1, true), (2, false), (0, true)]);
        let second = action(&[0, 2, 1], &[(2, false), (1, true), (0, true)]);
        let composed = second.compose(&first).unwrap();

        assert_eq!(composed.vertex_image(0), Ok(2));
        assert_eq!(
            composed.site_image(0),
            Ok(SignedSiteImage {
                target: 1,
                odd: false
            })
        );
        let phases = [true, false, false];
        assert_eq!(
            composed.apply_site_phases(&phases).unwrap(),
            second
                .apply_site_phases(&first.apply_site_phases(&phases).unwrap())
                .unwrap()
        );
        assert!(first.inverse().compose(&first).unwrap().is_identity());
        assert!(first.compose(&first.inverse()).unwrap().is_identity());
    }

    #[test]
    fn graphica_outer_entry_is_one_complete_generator() {
        let frames = vec![
            SignSiteFrame {
                owner: 0,
                layout_path: vec![0],
                members: vec![1, 2],
            },
            SignSiteFrame {
                owner: 3,
                layout_path: vec![0],
                members: vec![4, 5],
            },
        ];
        let group = SignedGroup::from_graphica(
            8,
            &frames,
            &[vec![vec![0, 3], vec![1, 5], vec![2, 4], vec![6, 7]]],
        )
        .unwrap();

        assert_eq!(group.generators().len(), 1);
        let generator = &group.generators()[0];
        assert_eq!(generator.vertex_image(0), Ok(3));
        assert_eq!(generator.vertex_image(6), Ok(7));
        assert_eq!(
            generator.site_image(0),
            Ok(SignedSiteImage {
                target: 1,
                odd: true
            })
        );
    }

    #[test]
    fn graphica_canonical_numbering_contract_matches_the_group_decoder() {
        // The two branches are interchangeable, while the differently colored
        // center and marker anchor the canonical numbering:
        //
        //     leaf_a -- branch_a -- center -- branch_b -- leaf_b
        //                              |
        //                            marker
        let mut graph = Graph::<u8, u8>::new();
        let leaf_a = graph.add_node(2);
        let marker = graph.add_node(3);
        let branch_b = graph.add_node(1);
        let center = graph.add_node(0);
        let leaf_b = graph.add_node(2);
        let branch_a = graph.add_node(1);
        graph.add_edge(leaf_a, branch_a, false, 0).unwrap();
        graph.add_edge(marker, center, false, 0).unwrap();
        graph.add_edge(branch_b, center, false, 0).unwrap();
        graph.add_edge(leaf_b, branch_b, false, 0).unwrap();
        graph.add_edge(branch_a, center, false, 0).unwrap();

        let canonical = graph.canonize();

        // Graphica maps input vertices to canonical vertices. Its generators
        // already use canonical numbering, but omit fixed one-cycles.
        assert_eq!(canonical.vertex_map, vec![3, 5, 1, 0, 4, 2]);
        assert_eq!(
            canonical.orbit_generators,
            vec![vec![vec![3, 4], vec![1, 2]]]
        );

        let group = SignedGroup::from_graphica(
            canonical.graph.nodes().len(),
            &[],
            &canonical.orbit_generators,
        )
        .unwrap();
        assert_eq!(group.generators().len(), 1);
        assert_eq!(group.generators()[0].vertex_map(), [0, 2, 1, 4, 3, 5]);
    }

    #[test]
    fn direct_frame_parity_obeys_composition_and_inverse() {
        let frames = vec![
            SignSiteFrame {
                owner: 0,
                layout_path: vec![0],
                members: vec![1, 2],
            },
            SignSiteFrame {
                owner: 3,
                layout_path: vec![0],
                members: vec![4, 5],
            },
            SignSiteFrame {
                owner: 6,
                layout_path: vec![0],
                members: vec![7, 8],
            },
        ];
        let first = SignedGroup::from_graphica(
            9,
            &frames,
            &[vec![vec![0, 3, 6], vec![1, 5, 7], vec![2, 4, 8]]],
        )
        .unwrap()
        .generators()[0]
            .clone();
        let second = SignedGroup::from_graphica(
            9,
            &frames,
            &[vec![vec![0, 6, 3], vec![1, 8, 4], vec![2, 7, 5]]],
        )
        .unwrap()
        .generators()[0]
            .clone();
        let composed = second.compose(&first).unwrap();
        let direct_sites =
            transport_site_frames(composed.vertices.map(), &frames, &frames).unwrap();

        assert_eq!(composed.sites, direct_sites);
        let inverse_sites =
            transport_site_frames(first.vertices.inverse().map(), &frames, &frames).unwrap();
        assert_eq!(first.inverse().sites, inverse_sites);
    }

    #[test]
    fn schreier_stabilizer_finds_composed_odd_element() {
        let swap_zero_one = action(&[1, 0, 2], &[(0, true)]);
        let swap_zero_two = action(&[2, 1, 0], &[(0, true)]);
        let group = SignedGroup::new(3, 1, vec![swap_zero_one, swap_zero_two]).unwrap();

        let stabilizer = group.pointwise_vertex_stabilizer(&[0]).unwrap();
        assert!(stabilizer.generators().iter().any(|generator| {
            generator.vertices.map() == [0, 2, 1] && generator.site_image(0).unwrap().odd
        }));
        assert!(stabilizer.has_odd_site_stabilizer(0).unwrap());
    }

    #[test]
    fn pointwise_stabilizer_uses_vertices_before_signed_sites() {
        let boundary_swap = action(&[0, 2, 1], &[(0, true)]);
        let group = SignedGroup::new(3, 1, vec![boundary_swap]).unwrap();

        assert!(group.has_odd_site_stabilizer(0).unwrap());
        let local = group.pointwise_vertex_stabilizer(&[1, 2]).unwrap();
        assert!(local.generators().is_empty());
        assert!(!local.has_odd_site_stabilizer(0).unwrap());
    }

    #[test]
    fn generated_group_is_independent_of_generators_and_order() {
        let a = action(&[1, 0, 2], &[(0, true)]);
        let b = action(&[0, 2, 1], &[(0, true)]);
        let ab = a.compose(&b).unwrap();
        let first = SignedGroup::new(3, 1, vec![a.clone(), b.clone()]).unwrap();
        let reordered = SignedGroup::new(3, 1, vec![b, ab, a]).unwrap();

        assert_eq!(generated_elements(&first), generated_elements(&reordered));
        let orbit = first.vertex_orbit_transversal(0).unwrap();
        assert_eq!(orbit.base(), 0);
        assert_eq!(orbit.orbit().collect::<Vec<_>>(), vec![0, 1, 2]);
        assert!(orbit.representative(2).is_some());
        assert_eq!(first.site_orbit_transversal(0).unwrap().orbit().count(), 1);
    }

    #[test]
    fn affine_decoration_minimum_is_map_and_generator_order_independent() {
        let swap = action(&[0], &[(1, true), (0, false)]);
        let first = SignedGroup::new(1, 2, vec![swap.clone()]).unwrap();
        let reordered = SignedGroup::new(1, 2, vec![swap.inverse(), swap]).unwrap();
        let input = vec![false, false];
        let mapped = first.generators()[0].apply_site_phases(&input).unwrap();

        let expected = first.canonical_decoration(&input).unwrap();
        assert_eq!(first.canonical_decoration(&mapped).unwrap(), expected);
        assert_eq!(reordered.canonical_decoration(&input).unwrap(), expected);
    }

    #[test]
    fn empty_decoration_orbit_skips_vertex_only_generators_but_keeps_the_limit() {
        let group = SignedGroup::new(2, 0, vec![action(&[1, 0], &[])]).unwrap();

        assert_eq!(
            group.decoration_orbit(&[], &[], None).unwrap(),
            vec![Vec::<bool>::new()]
        );
        assert_eq!(
            group.decoration_orbit(&[], &[], Some(0)),
            Err(SignedGroupError::DecorationOrbitLimit {
                observed_at_least: 1,
                limit: 0,
                sites: 0,
            })
        );
    }

    #[test]
    fn inactive_sites_are_quotiented_from_the_affine_decoration_orbit() {
        let site_count = 12;
        let generators = (0..site_count)
            .map(|toggle| {
                let sites = (0..site_count)
                    .map(|site| (site, site == toggle))
                    .collect::<Vec<_>>();
                action(&[], &sites)
            })
            .collect();
        let group = SignedGroup::new(0, site_count, generators).unwrap();
        let phases = vec![false; site_count];

        assert_eq!(
            group
                .decoration_orbit(&phases, &vec![false; site_count], None)
                .unwrap(),
            vec![phases.clone()]
        );

        let mut active_sites = vec![false; site_count];
        active_sites[..2].fill(true);
        let orbit = group
            .decoration_orbit(&phases, &active_sites, None)
            .unwrap();
        assert_eq!(orbit.len(), 4);
        assert!(
            orbit
                .iter()
                .all(|decoration| decoration[2..].iter().all(|phase| !phase))
        );
    }

    #[test]
    fn decoration_orbit_limit_accepts_the_limit_and_rejects_the_next_state() {
        let site_count = 12;
        let generators = (0..site_count - 1)
            .map(|left| {
                let mut sites = (0..site_count)
                    .map(|site| (site, false))
                    .collect::<Vec<_>>();
                sites.swap(left, left + 1);
                action(&[], &sites)
            })
            .collect();
        let group = SignedGroup::new(0, site_count, generators).unwrap();
        let mut phases = vec![false; site_count];
        phases[site_count / 2..].fill(true);
        let active_sites = vec![true; site_count];
        let exact = group
            .decoration_orbit(&phases, &active_sites, None)
            .unwrap();
        assert_eq!(exact.len(), 924);
        assert_eq!(
            group
                .decoration_orbit(&phases, &active_sites, Some(exact.len()))
                .unwrap(),
            exact
        );
        assert_eq!(
            group.decoration_orbit(&phases, &active_sites, Some(256)),
            Err(SignedGroupError::DecorationOrbitLimit {
                observed_at_least: 257,
                limit: 256,
                sites: site_count,
            })
        );
    }

    #[test]
    fn active_site_mask_must_be_invariant_under_site_transport() {
        let group = SignedGroup::new(
            0,
            3,
            vec![action(&[], &[(2, false), (1, false), (0, false)])],
        )
        .unwrap();

        assert_eq!(
            group.decoration_orbit(&[false; 3], &[true, false, false], None),
            Err(SignedGroupError::NonInvariantSiteMask {
                source_site: 0,
                target: 2,
            })
        );
    }

    #[test]
    fn small_signed_groups_match_exhaustive_decoration_and_stabilizer_reference() {
        let cases = [
            (
                "even site swap",
                SignedGroup::new(2, 2, vec![action(&[1, 0], &[(1, false), (0, false)])]).unwrap(),
            ),
            (
                "affine site swap",
                SignedGroup::new(2, 2, vec![action(&[1, 0], &[(1, true), (0, false)])]).unwrap(),
            ),
            (
                "coupled signed transpositions",
                SignedGroup::new(
                    3,
                    3,
                    vec![
                        action(&[1, 0, 2], &[(1, true), (0, false), (2, false)]),
                        action(&[0, 2, 1], &[(0, false), (2, true), (1, false)]),
                    ],
                )
                .unwrap(),
            ),
        ];

        for (case, group) in cases {
            let elements = generated_elements(&group);
            for bits in 0..(1usize << group.site_count) {
                let phases = (0..group.site_count)
                    .map(|site| bits & (1 << site) != 0)
                    .collect::<Vec<_>>();
                let exhaustive = elements
                    .iter()
                    .map(|element| element.apply_site_phases(&phases).unwrap())
                    .min()
                    .unwrap();

                assert_eq!(
                    group.canonical_decoration(&phases).unwrap(),
                    exhaustive,
                    "canonical decoration differs for {case} and input {phases:?}",
                );
            }

            for site in 0..group.site_count {
                let exhaustive = elements.iter().any(|element| {
                    let image = element.site_image(site).unwrap();
                    image.target == site && image.odd
                });
                assert_eq!(
                    group.has_odd_site_stabilizer(site).unwrap(),
                    exhaustive,
                    "odd-site stabilizer differs for {case} at site {site}",
                );
            }
        }
    }

    #[test]
    fn setwise_stabilizer_retains_composed_words() {
        let move_zero_out = action(&[0], &[(2, true), (1, false), (0, false)]);
        let move_one_out = action(&[0], &[(0, false), (2, false), (1, false)]);
        let group = SignedGroup::new(1, 3, vec![move_zero_out, move_one_out]).unwrap();
        let sites = BTreeSet::from([0, 1]);

        assert!(group.generators().iter().all(|generator| {
            sites
                .iter()
                .map(|site| generator.site_images()[*site].target)
                .collect::<BTreeSet<_>>()
                != sites
        }));
        let stabilizer = group.setwise_site_stabilizer(&sites).unwrap();
        assert!(stabilizer.generators().iter().any(|generator| {
            sites
                .iter()
                .fold(false, |odd, site| odd ^ generator.site_images()[*site].odd)
        }));
    }
}
