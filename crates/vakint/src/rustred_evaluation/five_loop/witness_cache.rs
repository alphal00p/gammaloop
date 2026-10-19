//! Bounded, package-owned memo of the unchanged matcher's selected witnesses.
//!
//! Bytes count packed Symbolica atoms, not RSS: graph/pattern/container metadata
//! and shared symbol storage are excluded. The separate entry limit bounds
//! metadata growth. Oversized entries and a poisoned optional cache are misses,
//! never reduction failures. Admission and scalar policy are checked on hits.

use std::{collections::VecDeque, sync::Arc};

use rustred::{family::IntegralKey, scalar_numerator::ScalarNumeratorLimits};
use symbolica::atom::{Atom, AtomView};

use crate::{
    ReplacementRules, TensorReductionMethod, VakintSettings,
    topologies::{MomentumMatchLimits, RegistryIdentity},
};

pub(super) type AcceptedKey = Vec<(IntegralKey, Atom)>;

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct Policy {
    matching: MomentumMatchLimits,
    scalar: ScalarNumeratorLimits,
    epsilon: String,
    verify_numerator: bool,
    allow_unknown: bool,
    dot_notation: bool,
    tensor_method: TensorReductionMethod,
}

impl From<&VakintSettings> for Policy {
    fn from(settings: &VakintSettings) -> Self {
        Self {
            matching: settings.rustred_momentum_match_limits,
            scalar: settings.rustred_scalar_numerator_limits,
            epsilon: settings.epsilon_symbol.clone(),
            verify_numerator: settings.verify_numerator_identification,
            allow_unknown: settings.allow_unknown_integrals,
            dot_notation: settings.use_dot_product_notation,
            tensor_method: settings.tensor_reduction_method,
        }
    }
}

/// Statistics for frontend witness reuse, separate from IBP coefficient caches.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct RustRedWitnessCacheStatistics {
    pub hits: u64,
    pub misses: u64,
    pub evictions: u64,
    pub oversized_skips: u64,
    pub invalidations: u64,
    pub retained_entries: usize,
    /// Encoded Atom payload only, not total allocator usage or shared symbols.
    pub tracked_atom_bytes: usize,
}

pub(super) struct Entry {
    registry: RegistryIdentity,
    policy: Policy,
    integral: Atom,
    numerator: Atom,
    pub witness: ReplacementRules,
    pub accepted: AcceptedKey,
    atom_bytes: usize,
}

pub(super) struct WitnessCache {
    entries: VecDeque<Arc<Entry>>,
    max_entries: usize,
    max_atom_bytes: usize,
    statistics: RustRedWitnessCacheStatistics,
}

impl Default for WitnessCache {
    fn default() -> Self {
        Self::new(32, 16 * 1024 * 1024)
    }
}

impl WitnessCache {
    pub fn new(max_entries: usize, max_atom_bytes: usize) -> Self {
        Self {
            entries: VecDeque::new(),
            max_entries,
            max_atom_bytes,
            statistics: RustRedWitnessCacheStatistics::default(),
        }
    }

    pub fn get(
        &mut self,
        registry: &RegistryIdentity,
        policy: &Policy,
        integral: AtomView,
        numerator: AtomView,
    ) -> Option<Arc<Entry>> {
        let position = self.entries.iter().position(|entry| {
            &entry.registry == registry
                && &entry.policy == policy
                && entry.integral.as_view() == integral
                && entry.numerator.as_view() == numerator
        });
        let Some(position) = position else {
            self.statistics.misses = self.statistics.misses.saturating_add(1);
            return None;
        };
        let entry = self.entries.remove(position).expect("observed entry");
        self.entries.push_back(Arc::clone(&entry));
        self.statistics.hits = self.statistics.hits.saturating_add(1);
        Some(entry)
    }

    pub fn invalidate(&mut self, old: &Arc<Entry>) {
        if let Some(position) = self
            .entries
            .iter()
            .position(|entry| Arc::ptr_eq(entry, old))
        {
            let entry = self.entries.remove(position).expect("observed entry");
            self.statistics.tracked_atom_bytes -= entry.atom_bytes;
            self.statistics.retained_entries = self.entries.len();
            self.statistics.invalidations = self.statistics.invalidations.saturating_add(1);
        }
    }

    pub fn insert(
        &mut self,
        registry: RegistryIdentity,
        policy: Policy,
        integral: AtomView,
        numerator: AtomView,
        witness: &ReplacementRules,
        accepted: AcceptedKey,
    ) {
        // Count borrowed payload before cloning any potentially large input.
        let atom_bytes = integral
            .get_byte_size()
            .checked_add(numerator.get_byte_size());
        let atom_bytes = atom_bytes.and_then(|bytes| witness_atom_bytes(witness, bytes));
        let atom_bytes = accepted.iter().fold(atom_bytes, |bytes, (_, coefficient)| {
            bytes.and_then(|bytes| bytes.checked_add(coefficient.as_view().get_byte_size()))
        });
        let Some(atom_bytes) = atom_bytes.filter(|bytes| *bytes <= self.max_atom_bytes) else {
            self.statistics.oversized_skips = self.statistics.oversized_skips.saturating_add(1);
            return;
        };
        if self.max_entries == 0 || policy.epsilon.len() > self.max_atom_bytes {
            self.statistics.oversized_skips = self.statistics.oversized_skips.saturating_add(1);
            return;
        }
        // A concurrent successful miss may already have inserted the same key.
        if self.entries.iter().any(|entry| {
            entry.registry == registry
                && entry.policy == policy
                && entry.integral.as_view() == integral
                && entry.numerator.as_view() == numerator
        }) {
            return;
        }
        while self.entries.len() >= self.max_entries
            || self.statistics.tracked_atom_bytes > self.max_atom_bytes - atom_bytes
        {
            let Some(old) = self.entries.pop_front() else {
                break;
            };
            self.statistics.tracked_atom_bytes -= old.atom_bytes;
            self.statistics.evictions = self.statistics.evictions.saturating_add(1);
        }
        self.entries.push_back(Arc::new(Entry {
            registry,
            policy,
            integral: integral.to_owned(),
            numerator: numerator.to_owned(),
            witness: witness.clone(),
            accepted,
            atom_bytes,
        }));
        self.statistics.retained_entries = self.entries.len();
        self.statistics.tracked_atom_bytes += atom_bytes;
    }

    pub fn clear(&mut self) {
        self.entries.clear();
        self.statistics.retained_entries = 0;
        self.statistics.tracked_atom_bytes = 0;
    }

    pub fn statistics(&self) -> RustRedWitnessCacheStatistics {
        self.statistics
    }
}

pub(super) fn witness_atom_bytes(witness: &ReplacementRules, mut bytes: usize) -> Option<usize> {
    let mut add = |value: &Atom| -> Option<()> {
        bytes = bytes.checked_add(value.as_view().get_byte_size())?;
        Some(())
    };
    for (source, target) in witness
        .canonical_expression_substitutions
        .iter()
        .chain(witness.numerator_substitutions.iter())
    {
        add(source)?;
        add(target)?;
    }
    let integral = witness.canonical_topology.get_integral();
    for value in [
        &integral.canonical_expression,
        &integral.short_expression,
        &integral.alphaloop_expression,
    ]
    .into_iter()
    .flatten()
    {
        add(value)?;
    }
    if let Some((parent, coordinates)) = &integral.parent_routing {
        add(parent)?;
        for coordinate in coordinates.iter() {
            add(coordinate)?;
        }
    }
    for edge in integral.graph.edges.values() {
        add(&edge.momentum)?;
        add(&edge.mass)?;
    }
    Some(bytes)
}

#[cfg(test)]
#[path = "witness_cache_tests.rs"]
mod tests;
