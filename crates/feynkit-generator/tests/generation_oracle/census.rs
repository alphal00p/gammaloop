//! Differential oracle for diagram generation.
//!
//! A [`Census`] records what one generation run produced, in output order. [`Census::diff`]
//! compares two runs and names the first difference, coarse to fine as in the discrepancy
//! triage: report counts, the diagram set (by canonical form), symmetry factors (|Aut|), overall
//! factors with their signs, the finalized frames (with routing and cuts), numerators, diagram
//! ids, numerator groups, names and finally the output order. [`Summary`] condenses a census into
//! the numbers the oracle table records.
//!
//! Digests are FNV-1a hashes of fixed byte encodings (std's hashers may change between Rust
//! releases, which would invalidate the recorded rows). The key digest hashes `canonical_key()`,
//! which normalizes particle orientations and identifies the isomorphism class, and the frame
//! digest covers everything the structural canonical form behind a [`DiagramId`] is computed
//! from. Ids hash that form with the names of the model's particles and rules only, so they do
//! not depend on the process (its symbol registration order or thread count) and the rows record
//! a digest of them too.

use std::{
    collections::{BTreeMap, BTreeSet},
    fmt::{self, Write as _},
    sync::Arc,
};

use feynkit_generator::{GenerationReport, GenerationResult, SourceIdentity};
use feynkit_graph::{
    CanonicalDiagramKey, DiagramId, FeynmanDiagram, LoopMomentumBasis,
    expressions::evaluate_overall_factor,
};
use serde_json::Value;
use symbolica::{
    atom::{Atom, AtomCore},
    domains::rational::Rational,
};

/// One diagram of a run.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Entry {
    /// Digest of the canonical key, the canonical form of the diagram's isomorphism class.
    pub key: u64,
    /// Full owner key pairs diagrams without relying on the recorded digest.
    full_key: Arc<CanonicalDiagramKey>,
    /// Every serialized field except the separately retained typed ID and optional name.
    /// Sharing this value keeps policy normalization and restricted censuses inexpensive.
    payload: Arc<Value>,
    pub id: DiagramId,
    pub name: Option<String>,
    pub symmetry_factor: u64,
    pub overall_factor: String,
    /// The overall factor with its bookkeeping heads evaluated: `1 / |Aut|` with the fermion
    /// signs and the number of diagrams that symmetrization or grouping identified with this one.
    pub weight: Rational,
    /// Digest of the diagram in its finalized frame: vertices, edges, slots, loop momentum
    /// basis, cuts and threshold candidates, without numerators.
    pub frame: u64,
    /// Exact routing, retained for comparisons across threads and parent builds.
    pub basis: LoopMomentumBasis,
    /// Digest of the numerator, its prefactor, the projector and the local numerator fragments.
    pub numerator: u64,
}

/// A numerator group: the master's position and its members.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Group {
    pub master: usize,
    pub members: Vec<Member>,
}

/// Group provenance retains physical indices/factors even when names and IDs were omitted.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Member {
    pub source_index: usize,
    pub identity: Option<SourceIdentity>,
    pub ratio: String,
    pub overall_factor: String,
}

#[derive(Clone, Debug)]
pub struct Census {
    pub report: GenerationReport,
    pub diagrams: Vec<Entry>,
    pub groups: Vec<Group>,
}

/// What [`Census::diff`] compares.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Scope {
    /// Everything, including report counts, ids, groups, names and order.
    Exact,
    /// The diagram set with symmetry and overall factors, frames and numerators, whatever the
    /// names, the order and the search that produced them.
    Content,
}

/// FNV-1a over explicit byte encodings.
struct Digest(u64);

impl Digest {
    fn new() -> Self {
        Digest(0xcbf2_9ce4_8422_2325)
    }

    fn bytes(mut self, bytes: &[u8]) -> Self {
        for byte in bytes {
            self.0 = (self.0 ^ u64::from(*byte)).wrapping_mul(0x0000_0100_0000_01b3);
        }
        self
    }

    /// A string and a terminator, so that consecutive strings cannot run together.
    fn text(self, text: &str) -> Self {
        self.bytes(text.as_bytes()).bytes(&[0xff])
    }

    /// Preserve the previous named encoding; absence has its own tag.
    fn name(self, name: Option<&str>) -> Self {
        match name {
            Some(name) => self.text(name),
            None => self.bytes(&[0xfe]),
        }
    }

    fn number(self, value: u64) -> Self {
        self.bytes(&value.to_le_bytes())
    }

    fn json(self, value: &Value) -> Self {
        self.bytes(&serde_json::to_vec(value).unwrap())
    }

    /// Digest of `values` as a multiset.
    fn multiset(mut values: Vec<u64>) -> u64 {
        values.sort_unstable();
        values.into_iter().fold(Digest::new(), Digest::number).0
    }
}

/// The numbers the oracle table records for a case.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Summary {
    pub diagrams: usize,
    pub topologies: usize,
    pub assignments: usize,
    /// Sum of `1 / |Aut|`.
    pub inverse_symmetry: Rational,
    /// Sum of the evaluated overall factors, the [`Entry::weight`]s. Without symmetrization or
    /// grouping it is `inverse_symmetry` with the fermion signs; a wrong orbit or group weight
    /// changes it.
    pub weights: Rational,
    /// Multiset of |Aut|: (symmetry factor, number of diagrams).
    pub symmetry_factors: Vec<(u64, usize)>,
    /// Order- and name-independent digest of canonical forms, symmetry and overall factors and
    /// frames.
    pub content: u64,
    /// Order-independent digest of the numerators by canonical form.
    pub numerators: u64,
    /// Digest of names and canonical forms in output order and of the numerator groups.
    pub named: u64,
    /// Order-independent digest of the diagram ids.
    pub ids: u64,
}

impl fmt::Display for Summary {
    /// The summary as a row of the oracle table.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let histogram = self
            .symmetry_factors
            .iter()
            .map(|(s, n)| format!("({s}, {n})"))
            .collect::<Vec<_>>()
            .join(", ");
        write!(
            f,
            "Row {{ diagrams: {}, topologies: {}, assignments: {}, inverse_symmetry: ({}, {}), \
             weights: ({}, {}), symmetry_factors: &[{histogram}], content: 0x{:016x}, \
             numerators: 0x{:016x}, named: 0x{:016x}, ids: 0x{:016x} }}",
            self.diagrams,
            self.topologies,
            self.assignments,
            self.inverse_symmetry.numerator(),
            self.inverse_symmetry.denominator(),
            self.weights.numerator(),
            self.weights.denominator(),
            self.content,
            self.numerators,
            self.named,
            self.ids
        )
    }
}

fn canonical(atom: &Atom) -> String {
    atom.to_canonical_string()
}

impl Entry {
    fn of(diagram: &FeynmanDiagram) -> Entry {
        let full_key = Arc::new(diagram.canonical_key().unwrap());
        let key = serde_json::to_value(full_key.as_ref()).unwrap();
        // Preserve the u128 identity through DiagramId below, never through Value's number
        // representation. Drop its lossy parsed copy and the separately compared name before
        // retaining the full payload; model identity and every physical field remain exact.
        // The existing oracle digests still split numerators from the remaining frame.
        let mut serialized: Value =
            serde_json::from_slice(&serde_json::to_vec(diagram).unwrap()).unwrap();
        let fields = serialized.as_object_mut().unwrap();
        fields.remove("id").unwrap();
        fields.remove("name").unwrap();
        let payload = Arc::new(serialized.clone());
        let fields = serialized.as_object_mut().unwrap();
        let mut numerator = Digest::new();
        for field in ["numerator", "numerator_prefactor", "projector"] {
            numerator = numerator.json(&fields.remove(field).unwrap());
        }
        let vertices = fields["vertices"].as_array_mut().unwrap().iter_mut();
        let mut fragments: Vec<_> = vertices.map(|vertex| vertex["numerator"].take()).collect();
        let edges = fields["edges"].as_array_mut().unwrap().iter_mut();
        fragments.extend(edges.map(|edge| edge[1]["numerator"].take()));
        for fragment in &fragments {
            numerator = numerator.json(fragment);
        }
        for field in ["model", "symmetry_factor", "overall_factor"] {
            fields.remove(field);
        }
        let weight = evaluate_overall_factor(diagram.overall_factor().as_view());
        Entry {
            key: Digest::new().json(&key).0,
            full_key,
            payload,
            id: diagram.id().unwrap(),
            name: diagram.name().map(str::to_owned),
            symmetry_factor: diagram.symmetry_factor(),
            overall_factor: canonical(diagram.overall_factor()),
            weight: Rational::try_from(&weight).unwrap_or_else(|_| {
                panic!(
                    "{:?}: overall factor {weight} is not a number",
                    diagram.name()
                )
            }),
            frame: Digest::new().json(&serialized).0,
            basis: diagram.loop_momentum_basis().clone(),
            numerator: numerator.0,
        }
    }
}

impl Census {
    pub fn of(result: &GenerationResult) -> Census {
        let groups = result
            .groups
            .iter()
            .map(|group| Group {
                master: group.master,
                members: group
                    .members
                    .iter()
                    .map(|member| Member {
                        source_index: member.source_diagram,
                        identity: member.source_identity.clone(),
                        ratio: canonical(&member.ratio),
                        overall_factor: canonical(&member.overall_factor),
                    })
                    .collect(),
            })
            .collect();
        Census {
            report: result.report.clone(),
            diagrams: result.diagrams.iter().map(Entry::of).collect(),
            groups,
        }
    }

    fn by_key(&self) -> BTreeMap<&CanonicalDiagramKey, &Entry> {
        let by_key: BTreeMap<_, _> = self
            .diagrams
            .iter()
            .map(|d| (d.full_key.as_ref(), d))
            .collect();
        assert_eq!(
            by_key.len(),
            self.diagrams.len(),
            "duplicate canonical forms"
        );
        by_key
    }

    /// The diagrams for which `keep` holds, renumbered, without groups.
    pub fn restricted(&self, keep: impl Fn(&Entry) -> bool) -> Census {
        Census {
            report: self.report.clone(),
            diagrams: self.diagrams.iter().filter(|d| keep(d)).cloned().collect(),
            groups: vec![],
        }
    }

    pub fn summary(&self) -> Summary {
        let mut inverse_symmetry = Rational::zero();
        let mut weights = Rational::zero();
        let mut histogram = BTreeMap::<u64, usize>::new();
        for d in &self.diagrams {
            inverse_symmetry += Rational::from((1, d.symmetry_factor));
            weights += &d.weight;
            *histogram.entry(d.symmetry_factor).or_default() += 1;
        }
        let content = self.diagrams.iter().map(|d| {
            Digest::new()
                .number(d.key)
                .number(d.symmetry_factor)
                .text(&d.overall_factor)
                .number(d.frame)
                .0
        });
        let numerators = self
            .diagrams
            .iter()
            .map(|d| Digest::new().number(d.key).number(d.numerator).0);
        let mut named = Digest::new();
        for d in &self.diagrams {
            named = named.name(d.name.as_deref()).number(d.key);
        }
        for group in &self.groups {
            named = named.number(group.master as u64);
            for member in &group.members {
                named = named
                    .name(
                        member
                            .identity
                            .as_ref()
                            .map(|identity| identity.name.as_str()),
                    )
                    .text(&member.ratio)
                    .text(&member.overall_factor);
            }
        }
        let ids: Vec<_> = self
            .diagrams
            .iter()
            .map(|d| Digest::new().bytes(&d.id.0.to_le_bytes()).0)
            .collect();
        let distinct: BTreeSet<_> = self.diagrams.iter().map(|d| d.id).collect();
        assert_eq!(distinct.len(), self.diagrams.len(), "duplicate diagram ids");
        Summary {
            diagrams: self.diagrams.len(),
            topologies: self.report.topology_count,
            assignments: self.report.interaction_assignment_count,
            inverse_symmetry,
            weights,
            symmetry_factors: histogram.into_iter().collect(),
            content: Digest::multiset(content.collect()),
            numerators: Digest::multiset(numerators.collect()),
            named: named.0,
            ids: Digest::multiset(ids),
        }
    }

    /// The first difference between two runs within `scope`, coarse to fine.
    pub fn diff(&self, other: &Census, scope: Scope) -> Option<String> {
        if scope == Scope::Exact && self.report != other.report {
            return Some(format!(
                "reports differ: {:?} vs {:?}",
                self.report, other.report
            ));
        }
        let (left, right) = (self.by_key(), other.by_key());
        let only = |a: &BTreeMap<&CanonicalDiagramKey, &Entry>,
                    b: &BTreeMap<&CanonicalDiagramKey, &Entry>| {
            a.values()
                .filter(|d| !b.contains_key(d.full_key.as_ref()))
                .take(5)
                .map(|d| d.name.clone())
                .collect::<Vec<_>>()
        };
        if left.keys().ne(right.keys()) {
            return Some(format!(
                "diagram sets differ: {} vs {} diagrams; only left (first 5): {:?}; only right \
                 (first 5): {:?}",
                left.len(),
                right.len(),
                only(&left, &right),
                only(&right, &left)
            ));
        }
        if let Some((l, r)) = left
            .values()
            .zip(right.values())
            .find(|(l, r)| l.payload != r.payload)
        {
            let left_payload = l.payload.as_object().unwrap();
            let right_payload = r.payload.as_object().unwrap();
            let field = left_payload
                .keys()
                .chain(right_payload.keys())
                .find(|field| left_payload.get(*field) != right_payload.get(*field))
                .expect("different payloads have a differing field");
            return Some(format!(
                "full {field} payload of {:?} vs {:?} differs",
                l.name, r.name
            ));
        }
        type Field = fn(&Entry) -> String;
        let mut fields: Vec<(&str, Field)> = vec![
            ("symmetry factor", |d| d.symmetry_factor.to_string()),
            ("overall factor", |d| d.overall_factor.clone()),
            ("frame digest", |d| format!("{:016x}", d.frame)),
            ("numerator digest", |d| format!("{:016x}", d.numerator)),
        ];
        if scope == Scope::Exact {
            fields.push(("id", |d| d.id.to_string()));
        }
        for (field, value) in fields {
            if let Some((l, r)) = left
                .values()
                .zip(right.values())
                .find(|(l, r)| value(l) != value(r))
            {
                return Some(format!(
                    "{field} of {:?} vs {:?} differs: {} vs {}",
                    l.name,
                    r.name,
                    value(l),
                    value(r)
                ));
            }
        }
        // A retained group's basis may intentionally come from its comparison
        // master before cut merging. Compare actual output bases, not a fresh
        // canonical routing computed after that metadata changed.
        if let Some((l, r)) = left
            .values()
            .zip(right.values())
            .find(|(l, r)| l.basis != r.basis)
        {
            return Some(format!(
                "momentum basis of {:?} vs {:?} differs: {:?} vs {:?}",
                l.name, r.name, l.basis, r.basis
            ));
        }
        if scope == Scope::Content {
            return None;
        }
        if self.groups.len() != other.groups.len() {
            return Some(format!(
                "group counts differ: {} vs {}",
                self.groups.len(),
                other.groups.len()
            ));
        }
        let group_sizes = |census: &Census| {
            census
                .groups
                .iter()
                .map(|g| (census.diagrams[g.master].full_key.clone(), g.members.len()))
                .collect::<BTreeMap<_, _>>()
        };
        if group_sizes(self) != group_sizes(other) {
            return Some(format!(
                "group member counts differ: {:?} vs {:?}",
                group_sizes(self),
                group_sizes(other)
            ));
        }
        if let Some((l, r)) = left
            .values()
            .zip(right.values())
            .find(|(l, r)| l.name != r.name)
        {
            return Some(format!("names differ: {:?} vs {:?}", l.name, r.name));
        }
        if let Some(position) = (0..self.diagrams.len())
            .find(|&i| self.diagrams[i].full_key != other.diagrams[i].full_key)
        {
            return Some(format!(
                "output order differs first at position {position}: {:?} vs {:?}",
                self.diagrams[position].name, other.diagrams[position].name
            ));
        }
        if let Some((l, r)) = self.groups.iter().zip(&other.groups).find(|(l, r)| l != r) {
            return Some(format!("groups differ: {l:?} vs {r:?}"));
        }
        None
    }

    pub fn assert_same(&self, other: &Census, scope: Scope, context: &str) {
        if let Some(difference) = self.diff(other, scope) {
            panic!("{context}: {difference}");
        }
    }

    /// Only requested metadata may differ between Named and Minimal output.
    pub fn assert_minimal(&self, minimal: &Census, context: &str) {
        let mut unnamed = self.clone();
        for entry in &mut unnamed.diagrams {
            assert!(
                entry.name.take().is_some(),
                "{context}: Named diagram lacks a name"
            );
        }
        for group in &mut unnamed.groups {
            for member in &mut group.members {
                assert!(
                    member.identity.take().is_some(),
                    "{context}: Named source lacks identity"
                );
            }
        }
        unnamed.assert_same(
            minimal,
            Scope::Exact,
            &format!("{context}: Named vs Minimal"),
        );
    }

    /// Resolve every Named member by identity against the full ungrouped inventory,
    /// including discarded sources. Source indices belong to each request: grouping's
    /// canonical encounter order can differ from the ungrouped search order.
    pub fn assert_group_sources(&self, ungrouped: &Census, context: &str) {
        assert_eq!(ungrouped.groups.len(), ungrouped.diagrams.len());
        let source_indices: BTreeSet<_> = ungrouped
            .groups
            .iter()
            .map(|group| {
                assert_eq!(group.members.len(), 1, "{context}: ungrouped reference");
                group.members[0].source_index
            })
            .collect();
        assert_eq!(
            source_indices.len(),
            ungrouped.diagrams.len(),
            "{context}: repeated ungrouped source index"
        );
        ungrouped.by_key(); // Require unique full canonical keys as well as unique IDs.
        let sources = ungrouped
            .diagrams
            .iter()
            .map(|source| (source.id, source))
            .collect::<BTreeMap<_, _>>();
        assert_eq!(
            sources.len(),
            ungrouped.diagrams.len(),
            "{context}: repeated source ID"
        );
        let width = sources.len().saturating_sub(1).to_string().len();
        let mut seen_indices = BTreeSet::new();
        let mut seen_sources = BTreeSet::new();
        for group in &self.groups {
            let master = group.members.first().expect("nonempty group");
            let identity = master.identity.as_ref().expect("Named source identity");
            assert_eq!(
                self.diagrams[group.master].name.as_ref(),
                Some(&identity.name),
                "{context}: master source name"
            );
            for member in &group.members {
                assert!(
                    member.source_index < sources.len(),
                    "{context}: grouped source index out of bounds"
                );
                assert!(
                    seen_indices.insert(member.source_index),
                    "{context}: repeated grouped source index"
                );
                let identity = member.identity.as_ref().expect("Named source identity");
                let source = sources.get(&identity.id).expect("group source exists");
                assert!(
                    seen_sources.insert(source.full_key.as_ref()),
                    "{context}: repeated grouped source key"
                );
                // check() verifies these fixtures use default naming and no selectors.
                assert_eq!(
                    identity.name,
                    format!("D{:0width$}", member.source_index),
                    "{context}: source name and request-local index"
                );
                assert_eq!(
                    member.overall_factor, source.overall_factor,
                    "{context}: source factor"
                );
            }
        }
    }

    /// One line per diagram and per group, for comparing the runs of two builds with `diff`.
    pub fn dump(&self) -> String {
        let mut out = format!("{:?}\n", self.report);
        for d in &self.diagrams {
            writeln!(
                out,
                "{} key={:016x} id={} S={} frame={:016x} numerator={:016x} factor={} basis={}",
                d.name.as_deref().unwrap_or("null"),
                d.key,
                d.id,
                d.symmetry_factor,
                d.frame,
                d.numerator,
                d.overall_factor,
                serde_json::to_string(&d.basis).unwrap()
            )
            .unwrap();
        }
        for group in &self.groups {
            writeln!(
                out,
                "group {} {:?}",
                self.diagrams[group.master]
                    .name
                    .as_deref()
                    .unwrap_or("null"),
                group.members
            )
            .unwrap();
        }
        out
    }
}
