//! Differential oracle for diagram generation.
//!
//! A [`Census`] records what one generation run produced, in output order. [`Census::diff`]
//! compares two runs and names the first difference, coarse to fine as in the discrepancy
//! triage: report counts, the diagram set (by canonical certificate), symmetry factors (|Aut|),
//! overall factors with their signs, the finalized frames (with routing and cuts), numerators,
//! diagram ids, numerator groups, names and finally the output order. [`Summary`] condenses a
//! census into the numbers the oracle table records.
//!
//! Digests are FNV-1a hashes of fixed byte encodings (std's hashers may change between Rust
//! releases, which would invalidate the recorded rows). They leave out the [`DiagramId`]s. An id
//! hashes the model fingerprint with the diagram's structural key, which keeps particle
//! orientations as they are. The fingerprint hashes the printed model definition, whose term order
//! follows Symbolica's symbol registration order, so the ids of one diagram differ between
//! processes that register symbols in different orders (for instance after other models were parsed
//! in the same process). The certificate digest instead hashes `canonical_key()`, which normalizes
//! particle orientations and identifies the isomorphism class, and the frame digest covers
//! everything the structural key is computed from. Ids are compared only between the runs of one
//! process, so the rows do not guard a change of how an id is hashed or of the fingerprint.

use std::{
    collections::BTreeMap,
    fmt::{self, Write as _},
};

use feynkit_generator::{GenerationReport, GenerationResult};
use feynkit_graph::{DiagramId, FeynmanDiagram, expressions::evaluate_overall_factor};
use serde_json::Value;
use symbolica::{
    atom::{Atom, AtomCore},
    domains::rational::Rational,
};

/// One diagram of a run.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Entry {
    /// Digest of the canonical key, the certificate of the diagram's isomorphism class.
    pub key: u64,
    pub id: DiagramId,
    pub name: String,
    pub symmetry_factor: u64,
    pub overall_factor: String,
    /// The overall factor with its bookkeeping heads evaluated: `1 / |Aut|` with the fermion
    /// signs and the number of diagrams that symmetrization or grouping identified with this one.
    pub weight: Rational,
    /// Digest of the diagram in its finalized frame: vertices, edges, slots, loop momentum
    /// basis, cuts and threshold candidates, without numerators.
    pub frame: u64,
    /// Digest of the numerator, its prefactor, the projector and the local numerator fragments.
    pub numerator: u64,
}

/// A numerator group: the master's position and its members.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Group {
    pub master: usize,
    /// Per member: source name, numerator ratio to the master and overall factor.
    pub members: Vec<[String; 3]>,
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
    /// Order- and name-independent digest of certificates, symmetry and overall factors and
    /// frames.
    pub content: u64,
    /// Order-independent digest of the numerators by certificate.
    pub numerators: u64,
    /// Digest of names and certificates in output order and of the numerator groups.
    pub named: u64,
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
             numerators: 0x{:016x}, named: 0x{:016x} }}",
            self.diagrams,
            self.topologies,
            self.assignments,
            self.inverse_symmetry.numerator(),
            self.inverse_symmetry.denominator(),
            self.weights.numerator(),
            self.weights.denominator(),
            self.content,
            self.numerators,
            self.named
        )
    }
}

fn canonical(atom: &Atom) -> String {
    atom.to_canonical_string()
}

impl Entry {
    fn of(diagram: &FeynmanDiagram) -> Entry {
        let key = serde_json::to_value(diagram.canonical_key().unwrap()).unwrap();
        // The serialized diagram without its id, model fingerprint and name, split into the
        // numerators and the rest. The id is a `u128`, which a `Value` holds only once parsed
        // (as a float), before it is dropped.
        let mut serialized: Value =
            serde_json::from_slice(&serde_json::to_vec(diagram).unwrap()).unwrap();
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
        for field in ["id", "model", "name", "symmetry_factor", "overall_factor"] {
            fields.remove(field);
        }
        let weight = evaluate_overall_factor(diagram.overall_factor().as_view());
        Entry {
            key: Digest::new().json(&key).0,
            id: diagram.id(),
            name: diagram.name().to_owned(),
            symmetry_factor: diagram.symmetry_factor(),
            overall_factor: canonical(diagram.overall_factor()),
            weight: Rational::try_from(&weight).unwrap_or_else(|_| {
                panic!(
                    "{}: overall factor {weight} is not a number",
                    diagram.name()
                )
            }),
            frame: Digest::new().json(&serialized).0,
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
                    .map(|member| {
                        [
                            member.source_name.clone(),
                            canonical(&member.ratio),
                            canonical(&member.overall_factor),
                        ]
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

    fn by_key(&self) -> BTreeMap<u64, &Entry> {
        let by_key: BTreeMap<_, _> = self.diagrams.iter().map(|d| (d.key, d)).collect();
        assert_eq!(by_key.len(), self.diagrams.len(), "duplicate certificates");
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
            named = named.text(&d.name).number(d.key);
        }
        for group in &self.groups {
            named = named.number(group.master as u64);
            for member in &group.members {
                named = member.iter().fold(named, |named, text| named.text(text));
            }
        }
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
        let only = |a: &BTreeMap<u64, &Entry>, b: &BTreeMap<u64, &Entry>| {
            a.values()
                .filter(|d| !b.contains_key(&d.key))
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
                    "{field} of {} vs {} differs: {} vs {}",
                    l.name,
                    r.name,
                    value(l),
                    value(r)
                ));
            }
        }
        if scope == Scope::Content {
            return None;
        }
        let group_sizes = |census: &Census| {
            census
                .groups
                .iter()
                .map(|g| (census.diagrams[g.master].key, g.members.len()))
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
            return Some(format!("names differ: {} vs {}", l.name, r.name));
        }
        if let Some(position) =
            (0..self.diagrams.len()).find(|&i| self.diagrams[i].key != other.diagrams[i].key)
        {
            return Some(format!(
                "output order differs first at position {position}: {} vs {}",
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

    /// One line per diagram and per group, for comparing the runs of two builds with `diff`.
    pub fn dump(&self) -> String {
        let mut out = format!("{:?}\n", self.report);
        for d in &self.diagrams {
            writeln!(
                out,
                "{} key={:016x} S={} frame={:016x} numerator={:016x} factor={}",
                d.name, d.key, d.symmetry_factor, d.frame, d.numerator, d.overall_factor
            )
            .unwrap();
        }
        for group in &self.groups {
            writeln!(
                out,
                "group {} {:?}",
                self.diagrams[group.master].name, group.members
            )
            .unwrap();
        }
        out
    }
}
