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
    fs::File,
    io::{Read as _, Seek as _, SeekFrom, Write as _},
    sync::{Arc, Mutex},
};

use feynkit_generator::{BuiltDiagramSet, GenerationReport, SourceIdentity};
use feynkit_graph::{
    CanonicalDiagramKey, DiagramId, FeynmanDiagram, LoopMomentumBasis,
    expressions::evaluate_overall_factor,
};
use serde::{Serialize, Serializer};
use serde_json::Value;
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    domains::rational::Rational,
    symbol,
};

/// One diagram of a run.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Entry {
    /// Digest of the canonical key, the canonical form of the diagram's isomorphism class.
    pub key: u64,
    /// Full owner key pairs diagrams without relying on the recorded digest.
    full_key: Arc<CanonicalDiagramKey>,
    /// Every serialized field except the separately retained typed ID and optional name.
    /// File-backed bytes keep complete payloads out of the retained census metadata.
    payload: Payload,
    pub id: DiagramId,
    pub name: Option<String>,
    pub symmetry_factor: u128,
    pub overall_factor: String,
    /// Original tree for binding grouping provenance without rewriting an expression.
    factor: Atom,
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

/// One exact serialized payload in an anonymous file shared by census clones.
/// Construction owns the only writer; comparisons begin after that writer is dropped.
#[derive(Clone, Debug)]
struct Payload {
    store: Arc<Mutex<File>>,
    offset: u64,
    length: usize,
}

impl Payload {
    fn write(value: &Value, writer: &mut File, store: &Arc<Mutex<File>>) -> Self {
        let bytes = serde_json::to_vec(value).unwrap();
        let offset = writer.stream_position().unwrap();
        writer.write_all(&bytes).unwrap();
        Self {
            store: Arc::clone(store),
            offset,
            length: bytes.len(),
        }
    }

    fn read(&self) -> Value {
        let mut bytes = vec![0; self.length];
        let mut file = self.store.lock().unwrap();
        file.seek(SeekFrom::Start(self.offset)).unwrap();
        file.read_exact(&mut bytes).unwrap();
        drop(file);
        serde_json::from_slice(&bytes).unwrap()
    }
}

impl PartialEq for Payload {
    fn eq(&self, other: &Self) -> bool {
        self.read() == other.read()
    }
}

impl Eq for Payload {}

#[test]
fn payload_storage_preserves_content_and_clone_lifetime() {
    use std::panic::{AssertUnwindSafe, catch_unwind};

    let value =
        serde_json::json!({"factor": "f(3,x+y,\\\"z\\\")", "nested": [null, 1, {"n": "-2/3"}]});
    let mut writer = tempfile::tempfile().unwrap();
    let store = Arc::new(Mutex::new(writer.try_clone().unwrap()));
    let first = Payload::write(&value, &mut writer, &store);
    let second = Payload::write(&value, &mut writer, &store);
    let different = Payload::write(
        &serde_json::json!({"factor": "changed"}),
        &mut writer,
        &store,
    );
    drop((writer, store));
    assert_eq!(first, second);
    assert_ne!(first, different);

    let mut writer = tempfile::tempfile().unwrap();
    let store = Arc::new(Mutex::new(writer.try_clone().unwrap()));
    let independent = Payload::write(&value, &mut writer, &store);
    drop((writer, store));
    assert_eq!(first, independent);
    let retained = first.clone();
    drop((first, second, different));
    assert_eq!(retained.read(), value);
    assert_eq!(retained, retained.clone());

    // A truncated record must fail instead of becoming an empty/default payload.
    retained
        .store
        .lock()
        .unwrap()
        .set_len(retained.length as u64 - 1)
        .unwrap();
    assert!(catch_unwind(AssertUnwindSafe(|| retained.read())).is_err());
    let lifetime = Arc::downgrade(&retained.store);
    drop(retained);
    assert!(lifetime.upgrade().is_none());
}

#[test]
fn streamed_exact_dump_matches_original_diagram_values() {
    use feynkit_generator::{Grouping, Naming, NumeratorOptions, Run};

    // This reference serializes the original diagrams, never the new stored payloads.
    let original_dump = |built: &BuiltDiagramSet| {
        let diagrams: Vec<_> = built
            .iter()
            .map(|diagram| {
                let mut payload: Value =
                    serde_json::from_slice(&serde_json::to_vec(diagram).unwrap()).unwrap();
                payload.as_object_mut().unwrap().remove("id").unwrap();
                payload.as_object_mut().unwrap().remove("name").unwrap();
                serde_json::json!({
                    "key": diagram.canonical_key().unwrap(),
                    "payload": payload,
                    "id": diagram.id().unwrap().to_string(),
                    "name": diagram.name(),
                    "aut": diagram.symmetry_factor().to_string(),
                    "overall_factor": canonical(diagram.overall_factor()),
                    "basis": diagram.loop_momentum_basis(),
                })
            })
            .collect();
        let groups: Vec<_> = built.groups().iter().map(|group| serde_json::json!({
            "master": group.master,
            "members": group.members.iter().map(|member| serde_json::json!({
                "source_index": member.source_diagram,
                "id": member.source_identity.as_ref().map(|identity| identity.id.to_string()),
                "name": member.source_identity.as_ref().map(|identity| &identity.name),
                "ratio": canonical(&member.ratio),
                "overall_factor": canonical(&member.overall_factor),
            })).collect::<Vec<_>>()
        })).collect();
        serde_json::to_vec(&serde_json::json!({"diagrams": diagrams, "groups": groups})).unwrap()
    };
    let case = super::cases::Case::named("qed_eeaa_scalar_L1");
    let model = Arc::new(case.model());
    let run = Run::default().threads(1);
    let prepared = case.request(&model).generate(model, &run).unwrap();
    for naming in [Naming::prefix("D"), Naming::None] {
        let built = prepared
            .build(
                &NumeratorOptions::default()
                    .with_naming(naming)
                    .with_zero_color_filter(false),
                &run,
            )
            .unwrap();
        for grouping in [Grouping::None, case.grouping()] {
            let grouped = built.group(&grouping, &run).unwrap();
            let census = Census::of(&grouped);
            for (entry, diagram) in census.diagrams.iter().zip(grouped.iter()) {
                let mut expected: Value =
                    serde_json::from_slice(&serde_json::to_vec(diagram).unwrap()).unwrap();
                expected.as_object_mut().unwrap().remove("id").unwrap();
                expected.as_object_mut().unwrap().remove("name").unwrap();
                assert_eq!(entry.payload.read(), expected);
            }
            if !matches!(grouping, Grouping::None) {
                assert!(census.groups.iter().any(|group| group.members.len() > 1));
            }
            let mut actual = Vec::new();
            census.write_exact_dump(&mut actual).unwrap();
            assert_eq!(actual, original_dump(&grouped));
            let restricted = census.restricted(|_| true);
            drop(census);
            assert_eq!(restricted.diagrams[0].payload.read(), {
                let mut expected: Value = serde_json::from_slice(
                    &serde_json::to_vec(grouped.iter().next().unwrap()).unwrap(),
                )
                .unwrap();
                expected.as_object_mut().unwrap().remove("id").unwrap();
                expected.as_object_mut().unwrap().remove("name").unwrap();
                expected
            });
        }
    }
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

impl Group {
    fn identity_factor(
        &self,
        factor: &Atom,
        sources: &BTreeMap<usize, (DiagramId, &str, &Atom)>,
    ) -> BTreeMap<DiagramId, (Atom, Atom)> {
        assert!(self.members.len() > 1);
        let members = self
            .members
            .iter()
            .map(|m| (m.source_index, m))
            .collect::<BTreeMap<_, _>>();
        assert_eq!(members.len(), self.members.len(), "duplicate source index");
        let AtomView::Add(sum) = factor.as_view() else {
            panic!("group factor must be the sum of its source terms");
        };
        assert_eq!(
            sum.get_nargs(),
            members.len(),
            "missing or extra source term"
        );
        let mut bound = BTreeMap::new();
        for term in sum.iter() {
            let AtomView::Fun(term) = term else {
                panic!("group source term has an extra coefficient or wrapper");
            };
            assert_eq!(
                term.get_symbol(),
                symbol!("feynkit_generator_factor::NumeratorDependentGrouping")
            );
            assert_eq!(term.get_nargs(), 3, "group source term arity");
            let mut args = term.iter();
            let index = usize::try_from(args.next().unwrap()).expect("natural source index");
            let member = members.get(&index).expect("group source index exists");
            let identity = member.identity.as_ref().expect("Named source identity");
            let source = sources.get(&index).expect("full source index exists");
            assert_eq!(identity.id, source.0, "source index identity");
            assert_eq!(identity.name, source.1, "source index name");
            let ratio = args.next().unwrap().to_owned();
            let source_factor = args.next().unwrap().to_owned();
            assert_eq!(canonical(&ratio), member.ratio, "group source ratio");
            assert_eq!(
                canonical(&source_factor),
                member.overall_factor,
                "group source factor"
            );
            assert_eq!(&source_factor, source.2);
            assert!(
                bound.insert(identity.id, (ratio, source_factor)).is_none(),
                "duplicate source identity"
            );
        }
        assert_eq!(
            bound.len(),
            members.len(),
            "missing or duplicate source reference"
        );
        bound
    }
}

#[test]
fn class_factor_bindings_reject_invalid_provenance_or_weights() {
    use std::panic::catch_unwind;
    use symbolica::function;

    let a = function!(symbol!("generation_oracle::source_a"), 1);
    let b = function!(symbol!("generation_oracle::source_b"), 2);
    let sources = BTreeMap::from([
        (3, (DiagramId(10), "D3", &a)),
        (7, (DiagramId(11), "D7", &b)),
    ]);
    let shifted = BTreeMap::from([
        (9, (DiagramId(10), "D9", &a)),
        (4, (DiagramId(11), "D4", &b)),
    ]);
    let group = |indices: [usize; 2]| Group {
        master: 0,
        members: [(DiagramId(10), &a), (DiagramId(11), &b)]
            .into_iter()
            .zip(indices)
            .map(|((id, factor), source_index)| Member {
                source_index,
                identity: Some(SourceIdentity {
                    id,
                    name: format!("D{source_index}"),
                }),
                ratio: "1".into(),
                overall_factor: canonical(factor),
            })
            .collect(),
    };
    let term = |index: i64, ratio: i64, factor: &Atom| {
        function!(
            symbol!("feynkit_generator_factor::NumeratorDependentGrouping"),
            index,
            ratio,
            factor
        )
    };
    let original = term(3, 1, &a) + term(7, 1, &b);
    let reordered = term(9, 1, &a) + term(4, 1, &b);
    assert_ne!(canonical(&original), canonical(&reordered));
    assert_eq!(
        group([3, 7]).identity_factor(&original, &sources),
        group([9, 4]).identity_factor(&reordered, &shifted)
    );
    for invalid in [
        term(42, 1, &a) + term(7, 1, &b),
        term(3, 2, &a) + term(7, 1, &b),
        term(3, 1, &b) + term(7, 1, &b),
        Atom::num(2) * term(3, 1, &a) + term(7, 1, &b),
        term(3, 1, &a),
        original.clone() + term(8, 1, &a),
        term(3, 1, &a) + term(3, 2, &b),
        function!(symbol!("generation_oracle::wrapped"), &original),
    ] {
        assert!(catch_unwind(|| group([3, 7]).identity_factor(&invalid, &sources)).is_err());
    }
    assert!(catch_unwind(|| group([3, 3]).identity_factor(&original, &sources)).is_err());
    let missing = BTreeMap::from([(3, (DiagramId(10), "D3", &a))]);
    assert!(catch_unwind(|| group([3, 7]).identity_factor(&original, &missing)).is_err());
    let changed = BTreeMap::from([
        (3, (DiagramId(10), "D3", &b)),
        (7, (DiagramId(11), "D7", &b)),
    ]);
    assert!(catch_unwind(|| group([3, 7]).identity_factor(&original, &changed)).is_err());
    // Equal ratios and factors must not hide swapped source identities.
    let mut equal = group([3, 7]);
    equal.members[1].overall_factor = canonical(&a);
    let equal_sources = BTreeMap::from([
        (3, (DiagramId(10), "D3", &a)),
        (7, (DiagramId(11), "D7", &a)),
    ]);
    let equal_factor = term(3, 1, &a) + term(7, 1, &a);
    assert_eq!(
        equal.identity_factor(&equal_factor, &equal_sources).len(),
        2
    );
    equal.members[0].identity.as_mut().unwrap().id = DiagramId(11);
    equal.members[1].identity.as_mut().unwrap().id = DiagramId(10);
    assert!(catch_unwind(|| equal.identity_factor(&equal_factor, &equal_sources)).is_err());
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
    fn of(diagram: &FeynmanDiagram, writer: &mut File, store: &Arc<Mutex<File>>) -> Entry {
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
        let payload = Payload::write(&serialized, writer, store);
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
            factor: diagram.overall_factor().clone(),
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
    pub fn of(result: &BuiltDiagramSet) -> Census {
        // Heavy runs put these anonymous files on the gate's disk-backed filesystem.
        // No extra output name survives, and the last cloned payload closes the file.
        let mut writer = if let Ok(directory) = std::env::var("FEYNKIT_ORACLE_DUMP") {
            std::fs::create_dir_all(&directory).unwrap();
            tempfile::tempfile_in(directory).unwrap()
        } else {
            tempfile::tempfile().unwrap()
        };
        let store = Arc::new(Mutex::new(writer.try_clone().unwrap()));
        let groups = result
            .groups()
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
            report: result.report().clone(),
            diagrams: result
                .iter()
                .map(|diagram| Entry::of(diagram, &mut writer, &store))
                .collect(),
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
            *histogram
                .entry(
                    u64::try_from(d.symmetry_factor)
                        .expect("recorded oracle symmetry factor fits u64"),
                )
                .or_default() += 1;
        }
        let content = self.diagrams.iter().map(|d| {
            Digest::new()
                .number(d.key)
                .number(
                    u64::try_from(d.symmetry_factor)
                        .expect("recorded oracle symmetry factor fits u64"),
                )
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
            topologies: usize::try_from(&self.report.topology_count)
                .expect("oracle topology count fits usize"),
            assignments: usize::try_from(&self.report.interaction_assignment_count)
                .expect("oracle assignment count fits usize"),
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
            let left_value = l.payload.read();
            let right_value = r.payload.read();
            let left_payload = left_value.as_object().unwrap();
            let right_payload = right_value.as_object().unwrap();
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
    /// including discarded sources. Grouped output positions may differ, but each
    /// source index refers to the original request's ungrouped inventory.
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

    /// Compare classes after joining by the full key. Search order and display names may
    /// change, but the retained frame, original rules, factors and group members may not.
    pub fn assert_class_equivalent(
        &self,
        other: &Census,
        sources: (&Census, &Census),
        context: &str,
    ) {
        self.assert_group_sources(sources.0, context);
        other.assert_group_sources(sources.1, context);
        // Only a grouping head's first argument is request-local provenance. Bind it
        // to the complete source identity; retain the ratio and nested factor trees.
        let factors = |census: &Census, sources: &Census| {
            let source_indices = sources
                .groups
                .iter()
                .map(|group| {
                    assert_eq!(group.members.len(), 1, "{context}: ungrouped source");
                    let member = &group.members[0];
                    let identity = member.identity.as_ref().expect("Named source identity");
                    let entry = &sources.diagrams[group.master];
                    assert_eq!(identity.id, entry.id, "{context}: source identity");
                    assert_eq!(
                        Some(&identity.name),
                        entry.name.as_ref(),
                        "{context}: source name"
                    );
                    (
                        member.source_index,
                        (entry.id, identity.name.as_str(), &entry.factor),
                    )
                })
                .collect::<BTreeMap<_, _>>();
            assert_eq!(
                source_indices.len(),
                sources.groups.len(),
                "{context}: unique source indices"
            );
            for entry in &census.diagrams {
                let serialized = canonical(&entry.factor);
                assert_eq!(
                    entry.overall_factor, serialized,
                    "{context}: factor serialization"
                );
                assert_eq!(
                    entry.payload.read()["overall_factor"].as_str(),
                    Some(serialized.as_str()),
                    "{context}: payload factor serialization"
                );
            }
            census
                .groups
                .iter()
                .filter(|group| group.members.len() > 1)
                .map(|group| {
                    let entry = &census.diagrams[group.master];
                    (
                        entry.id,
                        group.identity_factor(&entry.factor, &source_indices),
                    )
                })
                .collect::<BTreeMap<_, _>>()
        };
        let left_factors = factors(self, sources.0);
        let right_factors = factors(other, sources.1);
        assert_eq!(
            left_factors, right_factors,
            "{context}: source-bound group factors"
        );
        let left = self.by_key();
        let right = other.by_key();
        assert!(left.keys().eq(right.keys()), "{context}: diagram keys");
        for (key, entry) in &left {
            let other = right[key];
            let value = entry.payload.read();
            let other_value = other.payload.read();
            let payload = value.as_object().unwrap();
            let other_payload = other_value.as_object().unwrap();
            assert!(
                payload.keys().eq(other_payload.keys()),
                "{context}: payload fields"
            );
            for (field, value) in payload {
                if field != "overall_factor" || !left_factors.contains_key(&entry.id) {
                    assert_eq!(
                        value, &other_payload[field],
                        "{context}: full {field} payload"
                    );
                }
            }
            assert_eq!(
                (
                    entry.id,
                    entry.symmetry_factor,
                    entry.frame,
                    entry.numerator,
                    &entry.basis
                ),
                (
                    other.id,
                    other.symmetry_factor,
                    other.frame,
                    other.numerator,
                    &other.basis
                ),
                "{context}: identity, symmetry, frame, numerator or basis"
            );
        }
        let groups = |census: &Census| {
            census
                .groups
                .iter()
                .map(|group| {
                    let members = group
                        .members
                        .iter()
                        .map(|member| {
                            (
                                member.identity.as_ref().expect("Named identity").id,
                                (member.ratio.clone(), member.overall_factor.clone()),
                            )
                        })
                        .collect::<BTreeMap<_, _>>();
                    assert_eq!(
                        members.len(),
                        group.members.len(),
                        "{context}: duplicate member"
                    );
                    (census.diagrams[group.master].id, members)
                })
                .collect::<BTreeMap<_, _>>()
        };
        let left_groups = groups(self);
        let right_groups = groups(other);
        assert_eq!(
            left_groups.len(),
            self.groups.len(),
            "{context}: duplicate master"
        );
        assert_eq!(
            right_groups.len(),
            other.groups.len(),
            "{context}: duplicate master"
        );
        assert_eq!(
            left_groups, right_groups,
            "{context}: member identity, ratio or factor"
        );
    }

    /// Retain exact physical fields for independent replay, in addition to the compact census.
    pub fn write_exact_dump(&self, writer: impl std::io::Write) -> serde_json::Result<()> {
        // Keep the old Value serializer's sorted field order and literal values,
        // while materializing only one stored diagram payload at a time.
        #[derive(Serialize)]
        struct Diagram<'a> {
            aut: String,
            basis: Value,
            id: String,
            key: Value,
            name: &'a Option<String>,
            overall_factor: &'a str,
            payload: Value,
        }
        struct Diagrams<'a>(&'a [Entry]);
        impl Serialize for Diagrams<'_> {
            fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
                serializer.collect_seq(self.0.iter().map(|entry| Diagram {
                    aut: entry.symmetry_factor.to_string(),
                    basis: serde_json::to_value(&entry.basis).unwrap(),
                    id: entry.id.to_string(),
                    key: serde_json::to_value(entry.full_key.as_ref()).unwrap(),
                    name: &entry.name,
                    overall_factor: &entry.overall_factor,
                    payload: entry.payload.read(),
                }))
            }
        }
        #[derive(Serialize)]
        struct Source<'a> {
            id: Option<String>,
            name: Option<&'a str>,
            overall_factor: &'a str,
            ratio: &'a str,
            source_index: usize,
        }
        struct Members<'a>(&'a [Member]);
        impl Serialize for Members<'_> {
            fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
                serializer.collect_seq(self.0.iter().map(|member| {
                    Source {
                        id: member
                            .identity
                            .as_ref()
                            .map(|identity| identity.id.to_string()),
                        name: member
                            .identity
                            .as_ref()
                            .map(|identity| identity.name.as_str()),
                        overall_factor: &member.overall_factor,
                        ratio: &member.ratio,
                        source_index: member.source_index,
                    }
                }))
            }
        }
        #[derive(Serialize)]
        struct GroupDump<'a> {
            master: usize,
            members: Members<'a>,
        }
        struct Groups<'a>(&'a [Group]);
        impl Serialize for Groups<'_> {
            fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
                serializer.collect_seq(self.0.iter().map(|group| GroupDump {
                    master: group.master,
                    members: Members(&group.members),
                }))
            }
        }
        #[derive(Serialize)]
        struct Dump<'a> {
            diagrams: Diagrams<'a>,
            groups: Groups<'a>,
        }
        serde_json::to_writer(
            writer,
            &Dump {
                diagrams: Diagrams(&self.diagrams),
                groups: Groups(&self.groups),
            },
        )
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
