//! Open Lorentz slots (ports) on R*-remainder-shaped values: the open value,
//! closed afterwards with every independent tensor built from the metric and
//! the external vectors, equals the CONTRACTED route (the same value with that
//! tensor as its projector, through the same entry), order by order in eps,
//! exactly.
//!
//! The comparison runs once per value with the GENERIC closing tensor
//! `T = SUM_S z_S T_S` (`z_S` independent caller symbols, allowed in a
//! projector): both sides are linear in the `z_S`, so the identity of the
//! polynomials is the identity for every structure `T_S` separately.
//!
//! * VACUUM-SHAPED (the open IR hard part after its Taylor step): the real
//!   `kind hard` sector records of the beta's L3 run (tests/fixtures/wh_ports,
//!   1 / 3 / 4 open slots, three or four independent crown vectors, one or two
//!   loops, one massive closing line), read as vacuum values (the crown legs
//!   dropped: their momenta become constant external vectors, the edge
//!   momenta inside the factors are written out in loop and crown momenta,
//!   the denominators are those at zero crown momenta), through
//!   `evaluate_factorized_vacuum`.
//! * PROPAGATOR-TYPE (kind direct, `evaluate_factorized_ports`): background
//!   field graphs of two, three and four loops from the beta's prepared DOT
//!   files (`VAKINT_TEST_FACTORIZED_DOT_ROOT`, default
//!   /common/ben/beta/results/prepared; SKIP without them) with their two
//!   ports, one port (one leg's slot closed by its own momentum), three and
//!   four ports (extra edge momenta on boundary slots), each with the record's
//!   `p(1)` in both signs; and a hand-made two-loop value with a SOFT leg
//!   carrying a port.
//! * Refusals, the text round trip, the caller's notation.
#![cfg(feature = "feynkit-ingress")]
#[path = "support/open_indices.rs"]
mod open_indices;

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use feynkit_graph::FeynmanDiagram;
use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};
use vakint::p_integral::factorized::{
    Endpoint, FactorizedEdge, FactorizedVertex, FromFeynkitOptions, LineKind, Prefactor, Projector,
    ProjectorChoice, evaluate_factorized, evaluate_factorized_ports, evaluate_factorized_vacuum,
};
use vakint::p_integral::{FactorizedError, FactorizedNumerator, FactorizedParts};
use vakint::refpath;
use vakint::symbols::S;
use vakint::{
    EvaluationMethod, EvaluationOrder, LoopNormalizationFactor, RustRedEvaluationOptions,
    RustRedMasterMode, Vakint, VakintSettings,
};

const MINK: &str = "spenso::mink";
const METRIC: &str = "spenso::g";
const MOMENTUM: &str = "gammalooprs::Q";
const HEDGE: &str = "gammalooprs::hedge";
const DIM: &str = "gammalooprs::dim";

/// Register the contract's symbols with their attributes (as the DOT files
/// and the text form write them) before any plain lookup by name.
fn init() {
    Vakint::initialize_vakint_symbols();
    let _ = symbolica::try_parse!(
        "spenso::{spenso::representation,spenso::self_dual}::mink(gammalooprs::{}::dim,gammalooprs::{spenso::index}::hedge(1,1))*spenso::{spenso::representation,spenso::self_dual}::bis(4,gammalooprs::{spenso::index}::hedge(2,1))*spenso::{symmetric,linear,real,spenso::tensor}::g(x,y)*gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript,spenso::tensor-label:q}::Q(1,z)"
    )
    .unwrap();
}

fn sym(name: &str) -> Symbol {
    symbolica::symbol!(name)
}

fn engine() -> Vakint {
    init();
    let mut vakint = Vakint::new().unwrap();
    if let Some(library) = std::env::var_os("VAKINT_TEST_P_LIBRARY_1_4") {
        vakint
            .load_p_library_directory(Path::new(&library), false)
            .unwrap();
    }
    vakint
}

fn epsilon(settings: &VakintSettings) -> Symbol {
    vakint::vakint_symbol!(settings.epsilon_symbol.as_str())
}

/// The index `a` of a slot `mink(D, a)`.
fn index_of_slot(slot: &Atom) -> Atom {
    match slot.as_view() {
        AtomView::Fun(f) if f.get_nargs() == 2 => f.iter().nth(1).unwrap().to_owned(),
        _ => panic!("not a slot: {slot}"),
    }
}

/// The value with an empty schema-v3 direct record: a value WITHOUT any
/// new-feature marker takes exactly the base's rules in the old entries, so a
/// test that exercises a new feature there (a soft first leg, a partial or
/// gamma projector, open slots through the old entries) marks its value.
fn marked(value: &FactorizedNumerator) -> FactorizedNumerator {
    FactorizedNumerator {
        record: Some(vakint::p_integral::factorized::DirectRecord {
            original_cycles: None,
            depth: None,
            weights: Vec::new(),
        }),
        ..value.clone()
    }
}

fn hedge_slot(dimension: Symbol, h: i64) -> Atom {
    symbolica::function!(
        sym(MINK),
        Atom::var(dimension),
        symbolica::function!(sym(HEDGE), h, 1)
    )
}

fn q(edge: i64, slot: &Atom) -> Atom {
    symbolica::function!(sym(MOMENTUM), edge, slot)
}

/// Every independent closing structure over `slots` with vectors `vectors`:
/// a set of slots carried by vectors (each by one of them) and a metric
/// pairing of the others. Returned as (caller tensor, vakint tensor) builders
/// through the two closures.
fn structures<T: Clone>(
    slots: &[T],
    vectors: usize,
    vector: &dyn Fn(usize, &T) -> Atom,
    metric: &dyn Fn(&T, &T) -> Atom,
) -> Vec<Atom> {
    fn pairings<T: Clone>(points: &[T]) -> Vec<Vec<(T, T)>> {
        if points.is_empty() {
            return vec![vec![]];
        }
        let mut out = Vec::new();
        for i in 1..points.len() {
            let rest = points[1..i]
                .iter()
                .chain(&points[i + 1..])
                .cloned()
                .collect::<Vec<_>>();
            for mut pairing in pairings(&rest) {
                pairing.push((points[0].clone(), points[i].clone()));
                out.push(pairing);
            }
        }
        out
    }
    let n = slots.len();
    let mut out = Vec::new();
    for mask in 0..(1u32 << n) {
        let carried = (0..n).filter(|i| mask & (1 << i) != 0).collect::<Vec<_>>();
        let rest = (0..n)
            .filter(|i| mask & (1 << i) == 0)
            .map(|i| slots[i].clone())
            .collect::<Vec<_>>();
        if rest.len() % 2 == 1 {
            continue;
        }
        // every assignment of the carried slots to the vectors
        let assignments = vectors.pow(carried.len() as u32);
        for assignment in 0..assignments {
            let mut a = assignment;
            let mut carried_part = Atom::num(1);
            for &i in &carried {
                carried_part *= vector(a % vectors, &slots[i]);
                a /= vectors;
            }
            for pairing in pairings(&rest) {
                let mut tensor = carried_part.clone();
                for (x, y) in &pairing {
                    tensor *= metric(x, y);
                }
                out.push(tensor);
            }
        }
    }
    out
}

/// `SUM_S z_S T_S`: the generic closing tensor, in the caller's notation and
/// in vakint's, with the same symbols `z_S`.
fn generic_closing(
    slots: &[(Atom, Atom)],
    vectors: usize,
    caller_vector: &dyn Fn(usize, &Atom) -> Atom,
    vakint_vector: &dyn Fn(usize, &Atom) -> Atom,
) -> (Atom, Atom, usize) {
    let caller_slots = slots.iter().map(|(s, _)| s.clone()).collect::<Vec<_>>();
    let vakint_slots = slots.iter().map(|(_, v)| v.clone()).collect::<Vec<_>>();
    let caller = structures(&caller_slots, vectors, caller_vector, &|a, b| {
        symbolica::function!(sym(METRIC), a, b)
    });
    let vakint = structures(&vakint_slots, vectors, vakint_vector, &|a, b| {
        symbolica::function!(S.g, a, b)
    });
    assert_eq!(caller.len(), vakint.len());
    let mut c = Atom::Zero;
    let mut v = Atom::Zero;
    for (i, (ct, vt)) in caller.iter().zip(&vakint).enumerate() {
        let z = Atom::var(sym(&format!("wh_ports_test::z{i}")));
        c += &z * ct;
        v += &z * vt;
    }
    (c, v, caller.len())
}

/// Assert `open` closed with `closing` (vakint notation) equals `contracted`.
fn assert_closes(open: &Atom, closing: &Atom, contracted: &Atom, eps: Symbol, what: &str) {
    let closed = open_indices::close_indices(&(open * closing), eps);
    let last = open_indices::last_order(contracted, &closed, eps);
    open_indices::assert_series_through(&closed, contracted, last, what, eps);
    assert!(
        !open_indices::epsilon_coefficients(contracted, eps).is_empty(),
        "{what}: the contracted reference is zero (a vacuous identity)"
    );
}

// ---------------------------------------------------------------------------
// Vacuum-shaped values: the beta's `kind hard` sector records

fn fixtures() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/wh_ports")
}

/// A `kind hard` sector record (schema v2 of the R* caller) as a vacuum value
/// of the v1 input type: crown legs dropped (their `Q(x, ..)` become constant
/// external vectors), internal momenta at zero crown momenta, every internal
/// edge momentum inside a factor written out as `SUM eta Q(b_i) + SUM xi
/// Q(x_j)` (a linear substitution). The cut / routing / weights lines are
/// provenance and are dropped. The SHAPE of the open hard part (its slots, its
/// external vectors, its closing line, its numerator) is what these tests
/// need; whether the producer's form-B factors mean full or zero-crown
/// momenta does not matter for an identity of two routes on the same value.
fn hard_as_vacuum(text: &str) -> FactorizedNumerator {
    let drop = [
        "hard",
        "softs",
        "soft",
        "closure",
        "closure_momentum",
        "routing",
        "route",
        "cut_topology",
        "weights",
        "weight",
        "original_cycles",
        // v3 lines of the closure kind (the Taylor step is not applied here)
        "open",
        "open_slot",
        "taylor",
        "depth",
    ];
    let mut loops = 0;
    let mut externals = 0;
    let mut legs = BTreeMap::<usize, usize>::new(); // leg edge -> j
    let mut vectors = BTreeMap::<usize, Vec<i64>>::new(); // internal edge -> full vector
    for line in text.lines() {
        let fields = line.split('\t').collect::<Vec<_>>();
        match fields[0] {
            "momenta" => {
                loops = fields[1].trim_start_matches("loops=").parse().unwrap();
                externals = fields[2].trim_start_matches("externals=").parse().unwrap();
            }
            "external" => {
                let j: usize = fields[1].parse().unwrap();
                let x: usize = fields[2].trim_start_matches("edge=").parse().unwrap();
                legs.insert(x, j);
            }
            "edge" if !(fields[2] == "source=ext" || fields[3] == "target=ext") => {
                let id: usize = fields[1].parse().unwrap();
                let momentum = fields
                    .iter()
                    .find_map(|f| f.strip_prefix("momentum="))
                    .unwrap()
                    .split(',')
                    .map(|c| c.parse::<i64>().unwrap())
                    .collect::<Vec<_>>();
                vectors.insert(id, momentum);
            }
            _ => {}
        }
    }
    let width = loops + externals;
    let mut out = Vec::new();
    for line in text.lines() {
        let fields = line.split('\t').collect::<Vec<_>>();
        let key = fields[0];
        if drop.contains(&key) || key == "lines" || key == "end" {
            continue;
        }
        match key {
            "schema" => out.push("schema\tpint6.factorized-numerator.v1".to_string()),
            "kind" => out.push("kind\tdirect".to_string()),
            // the crown legs stay, at zero momentum: no edge carries p(j)
            "external" => out.push(format!("external\t{}\tedge=none", fields[1])),
            "edge" => {
                let id: usize = fields[1].parse().unwrap();
                let v = match vectors.get(&id) {
                    Some(v) => v[..loops]
                        .iter()
                        .copied()
                        .chain(std::iter::repeat_n(0, externals))
                        .collect::<Vec<_>>(),
                    None => vec![0; width],
                };
                let fields = fields
                    .iter()
                    .map(|f| {
                        if f.starts_with("momentum=") {
                            format!(
                                "momentum={}",
                                v.iter()
                                    .map(|c| c.to_string())
                                    .collect::<Vec<_>>()
                                    .join(",")
                            )
                        } else {
                            f.to_string()
                        }
                    })
                    .collect::<Vec<_>>();
                out.push(fields.join("\t"));
            }
            _ => out.push(line.to_string()),
        }
    }
    out.push(format!("lines\t{}", out.len() + 1));
    out.push("end".into());
    let mut v1 = out.join("\n");
    v1.push('\n');
    init();
    let mut value = FactorizedNumerator::parse_text(&v1).unwrap_or_else(|e| panic!("{e}\n{v1}"));
    // open fermion lines, gamma / partial projectors: a marked value (the old
    // entries keep the base's rules for unmarked ones)
    value.record = marked(&value).record;
    // the loop edges b_i; the crown vectors become constant vectors 10000 + x
    let loop_edges = value
        .loop_edges
        .iter()
        .map(|b| b.expect("every R* loop has an edge") as i64)
        .collect::<Vec<_>>();
    let constant_of = legs
        .iter()
        .map(|(x, j)| (*j, 10000 + *x as i64))
        .collect::<BTreeMap<_, _>>();
    let written_out = |factor: &Atom| {
        factor.replace_map(|part, _, out| {
            if let AtomView::Fun(f) = part
                && f.get_symbol().get_name() == MOMENTUM
                && f.get_nargs() == 2
            {
                let mut args = f.iter();
                let edge = i64::try_from(args.next().unwrap()).unwrap() as usize;
                let slot = args.next().unwrap().to_owned();
                if legs.contains_key(&edge) {
                    **out = q(10000 + edge as i64, &slot);
                    return;
                }
                let Some(vector) = vectors.get(&edge) else {
                    return;
                };
                if vector[loops..].iter().all(|c| *c == 0) {
                    return; // a pure loop combination: the edge's own momentum
                }
                let mut sum = Atom::Zero;
                for (i, c) in vector[..loops].iter().enumerate() {
                    if *c != 0 {
                        sum += Atom::num(*c) * q(loop_edges[i], &slot);
                    }
                }
                for (j, c) in vector[loops..].iter().enumerate() {
                    if *c != 0 {
                        sum += Atom::num(*c) * q(constant_of[&(j + 1)], &slot);
                    }
                }
                **out = sum;
            }
        })
    };
    for edge in &mut value.edges {
        edge.factor = written_out(&edge.factor);
    }
    for vertex in &mut value.vertices {
        vertex.factor = written_out(&vertex.factor);
    }
    if let Some(projector) = &mut value.projector {
        projector.tensor = written_out(&projector.tensor);
    }
    // the carried graph's weights for the v1 validator (the value does not use
    // the prefactor fields: the prefactor is passed explicitly)
    let closed = value
        .closed_loops()
        .unwrap_or_else(|e| panic!("{}: {e}", value.label));
    value.prefactor.loop_sign = if closed.len().is_multiple_of(2) {
        1
    } else {
        -1
    };
    value.prefactor.nf_power = closed
        .iter()
        .filter(|l| l.line == LineKind::Fermion)
        .count() as u32;
    value
        .validate()
        .unwrap_or_else(|e| panic!("{}: {e}", value.label));
    value
}

fn vacuum_settings() -> VakintSettings {
    VakintSettings {
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions {
                master_mode: RustRedMasterMode::Symbolic,
            },
        )]),
        form_exe_path: "/nonexistent/wh-ports/form".into(),
        integral_normalization_factor: LoopNormalizationFactor::MSbar,
        number_of_terms_in_epsilon_expansion: 3,
        use_dot_product_notation: true,
        n_cores: 1,
        ..Default::default()
    }
}

/// The open hard part (vacuum entry, nothing closes the slots) closed with the
/// generic tensor in the constant crown vectors and the metric equals the
/// contracted route; and the open value is in the caller's notation.
#[test]
fn open_hard_part_slots_close_to_the_contracted_route() {
    let vakint = engine();
    let settings = vacuum_settings();
    let eps = epsilon(&settings);
    let mut checked = 0;
    let mut zero = 0;
    for name in [
        "hard-L3-006078.txt", // 1 loop, 3 crown vectors, 3 slots, Taylor sector 1
        "hard-L3-008385.txt", // 2 loops, 3 crown vectors, 3 slots, Taylor sector 4
        "hard-L3-012034.txt", // 2 loops, 3 crown vectors, 4 slots
        "hard-L3-000993.txt", // 1 loop, 4 crown vectors, 4 slots
        "k3-hs.input.txt", // the beta's k3-hs record (L3 #36), its shape: 1 loop, 3 vectors, 3 slots
        "k3-dirac.input.txt", // k3-dirac (L4 #4966): bispinor pair closed by gamma, fresh nu open
    ] {
        let text = std::fs::read_to_string(fixtures().join(name)).unwrap();
        let value = hard_as_vacuum(&text);
        let constants = value.constant_vectors();
        let first_constant = value.externals(); // vakint's p(j) of the constants come after the legs
        let open = evaluate_factorized_vacuum(&vakint, &settings, &value, &Atom::num(1))
            .unwrap_or_else(|e| panic!("{name}: open: {e}"));
        // caller notation: nothing of vakint's g / p / k survives
        let mut leftover = false;
        open.value.as_view().visitor(&mut |part| {
            if let AtomView::Fun(f) = part
                && [S.g, S.p, S.k].contains(&f.get_symbol())
            {
                leftover = true;
            }
            !leftover
        });
        assert!(!leftover, "{name}: vakint notation in {}", open.value);
        // open slots: (caller slot, vakint raw index 9001 + position)
        let slots = value
            .open_slots()
            .unwrap()
            .iter()
            .enumerate()
            .map(|(i, port)| {
                (
                    port.clone(),
                    Atom::num(vakint::p_integral::factorized::vacuum::PORT + i as i64),
                )
            })
            .collect::<Vec<_>>();
        let (caller, closing, count) = generic_closing(
            &slots,
            constants.len(),
            &|j, slot| q(constants[j] as i64, slot),
            &|j, index| symbolica::function!(S.p, (first_constant + j + 1) as i64, index),
        );
        let mut contracted_value = value.clone();
        contracted_value.projector = Some(match &value.projector {
            Some(projector) => Projector {
                tensor: &projector.tensor * caller,
                p2_power: projector.p2_power,
            },
            None => Projector {
                tensor: caller,
                p2_power: 0,
            },
        });
        let contracted =
            evaluate_factorized_vacuum(&vakint, &settings, &contracted_value, &Atom::num(1))
                .unwrap_or_else(|e| panic!("{name}: contracted: {e}"));
        if contracted.vakint_value.expand().is_zero() {
            // both routes must then give exactly zero
            let closed = open_indices::close_indices(&(&open.vakint_value * &closing), eps);
            assert!(
                closed.expand().is_zero(),
                "{name}: contracted 0, open closed {closed}"
            );
            zero += 1;
        } else {
            assert_closes(
                &open.vakint_value,
                &closing,
                &contracted.vakint_value,
                eps,
                name,
            );
        }
        eprintln!(
            "WH_PORTS_HARD fixture={name} loops={} vectors={} slots={} structures={count} open_s={:.3} contracted_s={:.3} identical=1 zero={}",
            value.loops(),
            constants.len(),
            slots.len(),
            open.contract_seconds + open.evaluate_seconds,
            contracted.contract_seconds + contracted.evaluate_seconds,
            contracted.vakint_value.expand().is_zero(),
        );
        checked += 1;
    }
    assert_eq!(checked, 6);
    assert!(
        checked - zero >= 4,
        "{zero} of {checked} fixtures are exactly zero"
    );
}

// ---------------------------------------------------------------------------
// Propagator-type values (kind direct)

fn dot_root() -> Option<PathBuf> {
    let root = std::env::var("VAKINT_TEST_FACTORIZED_DOT_ROOT")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("/common/ben/beta/results/prepared"));
    if root.join("L2").join("model.json").is_file() {
        Some(root)
    } else {
        eprintln!("SKIP: no DOT data under {}", root.display());
        None
    }
}

/// The first colour branch of a background-field graph, ports open.
fn graph_value(root: &Path, level: &str, number: usize) -> FactorizedNumerator {
    let directory = root.join(level);
    let model = refpath::load_model(&directory).unwrap();
    let path = directory.join(format!("diagram-{number:06}.dot"));
    let text = std::fs::read_to_string(&path).unwrap();
    let diagram = FeynmanDiagram::from_dot(Arc::clone(&model), &text).unwrap();
    let values = FactorizedNumerator::from_feynkit(
        &diagram,
        &model,
        &FromFeynkitOptions {
            projector: ProjectorChoice::None,
            ..FromFeynkitOptions::default()
        },
    )
    .unwrap();
    let mut value = values.into_iter().next().expect("a non-zero colour branch");
    value.label = format!("{level}/{number}");
    value.prefactor = Prefactor {
        loop_sign: value.prefactor.loop_sign,
        nf_power: value.prefactor.nf_power,
        ..Prefactor::default()
    };
    assert_eq!(value.ports.len(), 2, "{}", value.label);
    value
}

/// The value with its record `p(1)` reversed (every external component of the
/// momentum vectors negated): the same integral, written in `p' = -p(1)`.
fn flipped(value: &FactorizedNumerator) -> FactorizedNumerator {
    let mut value = value.clone();
    let loops = value.loops();
    for edge in &mut value.edges {
        for c in &mut edge.momentum[loops..] {
            *c = -*c;
        }
    }
    value.external_edges = vec![None; value.externals()];
    value.label.push_str(" flipped");
    value.validate().unwrap();
    value
}

/// One port fewer: the slot of the external leg of smallest id closed by the
/// leg's own momentum (a factor on the leg).
fn one_port(value: &FactorizedNumerator) -> FactorizedNumerator {
    let mut value = value.clone();
    let first = value
        .edges
        .iter()
        .filter(|e| e.is_external())
        .map(|e| e.id)
        .min()
        .unwrap();
    let leg = value.edges.iter().position(|e| e.id == first).unwrap();
    let hedge = value.edges[leg]
        .source_hedge
        .or(value.edges[leg].target_hedge)
        .unwrap();
    let slot = hedge_slot(value.dimension, hedge as i64);
    let position = value
        .ports
        .iter()
        .position(|p| *p == slot)
        .unwrap_or_else(|| panic!("{}: the leg slot is not a port", value.label));
    let _ = value.ports.remove(position);
    value.edges[leg].factor = &value.edges[leg].factor * q(first as i64, &slot);
    value.label.push_str(" one-port");
    value.validate().unwrap();
    value
}

/// `extra` more ports: internal edge momenta on slots of absent (boundary)
/// half-edges `9000 + i`, on the internal edges carrying `p(1)`.
fn more_ports(value: &FactorizedNumerator, extra: usize) -> FactorizedNumerator {
    let mut value = value.clone();
    let loops = value.loops();
    let carriers = value
        .edges
        .iter()
        .enumerate()
        .filter(|(_, e)| !e.is_external() && e.momentum[loops] != 0)
        .map(|(i, _)| i)
        .collect::<Vec<_>>();
    for i in 0..extra {
        let at = carriers[i % carriers.len()];
        let slot = hedge_slot(value.dimension, 9000 + i as i64);
        let id = value.edges[at].id as i64;
        value.edges[at].factor = &value.edges[at].factor * q(id, &slot);
        value.ports.push(slot);
    }
    value.ports.sort_by_key(|p| p.to_canonical_string());
    value.label.push_str(&format!(" +{extra}ports"));
    value.validate().unwrap();
    value
}

/// The p-route identity for one direct value: the open value in the caller's
/// notation (`evaluate_factorized_ports`) read back with `Q(x, slot) -> xi_x
/// p(1, a)` (the record's p(1)) and closed with the generic tensor equals the
/// contracted route with that tensor as projector (`evaluate_factorized`).
fn check_direct(vakint: &Vakint, value: &FactorizedNumerator) -> usize {
    let settings = refpath::vakint_settings(value.loops() as i64 + 1);
    let eps = epsilon(&settings);
    let loops = value.loops();
    let open = evaluate_factorized_ports(vakint, &settings, value, &Atom::num(1))
        .unwrap_or_else(|e| panic!("{}: open: {e}", value.label));
    assert_eq!(open.open, value.open_slots().unwrap());
    // the leg named for p(1), else the first leg carrying it
    let x = value.external_edges[0].unwrap_or_else(|| {
        value
            .edges
            .iter()
            .filter(|e| e.is_external() && e.momentum.iter().any(|c| *c != 0))
            .map(|e| e.id)
            .min()
            .unwrap()
    });
    let xi = value.edge(x).unwrap().momentum[loops];
    // caller notation -> vakint notation in the RECORD's p(1)
    let head = symbolica::function!(S.p, 1);
    let back = open.value.replace_map(|part, _, out| {
        let AtomView::Fun(f) = part else { return };
        let args = f.iter().collect::<Vec<_>>();
        let symbol = f.get_symbol();
        let name = symbol.get_name();
        if name == MOMENTUM && args.len() == 2 {
            assert_eq!(i64::try_from(args[0]).unwrap(), x as i64);
            **out =
                Atom::num(xi) * symbolica::function!(S.p, 1, index_of_slot(&args[1].to_owned()));
        } else if name == "spenso::dot" && args.len() == 2 {
            **out = symbolica::function!(S.dot, &head, &head);
        } else if name == METRIC && args.len() == 2 {
            **out = symbolica::function!(
                S.g,
                index_of_slot(&args[0].to_owned()),
                index_of_slot(&args[1].to_owned())
            );
        }
    });
    let slots = open
        .open
        .iter()
        .map(|s| (s.clone(), index_of_slot(s)))
        .collect::<Vec<_>>();
    let (caller, closing, count) =
        generic_closing(&slots, 1, &|_, slot| q(x as i64, slot), &|_, index| {
            Atom::num(xi) * symbolica::function!(S.p, 1, index)
        });
    let mut contracted_value = value.clone();
    contracted_value.projector = Some(match &value.projector {
        // a partial projector: the generic tensor closes what it leaves open
        Some(projector) => Projector {
            tensor: &projector.tensor * caller,
            p2_power: projector.p2_power,
        },
        None => Projector {
            tensor: caller,
            p2_power: 0,
        },
    });
    let contracted = evaluate_factorized(vakint, &settings, &contracted_value, &Atom::num(1))
        .unwrap_or_else(|e| panic!("{}: contracted: {e}", value.label));
    assert_closes(&back, &closing, &contracted.value, eps, &value.label);
    eprintln!(
        "WH_PORTS_DIRECT value=\"{}\" loops={loops} ports={} structures={count} open_s={:.3} contracted_s={:.3} identical=1",
        value.label,
        value.ports.len(),
        open.contract_seconds + open.evaluate_seconds,
        contracted.contract_seconds + contracted.evaluate_seconds,
    );
    count
}

/// Two loops: 1, 2, 3, 4 ports on a ghost and a quark-loop graph, both signs
/// of the record's p(1).
#[test]
fn direct_ports_close_to_the_contracted_route_two_loops() {
    let Some(root) = dot_root() else { return };
    let vakint = engine();
    for number in [17, 26] {
        let two = graph_value(&root, "L2", number);
        for value in [
            one_port(&two),
            two.clone(),
            more_ports(&two, 1),
            more_ports(&two, 2),
        ] {
            check_direct(&vakint, &value);
            check_direct(&vakint, &flipped(&value));
        }
    }
}

/// Three loops: 1, 2, 3 ports (ghost and quark-loop graphs).
#[test]
fn direct_ports_close_to_the_contracted_route_three_loops() {
    let Some(root) = dot_root() else { return };
    let vakint = engine();
    for number in [354, 67] {
        let two = graph_value(&root, "L3", number);
        for value in [one_port(&two), two.clone(), more_ports(&two, 1)] {
            check_direct(&vakint, &value);
        }
        check_direct(&vakint, &flipped(&one_port(&two)));
    }
}

/// Four loops (needs the 1-4 loop p library): 1 and 2 ports on a ghost graph.
#[test]
fn direct_ports_close_to_the_contracted_route_four_loops() {
    let Some(root) = dot_root() else { return };
    if std::env::var_os("VAKINT_TEST_P_LIBRARY_1_4").is_none() {
        eprintln!("SKIP: no VAKINT_TEST_P_LIBRARY_1_4");
        return;
    }
    let vakint = engine();
    let two = graph_value(&root, "L4", 2519);
    check_direct(&vakint, &one_port(&two));
    check_direct(&vakint, &two);
    check_direct(&vakint, &flipped(&one_port(&two)));
}

/// A hand-made two-loop sunset `k1, k2, p - k1 - k2` with a SOFT leg (edge 0,
/// momentum 0, the smallest id) whose slot is open, coupled to the momentum
/// of edge 3 (`Q(3, slot)` at vertex 0), and a second port on the outgoing
/// leg (`Q(4, slot)` at vertex 1).
/// Marked (an empty v3 direct record): a soft first leg is a new feature.
fn soft_sunset(ports: usize) -> FactorizedNumerator {
    let d = sym(DIM);
    let edge = |id, source, target, hs: Option<usize>, ht: Option<usize>, momentum: Vec<i64>| {
        FactorizedEdge {
            id,
            source,
            target,
            source_hedge: hs,
            target_hedge: ht,
            particle: "phi".into(),
            line: LineKind::Scalar,
            power: 1,
            mass_squared: Atom::Zero,
            momentum,
            factor: Atom::num(1),
        }
    };
    use Endpoint::{External, Vertex};
    let edges = vec![
        edge(0, External, Vertex(0), None, Some(0), vec![0, 0, 0]),
        edge(1, External, Vertex(0), None, Some(1), vec![0, 0, 1]),
        edge(2, Vertex(1), External, Some(2), None, vec![0, 0, 1]),
        edge(3, Vertex(0), Vertex(1), Some(3), Some(4), vec![1, 0, 0]),
        edge(4, Vertex(0), Vertex(1), Some(5), Some(6), vec![0, 1, 0]),
        edge(5, Vertex(0), Vertex(1), Some(7), Some(8), vec![-1, -1, 1]),
    ];
    let soft = hedge_slot(d, 0);
    let out = hedge_slot(d, 2);
    let mut port_list = vec![soft.clone()];
    let mut v1 = Atom::num(1);
    if ports == 2 {
        port_list.push(out.clone());
        v1 = q(4, &out);
    }
    FactorizedNumerator::from_parts(FactorizedParts {
        label: format!("soft-sunset ports={ports}"),
        dimension: d,
        loop_edges: vec![Some(3), Some(4)],
        external_edges: vec![Some(1)],
        edges,
        vertices: vec![
            FactorizedVertex {
                id: 0,
                rule: "soft".into(),
                factor: q(3, &soft),
            },
            FactorizedVertex {
                id: 1,
                rule: "out".into(),
                factor: v1,
            },
        ],
        ports: port_list,
        projector: None,
        prefactor: Prefactor::default(),
    })
    .map(|value| marked(&value))
    .unwrap()
}

#[test]
fn a_soft_leg_carries_a_port() {
    let vakint = engine();
    for ports in [1, 2] {
        let value = soft_sunset(ports);
        check_direct(&vakint, &value);
        check_direct(&vakint, &flipped(&value));
    }
}

// ---------------------------------------------------------------------------
// Refusals, text, notation

#[test]
fn port_refusals_are_precise() {
    let vakint = engine();
    let settings = refpath::vakint_settings(3);
    let base = soft_sunset(2);
    let refused = |value: &FactorizedNumerator| {
        evaluate_factorized_ports(&vakint, &settings, value, &Atom::num(1))
            .err()
            .unwrap_or_else(|| panic!("{}: accepted", value.label))
    };
    // no ceiling on the open slots (a slot count is not the loop rank): seven
    // open slots build the vakint input
    let with_extra = |extra: usize| {
        let mut value = base.clone();
        for i in 0..extra {
            let slot = hedge_slot(value.dimension, 9100 + i as i64);
            let at = value.edges.iter().position(|e| e.id == 5).unwrap();
            value.edges[at].factor = &value.edges[at].factor * q(5, &slot);
            value.ports.push(slot);
        }
        value.ports.sort_by_key(|p| p.to_canonical_string());
        value.validate().unwrap();
        value
    };
    let seven = marked(&with_extra(5));
    assert_eq!(seven.ports.len(), 7);
    let _ = vakint::p_integral::factorized::vakint_input(&settings, &seven, &Atom::num(1)).unwrap();
    // an integer index on an open slot (a projector slot that is not a port)
    let mut integer = base.clone();
    let seven = symbolica::function!(sym(MINK), Atom::var(integer.dimension), 7);
    integer.projector = Some(Projector {
        tensor: symbolica::function!(sym(METRIC), &integer.ports[0], &seven),
        p2_power: 0,
    });
    integer.validate().unwrap();
    match refused(&integer) {
        FactorizedError::Ports { detail } => assert!(detail.contains("integer index"), "{detail}"),
        other => panic!("{other}"),
    }
    // an open spinor slot (a projector slot of the spinor representation)
    let mut spinor = base.clone();
    let bis = |h: i64| {
        symbolica::function!(
            sym("spenso::bis"),
            4,
            symbolica::function!(sym(HEDGE), h, 1)
        )
    };
    spinor.projector = Some(Projector {
        tensor: symbolica::function!(sym(METRIC), &spinor.ports[0], &spinor.ports[1])
            * symbolica::function!(sym(METRIC), bis(9500), bis(9501)),
        p2_power: 0,
    });
    let error = match spinor.validate() {
        Err(error) => error,
        Ok(()) => refused(&spinor),
    };
    eprintln!("open spinor slot refused: {error}");
    // every leg soft: no external momentum to return the tensor in
    let mut all_soft = base.clone();
    for edge in &mut all_soft.edges {
        let loops = 2;
        edge.momentum[loops] = 0;
    }
    all_soft.external_edges = vec![None];
    if let Some(at) = all_soft.edges.iter().position(|e| e.id == 5) {
        all_soft.edges[at].momentum = vec![-1, -1, 0];
    }
    all_soft.validate().unwrap();
    match refused(&all_soft) {
        FactorizedError::Evaluation { detail } => assert!(detail.contains("soft"), "{detail}"),
        other => panic!("{other}"),
    }
    // a projector slot that is not a port but occurs in the factors
    let mut stray = base.clone();
    stray.ports.retain(|p| *p != hedge_slot(stray.dimension, 2));
    stray.projector = Some(Projector {
        tensor: symbolica::function!(
            sym(METRIC),
            hedge_slot(stray.dimension, 0),
            hedge_slot(stray.dimension, 2)
        ),
        p2_power: 0,
    });
    assert!(matches!(
        stray.validate(),
        Err(FactorizedError::Ports { .. })
    ));
}

/// A projector closing SOME ports: the others stay open, and a fresh slot of
/// the projector is open too; the value equals the fully open one with the
/// projector applied afterwards (route identity).
#[test]
fn a_partial_projector_leaves_the_other_slots_open() {
    let vakint = engine();
    let base = soft_sunset(2);
    let d = base.dimension;
    let nu = symbolica::function!(sym(MINK), Atom::var(d), Atom::var(sym("wh_ports_test::nu")));
    let mut partial = base.clone();
    partial.projector = Some(Projector {
        tensor: symbolica::function!(sym(METRIC), hedge_slot(d, 0), &nu),
        p2_power: 0,
    });
    partial.validate().unwrap();
    let open = partial.open_slots().unwrap();
    let mut expected = vec![hedge_slot(d, 2), nu.clone()];
    expected.sort_by_key(|p| p.to_canonical_string());
    assert_eq!(open, expected);
    // the partial value: hedge(0,1) renamed nu, i.e. the open value with the
    // soft slot renamed
    let settings = refpath::vakint_settings(3);
    let got = evaluate_factorized_ports(&vakint, &settings, &partial, &Atom::num(1)).unwrap();
    let full = evaluate_factorized_ports(&vakint, &settings, &base, &Atom::num(1)).unwrap();
    let renamed = full.value.replace(hedge_slot(d, 0)).with(nu.clone());
    assert!(
        (&got.value - &renamed).expand().is_zero(),
        "partial {}\nrenamed open {}",
        got.value,
        renamed
    );
    check_direct(&vakint, &partial);
}

#[test]
fn ports_round_trip_through_the_text_form() {
    init();
    let mut values = vec![soft_sunset(1), soft_sunset(2)];
    for name in ["hard-L3-012034.txt", "hard-L3-000993.txt"] {
        values.push(hard_as_vacuum(
            &std::fs::read_to_string(fixtures().join(name)).unwrap(),
        ));
    }
    if let Some(root) = dot_root() {
        let two = graph_value(&root, "L2", 17);
        values.push(more_ports(&two, 2));
        values.push(flipped(&one_port(&two)));
    }
    for value in values {
        let text = value.to_text().unwrap();
        let back = FactorizedNumerator::from_text(&text).unwrap();
        assert_eq!(back, value, "{}", value.label);
        assert_eq!(back.to_text().unwrap(), text);
    }
}

/// Without open slots nothing changes: a closed value through the new entry
/// returns the adapter's value (dot(p(1),p(1)) in the caller's notation).
#[test]
fn a_closed_value_through_the_ports_entry_is_the_adapters() {
    let vakint = engine();
    let settings = refpath::vakint_settings(3);
    let mut closed = soft_sunset(2);
    closed.projector = Some(Projector {
        tensor: symbolica::function!(sym(METRIC), &closed.ports[0], &closed.ports[1]),
        p2_power: -1,
    });
    let adapter = evaluate_factorized(&vakint, &settings, &closed, &Atom::num(1)).unwrap();
    let ports = evaluate_factorized_ports(&vakint, &settings, &closed, &Atom::num(1)).unwrap();
    assert!(ports.open.is_empty());
    assert_eq!(ports.vakint_value, adapter.value);
    assert!(!adapter.value.is_zero());
}

/// OQ-23: the imaginary unit and UFO couplings INSIDE edge / vertex factors
/// (scalars free of loop momenta and indices; `prefactor_scalars = 1`) are
/// accepted, and the value is the same as with them in the prefactor.
#[test]
fn scalars_inside_factors_are_accepted() {
    let vakint = engine();
    let settings = vacuum_settings();
    let text = std::fs::read_to_string(fixtures().join("hard-L3-012034.txt")).unwrap();
    let value = hard_as_vacuum(&text);
    // zeta(3) and pi are legal scalars too (WH-OQ-ANSWERS addendum, OQ-23)
    let zeta3 = symbolica::try_parse!("symbolica::{}::zeta(3)*𝜋^(-2)").unwrap();
    let couplings =
        Atom::i() * Atom::var(sym("UFO::gs_minus")) * Atom::var(sym("UFO::gs2_i")) * &zeta3;
    let mut inside = value.clone();
    inside.vertices[0].factor =
        &inside.vertices[0].factor * Atom::i() * Atom::var(sym("UFO::gs_minus"));
    let last = inside.edges.len() - 1;
    inside.edges[last].factor = &inside.edges[last].factor * Atom::var(sym("UFO::gs2_i")) * &zeta3;
    inside.validate().unwrap();
    let got = evaluate_factorized_vacuum(&vakint, &settings, &inside, &Atom::num(1)).unwrap();
    let reference = evaluate_factorized_vacuum(&vakint, &settings, &value, &couplings).unwrap();
    assert!(
        (&got.value - &reference.value).expand().is_zero(),
        "inside {}\nprefactor {}",
        got.value,
        reference.value
    );
    assert!(!reference.value.is_zero());
}

// ---------------------------------------------------------------------------
// The beta's VERIFIED K6 records (insertion normalisations, kind direct)

/// The caller's settings of the R* sites (`true_settings`: eps =
/// `gammalooprs::ε`, mu_r^2 = `gammalooprs::μᵣ²`, MSbar), RustRed only,
/// `terms` = the record's `depth`.
fn rstar_site_settings(terms: i64) -> VakintSettings {
    VakintSettings {
        epsilon_symbol: "gammalooprs::ε".into(),
        mu_r_sq_symbol: "gammalooprs::μᵣ²".into(),
        number_of_terms_in_epsilon_expansion: terms,
        ..vacuum_settings()
    }
}

/// A v3 record of kind direct, read by the v3 reader (its `open` lines are
/// checked by the reader against `open_slots`), and its `depth`.
fn direct_v3(text: &str) -> (FactorizedNumerator, i64) {
    init();
    let value = FactorizedNumerator::from_text(text).unwrap_or_else(|e| panic!("{e}"));
    let depth = value
        .record
        .as_ref()
        .and_then(|r| r.depth)
        .expect("a v3 direct record with a depth") as i64;
    // the written form reads back to the same value, and is stable
    let written = value.to_text().unwrap();
    let again = FactorizedNumerator::from_text(&written).unwrap();
    assert_eq!(again, value, "{}", value.label);
    assert_eq!(again.to_text().unwrap(), written);
    (value, depth)
}

fn fixture_value(text: &str) -> (Atom, i64) {
    let value = text
        .lines()
        .find_map(|l| l.strip_prefix("value\t"))
        .expect("a value line");
    let through = text
        .lines()
        .find_map(|l| l.strip_prefix("valid_through\t"))
        .expect("a valid_through line")
        .trim_start_matches("eps^")
        .split_whitespace()
        .next()
        .unwrap()
        .parse()
        .unwrap();
    init();
    (symbolica::try_parse!(value).unwrap(), through)
}

/// Exact equality coefficient by coefficient in `gammalooprs::ε` through
/// `order` and no further (the fixtures' identity rule), algebraic.
fn assert_equal_through(got: &Atom, expected: &Atom, order: i64, what: &str) {
    let eps = sym("gammalooprs::ε");
    let got = open_indices::epsilon_coefficients(got, eps);
    let expected = open_indices::epsilon_coefficients(expected, eps);
    let first = got
        .keys()
        .chain(expected.keys())
        .copied()
        .min()
        .unwrap_or(0);
    for n in first..=order {
        let g = got.get(&n).cloned().unwrap_or(Atom::Zero);
        let e = expected.get(&n).cloned().unwrap_or(Atom::Zero);
        let difference = (&g - &e).expand();
        assert!(
            difference.is_zero(),
            "{what}: eps^{n}: entry {g}\nbeta  {e}\ndifference {difference}"
        );
    }
    assert!(
        expected.keys().any(|n| *n <= order),
        "{what}: the reference is zero through eps^{order}"
    );
}

/// K6 (f): each VERIFIED record (read as written: schema v3, kind direct)
/// through `evaluate_factorized_ports` equals the value the beta's contracted
/// path returns at that call (the return of `integrate_propagator_with`,
/// their value.txt), through
/// `valid_through`: the quark form factor (compact slash projector, open
/// quark legs closed by it), the gluon (both channel markers in one
/// projector), the gluon with a quark loop, the scalar bubble.
#[test]
fn the_beta_k6_records_equal_their_contracted_values() {
    let vakint = engine();
    for name in ["quark", "gluon", "gluon-qloop", "bubble"] {
        let input =
            std::fs::read_to_string(fixtures().join(format!("k6-{name}.input.txt"))).unwrap();
        let reference =
            std::fs::read_to_string(fixtures().join(format!("k6-{name}.value.txt"))).unwrap();
        let (value, depth) = direct_v3(&input);
        // the comparison point: the return of the site's
        // integrate_propagator_with (no 1/p^2, no 1/Tr[1], no channel solve:
        // p2_power 0, nothing else applied by the entry)
        assert_eq!(value.projector.as_ref().map_or(0, |p| p.p2_power), 0);
        let (expected, through) = fixture_value(&reference);
        let settings = rstar_site_settings(depth);
        let prefactor =
            &value.prefactor.scalars * &value.prefactor.colour * &value.prefactor.overall;
        let got = evaluate_factorized_ports(&vakint, &settings, &value, &prefactor)
            .unwrap_or_else(|e| panic!("{name}: {e}"));
        eprintln!(
            "K6 {name} entry {}\nK6 {name} beta  {}",
            got.value, expected
        );
        assert_equal_through(&got.value, &expected, through, name);
        eprintln!(
            "WH_PORTS_K6 record={name} identical_through=eps^{through} contract_s={:.3} evaluate_s={:.3}",
            got.contract_seconds, got.evaluate_seconds
        );
    }
}

/// R4: the compact slash in a projector is written out on read (a fresh
/// dummy per occurrence), idempotently; refused under a Taylor step.
#[test]
fn the_compact_slash_is_written_out_on_read() {
    init();
    let input = std::fs::read_to_string(fixtures().join("k6-quark.input.txt")).unwrap();
    let (value, _) = direct_v3(&input);
    let tensor = value
        .projector
        .as_ref()
        .unwrap()
        .tensor
        .to_canonical_string();
    assert!(
        tensor.contains("vakint::{}::slash(-1,0)") || tensor.contains("slash(-1,0)"),
        "{tensor}"
    );
    let mut again = value.clone();
    again.expand_compact_slashes().unwrap();
    assert_eq!(again, value);
    let compact = symbolica::try_parse!(
        "spenso::{linear,spenso::tensor}::gamma(spenso::{spenso::representation,spenso::self_dual}::bis(4,gammalooprs::{spenso::index}::hedge(0,1)),spenso::{spenso::representation,spenso::self_dual}::bis(4,gammalooprs::{spenso::index}::hedge(1,1)),gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript,spenso::tensor-label:q}::Q(0,spenso::{spenso::representation,spenso::self_dual}::mink(gammalooprs::{}::dim)))"
    )
    .unwrap();
    let other = compact
        .replace(symbolica::try_parse!("gammalooprs::{spenso::index}::hedge(0,1)").unwrap())
        .with(symbolica::try_parse!("gammalooprs::{spenso::index}::hedge(5,1)").unwrap());
    let written = vakint::p_integral::factorized::expand_compact_slash(&(&compact * &other), 7);
    let text = written.to_canonical_string();
    assert!(
        text.contains("slash(7,0)") && text.contains("slash(7,1)"),
        "{text}"
    );
    // under a Taylor step: refused (the beta's k3-hs closure record, degree 1)
    let k3 = std::fs::read_to_string(fixtures().join("k3-hs.input.txt")).unwrap();
    let mut taylor = FactorizedNumerator::parse_text(&k3).unwrap();
    taylor.projector = Some(Projector {
        tensor: compact.clone(),
        p2_power: 0,
    });
    assert!(matches!(
        taylor.expand_compact_slashes(),
        Err(FactorizedError::Ports { .. })
    ));
}
