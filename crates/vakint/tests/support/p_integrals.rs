//! Shared helpers of the p-integral (massless propagator-type) tests.
//!
//! Include it with `#[path = "support/p_integrals.rs"] mod p_integrals;`.
//! Integrals are written as strings in vakint's `topo(prop(...))` notation,
//! with an optional numerator factor in front.

// The frontend and Forcer test binaries each use a part of these helpers.
#![allow(dead_code)]

use std::sync::LazyLock;

use feynkit_graph::FeynmanDiagram;
use feynkit_kinematics::Kinematics;
use feynkit_model::Model;
use symbolica::atom::{Atom, AtomCore, AtomView};
use symbolica::id::Replacement;
use vakint::diagram_integral::DiagramIntegralOptions;
use vakint::symbols::S;
use vakint::{
    EvaluationMethod, EvaluationOrder, ForcerOptions, LoopNormalizationFactor, RustRedMasterMode,
    Vakint, VakintDependency, VakintExpression, VakintSettings, p_integral, vakint_parse,
};

pub static VAKINT: LazyLock<Vakint> = LazyLock::new(|| Vakint::new().unwrap());

/// The scalar `phi` of the DOT tests: massless, vertices of valence 3 to 10.
pub static PHI: LazyLock<Model> = LazyLock::new(|| {
    Model::from_json(include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/phi_p_graphs.json"
    )))
    .unwrap()
});

pub fn parse(input: &str) -> Atom {
    vakint_parse!(input).unwrap_or_else(|error| panic!("{input}: {error}"))
}

/// Vakint's canonical long form (family, matched powers, `σ·p(i)`) of `input`.
pub fn canonical(input: &Atom) -> Atom {
    VAKINT
        .to_canonical(&VakintSettings::default(), input.as_view(), false)
        .unwrap_or_else(|error| panic!("{input}: {error}"))
}

/// The canonical short form `topo(NAME(σ·p(i), a_1, ..., a_N))` of `input`.
pub fn short_form(input: &str) -> Atom {
    VAKINT
        .to_canonical(&VakintSettings::default(), parse(input).as_view(), true)
        .unwrap_or_else(|error| panic!("{input}: {error}"))
}

/// `input` is scaleless: it canonicalizes and evaluates to an exact zero.
pub fn assert_p_scaleless(input: &str) {
    let atom = parse(input);
    assert!(canonical(&atom).is_zero(), "{input} must vanish");
    assert!(
        VAKINT
            .evaluate(&VakintSettings::default(), atom.as_view())
            .unwrap_or_else(|error| panic!("{input}: {error}"))
            .is_zero(),
        "{input} must evaluate to 0"
    );
}

/// `input` is rejected with an error mentioning `fragment` when unknown
/// integrals are not allowed (an out-of-domain p-integral is otherwise an
/// unknown integral, see [`assert_p_unknown`]).
pub fn assert_p_error(input: &str, fragment: &str) {
    let settings = VakintSettings {
        allow_unknown_integrals: false,
        ..VakintSettings::default()
    };
    assert_canonical_error(input, fragment, &settings);
}

fn assert_canonical_error(input: &str, fragment: &str, settings: &VakintSettings) {
    match VAKINT.to_canonical(settings, parse(input).as_view(), false) {
        Ok(result) => panic!("{input} was accepted as {result}"),
        Err(error) => assert!(
            error.to_string().contains(fragment),
            "{input}: '{error}' does not mention '{fragment}'"
        ),
    }
}

/// The short form of `input` is a sector of `family` with `pinched` lines.
pub fn assert_p_family(input: &str, family: &str, pinched: usize) {
    let short = short_form(input);
    let factors = match short.as_view() {
        AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
        view => vec![view],
    };
    let head = factors
        .into_iter()
        .find_map(|factor| match factor {
            AtomView::Fun(topo) if topo.get_symbol() == S.topo => topo.iter().next(),
            _ => None,
        })
        .unwrap_or_else(|| panic!("{input}: {short} has no topo(...)"));
    let AtomView::Fun(head) = head else {
        panic!("{input}: {short} is not a short form");
    };
    assert_eq!(
        head.get_symbol().get_stripped_name(),
        family,
        "{input}: {short}"
    );
    let zeros = head.iter().skip(1).filter(|power| power.is_zero()).count();
    assert_eq!(zeros, pinched, "{input}: {short}");
}

/// The top of `family` with line powers `1, 2, ..., N` (no graph symmetry
/// survives distinct powers) in the family routing, times `numerator`.
pub fn distinct_power_top(family: &str, lines: usize, numerator: &str) -> Atom {
    let powers = (1..=lines)
        .map(|n| n.to_string())
        .collect::<Vec<_>>()
        .join(",");
    parse(numerator) * canonical(&parse(&format!("topo({family}(p(1),{powers}))")))
}

/// Two inputs have the same canonical form.
pub fn assert_p_same_canonical(left: &str, right: &str) {
    assert_eq!(
        canonical(&parse(left)),
        canonical(&parse(right)),
        "{left} vs {right}"
    );
}

/// A deterministic xorshift generator for reproducible random routings.
struct Random(u64);

impl Random {
    fn next(&mut self, bound: usize) -> usize {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 % bound as u64) as usize
    }

    fn sign(&mut self) -> i64 {
        if self.next(2) == 0 { 1 } else { -1 }
    }
}

/// Number of distinct loop momenta `k(j)` of `input`, which uses `k(1..=L)`.
fn loops(input: &Atom) -> usize {
    input
        .pattern_match(&parse("k(i_)").to_pattern(), None, None)
        .map(|matched| matched[&wildcard("i_")].clone())
        .collect::<std::collections::BTreeSet<_>>()
        .len()
}

fn wildcard(name: &str) -> symbolica::atom::Symbol {
    parse(name).get_symbol().unwrap()
}

/// Rewrite every `prop(...)` of `input` with `rewrite(id, from, to, momentum)`.
fn map_props(input: &Atom, rewrite: impl Fn(i64, i64, i64, Atom) -> (i64, i64, Atom)) -> Atom {
    input.replace_map(|view, _, out| {
        let AtomView::Fun(prop) = view else {
            return;
        };
        if prop.get_symbol() != S.prop || prop.get_nargs() != 5 {
            return;
        }
        let arguments = prop.iter().collect::<Vec<_>>();
        let integer = |view: AtomView| i64::try_from(view).unwrap();
        let AtomView::Fun(edge) = arguments[1] else {
            panic!("{view} has no edge(a,b)");
        };
        let mut ends = edge.iter().map(integer);
        let id = integer(arguments[0]);
        let (from, to, momentum) = rewrite(
            id,
            ends.next().unwrap(),
            ends.next().unwrap(),
            arguments[2].to_owned(),
        );
        **out = S.prop.call((
            id,
            S.edge.call((from, to)),
            momentum.expand(),
            arguments[3].to_owned(),
            arguments[4].to_owned(),
        ));
    })
}

/// The whole integrand (numerator and props) with `rules` applied at once.
fn substitute(input: &Atom, rules: &[(Atom, Atom)]) -> Atom {
    let rules = rules
        .iter()
        .map(|(source, target)| Replacement::new(source.to_pattern(), target.to_pattern()))
        .collect::<Vec<_>>();
    map_props(&input.replace_multiple(&rules), |_, a, b, q| (a, b, q))
}

/// Routing invariance (section 7): for every seed, a random unimodular change
/// of loop momenta `k -> M k + c·p(1)`, random edge flips and a random vertex
/// relabeling, then the reversal `p(1) -> -p(1)` and the renaming
/// `p(1) -> p(3)` of the whole integrand, all give exactly the canonical form
/// of `input` (with the reversal and the renaming undone). Choose powers that
/// break the graph symmetries: a vertex relabeling is only exact up to the
/// automorphisms that preserve the powers.
pub fn assert_p_routing_invariant(input: &str, seeds: std::ops::Range<u64>) {
    routing_invariance(&parse(input), seeds, true, canonical);
}

/// [`assert_p_routing_invariant`] for a parsed integral.
pub fn assert_p_atom_routing_invariant(original: &Atom, seeds: std::ops::Range<u64>) {
    routing_invariance(original, seeds, true, canonical);
}

/// [`assert_p_routing_invariant`] without vertex relabelings, for graphs with
/// two vertices: exchanging them maps every line to itself, so the reversal
/// `q -> -q` is a symmetry and a relabeling may pick either orientation.
pub fn assert_p_reroutings_invariant(input: &str, seeds: std::ops::Range<u64>) {
    routing_invariance(&parse(input), seeds, false, canonical);
}

/// `observe` (the canonical form, or an evaluation) of every rerouting,
/// reversal and renaming of `original` equals that of `original`.
fn routing_invariance(
    original: &Atom,
    seeds: std::ops::Range<u64>,
    relabel: bool,
    observe: impl Fn(&Atom) -> Atom,
) {
    let original = original.clone();
    let expected = observe(&original);
    let loops = loops(&original);
    let (k, p, p3) = (
        |axis: usize| vakint_parse!(format!("k({axis})")).unwrap(),
        parse("p(1)"),
        parse("p(3)"),
    );
    for seed in seeds {
        let mut random =
            Random(0x9E37_79B9_7F4A_7C15 ^ (seed + 1).wrapping_mul(0x2545_F491_4F6C_DD1D));
        // A unimodular integer matrix from random elementary operations.
        let mut matrix = (0..loops)
            .map(|i| (0..loops).map(|j| i64::from(i == j)).collect::<Vec<_>>())
            .collect::<Vec<_>>();
        for _ in 0..3 * loops {
            let (i, j) = (random.next(loops), random.next(loops));
            if i != j {
                let sign = random.sign();
                let row = matrix[j].clone();
                matrix[i]
                    .iter_mut()
                    .zip(row)
                    .for_each(|(a, b)| *a += sign * b);
            } else {
                matrix[i].iter_mut().for_each(|a| *a = -*a);
            }
        }
        let images = (0..loops)
            .map(|i| {
                let shift = random.next(3) as i64 - 1;
                let image = (0..loops).fold(Atom::num(shift) * &p, |sum, j| {
                    sum + Atom::num(matrix[i][j]) * k(j + 1)
                });
                (k(i + 1), image)
            })
            .collect::<Vec<_>>();
        let rerouted = substitute(&original, &images);
        let flips = (0..64).map(|_| random.next(2) == 1).collect::<Vec<_>>();
        let shift = 10 + random.next(5) as i64;
        let mut order = (0..32).collect::<Vec<i64>>();
        for i in (1..order.len()).rev() {
            order.swap(i, random.next(i + 1));
        }
        let relabeled = map_props(&rerouted, |id, a, b, q| {
            let (a, b) = if relabel {
                (order[a as usize] + shift, order[b as usize] + shift)
            } else {
                (a, b)
            };
            if flips[id as usize % flips.len()] {
                (b, a, -q)
            } else {
                (a, b, q)
            }
        });
        let found = observe(&relabeled);
        assert_eq!(found, expected, "seed {seed}: {relabeled}");
        let reversed = substitute(&relabeled, &[(p.clone(), -p.clone())]);
        let found = substitute(&observe(&reversed), &[(p.clone(), -p.clone())]);
        assert_eq!(found, expected, "seed {seed}, p(1) -> -p(1): {reversed}");
        let renamed = substitute(&relabeled, &[(p.clone(), p3.clone())]);
        let found = substitute(&observe(&renamed), &[(p3.clone(), p.clone())]);
        assert_eq!(found, expected, "seed {seed}, p(1) -> p(3): {renamed}");
    }
}

/// Read a compact DOT p-graph back into vakint (`from_diagram`).
pub fn from_dot(dot: &str) -> Atom {
    from_dot_with(dot, &DiagramIntegralOptions::default())
        .unwrap_or_else(|error| panic!("{dot}: {error}"))
}

/// [`from_dot`] with explicit options.
pub fn from_dot_with(
    dot: &str,
    options: &DiagramIntegralOptions,
) -> Result<Atom, vakint::VakintError> {
    let diagram =
        FeynmanDiagram::from_dot(PHI.clone(), dot).unwrap_or_else(|error| panic!("{dot}: {error}"));
    let family = diagram.propagator_family(&Kinematics::new()).unwrap();
    VakintExpression::from_diagram(&diagram, &family, &Atom::one(), options).map(Atom::from)
}

/// DOT round trip (section 3): `atom -> to_dot -> from_diagram` keeps the
/// canonical form, and `dot -> atom -> dot -> atom` is a fixed point after
/// the first export.
pub fn assert_p_dot_round_trip(input: &str) {
    let atom = parse(input);
    let dot = p_integral::to_dot(atom.as_view(), "phi").unwrap();
    let imported = from_dot(&dot);
    assert_eq!(canonical(&imported), canonical(&atom), "{input}\n{dot}");
    let again = p_integral::to_dot(imported.as_view(), "phi").unwrap();
    assert_eq!(from_dot(&again), imported, "{dot}\n{again}");
}

/// FORM 5 for live Forcer runs: `VAKINT_TEST_FORM5`, otherwise the default
/// FORM (`FORM_PATH` or `form`) when it is FORM 5, as in the dev shell and CI.
pub fn forcer_form() -> Option<String> {
    std::env::var("VAKINT_TEST_FORM5").ok().or_else(|| {
        let form = std::env::var("FORM_PATH").unwrap_or_else(|_| "form".into());
        VakintDependency::Forcer
            .version(&form, "python3")
            .is_ok()
            .then_some(form)
    })
}

/// Settings that evaluate p-integrals with Forcer alone (FORM at `form`), in
/// `master_mode` and `normalization`, through `terms` Laurent terms.
pub fn forcer_settings(
    form: &str,
    master_mode: RustRedMasterMode,
    normalization: LoopNormalizationFactor,
    terms: i64,
) -> VakintSettings {
    VakintSettings {
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::Forcer(ForcerOptions {
            master_mode,
        })]),
        form_exe_path: form.into(),
        integral_normalization_factor: normalization,
        number_of_terms_in_epsilon_expansion: terms,
        use_dot_product_notation: true,
        n_cores: 1,
        ..VakintSettings::default()
    }
}

/// The Forcer evaluation of `input` with `settings`.
pub fn forcer_evaluate(input: &Atom, settings: &VakintSettings) -> Atom {
    VAKINT
        .evaluate(settings, input.as_view())
        .unwrap_or_else(|error| panic!("{input}: {error}"))
}

/// Laurent coefficients `ε^first..=ε^last` of `value` at `Q² = μ² = 1`:
/// `dot(p(1),p(1)) = -1` and every logarithm zero.
pub fn at_unit_scale(value: &Atom, first: i64, last: i64) -> Vec<Atom> {
    let x = Atom::var(vakint::vakint_symbol!("x_"));
    let p = symbolica::function!(S.p, &x);
    let unit = value
        .replace(symbolica::function!(S.dot, &p, &p).to_pattern())
        .with(Atom::num(-1).to_pattern())
        .replace(symbolica::function!(symbolica::atom::Symbol::LOG, &x).to_pattern())
        .with(Atom::Zero.to_pattern())
        .expand();
    let epsilon = Atom::var(vakint::vakint_symbol!("ε"));
    let shifted = (unit * epsilon.pow(Atom::num(-first))).expand();
    let mut coefficients = vec![Atom::Zero; (last - first + 1) as usize];
    for (monomial, coefficient) in shifted.coefficient_list::<i8>(std::slice::from_ref(&epsilon)) {
        let power = match monomial.as_view() {
            AtomView::Num(_) => 0,
            AtomView::Pow(power) => i64::try_from(power.get_base_exp().1).unwrap(),
            _ => 1,
        };
        if let Some(slot) = coefficients.get_mut(power as usize) {
            *slot = coefficient.expand();
        }
    }
    coefficients
}

/// The G-scheme Laurent coefficients of a Forcer evaluation of `input` at
/// `Q² = μ² = 1`, from `ε^-L` through `ε^last`.
pub fn assert_forcer_g_scheme(input: &str, loops: i64, last: i64, expected: &[Atom]) {
    let Some(form) = forcer_form() else {
        eprintln!("FORM 5 is unavailable: skipped the live Forcer evaluation of {input}");
        return;
    };
    let settings = forcer_settings(
        &form,
        RustRedMasterMode::Symbolic,
        LoopNormalizationFactor::GScheme,
        last + loops + 1,
    );
    let value = forcer_evaluate(&parse(input), &settings);
    let found = at_unit_scale(&value, -loops, last);
    let expected = expected.iter().map(|e| e.expand()).collect::<Vec<_>>();
    assert_eq!(found, expected, "{input}: {value}");
}

/// [`assert_p_routing_invariant`] on exact (`Retain`) Forcer evaluations.
pub fn assert_forcer_routing_invariant(input: &str, seeds: std::ops::Range<u64>) {
    assert_forcer_atom_routing_invariant(&parse(input), seeds);
}

/// [`assert_forcer_routing_invariant`] for a parsed integral.
pub fn assert_forcer_atom_routing_invariant(original: &Atom, seeds: std::ops::Range<u64>) {
    let Some(settings) = live_forcer_retain_settings(original) else {
        return;
    };
    routing_invariance(original, seeds, true, |integral| {
        forcer_evaluate(integral, &settings)
    });
}

/// The DOT export of `input` (a normalized p-integral with positive powers),
/// read back through `from_diagram`, evaluates exactly like `input`.
pub fn assert_forcer_dot_like_atom(input: &Atom) {
    let Some(settings) = live_forcer_retain_settings(input) else {
        return;
    };
    let dot = p_integral::to_dot(input.as_view(), "phi").unwrap();
    assert_eq!(
        forcer_evaluate(&from_dot(&dot), &settings),
        forcer_evaluate(input, &settings),
        "{input}\n{dot}"
    );
}

/// The default evaluation order (RustRed first) hands the p-integral terms of
/// `vacuum + p_integral` to Forcer and the vacuum term to its own backend: the
/// sum evaluates like its two terms evaluated apart, and the p-integral alone
/// like Forcer alone.
pub fn assert_default_order_dispatches_to_forcer(vacuum: &str, p_integral: &str) {
    let Some(form) = forcer_form() else {
        eprintln!("FORM 5 is unavailable: skipped the default-order evaluation of {p_integral}");
        return;
    };
    let forcer = forcer_settings(
        &form,
        RustRedMasterMode::Symbolic,
        LoopNormalizationFactor::MSbar,
        3,
    );
    // The default order with Forcer in the master mode of `forcer`.
    let mut order = EvaluationOrder::available(&form, "python3");
    assert!(matches!(order.0[0], EvaluationMethod::RustRed(_)));
    for method in order.0.iter_mut() {
        if let EvaluationMethod::Forcer(options) = method {
            options.master_mode = RustRedMasterMode::Symbolic;
        }
    }
    let default = VakintSettings {
        evaluation_order: order,
        ..forcer.clone()
    };
    let (vacuum, p_integral) = (parse(vacuum), parse(p_integral));
    let p_value = forcer_evaluate(&p_integral, &default);
    assert_eq!(p_value, forcer_evaluate(&p_integral, &forcer));
    assert_eq!(
        forcer_evaluate(&(&vacuum + &p_integral), &default).expand(),
        (forcer_evaluate(&vacuum, &default) + p_value).expand()
    );
}

/// Exact (`Retain`) G-scheme Forcer settings, or `None` (with a skip
/// message) without FORM 5.
fn live_forcer_retain_settings(input: &Atom) -> Option<VakintSettings> {
    let Some(form) = forcer_form() else {
        eprintln!("FORM 5 is unavailable: skipped the live Forcer evaluation of {input}");
        return None;
    };
    Some(forcer_settings(
        &form,
        RustRedMasterMode::Retain,
        LoopNormalizationFactor::GScheme,
        4,
    ))
}

/// `ζ(n)`, as the finalized results write it.
pub fn zeta(n: i64) -> Atom {
    symbolica::transcendental::zeta().call(Atom::num(n))
}

/// `input` is malformed: rejected with an error mentioning `fragment` with
/// and without `allow_unknown_integrals` (section 1, "Out of domain").
pub fn assert_p_malformed(input: &str, fragment: &str) {
    assert_p_error(input, fragment);
    assert_canonical_error(input, fragment, &VakintSettings::default());
}

/// `input` mentions an external momentum but is no p-integral vakint covers
/// (section 1, "Out of domain"): without `allow_unknown_integrals` it is an
/// error mentioning `fragment`; with it (the default) it canonicalizes to the
/// unknown integral `topo(UNKNOWN(...))` of its long form, which is `input`
/// itself, and evaluating it is the out-of-domain error, not a value.
pub fn assert_p_unknown(input: &str, fragment: &str) {
    assert_p_error(input, fragment);
    let settings = VakintSettings::default();
    assert!(settings.allow_unknown_integrals);
    let atom = parse(input);
    let unknown = VAKINT
        .to_canonical(&settings, atom.as_view(), false)
        .unwrap_or_else(|error| panic!("{input}: {error}"));
    let expected = atom
        .replace(parse("topo(x_)").to_pattern())
        .with(parse("topo(UNKNOWN(x_))").to_pattern());
    assert_eq!(unknown, expected, "{input}");
    match VAKINT.evaluate(&settings, atom.as_view()) {
        Ok(value) => panic!("{input} is no p-integral, yet it evaluated to {value}"),
        Err(error) => assert!(
            error.to_string().contains("out-of-domain p-integral"),
            "{input}: evaluating gave '{error}', not the out-of-domain error"
        ),
    }
}

/// `lhs` and `rhs` evaluate to exactly the same value with Forcer (symbolic
/// masters, G-scheme, six Laurent terms: through `ε^(5-L)`, so `ε^1` at four
/// loops), and that value is not zero: an analytic identity between
/// p-integrals. Skipped without FORM 5.
pub fn assert_p_identity(lhs: &str, rhs: &str) {
    assert_p_identity_with(lhs, rhs, LoopNormalizationFactor::GScheme);
}

/// [`assert_p_identity`] in the loop normalization `normalization`.
pub fn assert_p_identity_with(lhs: &str, rhs: &str, normalization: LoopNormalizationFactor) {
    let Some((left, right)) = p_identity_sides(lhs, rhs, normalization) else {
        eprintln!("FORM 5 is unavailable: skipped the live Forcer identity {lhs} = {rhs}");
        return;
    };
    assert!(
        !left.is_zero(),
        "{lhs} vanishes at this truncation: the identity tests nothing"
    );
    let difference = (left - right).expand();
    assert!(difference.is_zero(), "{lhs} - ({rhs}) = {difference}");
}

/// Both sides of an identity at the truncation of [`assert_p_identity`], which
/// reaches every order of a finite four-loop integral up to `ε^1`.
pub fn p_identity_sides(
    lhs: &str,
    rhs: &str,
    normalization: LoopNormalizationFactor,
) -> Option<(Atom, Atom)> {
    let form = forcer_form()?;
    let settings = forcer_settings(&form, RustRedMasterMode::Symbolic, normalization, 6);
    Some((
        forcer_evaluate(&parse(lhs), &settings),
        forcer_evaluate(&parse(rhs), &settings),
    ))
}

/// A fake executable that prints one line, whatever its arguments, and does
/// nothing else; it lives as long as this value.
pub struct FakeExecutable {
    _directory: tempfile::TempDir,
    pub path: String,
}

fn fake_executable(name: &str, output: &str) -> FakeExecutable {
    use std::io::Write;
    use std::os::unix::fs::PermissionsExt;
    let directory = tempfile::tempdir().unwrap();
    let path = directory.path().join(name);
    {
        let mut script = std::fs::File::create(&path).unwrap();
        writeln!(script, "#!/bin/sh\necho \"{output}\"").unwrap();
    }
    std::fs::set_permissions(&path, std::fs::Permissions::from_mode(0o755)).unwrap();
    FakeExecutable {
        path: path.to_str().unwrap().to_owned(),
        _directory: directory,
    }
}

/// A fake FORM of `version`, e.g. "4.3.1" for a host without FORM 5.
pub fn fake_form(version: &str) -> FakeExecutable {
    fake_executable(
        "form",
        &format!("FORM {version} (vakint test fake) 64-bits"),
    )
}

/// A fake python whose every run reports the pySecDec `version`.
pub fn fake_python(version: &str) -> FakeExecutable {
    fake_executable("python3", version)
}
