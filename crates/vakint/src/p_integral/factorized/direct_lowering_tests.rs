//! The open tensor of a direct value reduced by the adapter
//! (`reduce_open_tensor`, then vakint's integral entry) against the former
//! route (vakint's tensor prepass inside `Vakint::evaluate`): the same value,
//! exactly (canonical strings), on the repository's direct records with open
//! ports; one test per fallback reason. And a replay of dumped direct records
//! (`GAMMALOOP_UNCONTRACTED_DUMP`'s `*.record.txt`) through
//! [`evaluate_factorized_ports`], one at a time: the wall seconds per record
//! and, with `DIRECT_REPLAY_OUT`, the value of each record written as its
//! canonical string (for exact comparisons of two builds);
//! `DIRECT_REPLAY_COMPARE=1` also evaluates the former route and asserts the
//! same value. `#[ignore]`: the records live outside the repository
//! (`DIRECT_REPLAY_RECORDS`: a file with one record path per line); four-loop
//! records need the 1-4 loop package of `VAKINT_TEST_P_LIBRARY_1_4`.
use super::evaluate::{
    FORCED_FALLBACK, FORMER_ROUTE, LAST_ROUTE, OpenTensorFallback,
    evaluate_factorized_prepass_route, has_loop_index, reduce_open_tensor,
};
use super::*;
use crate::{
    EvaluationMethod, EvaluationOrder, LoopNormalizationFactor, RustRedEvaluationOptions,
    RustRedMasterMode, Vakint, VakintSettings,
};

fn site_settings(terms: i64) -> VakintSettings {
    VakintSettings {
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions {
                master_mode: RustRedMasterMode::Symbolic,
            },
        )]),
        form_exe_path: "/nonexistent/form".into(),
        integral_normalization_factor: LoopNormalizationFactor::MSbar,
        number_of_terms_in_epsilon_expansion: terms,
        use_dot_product_notation: true,
        n_cores: 1,
        epsilon_symbol: "gammalooprs::ε".into(),
        mu_r_sq_symbol: "gammalooprs::μᵣ²".into(),
        ..Default::default()
    }
}

fn engine() -> Vakint {
    let mut vakint = Vakint::new().unwrap();
    if let Some(directory) = std::env::var_os("VAKINT_TEST_P_LIBRARY_1_4") {
        let library = rustred_app::propagator::load_p_library(
            rustred_app::propagator::PLibrarySource::Directory(directory.into()),
            &rustred_app::propagator::PLibraryOptions {
                verify_descent: false,
            },
        )
        .unwrap();
        vakint.inject_p_library(library).unwrap();
    }
    vakint
}

/// `value` marked as the ports entry marks a value without new features (an
/// empty v3 direct record), so that both routes see the same value.
fn marked(value: FactorizedNumerator) -> FactorizedNumerator {
    if value.has_new_features() {
        return value;
    }
    FactorizedNumerator {
        record: Some(DirectRecord {
            original_cycles: None,
            depth: None,
            weights: Vec::new(),
        }),
        ..value
    }
}

/// Both routes, exactly the same value (vakint notation), and nothing of the
/// loop tensor in the reduced route's integrand.
fn assert_routes_agree(vakint: &Vakint, value: &FactorizedNumerator, terms: i64, what: &str) {
    let settings = site_settings(terms);
    let prefactor = value.prefactor.total();
    let got = evaluate_factorized_ports(vakint, &settings, value, &prefactor)
        .unwrap_or_else(|e| panic!("{what}: {e}"));
    let former = evaluate_factorized_prepass_route(vakint, &settings, value, &prefactor)
        .unwrap_or_else(|e| panic!("{what}: {e}"))
        .unwrap_or_else(|| panic!("{what}: the adapter did not reduce the open tensor"));
    assert!(
        !got.vakint_value.is_zero(),
        "{what}: a zero value compares nothing"
    );
    assert_eq!(
        got.vakint_value.to_canonical_string(),
        former.to_canonical_string(),
        "{what}: the reduced route changed the value"
    );
    assert!(
        !has_loop_index(got.input.as_view()),
        "{what}: a loop momentum kept an open index"
    );
}

/// The beta's real two-loop K5 direct record (two open ports, one insertion)
/// and the one-loop gluon and quark bubbles with their projector removed
/// (two open ports each): the reduced route gives the former route's value.
#[test]
fn open_tensor_reduction_keeps_the_prepass_value() {
    let vakint = engine();
    assert_routes_agree(&vakint, &k5_direct(), 2, "k5-direct");
    for (what, text) in [
        ("gluon bubble", include_str!("fixtures/l1-gluon-bubble.txt")),
        ("quark bubble", include_str!("fixtures/l1-quark-bubble.txt")),
    ] {
        let mut value = FactorizedNumerator::from_text(text).unwrap();
        value.projector = None;
        assert_routes_agree(&vakint, &marked(value), 2, what);
    }
}

/// A closed value (a projector closing every port) is not reduced by the
/// adapter: the former route is its only one.
#[test]
fn a_closed_value_keeps_the_former_route() {
    let vakint = engine();
    let value = marked(
        FactorizedNumerator::from_text(include_str!("fixtures/l1-gluon-bubble.txt")).unwrap(),
    );
    let former = evaluate_factorized_prepass_route(
        &vakint,
        &site_settings(2),
        &value,
        &value.prefactor.total(),
    )
    .unwrap();
    assert!(
        former.is_none(),
        "a closed value was reduced by the adapter"
    );
}

/// Parsing outside this module's `Result` alias (`vk_parse!` names `Result`).
mod parse {
    pub(super) fn atom(text: &str) -> symbolica::atom::Atom {
        crate::utils::vakint_macros::vk_parse!(text).unwrap()
    }
}

fn k5_direct() -> FactorizedNumerator {
    FactorizedNumerator::from_text(include_str!("fixtures/wh-k5-direct.record.txt")).unwrap()
}

/// The fixture is reduced by the adapter (the reference of the fallbacks
/// below).
#[test]
fn the_k5_direct_record_is_reduced_by_the_adapter() {
    assert!(reduce_open_tensor(&k5_direct()).is_ok());
}

/// The ports entry's outcome on `value` (its value in vakint notation, or
/// its error) and the last open-tensor route it logged; `former`: the former
/// route everywhere (the code before the adapter's reduction); `forced`: the
/// route takes that fallback at its own site.
fn outcome(
    vakint: &Vakint,
    value: &FactorizedNumerator,
    former: bool,
    forced: Option<OpenTensorFallback>,
) -> (std::result::Result<String, String>, Option<String>) {
    FORMER_ROUTE.with(|f| f.set(former));
    FORCED_FALLBACK.with(|f| f.set(forced));
    LAST_ROUTE.with(|last| last.take());
    let got = evaluate_factorized_ports(vakint, &site_settings(2), value, &value.prefactor.total())
        .map(|got| got.vakint_value.to_canonical_string())
        .map_err(|error| error.to_string());
    FORMER_ROUTE.with(|f| f.set(false));
    FORCED_FALLBACK.with(|f| f.set(None));
    (got, LAST_ROUTE.with(|last| last.take()))
}

/// EVERY fallback reason, forced at its own site on the beta's k5-direct
/// record: the route logs that reason and the entry returns exactly the
/// former route's value (which is also the reduced route's).
#[test]
fn every_fallback_reason_returns_the_former_routes_value() {
    let vakint = engine();
    let value = k5_direct();
    let (former, route) = outcome(&vakint, &value, true, None);
    assert_eq!(route, None, "the former route reduces nothing");
    assert!(former.is_ok(), "{former:?}");
    let (reduced, route) = outcome(&vakint, &value, false, None);
    assert_eq!(route.as_deref(), Some("reduced"));
    assert_eq!(reduced, former, "the reduced route's value");
    for reason in OpenTensorFallback::ALL {
        let (got, route) = outcome(&vakint, &value, false, Some(reason));
        assert_eq!(
            route.as_deref(),
            Some(format!("fallback:{}", reason.as_str()).as_str()),
            "{reason:?}"
        );
        assert_eq!(
            got, former,
            "fallback {reason:?}: not the former route's value"
        );
    }
}

/// `value`'s outcome through the entry equals the former route's (value or
/// error); returns the route the entry logged.
fn assert_same_outcome(vakint: &Vakint, value: &FactorizedNumerator, what: &str) -> Option<String> {
    let (former, route) = outcome(vakint, value, true, None);
    assert_eq!(route, None, "{what}: the former route reduces nothing");
    let (got, route) = outcome(vakint, value, false, None);
    assert_eq!(got, former, "{what}: not the former route's outcome");
    eprintln!(
        "FALLBACK_CASE {what}: route={route:?} outcome={}",
        if got.is_ok() { "value" } else { "error" }
    );
    route
}

/// Fallback `loop_basis`, constructed: a loop momentum without an edge of its
/// own, or a `loop` line naming an edge that carries another momentum; the
/// reason, and the entry's outcome = the former route's.
#[test]
fn a_loop_without_its_own_edge_falls_back() {
    let vakint = engine();
    let mut missing = k5_direct();
    missing.loop_edges[0] = None;
    assert_eq!(
        reduce_open_tensor(&missing),
        Err(OpenTensorFallback::LoopBasis)
    );
    let route = assert_same_outcome(&vakint, &missing, "missing loop edge");
    assert!(
        matches!(route.as_deref(), None | Some("fallback:loop_basis")),
        "{route:?}"
    );
    let value = k5_direct();
    let loops = value.loops();
    let other = value
        .edges
        .iter()
        .find(|edge| {
            !edge.is_external() && edge.momentum[..loops].iter().filter(|c| **c != 0).count() > 1
        })
        .expect("an internal edge carrying a sum of loop momenta")
        .id;
    let mut misnamed = value;
    misnamed.loop_edges[0] = Some(other);
    assert_eq!(
        reduce_open_tensor(&misnamed),
        Err(OpenTensorFallback::LoopBasis)
    );
    let route = assert_same_outcome(&vakint, &misnamed, "misnamed loop edge");
    assert!(
        matches!(route.as_deref(), None | Some("fallback:loop_basis")),
        "{route:?}"
    );
}

/// Fallback `orientation`: every external leg soft. Through the entries the
/// adapter is not reached (vakint's momenta refuse the value first, on both
/// routes, with the same error): asserted here.
#[test]
fn an_unorientable_external_momentum_falls_back() {
    let vakint = engine();
    let mut value = k5_direct();
    for edge in value.edges.iter_mut().filter(|edge| edge.is_external()) {
        edge.momentum.iter_mut().for_each(|c| *c = 0);
    }
    assert_eq!(
        reduce_open_tensor(&value),
        Err(OpenTensorFallback::Orientation)
    );
    let route = assert_same_outcome(&vakint, &value, "soft external legs");
    assert_eq!(route, None, "the adapter is not reached");
}

/// Fallback `contraction`: idenso refuses the routed product (one Lorentz
/// index three times in a factor); the entry's outcome = the former route's.
#[test]
fn a_product_idenso_refuses_falls_back() {
    let vakint = engine();
    let mut value = k5_direct();
    let ambiguous = parse::atom(
        "gammalooprs::Q(4,spenso::mink(gammalooprs::dim,a))*gammalooprs::Q(8,spenso::mink(gammalooprs::dim,a))*gammalooprs::Q(10,spenso::mink(gammalooprs::dim,a))",
    );
    let edge = value
        .edges
        .iter_mut()
        .find(|edge| !edge.is_external())
        .unwrap();
    edge.factor = &edge.factor * &ambiguous;
    assert_eq!(
        reduce_open_tensor(&value),
        Err(OpenTensorFallback::Contraction)
    );
    let route = assert_same_outcome(&vakint, &value, "ambiguous index");
    assert!(
        matches!(route.as_deref(), None | Some("fallback:contraction")),
        "{route:?}"
    );
}

/// Fallback `five_loop`: a value of five loops (here the fixture's basis
/// widened by three unused loop momenta) falls back before anything else.
/// Its value through the entry is not constructible in a unit test (no
/// producer emits a five-loop direct open-slot record, and the unit tests
/// have no five-loop package); the forced `five_loop` fallback above returns
/// the former route's value.
#[test]
fn a_five_loop_value_falls_back_first() {
    let mut value = k5_direct();
    let loops = value.loops();
    for edge in &mut value.edges {
        for _ in 0..3 {
            edge.momentum.insert(loops, 0);
        }
    }
    value.loop_edges.extend([None, None, None]);
    assert_eq!(value.loops(), 5);
    assert_eq!(
        reduce_open_tensor(&value),
        Err(OpenTensorFallback::FiveLoop)
    );
}

/// Fallback `loop_index`: the check that sends a reduced integrand whose loop
/// momentum still carries an index back to the former route.
#[test]
fn a_loop_index_left_is_detected() {
    crate::Vakint::initialize_vakint_symbols();
    let k = parse::atom;
    assert!(has_loop_index(k("k(1,a)*p(1,b)*dot(k(1),k(2))").as_view()));
    assert!(!has_loop_index(
        k("g(a,b)*dot(k(1),k(2))+p(1,a)*p(1,b)*dot(k(1),p(1))").as_view()
    ));
}

#[test]
#[ignore = "needs DIRECT_REPLAY_RECORDS (dumped records outside the repository)"]
fn direct_record_replay() {
    let list = std::env::var("DIRECT_REPLAY_RECORDS").expect("DIRECT_REPLAY_RECORDS");
    let out = std::env::var_os("DIRECT_REPLAY_OUT").map(std::path::PathBuf::from);
    let compare = std::env::var("DIRECT_REPLAY_COMPARE").is_ok_and(|v| v == "1");
    let vakint = engine();
    let paths = std::fs::read_to_string(&list).unwrap();
    let mut total = 0.0;
    for path in paths.lines().filter(|line| !line.trim().is_empty()) {
        let text = std::fs::read_to_string(path).unwrap();
        let value = FactorizedNumerator::from_text(&text).unwrap_or_else(|e| panic!("{path}: {e}"));
        let settings = site_settings(value.loops() as i64);
        let prefactor = value.prefactor.total();
        let started = std::time::Instant::now();
        let got = evaluate_factorized_ports(&vakint, &settings, &value, &prefactor)
            .unwrap_or_else(|e| panic!("{path}: {e}"));
        let seconds = started.elapsed().as_secs_f64();
        total += seconds;
        if compare {
            let former = evaluate_factorized_prepass_route(&vakint, &settings, &value, &prefactor)
                .unwrap_or_else(|e| panic!("{path}: {e}"));
            if let Some(former) = &former {
                assert_eq!(
                    got.vakint_value.to_canonical_string(),
                    former.to_canonical_string(),
                    "{path}: the reduced route changed the value"
                );
            }
            eprintln!(
                "DIRECT_REPLAY_COMPARE path={path} reduced={} equal=true",
                former.is_some()
            );
        }
        eprintln!(
            "DIRECT_REPLAY path={path} seconds={seconds:.3} input_bytes={} value_bytes={}",
            got.input.as_view().get_byte_size(),
            got.value.as_view().get_byte_size()
        );
        if let Some(out) = &out {
            std::fs::create_dir_all(out).unwrap();
            let name = std::path::Path::new(path)
                .file_name()
                .unwrap()
                .to_string_lossy()
                .into_owned();
            std::fs::write(
                out.join(format!("{name}.value")),
                format!(
                    "{}\n{}\n",
                    got.vakint_value.to_canonical_string(),
                    got.value.to_canonical_string()
                ),
            )
            .unwrap();
        }
    }
    eprintln!("DIRECT_REPLAY_TOTAL seconds={total:.3}");
}
