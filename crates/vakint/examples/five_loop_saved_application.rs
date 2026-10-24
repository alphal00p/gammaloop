//! Apply an externally hosted package to a native five-loop DOT fixture.
//! The package's exact entry union decides whether this example is admitted.

#[allow(dead_code)] // The example needs parent definitions, not the full matcher census.
#[path = "../src/topologies/five_loop_data.rs"]
mod five_loop_data;

#[path = "five_loop_cancellation_fixtures/mod.rs"]
mod cancellation_fixtures;
#[path = "saved_application/dot.rs"]
mod saved_application_dot;

use std::collections::{BTreeMap, BTreeSet};
use std::{
    fs::{self, File, OpenOptions},
    io::{Read, Write},
    path::{Path, PathBuf},
    time::Instant,
};

use feynkit_graph::{FeynmanDiagram, IntegralFamily};
use feynkit_kinematics::Kinematics;
use feynkit_model::Model;
use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use symbolica::atom::{Atom, AtomCore, AtomView, FunctionBuilder};
use symbolica::id::Replacement;
use symbolica::prelude::{Integer, Rational};
use vakint::{
    EvaluationMethod, EvaluationOrder, LoopNormalizationFactor, RustRedEvaluationOptions, Vakint,
    VakintExpression, VakintSettings, diagram_integral::DiagramIntegralOptions, vakint_parse,
    vakint_symbol,
};

struct ApplicationLogger;
impl log::Log for ApplicationLogger {
    fn enabled(&self, metadata: &log::Metadata<'_>) -> bool {
        matches!(
            metadata.target(),
            "vakint::saved_application" | "vakint::five_loop_match"
        )
    }
    fn log(&self, record: &log::Record<'_>) {
        if self.enabled(record.metadata()) {
            eprintln!("{}", record.args());
        }
    }
    fn flush(&self) {}
}
static LOGGER: ApplicationLogger = ApplicationLogger;

/// A supplementary cold/warm API measurement, without the demo's later
/// expansion, routing and guard-factorization checks. Output is append-only;
/// the completed cold record is synced before the warm evaluation starts.
struct TimedPair {
    directory: PathBuf,
}

impl TimedPair {
    fn take_argument(
        args: &mut Vec<PathBuf>,
    ) -> Result<Option<PathBuf>, Box<dyn std::error::Error>> {
        let mut directory = None;
        let mut remaining = Vec::new();
        let mut arguments = args.iter();
        while let Some(argument) = arguments.next() {
            if argument.as_os_str() == "--timed-pair-only" {
                if directory.is_some() {
                    return Err("--timed-pair-only must appear exactly once".into());
                }
                directory = Some(
                    arguments
                        .next()
                        .filter(|value| !value.as_os_str().is_empty())
                        .ok_or("--timed-pair-only requires a new OUTPUT_DIR")?
                        .clone(),
                );
            } else {
                remaining.push(argument.clone());
            }
        }
        *args = remaining;
        Ok(directory)
    }

    fn create(directory: PathBuf) -> Result<Self, Box<dyn std::error::Error>> {
        // create_dir, rather than create_dir_all, reserves a new run and rejects
        // accidental reuse. The caller supplies an existing parent directory.
        fs::create_dir(&directory)?;
        let parent = directory
            .parent()
            .filter(|path| !path.as_os_str().is_empty())
            .unwrap_or(Path::new("."));
        File::open(parent)?.sync_all()?;
        Ok(Self { directory })
    }

    fn write(&self, name: &str, bytes: &[u8]) -> Result<Value, Box<dyn std::error::Error>> {
        let path = self.directory.join(name);
        let mut file = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&path)?;
        file.write_all(bytes)?;
        file.sync_all()?;
        File::open(&self.directory)?.sync_all()?;
        Ok(json!({"file":name,"bytes":bytes.len(),"sha256":format!("{:x}",Sha256::digest(bytes))}))
    }

    fn write_json(&self, name: &str, value: &Value) -> Result<Value, Box<dyn std::error::Error>> {
        self.write(name, &serde_json::to_vec_pretty(value)?)
    }

    fn file_identity(path: &Path) -> Result<Value, Box<dyn std::error::Error>> {
        let path = path.canonicalize()?;
        let mut file = File::open(&path)?;
        let mut digest = Sha256::new();
        let mut buffer = [0; 65536];
        loop {
            let count = file.read(&mut buffer)?;
            if count == 0 {
                break;
            }
            digest.update(&buffer[..count]);
        }
        Ok(
            json!({"path":path,"bytes":file.metadata()?.len(),"sha256":format!("{:x}",digest.finalize())}),
        )
    }

    fn cache_statistics(vakint: &Vakint) -> Result<Value, Box<dyn std::error::Error>> {
        let mut loops = BTreeMap::new();
        for count in 1..=5 {
            let statistics = vakint.cache_statistics(count)?;
            loops.insert(count, json!({
                "attached":statistics.attached,
                "statistics_debug":format!("{statistics:?}"),
                "saved_application":statistics.saved_application.as_ref().map(|value|value.to_document()),
            }));
        }
        Ok(json!({
            "persistent":loops,
            "application_debug":format!("{:?}",vakint.rustred_application_cache_statistics()?),
            "application_work_debug":format!("{:?}",vakint.rustred_application_cache_work_statistics()?),
            "owner_preparation_debug":format!("{:?}",vakint.rustred_application_owner_preparation_statistics()?),
            "factor_debug":format!("{:?}",vakint.rustred_application_factor_cache_statistics()?),
            "witness_debug":format!("{:?}",vakint.rustred_application_witness_cache_statistics()?),
        }))
    }

    fn prepare(
        &self,
        application: &Path,
        scalar: &saved_application_dot::ScalarDotInput,
        input: &Atom,
        settings: &VakintSettings,
        vakint: &Vakint,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let application = application.canonicalize()?;
        let root = application
            .ancestors()
            .find(|root| {
                root.join("manifest.json").is_file() && root.join("5L/manifest.json").is_file()
            })
            .ok_or("timed pairs require a local manifested external-data checkout")?;
        let mut manifests = BTreeMap::new();
        for name in [
            "manifest.json",
            "1L/manifest.json",
            "2L/manifest.json",
            "3L/manifest.json",
            "4L/manifest.json",
            "5L/manifest.json",
        ] {
            manifests.insert(name, Self::file_identity(&root.join(name))?);
        }
        let dot = self.write("input.dot", &fs::read(&scalar.path)?)?;
        let input = self.write("input.symbolica", input.to_canonical_string().as_bytes())?;
        self.write_json("input.json", &json!({
            "schema":"vakint.five-loop-connected-timed-pair.v1",
            "scope":"supplementary same-process cold/warm whole public evaluate_rustred_symbolic calls; not a persistent-cache comparison or pure IBP timer",
            "timing_excludes":"pre-evaluation data/application loading, input preparation, result serialization, cache-statistics snapshots and parity checks",
            "owner_preparation_scope":"on-demand owner preparation inside either timed public API call is included",
            "warm_policy":"same engine, settings and input; no cache clear or extra evaluation between the two calls",
            "cache_policy":"with_rustred_application default local-checkout loading, including any available prepared lower-loop caches",
            "dot":dot,"dot_source":Self::file_identity(&scalar.path)?,"powers_native_order":scalar.powers,
            "input":input,"mass_squared":"1","numerator":"1","settings_debug":format!("{settings:?}"),
            "application_directory":application,"data_manifests":manifests,
            "application":Self::file_identity(&application.join("application.json"))?,
            "model_sha256":format!("{:x}",Sha256::digest(include_bytes!("../tests/fixtures/phi_3_4.json"))),
            "compiled_example_sha256":format!("{:x}",Sha256::digest(include_bytes!("five_loop_saved_application.rs"))),
            "compiled_dot_adapter_sha256":format!("{:x}",Sha256::digest(include_bytes!("saved_application/dot.rs"))),
            "compiled_cargo_lock_sha256":format!("{:x}",Sha256::digest(include_bytes!("../../../Cargo.lock"))),
            "executable":Self::file_identity(&std::env::current_exe()?)?,
            "build_rustc":option_env!("VAKINT_CACHE_BUILD_RUSTC"),
            "build_profile":option_env!("VAKINT_CACHE_BUILD_PROFILE"),
            "build_git_revision":option_env!("VAKINT_CACHE_BUILD_GIT_REVISION"),
            "build_rustred_revision":option_env!("VAKINT_CACHE_BUILD_RUSTRED_REVISION"),
            "build_symbolica_revision":option_env!("VAKINT_CACHE_BUILD_SYMBOLICA_REVISION"),
            "build_main_source_snapshot_sha256":option_env!("VAKINT_CACHE_BUILD_MAIN_SOURCE_SNAPSHOT_SHA256"),
            "build_native_source_snapshot_sha256":option_env!("VAKINT_CACHE_BUILD_NATIVE_SOURCE_SNAPSHOT_SHA256"),
            "cache_statistics_before_cold":Self::cache_statistics(vakint)?,
        }))?;
        Ok(())
    }

    fn record(
        &self,
        phase: &str,
        seconds: f64,
        result: &vakint::RustRedSymbolicEvaluation,
        vakint: &Vakint,
    ) -> Result<Value, Box<dyn std::error::Error>> {
        let value = self.write(
            &format!("{phase}.symbolica"),
            result.value.to_canonical_string().as_bytes(),
        )?;
        let guards = result
            .nonzero_conditions
            .iter()
            .map(AtomCore::to_canonical_string)
            .collect::<Vec<_>>();
        let guards = self.write_json(&format!("{phase}-guards.json"), &json!(guards))?;
        let families = result
            .families
            .iter()
            .map(|(symbol, family)| (symbol.to_canonical_string(), family))
            .collect::<Vec<_>>();
        let families = self.write_json(&format!("{phase}-families.json"), &json!(families))?;
        self.write_json(&format!("{phase}.json"), &json!({
            "phase":phase,"completed":true,"whole_public_application_seconds":seconds,
            "value":value,"ordered_nonzero_conditions":guards,"families":families,
            "nonzero_condition_count":result.nonzero_conditions.len(),"family_count":result.families.len(),
            "cache_statistics_after":Self::cache_statistics(vakint)?,
        }))
    }
}

#[cfg(test)]
mod timed_pair_tests {
    use super::TimedPair;
    use serde_json::json;
    use std::{fs, path::PathBuf};

    #[test]
    fn optional_pair_arguments_preserve_the_existing_dot_parser() {
        let original = ["application", "--dot", "fixture.dot", "--powers", "1,2,1"]
            .map(PathBuf::from)
            .to_vec();
        let mut plain = original.clone();
        assert!(TimedPair::take_argument(&mut plain).unwrap().is_none());
        assert_eq!(plain, original);
        for position in 0..=original.len() {
            let mut args = original.clone();
            args.splice(
                position..position,
                ["--timed-pair-only", "output"].map(PathBuf::from),
            );
            assert_eq!(
                TimedPair::take_argument(&mut args).unwrap(),
                Some(PathBuf::from("output"))
            );
            assert_eq!(args, original);
        }
    }

    #[test]
    fn missing_and_duplicate_pair_options_are_errors() {
        for args in [
            vec!["--timed-pair-only"],
            vec!["--timed-pair-only", "one", "--timed-pair-only", "two"],
        ] {
            assert!(
                TimedPair::take_argument(&mut args.into_iter().map(PathBuf::from).collect())
                    .is_err()
            );
        }
    }

    #[test]
    fn completed_records_cannot_be_overwritten() {
        let parent = tempfile::tempdir().unwrap();
        let directory = parent.path().join("pair");
        let pair = TimedPair::create(directory.clone()).unwrap();
        assert!(TimedPair::create(directory.clone()).is_err());
        let cold = json!({"completed":true});
        let identity = pair.write_json("cold.json", &cold).unwrap();
        assert_eq!(
            serde_json::from_slice::<serde_json::Value>(
                &fs::read(directory.join("cold.json")).unwrap()
            )
            .unwrap(),
            cold
        );
        assert_eq!(
            TimedPair::file_identity(&directory.join("cold.json")).unwrap()["sha256"],
            identity["sha256"]
        );
        assert!(
            pair.write_json("cold.json", &json!({"completed":false}))
                .is_err()
        );
        assert!(!directory.join("warm.json").exists());
        assert!(!directory.join("pair.json").exists());
    }
}

fn rerouted_fixture(input: &Atom) -> Result<Atom, Box<dyn std::error::Error>> {
    let routing = [
        ("k(1,idx___)", "k(2,idx___)+k(3,idx___)"),
        ("k(2,idx___)", "-k(1,idx___)"),
        ("k(3,idx___)", "k(3,idx___)"),
        ("k(4,idx___)", "k(5,idx___)"),
        ("k(5,idx___)", "k(4,idx___)"),
    ]
    .iter()
    .map(|(source, target)| {
        Ok(Replacement::new(
            vakint_parse!(source)?.to_pattern(),
            vakint_parse!(target)?.to_pattern(),
        ))
    })
    .collect::<Result<Vec<_>, Box<dyn std::error::Error>>>()?;
    Ok(input
        .replace_multiple(&routing)
        .replace(vakint_parse!("prop(id_,edge(a_,b_),q_,mass_,power_)")?.to_pattern())
        .with(vakint_parse!("prop(100-id_,edge(100-b_,100-a_),-q_,mass_,power_)")?.to_pattern()))
}

fn expanded_propagator_momenta(input: &Atom) -> Atom {
    input.replace_map(|term, _context, out| {
        if let AtomView::Fun(function) = term
            && function.get_symbol() == vakint_symbol!("prop")
        {
            let arguments = function.iter().collect::<Vec<_>>();
            **out = FunctionBuilder::new(function.get_symbol())
                .add_arg(arguments[0])
                .add_arg(arguments[1])
                .add_arg(arguments[2].expand())
                .add_arg(arguments[3])
                .add_arg(arguments[4])
                .finish();
        }
    })
}

fn verify_routing_inverse(
    original: &Atom,
    transformed: &Atom,
) -> Result<(), Box<dyn std::error::Error>> {
    let inverse = [
        ("k(1,idx___)", "-k(2,idx___)"),
        ("k(2,idx___)", "k(1,idx___)-k(3,idx___)"),
        ("k(3,idx___)", "k(3,idx___)"),
        ("k(4,idx___)", "k(5,idx___)"),
        ("k(5,idx___)", "k(4,idx___)"),
    ]
    .iter()
    .map(|(source, target)| {
        Ok(Replacement::new(
            vakint_parse!(source)?.to_pattern(),
            vakint_parse!(target)?.to_pattern(),
        ))
    })
    .collect::<Result<Vec<_>, Box<dyn std::error::Error>>>()?;
    let restored = transformed
        .replace(vakint_parse!("prop(id_,edge(a_,b_),q_,mass_,power_)")?.to_pattern())
        .with(vakint_parse!("prop(100-id_,edge(100-b_,100-a_),-q_,mass_,power_)")?.to_pattern())
        .replace_multiple(&inverse);
    if !(expanded_propagator_momenta(&restored) - expanded_propagator_momenta(original))
        .expand()
        .together()
        .is_zero()
    {
        return Err("the alternate-basis fixture failed exact inverse reconstruction".into());
    }
    println!("inverse_fixture_reconstruction=passed");
    Ok(())
}

fn admitted_denominator_numerator(
    vakint: &Vakint,
    settings: &VakintSettings,
    degree: usize,
) -> Result<Atom, Box<dyn std::error::Error>> {
    let topology = vakint_parse!(format!(
        "topo(I5L_32745(1,0,0,0,0,1,0,0,1,1,0,1,{}))",
        5 + degree,
    ))?;
    let canonical = vakint.to_canonical(settings, topology.as_view(), false)?;
    // Short-form loop names denote Vakint's canonical LMB, not the external
    // package's original k1..k5. Resolve the parent's physical q1 line through
    // its retained propagator ID before constructing the cancelling numerator.
    let square = physical_square(&canonical, 1)?;
    Ok((square - Atom::one()).pow(Atom::num(degree)) * canonical)
}

fn physical_square(canonical: &Atom, id: usize) -> Result<Atom, Box<dyn std::error::Error>> {
    let pattern = vakint_parse!(format!("prop({id},edge(a_,b_),q_,mass_,power_)"))?;
    let matched = canonical
        .pattern_match(&pattern.to_pattern(), None, None)
        .next()
        .ok_or("canonical parent is missing the requested physical line")?;
    let momentum = matched[&vakint_symbol!("q_")].clone();
    Ok(vakint_symbol!("dot").call((&momentum, &momentum)))
}

fn homogeneous_tadpole_control(
    vakint: &Vakint,
    settings: &VakintSettings,
) -> Result<Atom, Box<dyn std::error::Error>> {
    let parent = vakint_parse!("topo(I5L_32745(1,0,0,0,0,1,0,0,1,1,0,1,10))")?;
    let parent = vakint.to_canonical(settings, parent.as_view(), false)?;
    // Both physical rows are absent from the positive support. The common
    // masses cancel in D1-D2, leaving a homogeneous momentum polynomial.
    let numerator = (physical_square(&parent, 1)? - physical_square(&parent, 2)?).pow(5);
    Ok(numerator * parent)
}

type GuardFactorCache = BTreeMap<Atom, BTreeSet<Atom>>;

/// Test-only domain comparison. Native factorization supplies every factor;
/// this only interprets products/powers as nonzero obligations. Returned
/// application guards are never replaced or discarded.
fn guard_domain(
    values: &[Atom],
    cache: &mut GuardFactorCache,
) -> Result<BTreeSet<Atom>, Box<dyn std::error::Error>> {
    fn collect(
        value: AtomView,
        out: &mut BTreeSet<Atom>,
    ) -> Result<(), Box<dyn std::error::Error>> {
        if value.is_zero() {
            return Err("a retained nonzero condition vanished".into());
        }
        match value {
            AtomView::Num(_) => {
                Rational::try_from(value)
                    .map_err(|_| "only proved nonzero rational units may be omitted")?;
            }
            AtomView::Mul(product) => {
                for factor in product.iter() {
                    collect(factor, out)?;
                }
            }
            AtomView::Pow(power) if Integer::try_from(power.get_exp()).is_ok() => {
                if Integer::try_from(power.get_exp())? != 0 {
                    collect(power.get_base(), out)?;
                }
            }
            _ => {
                out.insert(value.to_owned());
            }
        }
        Ok(())
    }
    let mut result = BTreeSet::new();
    for condition in values {
        let condition = condition.expand().together();
        if !cache.contains_key(&condition) {
            let mut factors = BTreeSet::new();
            collect(condition.factor().as_view(), &mut factors)?;
            cache.insert(condition.clone(), factors);
        }
        result.extend(cache[&condition].iter().cloned());
    }
    Ok(result)
}

fn require_same_value(
    left: &vakint::RustRedSymbolicEvaluation,
    right: &vakint::RustRedSymbolicEvaluation,
    cache: &mut GuardFactorCache,
    guard_discrepancies: &mut Vec<String>,
) -> Result<(), Box<dyn std::error::Error>> {
    if !(&left.value - &right.value).expand().together().is_zero()
        || left.families != right.families
    {
        return Err("equivalent input changed the exact qualified master expression".into());
    }
    let left_guards = guard_domain(&left.nonzero_conditions, cache)?;
    let right_guards = guard_domain(&right.nonzero_conditions, cache)?;
    if left_guards != right_guards {
        let bounded = |values: Vec<&Atom>| {
            values
                .iter()
                .take(8)
                .map(|value| {
                    value
                        .to_canonical_string()
                        .chars()
                        .take(512)
                        .collect::<String>()
                })
                .collect::<Vec<_>>()
        };
        let left_only = left_guards.difference(&right_guards).collect::<Vec<_>>();
        let right_only = right_guards.difference(&left_guards).collect::<Vec<_>>();
        let diagnostic = format!(
            "equivalent value has distinct native-factored guard domains: left_only_count={} right_only_count={} left_only_prefix={:?} right_only_prefix={:?}",
            left_only.len(),
            right_only.len(),
            bounded(left_only),
            bounded(right_only),
        );
        eprintln!("guard_domain_discrepancy={diagnostic}");
        guard_discrepancies.push(diagnostic);
        return Ok(());
    }
    println!(
        "native_factored_guard_domain_parity=passed factors={}",
        left_guards.len()
    );
    Ok(())
}

fn check_independent_tadpoles(
    result: &vakint::RustRedSymbolicEvaluation,
) -> Result<(), Box<dyn std::error::Error>> {
    use rustred::{
        family::IntegralKey,
        foundry::artifact::ClosedArtifact,
        reduction::{Reducer, ReductionLimits},
    };
    let artifact = ClosedArtifact::decode_durable(include_bytes!("../data/1L/ibps/rules.rrbin"))?;
    let label = result
        .families
        .iter()
        .find_map(|(label, fingerprint)| {
            (fingerprint == artifact.family_fingerprint()).then_some(label)
        })
        .ok_or("five-tadpole result did not use the independently shipped one-loop family")?;
    let dimension = artifact.family().dimension().to_expression();
    let physical_dimension = vakint_parse!("4-2*eps")?;
    let mut expected = Atom::one();
    for power in [1, 1, 1, 1, 5] {
        let (reduced, conditions, _) = Reducer::reduce_artifact_unit_mass_guarded(
            &artifact,
            &IntegralKey::try_new([power])?,
            ReductionLimits::default(),
        )?;
        let mut factor = Atom::Zero;
        for (integral, coefficient) in reduced.terms() {
            let master = vakint_parse!(format!(
                "RustRedMaster({},{})",
                label.to_canonical_string(),
                integral.powers()[0]
            ))?;
            factor += coefficient
                .to_expression()
                .replace(dimension.to_pattern())
                .with(physical_dimension.to_pattern())
                * master;
        }
        for condition in conditions {
            let condition = condition
                .to_expression()
                .replace(dimension.to_pattern())
                .with(physical_dimension.to_pattern());
            if !result.nonzero_conditions.contains(&condition) {
                return Err(format!("independent factor condition was lost: {condition}").into());
            }
        }
        expected *= factor;
    }
    if !(&result.value - expected).expand().together().is_zero() {
        return Err(
            "five-loop result differs from five independent guarded one-loop reductions".into(),
        );
    }
    println!("independent_five_tadpole_parity=passed");
    Ok(())
}

/// Native five-tadpole chain, with four zero-flow bridges authenticated from
/// its stored routing. Only one bridge is restored in the second input.
type MassiveBridgeFixture = (String, Atom, Atom, Atom);

fn massive_bridge_fixtures() -> Result<Vec<MassiveBridgeFixture>, Box<dyn std::error::Error>> {
    let diagram = FeynmanDiagram::from_dot(
        Model::from_json(include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/fixtures/phi_3_4.json"
        )))
        .expect("valid phi_3_4 fixture"),
        "digraph five_tadpoles {\
         a -> a [particle=\"phi\"]; a -> b [particle=\"phi\"];\
         b -> b [particle=\"phi\"]; b -> c [particle=\"phi\"];\
         c -> c [particle=\"phi\"]; c -> d [particle=\"phi\"];\
         d -> d [particle=\"phi\"]; d -> e [particle=\"phi\"];\
         e -> e [particle=\"phi\"]; }",
    )?;
    assert_eq!(diagram.loop_count(), 5);
    let original = diagram.propagator_family(&Kinematics::new())?;
    let mut edges = diagram.edges().collect::<Vec<_>>();
    edges.sort_by_key(|(id, _, _)| *id);
    let mut powers = vec![Integer::from(0); original.denominators().len()];
    let mut bridges = Vec::new();
    let mut tadpoles = Vec::new();
    for (slot, (id, _, _)) in edges.iter().enumerate() {
        if diagram.loop_momentum_basis().edge_signatures[id]
            .loops
            .integer_coefficients()
            .iter()
            .all(|value| *value == 0)
        {
            bridges.push(slot);
        } else {
            tadpoles.push(slot);
            powers[slot] = Integer::from(1);
        }
    }
    assert_eq!(bridges.len(), 4);
    assert_eq!(tadpoles.len(), 5);
    powers[tadpoles[4]] = Integer::from(5);
    let symbols = original.denominators()[bridges[0]].get_all_symbols(true);
    assert_eq!(symbols.len(), 1);
    let mass = Atom::var(*symbols.iter().next().unwrap());
    let mut fixtures = Vec::new();
    for mass_value in [Some(1), Some(2), None] {
        let substitutions = mass_value
            .map(|value| (mass.clone(), Atom::num(value)))
            .into_iter()
            .collect::<Vec<_>>();
        let replacements = substitutions
            .iter()
            .map(|(source, target)| Replacement::new(source.clone(), target.clone()))
            .collect::<Vec<_>>();
        let family = IntegralFamily::new(
            original.loop_momenta().to_vec(),
            vec![],
            original
                .denominators()
                .iter()
                .map(|denominator| denominator.replace_multiple(&replacements))
                .collect(),
            original.kinematics(),
        )?;
        let convert = |bridge_power: i64| -> Result<Atom, Box<dyn std::error::Error>> {
            let mut powers = powers.clone();
            powers[bridges[0]] = Integer::from(bridge_power);
            Ok(Atom::from(VakintExpression::from_diagram(
                &diagram,
                &family,
                &Atom::one(),
                &DiagramIntegralOptions {
                    powers: Some(powers),
                    parameter_substitutions: substitutions.clone(),
                    ..Default::default()
                },
            )?))
        };
        let baseline = convert(0)?;
        let bridge = convert(1)?;
        let denominator = family.denominators()[bridges[0]].clone();
        if !(&bridge - &baseline * denominator.clone().pow(-1))
            .expand()
            .together()
            .is_zero()
        {
            return Err("native massive bridge extraction changed the five-loop input".into());
        }
        fixtures.push((
            format!("mass_{mass_value:?}"),
            baseline,
            bridge,
            denominator,
        ));
    }
    Ok(fixtures)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    log::set_logger(&LOGGER)?;
    log::set_max_level(log::LevelFilter::Debug);
    let mut args = std::env::args_os()
        .skip(1)
        .map(PathBuf::from)
        .collect::<Vec<_>>();
    let pair_directory = TimedPair::take_argument(&mut args)?;
    let scalar_dot = saved_application_dot::ScalarDotInput::from_arguments(&args)?;
    if pair_directory.is_some() && scalar_dot.is_none() {
        return Err("--timed-pair-only requires explicit --dot PATH --powers CSV input".into());
    }
    if args.is_empty() || (scalar_dot.is_none() && !(1..=2).contains(&args.len())) {
        return Err("usage: five_loop_saved_application APPLICATION_DIRECTORY [VACUUM_DOT|--admitted-tadpoles|--admitted-rank11|--admitted-suite|--dot PATH --powers CSV] [--timed-pair-only OUTPUT_DIR]".into());
    }
    let timed_pair = pair_directory.map(TimedPair::create).transpose()?;
    let mut settings = VakintSettings {
        epsilon_symbol: "eps".into(),
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions {
                master_mode: vakint::RustRedMasterMode::Retain,
            },
        )]),
        form_exe_path: "/five-loop-example-must-not-use-form".into(),
        integral_normalization_factor: LoopNormalizationFactor::Custom("1".into()),
        ..Default::default()
    };
    let catalog_started = Instant::now();
    let vakint = Vakint::new()?;
    println!(
        "vakint_catalogue_construction_seconds={:.6}",
        catalog_started.elapsed().as_secs_f64()
    );
    println!(
        "five_loop_catalogue_before_load={:?}",
        vakint.topologies.five_loop_catalogue_statistics()
    );
    if args.len() == 1 && args[0].as_os_str() == "--catalogue-only" {
        return Ok(());
    }
    if args.len() == 1 && args[0].as_os_str() == "--native-fixtures" {
        let fixtures = cancellation_fixtures::fixtures(&vakint)?;
        cancellation_fixtures::verify_fixture_shapes(&fixtures)?;
        println!(
            "native_positive_line_cancellation_fixtures={}",
            fixtures.len()
        );
        println!(
            "native_massive_bridge_fixtures={}",
            massive_bridge_fixtures()?.len()
        );
        return Ok(());
    }
    if args.len() == 1 && args[0].as_os_str() == "--routing-fixtures" {
        return Err("routing fixtures require APPLICATION_DIRECTORY --routing-fixtures; --native-fixtures needs no package".into());
    }
    let started = Instant::now();
    let vakint = vakint.with_rustred_application(&mut settings, &args[0])?;
    println!(
        "load_and_prepare_seconds={:.6}",
        started.elapsed().as_secs_f64()
    );
    println!(
        "five_loop_catalogue_after_load={:?}",
        vakint.topologies.five_loop_catalogue_statistics()
    );
    println!(
        "products={:?}",
        vakint.rustred_application_product_statistics()?
    );
    println!(
        "product_charts={:?}",
        vakint.rustred_application_product_chart_statistics()?
    );
    println!(
        "owner_preparation={:?}",
        vakint.rustred_application_owner_preparation_statistics()?
    );
    if args.len() == 1 {
        return Ok(());
    }
    let fixtures_started = Instant::now();
    let cancellation_controls = if args.get(1).is_some_and(|arg| {
        arg.as_os_str() == "--admitted-suite" || arg.as_os_str() == "--routing-fixtures"
    }) {
        let fixtures = cancellation_fixtures::fixtures(&vakint)?;
        cancellation_fixtures::verify_fixture_shapes(&fixtures)?;
        println!(
            "native_positive_line_cancellation_fixtures={}",
            fixtures.len()
        );
        fixtures
    } else {
        Vec::new()
    };
    println!(
        "native_fixture_preparation_seconds={:.6}",
        fixtures_started.elapsed().as_secs_f64()
    );
    if args[1].as_os_str() == "--routing-fixtures" {
        let bridges = massive_bridge_fixtures()?;
        println!("native_massive_bridge_fixtures={}", bridges.len());
        let mut scalar_settings = settings.clone();
        scalar_settings.use_dot_product_notation = true;
        for degree in [5, 11] {
            let original = admitted_denominator_numerator(&vakint, &settings, degree)?;
            let transformed = rerouted_fixture(&original)?;
            let indexed = Vakint::convert_from_dot_notation(transformed.as_view());
            println!("scope_degree={degree} stage=first_canonicalization");
            let direct = vakint.to_canonical(&settings, transformed.as_view(), false)?;
            let indexed = vakint.to_canonical(&settings, indexed.as_view(), false)?;
            println!("scope_degree={degree} stage=direct_dot_tensor_prepass");
            let direct = vakint.tensor_reduce(&scalar_settings, direct.as_view())?;
            println!("scope_degree={degree} stage=raw_indexed_tensor_prepass");
            let indexed = vakint.tensor_reduce(&scalar_settings, indexed.as_view())?;
            if !(&direct - indexed).expand().together().is_zero() {
                return Err(format!(
                    "initial canonicalization broke indexed scalar scope at degree{degree}"
                )
                .into());
            }
            let expanded = original.expand();
            println!("scope_degree={degree} stage=expanded_internal_scalar_tensor_prepass");
            let first = vakint.tensor_reduce(&scalar_settings, expanded.as_view())?;
            println!("scope_degree={degree} stage=internal_scalar_output_restore");
            let restored = vakint.tensor_reduce(&scalar_settings, first.as_view())?;
            if !(restored - expanded).expand().together().is_zero() {
                return Err(format!(
                    "internal tensor output lost expanded scalar scope at degree{degree}"
                )
                .into());
            }
            println!("full_indexed_and_expanded_scalar_scope_degree={degree}:passed");
        }
        let input = admitted_denominator_numerator(&vakint, &settings, 5)?;
        let canonical = vakint.to_canonical(&settings, input.as_view(), false)?;
        let transformed = rerouted_fixture(&canonical)?;
        verify_routing_inverse(&canonical, &transformed)?;
        let rematched = vakint.to_canonical(&settings, transformed.as_view(), false)?;
        let tensor_reduced = vakint.tensor_reduce(&settings, rematched.as_view())?;
        let scalar = vakint.tensor_reduce(&scalar_settings, tensor_reduced.as_view())?;
        if !(&scalar - &rematched).expand().together().is_zero() {
            return Err(
                "the native closed-scalar bridge changed the routed five-loop numerator".into(),
            );
        }
        println!("native_closed_scalar_scope_restoration=passed");
        for (stage, expression) in [
            ("constructed", input),
            ("canonical", canonical),
            ("transformed", transformed),
            ("rematched", rematched),
            ("tensor_reduced", tensor_reduced),
            ("scalar_restored", scalar),
        ] {
            println!(
                "stage={stage} expression={}",
                expression.to_canonical_string()
            );
        }
        return Ok(());
    }
    let input = if let Some(input) = &scalar_dot {
        let prepared = input.prepare()?;
        if timed_pair.is_some() && prepared.diagram.loop_count() != 5 {
            return Err("timed pairs require a five-loop DOT fixture".into());
        }
        println!(
            "scalar_dot_loops={} physical_denominators={} powers={:?} mass_squared=1 numerator=1",
            prepared.diagram.loop_count(),
            prepared.family.denominators().len(),
            input.powers
        );
        prepared.expression
    } else if args[1].as_os_str() == "--admitted-tadpoles"
        || args[1].as_os_str() == "--admitted-suite"
    {
        // Generation54 point00162: 15-slot powers
        // [0,0,0,0,1,0,0,1,1,0,0,1,0,0,5]. The explicit parent physical
        // slots place q12 and q15 at positions 11 and 12 of this short form.
        vakint_parse!("topo(I5L_32745(1,0,0,0,0,1,0,0,1,1,0,1,5))")?
    } else if args[1].as_os_str() == "--admitted-rank11" {
        // Generation54 point00000. This is an inverse-denominator numerator,
        // not a claim that all homogeneous momentum polynomials are admitted.
        admitted_denominator_numerator(&vakint, &settings, 11)?
    } else {
        let diagram = FeynmanDiagram::from_dot(
            Model::from_json(include_str!(concat!(
                env!("CARGO_MANIFEST_DIR"),
                "/tests/fixtures/phi_3_4.json"
            )))
            .expect("valid phi_3_4 fixture"),
            &std::fs::read_to_string(&args[1])?,
        )?;
        if diagram.loop_count() != 5 {
            return Err("the DOT fixture must have five loops".into());
        }
        let family = diagram.propagator_family(&Kinematics::new())?;
        if family.denominators().len() != 12 {
            return Err("this degree-ten example requires a maximal twelve-line fixture".into());
        }
        let mut powers = vec![1.into(); family.denominators().len()];
        powers[0] = 2.into();
        powers[1] = 2.into();
        // The common masses cancel: momentum degree ten, scalar signed difference 9.
        // Each lowered/routed endpoint must still pass the original package union.
        let numerator = (&family.denominators()[0] - &family.denominators()[1]).pow(Atom::num(5));
        let input = VakintExpression::from_diagram(
            &diagram,
            &family,
            &numerator,
            &DiagramIntegralOptions {
                powers: Some(powers),
                ..Default::default()
            },
        )?;
        Atom::from(input)
    };
    if let Some(pair) = &timed_pair {
        pair.prepare(
            &args[0],
            scalar_dot.as_ref().unwrap(),
            &input,
            &settings,
            &vakint,
        )?;
    }
    let mut inputs = vec![input];
    let mut guard_factors = GuardFactorCache::new();
    let mut guard_discrepancies = Vec::new();
    if args[1].as_os_str() == "--admitted-suite" {
        inputs.push(admitted_denominator_numerator(&vakint, &settings, 5)?);
        inputs.push(admitted_denominator_numerator(&vakint, &settings, 11)?);
    }
    for (case, input) in inputs.into_iter().enumerate() {
        println!("case={case}");
        let started = Instant::now();
        let first = vakint.evaluate_rustred_symbolic(&settings, input.as_view())?;
        let first_seconds = started.elapsed().as_secs_f64();
        println!("first_application_seconds={:.6}", first_seconds);
        let cold_record = timed_pair
            .as_ref()
            .map(|pair| pair.record("cold", first_seconds, &first, &vakint))
            .transpose()?;
        let started = Instant::now();
        let warm = match vakint.evaluate_rustred_symbolic(&settings, input.as_view()) {
            Ok(value) => value,
            Err(error) => {
                if let Some(pair) = &timed_pair {
                    pair.write_json("warm-error.json", &json!({"completed":false,"error":error.to_string(),"whole_public_application_seconds":started.elapsed().as_secs_f64()}))?;
                }
                return Err(error.into());
            }
        };
        let warm_seconds = started.elapsed().as_secs_f64();
        println!("warm_application_seconds={:.6}", warm_seconds);
        let exact_parity = first.value == warm.value
            && first.nonzero_conditions == warm.nonzero_conditions
            && first.families == warm.families;
        if let Some(pair) = &timed_pair {
            let warm_record = pair.record("warm", warm_seconds, &warm, &vakint)?;
            pair.write_json(
                "pair.json",
                &json!({
                    "schema":"vakint.five-loop-connected-timed-pair.v1","completed":exact_parity,
                    "cold_record":cold_record,"warm_record":warm_record,
                    "exact_value_parity":first.value == warm.value,
                    "ordered_full_guard_parity":first.nonzero_conditions == warm.nonzero_conditions,
                    "family_identity_parity":first.families == warm.families,
                    "additional_demo_checks_run":false,
                }),
            )?;
        }
        if !exact_parity {
            return Err("warm application changed the exact result".into());
        }
        if timed_pair.is_some() {
            return Ok(());
        }
        println!("phase=expanded_input");
        let expanded = vakint.evaluate_rustred_symbolic(&settings, input.expand().as_view())?;
        require_same_value(
            &first,
            &expanded,
            &mut guard_factors,
            &mut guard_discrepancies,
        )?;
        let standard = vakint.evaluate(&settings, input.expand().as_view())?;
        if !(&first.value - standard).expand().together().is_zero() {
            return Err("standard evaluation and detailed scalar handoff disagree".into());
        }
        println!("standard_and_detailed_expanded_value_parity=passed");
        let long_form = vakint.to_canonical(&settings, input.as_view(), false)?;
        let transformed = rerouted_fixture(&long_form)?;
        println!("canonical_fixture={}", long_form.to_canonical_string());
        println!("rerouted_fixture={}", transformed.to_canonical_string());
        verify_routing_inverse(&long_form, &transformed)?;
        println!("phase=rerouted_input");
        let rerouted = vakint.evaluate_rustred_symbolic(&settings, transformed.as_view())?;
        require_same_value(
            &first,
            &rerouted,
            &mut guard_factors,
            &mut guard_discrepancies,
        )?;
        println!("expanded_numerator_and_simultaneous_routing_value_parity=passed");
        let indexed_transformed = Vakint::convert_from_dot_notation(transformed.as_view());
        let indexed = vakint.evaluate_rustred_symbolic(&settings, indexed_transformed.as_view())?;
        require_same_value(
            &first,
            &indexed,
            &mut guard_factors,
            &mut guard_discrepancies,
        )?;
        println!("indexed_simultaneous_routing_value_parity=passed");
        println!("cache={:?}", vakint.rustred_application_cache_statistics()?);
        println!(
            "cache_work={:?}",
            vakint.rustred_application_cache_work_statistics()?
        );
        println!(
            "owner_preparation={:?}",
            vakint.rustred_application_owner_preparation_statistics()?
        );
        println!(
            "factor_cache={:?}",
            vakint.rustred_application_factor_cache_statistics()?
        );
        println!(
            "witness_cache={:?}",
            vakint.rustred_application_witness_cache_statistics()?
        );
        println!("nonzero_condition_count={}", first.nonzero_conditions.len());
        let output = first.value.to_canonical_string();
        println!("output_bytes={}", output.len());
        println!("{output}");
        if case == 0
            && (args[1].as_os_str() == "--admitted-tadpoles"
                || args[1].as_os_str() == "--admitted-suite")
        {
            check_independent_tadpoles(&first)?;
            let spectator = vakint_parse!("dot(p(1),p(2))*g(71,72)")?;
            let decorated_input = &input * &spectator;
            let decorated =
                vakint.evaluate_rustred_symbolic(&settings, decorated_input.as_view())?;
            let mut dot_settings = settings.clone();
            dot_settings.use_dot_product_notation = true;
            let mut decorated_dots =
                vakint.evaluate_rustred_symbolic(&dot_settings, decorated_input.as_view())?;
            decorated_dots.value =
                Vakint::convert_from_dot_notation(decorated_dots.value.as_view());
            require_same_value(
                &decorated,
                &decorated_dots,
                &mut guard_factors,
                &mut guard_discrepancies,
            )?;
            let expected = Vakint::convert_from_dot_notation((&first.value * spectator).as_view());
            if !(&decorated.value - expected).expand().together().is_zero() {
                return Err(
                    "internal compact handoff changed external scalar/tensor spectators".into(),
                );
            }
            println!("external_spectator_presentation_parity=passed");
            let massive = vakint_parse!("topo(I5L_32745(4,0,0,0,0,1,0,0,1,1,0,1,5))")?;
            let massive = vakint.evaluate_rustred_symbolic(&settings, massive.as_view())?;
            let expected = &first.value * vakint_parse!("4^(1-5*eps)")?;
            if !(&massive.value - expected).expand().together().is_zero() {
                return Err("nonunit common mass changed the five-loop scaling law".into());
            }
            println!("nonunit_mass_scaling=passed");
        }
        // Settings remain the sole enable path, including on a clone with a warm cache.
        let clone = vakint.clone();
        let mut disabled = settings.clone();
        disabled.rustred_artifact_path = None;
        if clone
            .evaluate_rustred_symbolic(&disabled, input.as_view())
            .is_ok()
        {
            return Err("removing the settings path did not disable a cached application".into());
        }
        let mut tightened = settings.clone();
        tightened
            .rustred_application_options
            .binary_limits
            .max_program_bytes = 0;
        if clone
            .evaluate_rustred_symbolic(&tightened, input.as_view())
            .is_ok()
        {
            return Err("a warm cache bypassed changed application options".into());
        }
        vakint.clear_rustred_application_cache()?;
        if vakint
            .rustred_application_witness_cache_statistics()?
            .retained_entries
            != 0
        {
            return Err("clearing the application left frontend witnesses cached".into());
        }
        if vakint
            .rustred_application_factor_cache_statistics()?
            .cached_factors
            != 0
        {
            return Err("clearing the application left lower-factor entries cached".into());
        }
        println!("cached_path_and_options_guards=passed");
    }
    if args[1].as_os_str() == "--admitted-suite" {
        for fixture in cancellation_controls {
            let cancelled =
                vakint.evaluate_rustred_symbolic(&settings, fixture.cancelled.as_view())?;
            let pinched = vakint.evaluate_rustred_symbolic(&settings, fixture.pinched.as_view())?;
            require_same_value(
                &cancelled,
                &pinched,
                &mut guard_factors,
                &mut guard_discrepancies,
            )?;
            let expanded = vakint
                .evaluate_rustred_symbolic(&settings, fixture.cancelled.expand().as_view())?;
            require_same_value(
                &cancelled,
                &expanded,
                &mut guard_factors,
                &mut guard_discrepancies,
            )?;
            println!(
                "positive_line_cancellation_value_parity={}:passed expected_key={:?} mass_squared={} positive_lines={:?} independent_dot_bytes={}",
                fixture.label,
                fixture.expected_key.powers(),
                fixture.mass_squared,
                fixture.positive_cancellations,
                fixture.pinched_dot.len(),
            );
        }
        for (label, baseline, bridge, denominator) in massive_bridge_fixtures()? {
            let baseline = vakint.evaluate_rustred_symbolic(&settings, baseline.as_view())?;
            let actual = vakint.evaluate_rustred_symbolic(&settings, bridge.as_view())?;
            let expected = &baseline.value * denominator.clone().pow(-1);
            if !(&actual.value - expected).expand().together().is_zero()
                || actual.families != baseline.families
            {
                return Err(format!("five-loop massive bridge parity failed: {label}").into());
            }
            if !denominator.get_all_symbols(true).is_empty()
                && !actual.nonzero_conditions.iter().any(|condition| {
                    (condition - &denominator).expand().together().is_zero()
                        || (condition + &denominator).expand().together().is_zero()
                })
            {
                return Err("the symbolic massive bridge lost its original mass pole guard".into());
            }
            println!("massive_bridge_parity={label}:passed");
        }
        // This particular scalar twelve-line corner (D=12) is outside the
        // pinned generation-54 entry union; this is not a global D cutoff.
        // Loading owners or cached terminals does not upgrade that coverage.
        // Vakint now warns and preserves any residual exact integrals.
        let unsupported = vakint_parse!("topo(I5L_32745(1,1,1,1,1,1,1,1,1,1,1,1,1))")?;
        let partial = vakint.evaluate_rustred_symbolic(&settings, unsupported.as_view())?;
        if partial.value.is_zero() || partial.families.is_empty() {
            return Err("the out-of-coverage D=12 control lost its integral identity".into());
        }
        println!("unsupported_entry_after_cache_reuse=retained_or_safely_reduced");
        println!("phase=homogeneous_degree10_tadpole_control");
        let homogeneous = homogeneous_tadpole_control(&vakint, &settings)?;
        let first = vakint.evaluate_rustred_symbolic(&settings, homogeneous.as_view())?;
        let warm = vakint.evaluate_rustred_symbolic(&settings, homogeneous.as_view())?;
        require_same_value(&first, &warm, &mut guard_factors, &mut guard_discrepancies)?;
        let expanded =
            vakint.evaluate_rustred_symbolic(&settings, homogeneous.expand().as_view())?;
        require_same_value(
            &first,
            &expanded,
            &mut guard_factors,
            &mut guard_discrepancies,
        )?;
        println!("homogeneous_degree10_tadpole_value_parity=passed");
    }
    if !guard_discrepancies.is_empty() {
        return Err(format!(
            "{} native-factored guard-domain discrepancies require provenance audit",
            guard_discrepancies.len()
        )
        .into());
    }
    println!("guard_domain_parity_all_comparisons=passed");
    Ok(())
}
