use std::{
    collections::BTreeMap,
    error::Error,
    fs::{self, File},
    io::{BufReader, BufWriter, Read, Seek, Write},
    path::{Path, PathBuf},
    time::Instant,
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Indeterminate, Symbol},
    evaluate::{FunctionMap, OptimizationSettings},
    function,
    state::State,
    symbol,
};

// Seed the captured prefix before Symbolica's special-function initializer
// allocates dynamic IDs. Let the native initializer restore its own callbacks;
// main then imports the complete registry and requires an identity state map.
// The fixed sixteen builtins are already initialized.
// This uses the public, re-entry-safe initializer protocol; no State reset or
// Atom remapping is needed, and the dependency's evaluator code stays unchanged.
symbolica::_inventory::submit! {
    symbolica::state::StateInitializer::new("exact_builder_registry", || {
        let args = std::env::args().skip(1).collect::<Vec<_>>();
        if args.first().is_some_and(|a| a == "--exact-builder") {
            assert_eq!(args.len(), 2);
            let path = Path::new(&args[1]).join("state.bin");
            let mut input = BufReader::new(File::open(path).expect("Captured registry"));
            let mut header = [0; 14];
            input.read_exact(&mut header).unwrap();
            let count = u64::from_le_bytes(header[6..].try_into().unwrap());
            for index in 0..count {
                let position = input.stream_position().unwrap();
                let mut size = [0; 4];
                input.read_exact(&mut size).unwrap();
                let mut name = vec![0; u32::from_le_bytes(size) as usize];
                input.read_exact(&mut name).unwrap();
                input.seek(std::io::SeekFrom::Start(position)).unwrap();
                // First symbol of pinned Symbolica's GeometricSymbols group.
                if name == b"symbolica::tan" {
                    return;
                }
                let symbol = Symbol::import(&mut input).unwrap();
                assert_eq!(symbol.get_id() as u64, index, "Captured prefix changed IDs");
            }
            panic!("Captured registry lacks the native special-function boundary");
        }
    }, &["symbolica"])
}

// Audit the pinned FunctionMap binary format without normalizing its Atoms.
// Sorting only the two serialized hash maps preserves definition IDs, argument
// order and every raw byte; it does not sort or rewrite symbolic expressions.
#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
struct RawAtom(Vec<u8>);
bincode::impl_borrow_decode!(RawAtom);

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
struct RawSymbol(Vec<u8>);
bincode::impl_borrow_decode!(RawSymbol);

impl<C> bincode::Decode<C> for RawSymbol {
    fn decode<D: bincode::de::Decoder<Context = C>>(
        decoder: &mut D,
    ) -> Result<Self, bincode::error::DecodeError> {
        let mut bytes = Vec::<u8>::decode(decoder)?;
        // Symbol::export ends with callback presence, which import explicitly
        // cannot restore. This flag includes display-only callbacks. Preserve
        // all names, namespaces, attributes, tags, aliases and user data; the
        // manifest separately records missing application callbacks.
        assert!(matches!(bytes.pop(), Some(0 | 1)));
        Ok(Self(bytes))
    }
}

impl<C> bincode::Decode<C> for RawAtom {
    fn decode<D: bincode::de::Decoder<Context = C>>(
        decoder: &mut D,
    ) -> Result<Self, bincode::error::DecodeError> {
        use bincode::de::read::Reader;
        let reader = decoder.reader();
        let mut header = [0; 9];
        reader.read(&mut header)?;
        assert_eq!(header[0], 0);
        let size = u64::from_le_bytes(header[1..].try_into().unwrap()) as usize;
        let mut bytes = vec![0; size];
        reader.read(&mut bytes)?;
        Ok(Self(bytes))
    }
}

#[derive(bincode::Decode, Debug, PartialEq, Eq)]
enum RawIndeterminate {
    Symbol(RawSymbol, ([u8; 16], u8)),
    Function(RawSymbol, RawAtom),
}

type RawFunctionMap = (
    BTreeMap<(RawSymbol, Vec<RawAtom>), (usize, usize, Vec<RawIndeterminate>, RawAtom)>,
    BTreeMap<RawSymbol, usize>,
);

fn stage(started: Instant, name: &str) {
    let status = fs::read_to_string("/proc/self/status").unwrap_or_default();
    let memory = |key: &str| {
        status
            .lines()
            .find_map(|line| {
                line.strip_prefix(key)?
                    .split_whitespace()
                    .next()?
                    .parse::<u64>()
                    .ok()
            })
            .map_or_else(|| "null".to_owned(), |kib| (kib * 1024).to_string())
    };
    println!(
        "{{\"stage\":\"{name}\",\"elapsed_seconds\":{},\"pid\":{},\"rss_bytes\":{},\"vmhwm_bytes\":{}}}",
        started.elapsed().as_secs_f64(),
        std::process::id(),
        memory("VmRSS:"),
        memory("VmHWM:")
    );
    std::io::stdout().flush().unwrap();
}

fn import(path: &Path) -> Result<Atom, std::io::Error> {
    Atom::import(&mut BufReader::new(File::open(path)?), None)
}

fn unpack(atom: AtomView<'_>, wrapper: Symbol) -> Vec<Atom> {
    let AtomView::Fun(f) = atom else {
        panic!("Expected a portable-context function wrapper");
    };
    assert_eq!(f.get_symbol(), wrapper, "Unexpected context wrapper");
    f.iter().map(|a| a.to_owned()).collect()
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    let started = Instant::now();
    stage(started, "start");
    if args.first().is_some_and(|a| a == "--exact-builder") {
        assert_eq!(args.len(), 2, "usage: --exact-builder CAPTURE_DIRECTORY");
        let directory = Path::new(&args[1]);
        let manifest: serde_json::Value =
            serde_json::from_reader(File::open(directory.join("manifest.json"))?)?;
        assert_eq!(manifest["format"], "symbolica-exact-builder-v1");
        assert_eq!(manifest["compile"], false);
        // Import the complete registry first, before registering any user symbol.
        // Refuse remapping: the raw expression must retain its original IDs.
        println!("registered_symbols={}", State::symbol_iter().count());
        let state_map = State::import(
            &mut BufReader::new(File::open(directory.join("state.bin"))?),
            None,
        )?;
        assert!(state_map.is_empty(), "Symbol registry is not identical");
        stage(started, "exact_registry_ready");
        let config = bincode::config::standard();
        let parameter_bytes = fs::read(directory.join("params.bin"))?;
        let (raw_params, consumed): (Vec<Vec<u8>>, _) =
            bincode::decode_from_slice(&parameter_bytes, config)?;
        assert_eq!(consumed, parameter_bytes.len());
        assert_eq!(raw_params.len() as u64, manifest["parameter_count"]);
        let params = raw_params
            .iter()
            .map(|p| AtomView::from(p))
            .collect::<Vec<_>>();
        let map_bytes = fs::read(directory.join("function_map.bin"))?;
        let (fn_map, consumed): (FunctionMap, _) =
            bincode::decode_from_slice_with_context(&map_bytes, config, state_map)?;
        assert_eq!(consumed, map_bytes.len());
        let (before, consumed): (RawFunctionMap, _) =
            bincode::decode_from_slice(&map_bytes, config)?;
        assert_eq!(consumed, map_bytes.len());
        let restored_map = bincode::encode_to_vec(&fn_map, config)?;
        let (after, consumed): (RawFunctionMap, _) =
            bincode::decode_from_slice(&restored_map, config)?;
        assert_eq!(consumed, restored_map.len());
        assert_eq!(before, after, "FunctionMap changed during loading");
        stage(started, "exact_function_map_verified");
        let settings_value: serde_json::Value =
            serde_json::from_reader(File::open(directory.join("optimization.json"))?)?;
        let settings: OptimizationSettings = serde_json::from_value(settings_value.clone())?;
        assert_eq!(serde_json::to_value(&settings)?, settings_value);
        assert!(
            settings_value["hot_start"].is_null(),
            "Replay expects no hot start"
        );
        // The live interrupt callback returns false throughout uninterrupted
        // construction. No generated evaluator compilation is performed here.
        let settings = settings.abort_check(Some(Box::new(|| false)));
        let bytes = fs::read(directory.join("expression.raw"))?;
        assert_eq!(bytes.len() as u64, manifest["atom_bytes"]);
        let expression = AtomView::from(&bytes);
        stage(started, "exact_expression_loaded_without_normalization");
        println!("settings={settings_value}");
        stage(started, "symbolica_build_start");
        let evaluator = expression
            .evaluator(&params)
            .function_map(fn_map)
            .optimization_settings(settings)
            .build()?;
        stage(started, "symbolica_build_done");
        let operations = format!("{:?}", evaluator.count_operations());
        println!("operations={operations}");
        let success = directory.join("build_success.json");
        if success.exists() {
            let original: serde_json::Value = serde_json::from_reader(File::open(success)?)?;
            assert_eq!(
                operations, original["operations"],
                "Live/replayed operation counts differ"
            );
            stage(started, "live_operation_counts_match");
        }
        return Ok(());
    }
    if args.first().is_some_and(|a| a == "--combine") {
        assert!(
            args.len() >= 3,
            "usage: --combine OUTPUT INPUT...|DIRECTORY"
        );
        let output = Path::new(&args[1]);
        assert!(!output.exists(), "Refusing to overwrite the combined dump");
        let paths = if args.len() == 3 && Path::new(&args[2]).is_dir() {
            let mut paths = fs::read_dir(&args[2])?
                .map(|entry| entry.map(|e| e.path()))
                .collect::<Result<Vec<_>, _>>()?
                .into_iter()
                .filter_map(|path| {
                    let index = path
                        .file_name()?
                        .to_str()?
                        .strip_prefix("scalar_")?
                        .strip_suffix(".symbolica")?
                        .parse::<usize>()
                        .ok()?;
                    Some((index, path))
                })
                .collect::<Vec<_>>();
            paths.sort_by_key(|(index, _)| *index);
            paths.into_iter().map(|(_, path)| path).collect::<Vec<_>>()
        } else {
            args[2..].iter().map(PathBuf::from).collect()
        };
        assert!(!paths.is_empty(), "No scalar dumps supplied");
        stage(started, "component_import_start");
        let mut terms = Vec::with_capacity(paths.len());
        for (index, path) in paths.iter().enumerate() {
            let term = import(path)?;
            println!(
                "component={index} bytes={} path={}",
                term.as_view().get_byte_size(),
                path.display()
            );
            terms.push(term);
        }
        stage(started, "component_import_done");
        // Native n-way addition preserves factorized subexpressions; no expansion.
        let combined = Atom::add_many(terms);
        stage(started, "combine_done");
        println!(
            "atom_bytes={} terms={}",
            combined.as_view().get_byte_size(),
            combined.nterms()
        );
        let mut partial = output.as_os_str().to_os_string();
        partial.push(".partial");
        let partial = PathBuf::from(partial);
        let mut writer = BufWriter::new(File::create_new(&partial)?);
        combined.as_view().export(&mut writer)?;
        writer.flush()?;
        writer.get_ref().sync_all()?;
        drop(writer);
        fs::rename(partial, output)?;
        stage(started, "combined_export_done");
        return Ok(());
    }

    let abstract_parameters = args.first().is_some_and(|a| a == "--abstract-parameters");
    let control = args.as_slice() == ["--control"]
        || args.as_slice() == ["--abstract-parameters", "--control"];
    let (params, fn_map, expression) = if abstract_parameters {
        assert_eq!(
            args.len(),
            2,
            "usage: --abstract-parameters EXPRESSION|--control"
        );
        stage(started, "scalar_import_start");
        let expression = if control {
            let x = Atom::var(symbol!("abstract_mre::x"));
            let y = Atom::var(symbol!("abstract_mre::y"));
            Symbol::SIN.call(&x)
                + function!(symbol!("abstract_mre::f"), &x).pow(2) / Symbol::COS.call(&y)
        } else {
            import(Path::new(&args[1]))?
        };
        stage(started, "scalar_import_done");
        stage(started, "abstract_parameters_start");
        // An explicit abstraction of the imported math, not a physical mapping:
        // exact Var/Fun parameter lookup precedes builtin or function-map dispatch.
        let mut leaves = expression
            .get_all_indeterminates(false)
            .into_iter()
            .collect::<Vec<_>>();
        leaves.sort();
        let params = leaves.into_iter().map(|a| a.to_owned()).collect::<Vec<_>>();
        println!(
            "abstract_parameters={} atom_bytes={} terms={}",
            params.len(),
            expression.as_view().get_byte_size(),
            expression.nterms()
        );
        stage(started, "abstract_parameters_ready");
        (params, FunctionMap::new(), expression)
    } else {
        assert!(
            control || args.len() == 3,
            "usage: --control | EXPRESSION PARAMS DEFINITIONS | --inspect-context PARAMS DEFINITIONS | --abstract-parameters EXPRESSION|--control | --combine OUTPUT INPUT...|DIRECTORY"
        );
        stage(started, "context_import_start");
        let (params_atom, definitions_atom) = if control {
            let x = Atom::var(symbol!("mre::x"));
            let t = Atom::var(symbol!("mre::t"));
            let base = &t + Atom::one();
            (
                function!(symbol!("evaluator_probe::params"), &x),
                function!(
                    symbol!("evaluator_probe::definitions"),
                    function!(
                        symbol!("evaluator_probe::definition"),
                        function!(symbol!("mre::f"), 7, &t),
                        base.pow(2) + base.pow(3),
                        function!(symbol!("evaluator_probe::tags"), 7),
                        function!(symbol!("evaluator_probe::args"), &t)
                    ),
                    function!(
                        symbol!("evaluator_probe::definition"),
                        Atom::var(symbol!("mre::c")),
                        3,
                        function!(symbol!("evaluator_probe::tags")),
                        function!(symbol!("evaluator_probe::args"))
                    )
                ),
            )
        } else {
            (import(Path::new(&args[1]))?, import(Path::new(&args[2]))?)
        };
        let params = unpack(params_atom.as_view(), symbol!("evaluator_probe::params"));
        let entries = unpack(
            definitions_atom.as_view(),
            symbol!("evaluator_probe::definitions"),
        );
        let entry_count = entries.len();
        let mut fn_map = FunctionMap::new();
        // Replay ParamBuilder.reps in insertion order, including exact tags and args.
        for entry in entries {
            let [lhs, rhs, tags, args]: [Atom; 4] =
                unpack(entry.as_view(), symbol!("evaluator_probe::definition"))
                    .try_into()
                    .expect("Expected lhs, rhs, tags and args");
            let tags = unpack(tags.as_view(), symbol!("evaluator_probe::tags"));
            let args = unpack(args.as_view(), symbol!("evaluator_probe::args"));
            match lhs.as_view() {
                AtomView::Var(_) => {
                    assert!(tags.is_empty() && args.is_empty());
                    fn_map.add_aliases([(lhs, rhs)])?;
                }
                AtomView::Fun(f) => {
                    // A constant function-valued key has its arguments as fixed tags.
                    if tags.is_empty() && args.is_empty() && f.get_nargs() != 0 {
                        assert!(matches!(rhs.as_view(), AtomView::Num(_)));
                        fn_map.add_aliases([(lhs, rhs)])?;
                    } else {
                        assert!(f.iter().eq(tags.iter().chain(&args).map(Atom::as_view)));
                        let args = args
                            .into_iter()
                            .map(Indeterminate::try_from)
                            .collect::<Result<Vec<_>, _>>()?;
                        fn_map.add_tagged_function(f.get_symbol(), tags, args, rhs)?;
                    }
                }
                _ => panic!("Function-map key is not an indeterminate"),
            }
        }
        // Same two aliases appended by the measured GenericEvaluator owner.
        fn_map.add_aliases([
            (Atom::var(symbol!("vakint::𝑖")), Atom::i()),
            (Atom::var(symbol!("symbolica::𝑖")), Atom::i()),
        ])?;
        println!(
            "ordered_parameters={} definitions={entry_count}",
            params.len()
        );
        stage(started, "context_ready");
        if args.first().is_some_and(|a| a == "--inspect-context") {
            for (index, param) in params.iter().enumerate() {
                println!("parameter[{index}]={param}");
                if let AtomView::Var(v) = param.as_view() {
                    println!(
                        "symbol={} id={} attributes={:?}",
                        v.get_symbol().get_name(),
                        v.get_symbol().get_id(),
                        v.get_symbol().get_attributes()
                    );
                }
            }
            println!("definitions={definitions_atom}");
            return Ok(());
        }
        stage(started, "scalar_import_start");
        let expression = if control {
            function!(symbol!("mre::f"), 7, &params[0]) + Atom::var(symbol!("mre::c"))
        } else {
            import(Path::new(&args[0]))?
        };
        println!(
            "atom_bytes={} terms={}",
            expression.as_view().get_byte_size(),
            expression.nterms()
        );
        stage(started, "scalar_import_done");
        (params, fn_map, expression)
    };
    let settings = OptimizationSettings::new()
        .horner_iterations(1)
        .cpe_iterations(Some(5))
        .cores(1)
        .direct_translation(true)
        .abort_check(Some(Box::new(|| false)))
        .abort_level(0)
        .max_horner_scheme_variables(500)
        .max_common_pair_cache_entries(1_000_000)
        .max_common_pair_distance(1000)
        .verbose(false);
    println!("H1 CPE5 cores=1 direct_translation=true; evaluator_code_compilation=false");
    stage(started, "symbolica_build_start");
    // Constructs Symbolica's interpreted instruction representation only.
    let evaluator = expression
        .evaluator(&params)
        .function_map(fn_map)
        .optimization_settings(settings)
        .build()?;
    stage(started, "symbolica_build_done");
    println!("operations={:?}", evaluator.count_operations());
    if control {
        let mut numeric = evaluator.map_coeff(&|r| r.re.to_f64());
        assert_eq!(
            numeric.evaluate_single(&vec![2.0; params.len()]),
            if abstract_parameters { 4.0 } else { 39.0 }
        );
        stage(started, "control_passed");
    }
    Ok(())
}
