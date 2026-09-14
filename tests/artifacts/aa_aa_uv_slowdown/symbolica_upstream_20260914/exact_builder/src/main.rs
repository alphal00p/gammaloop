use std::{
    collections::BTreeMap,
    error::Error,
    fs::{self, File},
    io::{BufReader, Read, Seek, Write},
    path::Path,
    time::Instant,
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Indeterminate, Symbol},
    evaluate::{FunctionMap, FunctionRegistrationOptions, OptimizationSettings},
    state::State,
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
        // Upstream added an inlining policy to each serialized definition.
        // Decode the captured layout and add its equivalent default (Always).
        // Preserve internal IDs, ordered arguments, aliases and every Atom;
        // the raw audit below rejects any expression change during conversion.
        type CapturedMap = (
            BTreeMap<(Symbol, Vec<Atom>), (usize, usize, Vec<Indeterminate>, Atom)>,
            BTreeMap<Symbol, usize>,
        );
        let ((definitions, tags), consumed): (CapturedMap, _) =
            bincode::decode_from_slice_with_context(&map_bytes, config, state_map)?;
        assert_eq!(consumed, map_bytes.len());
        let definitions = definitions
            .into_iter()
            .map(|(key, (id, tags, args, body))| {
                (
                    key,
                    (id, tags, args, body, FunctionRegistrationOptions::default()),
                )
            })
            .collect::<BTreeMap<_, _>>();
        let upgraded = bincode::encode_to_vec((definitions, tags), config)?;
        let (fn_map, consumed): (FunctionMap, _) = bincode::decode_from_slice_with_context(
            &upgraded,
            config,
            symbolica::state::StateMap::default(),
        )?;
        assert_eq!(consumed, upgraded.len());
        let (before, consumed): (RawFunctionMap, _) =
            bincode::decode_from_slice(&map_bytes, config)?;
        assert_eq!(consumed, map_bytes.len());
        let restored_map = bincode::encode_to_vec(&fn_map, config)?;
        type CurrentRawMap = (
            BTreeMap<
                (RawSymbol, Vec<RawAtom>),
                (
                    usize,
                    usize,
                    Vec<RawIndeterminate>,
                    RawAtom,
                    FunctionRegistrationOptions,
                ),
            >,
            BTreeMap<RawSymbol, usize>,
        );
        let ((definitions, tags), consumed): (CurrentRawMap, _) =
            bincode::decode_from_slice(&restored_map, config)?;
        assert_eq!(consumed, restored_map.len());
        let after = (
            definitions
                .into_iter()
                .map(|(key, (id, tags, args, body, options))| {
                    assert_eq!(options, FunctionRegistrationOptions::default());
                    (key, (id, tags, args, body))
                })
                .collect(),
            tags,
        );
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
            println!("pinned_live_operations={}", original["operations"]);
            println!(
                "operation_counts_match={}",
                operations == original["operations"]
            );
            // Upstream may legitimately produce a different instruction count.
            // Raw input/context identity is certified before construction.
        }
        return Ok(());
    }
    panic!("usage: --exact-builder CAPTURE_DIRECTORY");
}
