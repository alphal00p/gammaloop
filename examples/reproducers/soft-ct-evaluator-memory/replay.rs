#!/usr/bin/env rust-script
//! ```cargo
//! [dependencies]
//! symbolica = { git = "https://github.com/symbolica-dev/symbolica", rev = "578dfcb55fb0662d871456ea8ecc57c1390caa77", default-features = false, features = ["float-mpfr", "integer-gmp", "serde"] }
//! serde_json = "1.0.149"
//! [patch.crates-io]
//! numerica = { git = "https://github.com/symbolica-dev/symbolica", rev = "578dfcb55fb0662d871456ea8ecc57c1390caa77" }
//! graphica = { git = "https://github.com/symbolica-dev/symbolica", rev = "578dfcb55fb0662d871456ea8ecc57c1390caa77" }
//! ```
// Usage: replay INPUT.json default|off|raw [--check-raw-numeric]
// Accepts the evaluator_input payload or its enclosing JSON log event. Run each
// mode in a fresh process under an external resource guard. Between default and
// off, only cpe_iterations changes: None versus Some(0). No function or
// alias is substituted algebraically; all recorded inlining policies are kept.
// Numeric checks use three deterministic complex parameter vectors, NOT physical
// momentum points. Compare the full printed complex outputs across the two runs.
// Raw mode requires the captured single constant-argument function root and
// Horner iteration count 1. It preserves that mode's preprocessing, then returns
// before CSE, CPE and stack compaction. --check-raw-numeric compacts the stack and
// checks small raw-mode fixtures; omit it for the large memory diagnostic.
use serde_json::{Value, json};
use std::{error::Error, io::Write, time::Instant};
use symbolica::{
    evaluate::{FunctionMap, FunctionRegistrationOptions, InliningPolicy, OptimizationSettings},
    prelude::*,
};

const REVISION: &str = "578dfcb55fb0662d871456ea8ecc57c1390caa77";
const PRECISION: u32 = 1000;

fn emit(value: Value) -> Result<(), Box<dyn Error>> {
    let mut out = std::io::stdout().lock();
    serde_json::to_writer(&mut out, &value)?;
    writeln!(out)?;
    out.flush()?;
    Ok(())
}

fn memory() -> Value {
    let status = std::fs::read_to_string("/proc/self/status").unwrap_or_default();
    let field = |key: &str| {
        status.lines().find_map(|line| {
            line.strip_prefix(key)?
                .split_whitespace()
                .next()?
                .parse::<u64>()
                .ok()
        })
    };
    json!({"rss_kib": field("VmRSS:"), "hwm_kib": field("VmHWM:")})
}

fn atom(value: &Value) -> Result<Atom, Box<dyn Error>> {
    let text = value.as_str().ok_or("expected an atom string")?;
    Ok(try_parse!(text, default_namespace = "gammalooprs")?)
}

fn atoms(value: &Value) -> Result<Vec<Atom>, Box<dyn Error>> {
    value
        .as_array()
        .ok_or("expected an atom array")?
        .iter()
        .map(atom)
        .collect()
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = std::env::args().collect::<Vec<_>>();
    if !(3..=4).contains(&args.len())
        || !["default", "off", "raw"].contains(&args[2].as_str())
        || (args.len() == 4 && (args[2] != "raw" || args[3] != "--check-raw-numeric"))
    {
        return Err("usage: replay INPUT.json default|off|raw [--check-raw-numeric]".into());
    }
    let mode = &args[2];
    let preparation = Instant::now();
    emit(json!({"stage":"parse_start", "mode":mode, "symbolica_revision":REVISION}))?;
    let text = std::fs::read_to_string(&args[1])?;
    let mut input: Value = serde_json::from_str(&text)?;
    drop(text);
    if let Some(payload) = input
        .get("evaluator_input")
        .or_else(|| input.get("file.evaluator_input"))
    {
        input = match payload {
            Value::String(text) => serde_json::from_str(text)?,
            Value::Object(_) => payload.clone(),
            _ => return Err("invalid evaluator_input payload".into()),
        };
    }
    if !input
        .get("dual_shape")
        .ok_or("missing dual_shape")?
        .is_null()
    {
        return Err("this replay requires a scalar evaluator, without dual vectorization".into());
    }
    let recorded_settings = input
        .get("optimization_settings")
        .ok_or("missing recorded optimization settings")?
        .clone();
    if !recorded_settings
        .get("cpe_iterations")
        .ok_or("missing cpe_iterations")?
        .is_null()
    {
        return Err("expected the recorded cpe_iterations to be None".into());
    }
    let mut settings: OptimizationSettings = serde_json::from_value(recorded_settings.clone())?;
    if mode == "off" {
        settings = settings.cpe_iterations(Some(0));
    }

    let mut root = atom(&input["root"])?;
    if mode == "raw" {
        let empty_scheme = match root.as_view() {
            AtomView::Fun(call) => call.iter().all(|arg| matches!(arg, AtomView::Num(_))),
            _ => false,
        };
        if !empty_scheme
            || recorded_settings["horner_iterations"] != 1
            || recorded_settings["direct_translation"] != true
            || !recorded_settings["hot_start"].is_null()
        {
            return Err("raw mode requires a single numeric-argument function root, direct translation, one Horner iteration, and no hot start".into());
        }
        // Pinned tree.rs counts indeterminates in the ROOT, keeping counts >1.
        // This root has exactly one function occurrence and no nonnumeric args,
        // hence an empty scheme. collect.rs returns the already-normalized atom
        // for that empty scheme, then calls collect_by_coefficient(). Applying
        // this public operation to every body reproduces the actual transform;
        // collect_horner(None) would choose a DIFFERENT scheme for each body.
        settings = settings.horner_iterations(0);
    }
    let mut preprocessing_seconds = 0.0;
    let mut preprocessing_count = 0usize;
    let mut preprocess = |body: Atom| {
        if mode == "raw" {
            let started = Instant::now();
            let result = body.collect_by_coefficient();
            preprocessing_seconds += started.elapsed().as_secs_f64();
            preprocessing_count += 1;
            result
        } else {
            body
        }
    };
    let parameters = atoms(&input["parameters"])?;
    // Canonical alias bodies carry symbol attributes. Parse them before archived
    // plain function LHSs, so first registration does not erase those attributes.
    let mut aliases = input["aliases"]
        .as_array()
        .ok_or("missing aliases")?
        .iter()
        .map(|pair| -> Result<_, Box<dyn Error>> {
            let fields = pair.as_array().ok_or("expected an alias pair")?;
            if fields.len() != 2 {
                return Err("an alias requires a key and a body".into());
            }
            Ok((atom(&fields[0])?, atom(&fields[1])?))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let source_bytes = root.as_view().get_byte_size()
        + aliases
            .iter()
            .map(|(key, rhs)| key.as_view().get_byte_size() + rhs.as_view().get_byte_size())
            .sum::<usize>();
    root = preprocess(root);
    for (_, rhs) in &mut aliases {
        *rhs = preprocess(std::mem::take(rhs));
    }
    let mut function_map = FunctionMap::new();
    let entries = input["functions"].as_array().ok_or("missing functions")?;
    let mut function_bytes = 0usize;
    let mut policy_counts = [0usize; 3];
    let mut function_aliases = 0usize;
    for entry in entries {
        let fields = entry
            .as_array()
            .ok_or("expected a function archive tuple")?;
        if fields.len() != 6 {
            return Err("expected lhs,rhs,tags,args,inlining,is_alias".into());
        }
        let lhs = atom(&fields[0])?;
        let rhs = atom(&fields[1])?;
        let tags = atoms(&fields[2])?;
        let formals = atoms(&fields[3])?
            .into_iter()
            .map(Indeterminate::try_from)
            .collect::<Result<Vec<_>, _>>()?;
        let policy: InliningPolicy = serde_json::from_value(fields[4].clone())?;
        let is_alias = fields[5].as_bool().ok_or("is_alias must be a boolean")?;
        function_bytes += rhs.as_view().get_byte_size();
        let rhs = preprocess(rhs);
        if is_alias {
            if !formals.is_empty() || policy != InliningPolicy::Always {
                return Err("caller-scope aliases require no formals and Always inlining".into());
            }
            function_map.add_aliases([(lhs, rhs)])?;
            function_aliases += 1;
        } else {
            let name = match lhs.as_view() {
                AtomView::Fun(call) => call.get_symbol(),
                AtomView::Var(variable) if tags.is_empty() => variable.get_symbol(),
                _ => return Err("invalid ordinary-function LHS".into()),
            };
            policy_counts[match policy {
                InliningPolicy::Always => 0,
                InliningPolicy::Never => 1,
                InliningPolicy::Auto => 2,
            }] += 1;
            function_map.add_tagged_function_with_options(
                name,
                tags,
                formals,
                rhs,
                FunctionRegistrationOptions::new().inlining(policy),
            )?;
        }
    }
    // Match GenericEvaluator's two legacy imaginary-unit aliases, which are
    // registered immediately before the captured build but absent from entries.
    function_map.add_aliases([
        (Atom::var(symbol!("vakint::𝑖")), preprocess(Atom::i())),
        (Atom::var(symbol!("symbolica::𝑖")), preprocess(Atom::i())),
    ])?;
    let source = json!({
        "root_bytes":root.as_view().get_byte_size(), "root_and_alias_bytes":source_bytes,
        "root_alias_count":aliases.len(), "function_count":entries.len(),
        "function_body_bytes":function_bytes, "function_alias_count":function_aliases,
        "policy_counts":{"Always":policy_counts[0],"Never":policy_counts[1],"Auto":policy_counts[2]},
        "parameter_count":parameters.len(), "recorded_settings":recorded_settings,
        "effective_settings":settings,
        "raw_preprocessing_seconds":preprocessing_seconds,
        "raw_preprocessed_atom_count":preprocessing_count,
    });
    drop(input);
    emit(json!({"stage":"build_start", "mode":mode, "source":source,
        "preparation_seconds":preparation.elapsed().as_secs_f64(), "memory":memory()}))?;
    let started = Instant::now();
    let mut evaluator = root
        .evaluator(&parameters)
        .function_map(function_map)
        .add_aliases(aliases)?
        .optimization_settings(settings)
        .build()?;
    let build_seconds = started.elapsed().as_secs_f64();
    let operations = evaluator.count_operations();
    let input_count = evaluator.get_input_len();
    let output_count = evaluator.get_output_len();
    emit(
        json!({"stage":"build_done", "mode":mode, "build_seconds":build_seconds,
        "memory":memory(), "input_count":input_count,
        "output_count":output_count, "operations":{
            "additions":operations.additions,"multiplications":operations.multiplications,
            "inversions":operations.inversions,"function_calls":operations.function_calls}}),
    )?;
    if mode == "raw" {
        emit(json!({"stage":"diagnostic_raw_complete", "mode":mode,
            "build_seconds":build_seconds, "preprocessing_seconds":preprocessing_seconds,
            "cse_applied":false, "cpe_applied":false, "stack_compacted":false,
            "memory":memory()}))?;
        if args.len() == 3 {
            return Ok(());
        }
        let started = Instant::now();
        evaluator.optimize_stack();
        emit(json!({"stage":"raw_numeric_stack_compacted",
            "elapsed_seconds":started.elapsed().as_secs_f64(), "memory":memory()}))?;
    }
    emit(json!({"stage":"numeric_start", "precision_bits":PRECISION,
        "point_kind":"synthetic complex parameter vectors, not physical momenta",
        "point_formula":"parameter i at point p: Re=(17+(73*i+137*p)%997)/257; Im=(11+(37*i+193*p)%991)/263, zero-based i,p"}))?;
    let mut evaluator = evaluator.map_coeff_with_prec(
        &|c| {
            Complex::new(
                c.re.to_multi_prec_float(PRECISION),
                c.im.to_multi_prec_float(PRECISION),
            )
        },
        PRECISION,
    );
    for point_index in 0..3u64 {
        let values = (0..input_count)
            .map(|i| {
                let i = i as u64;
                Complex::new(
                    Float::with_val(PRECISION, 17 + (73 * i + 137 * point_index) % 997)
                        / Float::with_val(PRECISION, 257),
                    Float::with_val(PRECISION, 11 + (37 * i + 193 * point_index) % 991)
                        / Float::with_val(PRECISION, 263),
                )
            })
            .collect::<Vec<_>>();
        let mut outputs =
            vec![
                Complex::new(Float::with_val(PRECISION, 0), Float::with_val(PRECISION, 0));
                output_count
            ];
        let started = Instant::now();
        evaluator.evaluate(&values, &mut outputs);
        let finite = outputs.iter().all(|v| v.re.is_finite() && v.im.is_finite());
        let zero = Float::with_val(PRECISION, 0);
        let nonzero = outputs.iter().any(|v| v.re != zero || v.im != zero);
        emit(
            json!({"stage":"numeric_point", "mode":mode, "point_index":point_index,
            "elapsed_seconds":started.elapsed().as_secs_f64(), "finite":finite, "nonzero":nonzero,
            "outputs":outputs.iter().map(|v| [format!("{:.310e}",v.re),format!("{:.310e}",v.im)]).collect::<Vec<_>>()}),
        )?;
        if !finite {
            return Err("nonfinite synthetic-point result; comparison is inconclusive".into());
        }
        if !nonzero {
            return Err("all synthetic-point outputs are zero; comparison is inconclusive".into());
        }
    }
    emit(json!({"stage":"complete", "mode":mode, "memory":memory()}))
}
