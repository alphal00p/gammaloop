//! Replay of dumped closure requests (`BETA_RUSTRED_DUMP_DIR`): the factored
//! finalization against the expansion of the product, the path before the
//! factored values and still the path of values without a spectator.
//!
//! `VAKINT_CLOSURE_DUMPS=<dir of closure-*.txt> VAKINT_P_LIBRARY=<p library>
//! cargo test --release -p vakint --lib closure_replay -- --ignored --nocapture`
//!
//! `VAKINT_REPLAY_TARGET` (default -1, the production closure depth: poles
//! only), `VAKINT_REPLAY_REPEATS` (default 2) and `VAKINT_REPLAY_ONLY` (a
//! file-name filter) refine it.
use std::time::Instant;

use super::*;
use crate::{EvaluationOrder, LoopNormalizationFactor, RustRedEvaluationOptions, VakintSettings};

fn environment<T: std::str::FromStr>(name: &str, default: T) -> T {
    std::env::var(name)
        .ok()
        .and_then(|value| value.parse().ok())
        .unwrap_or(default)
}

#[test]
#[ignore = "needs VAKINT_CLOSURE_DUMPS and VAKINT_P_LIBRARY"]
fn closure_replay() {
    let dumps = std::env::var("VAKINT_CLOSURE_DUMPS").expect("VAKINT_CLOSURE_DUMPS");
    let library = std::env::var("VAKINT_P_LIBRARY").expect("VAKINT_P_LIBRARY");
    let repeats: usize = environment("VAKINT_REPLAY_REPEATS", 2);
    let target: i64 = environment("VAKINT_REPLAY_TARGET", -1);
    let only = std::env::var("VAKINT_REPLAY_ONLY").ok();
    Vakint::initialize_vakint_symbols();
    let package = rustred_app::propagator::load_p_library(
        rustred_app::propagator::PLibrarySource::Directory(library.into()),
        &rustred_app::propagator::PLibraryOptions {
            verify_descent: false,
        },
    )
    .unwrap();
    let mut vakint = Vakint::new().unwrap();
    vakint
        .inject_p_library(std::sync::Arc::clone(&package))
        .unwrap();
    let base = VakintSettings {
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions {
                master_mode: RustRedMasterMode::Symbolic,
            },
        )]),
        integral_normalization_factor: LoopNormalizationFactor::MSbar,
        epsilon_symbol: "gammalooprs::ε".into(),
        mu_r_sq_symbol: "gammalooprs::mursq".into(),
        use_dot_product_notation: true,
        n_cores: 1,
        number_of_terms_in_epsilon_expansion: 4,
        ..VakintSettings::default()
    };
    let mut files = std::fs::read_dir(&dumps)
        .unwrap()
        .map(|entry| entry.unwrap().path())
        .filter(|path| path.extension().is_some_and(|e| e == "txt"))
        .collect::<Vec<_>>();
    files.sort();
    let (mut total_old, mut total_new, mut total_reduce) = (0.0, 0.0, 0.0);
    let mut all_equal = true;
    let mut closures = 0;
    println!("REPLAY file cut_loops values reduce_s old_s new_s ratio equal out_bytes");
    for path in files {
        let name = path.file_name().unwrap().to_string_lossy().to_string();
        if only
            .as_ref()
            .is_some_and(|only| !name.contains(only.as_str()))
        {
            continue;
        }
        let text = std::fs::read_to_string(&path).unwrap();
        let mut input = Atom::Zero;
        for line in text
            .lines()
            .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        {
            input += vk_parse!(line).unwrap();
        }
        let mut topologies = Vec::new();
        input.visitor(&mut |part| {
            if let AtomView::Fun(function) = part
                && function.get_symbol() == S.topo
            {
                topologies.push(part.to_owned());
                return false;
            }
            true
        });
        topologies.sort();
        topologies.dedup();
        let mut jobs = Vec::new();
        let mut head = None;
        for (topology, numerator) in input.coefficient_list::<i16>(&topologies) {
            let mut expression = VakintExpression::try_from(topology * numerator).unwrap();
            for term in &mut expression.0 {
                if !term.normalize_propagator().unwrap() {
                    continue;
                }
                let external = term.propagator_external().unwrap().unwrap();
                term.numerator = crate::p_integral::scalar_numerator(term.numerator.as_view(), &external).unwrap();
                let matched = term.prepare_integral_specs(&vakint, &base, None).unwrap();
                head = Some(external);
                jobs.push((matched, term.numerator.clone()));
            }
        }
        let head = head.unwrap();
        let requests = jobs
            .iter()
            .map(|(matched, numerator)| matched.propagator_request(numerator.as_view()).unwrap())
            .collect::<Vec<_>>();
        let started = Instant::now();
        let reduced = crate::rustred_evaluation::p_integral::reduce_requests_factored(
            &vakint, &base, &requests, None,
        )
        .unwrap();
        let reduce_seconds = started.elapsed().as_secs_f64();
        total_reduce += reduce_seconds;
        for ((cut_loops, _), values) in reduced {
            let mut settings = base.clone();
            settings.number_of_terms_in_epsilon_expansion = target + cut_loops as i64 + 2;
            let products = values
                .iter()
                .map(|value| ReducedPropagatorValue {
                    scale_power: value.scale_power,
                    reduced: value.value().into_owned(),
                    spectator: None,
                })
                .collect::<Vec<_>>();
            // A logarithmically divergent closure of the first channel.
            let radial = MassiveClosure::radial_terms(
                &values[0].spectator.as_ref().unwrap().expand(),
                &head,
            )
            .unwrap()
            .first()
            .map_or(0, |(power, _)| *power);
            let closure = MassiveClosure {
                momentum: head.clone(),
                mass_squared: vk_parse!("mUVsq").unwrap(),
                power: (2 + values[0].scale_power + radial).max(1),
            };
            let run = |values: &[ReducedPropagatorValue]| {
                let started = Instant::now();
                let value = closure
                    .finalize(
                        values,
                        cut_loops,
                        &head,
                        &settings,
                        RustRedMasterMode::Symbolic,
                    )
                    .unwrap();
                (started.elapsed().as_secs_f64(), value)
            };
            let (mut old_best, mut new_best) = (f64::INFINITY, f64::INFINITY);
            let (mut old, mut new) = (Atom::Zero, Atom::Zero);
            // Alternate the order to pair both on a loaded host.
            for repeat in 0..repeats.max(1) {
                for factored in [repeat % 2 == 0, repeat % 2 != 0] {
                    if factored {
                        let (seconds, value) = run(&values);
                        new_best = new_best.min(seconds);
                        new = value;
                    } else {
                        let (seconds, value) = run(&products);
                        old_best = old_best.min(seconds);
                        old = value;
                    }
                }
            }
            let equal = old.to_canonical_string() == new.to_canonical_string();
            all_equal &= equal;
            closures += 1;
            total_old += old_best;
            total_new += new_best;
            println!(
                "REPLAY {name} {cut_loops} {} {reduce_seconds:.3} {old_best:.4} {new_best:.4} {:.2} {equal} {}",
                values.len(),
                old_best / new_best,
                old.as_view().get_byte_size(),
            );
        }
    }
    println!(
        "REPLAY_TOTAL closures {closures} reduce_s {total_reduce:.2} old_s {total_old:.3} new_s {total_new:.3} ratio {:.2} all_canonical_equal {all_equal}",
        total_old / total_new
    );
    assert!(all_equal);
}
