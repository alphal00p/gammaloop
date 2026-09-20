//! Shipped, finite-tested four-loop candidate programs.
//!
//! These are not `ClosedArtifact`s: a requested point succeeds only when
//! RustRed's candidate applier finds guarded, descending rules ending in the
//! declared terminal catalog. Unknown leaves remain errors. Discovery and
//! FMFT terminal preparation are offline tasks, never runtime fallbacks.
//! Native Symbolica programs are trusted, build-time embedded data from the
//! pinned RustRed/Symbolica stack, not a public untrusted-file input boundary.

use std::io::Read;
use std::sync::{Arc, LazyLock};

use flate2::read::GzDecoder;
use rustred::input::{Compiler, Limits, LoweringLimits, TextProject, TextPropagator};
use rustred::reduction::ReductionLimits;
use rustred_app::{
    CandidateBundleLimits, MAX_CANDIDATE_BUNDLE_BYTES, load_generated_candidate_bundle,
};
use symbolica::atom::{Atom, AtomView};
use symbolica::function;

use super::experimental::ExperimentalRustRed;
use super::experimental::catalog::OfflineTerminalCatalog;
use super::experimental::native::NativeCandidate;
use super::{RustRedEvaluationError, RustRedEvaluationOptions};
use crate::symbols::S;
use crate::utils::vakint_macros::vk_symbol;
use crate::{
    ReplacementRules, Topology, Vakint, VakintError, VakintSettings, get_integer_from_atom,
    get_prop_with_id,
};

// Explicit shipped-data policy. These caps bound decompressed/framed data;
// Symbolica's native readers still require trusted generated provenance.
const MAX_PROGRAM_BYTES: usize = MAX_CANDIDATE_BUNDLE_BYTES;
const MAX_COEFFICIENT_BYTES: usize = 512 * 1024 * 1024;

struct Input {
    label: &'static str,
    csv: &'static str,
    program: &'static [u8],
    catalog: &'static [u8],
    normalization: &'static [u8],
}

macro_rules! input {
    ($name:literal) => {
        Input {
            label: concat!("four-loop-", $name),
            csv: include_str!(concat!("../../data/rustred/four_loop/", $name, ".csv")),
            program: include_bytes!(concat!(
                "../../data/rustred/four_loop/",
                $name,
                ".candidates.rrbin.gz"
            )),
            catalog: include_bytes!(concat!(
                "../../data/rustred/four_loop/",
                $name,
                ".rrcat.bin"
            )),
            normalization: include_bytes!(concat!(
                "../../data/rustred/four_loop/",
                $name,
                ".rrnorm.bin"
            )),
        }
    };
}

// Labels are asset identifiers only. Selection below uses the retained
// simultaneous momentum witness, never a topology/integral name.
const INPUTS: [Input; 4] = [input!("h"), input!("fg"), input!("bmw"), input!("x")];

struct Descriptor {
    physical_momenta: Vec<Atom>,
    family_fingerprint: String,
}

static DESCRIPTORS: LazyLock<Result<Vec<Descriptor>, String>> = LazyLock::new(|| {
    Vakint::initialize_vakint_symbols();
    INPUTS.iter().map(|input| descriptor(input.csv)).collect()
});

// Independent lazy cells avoid loading all large programs when a calculation
// only visits one parent. NativeCandidate owns one shared memoized applier.
static H: LazyLock<Result<ExperimentalRustRed, String>> = LazyLock::new(|| load(0));
static FG: LazyLock<Result<ExperimentalRustRed, String>> = LazyLock::new(|| load(1));
static BMW: LazyLock<Result<ExperimentalRustRed, String>> = LazyLock::new(|| load(2));
static X: LazyLock<Result<ExperimentalRustRed, String>> = LazyLock::new(|| load(3));

fn descriptor(source: &str) -> Result<Descriptor, String> {
    let mut physical_momenta = Vec::new();
    let mut propagators = Vec::new();
    let mut auxiliary_seen = false;
    for line in source
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
    {
        let row = line
            .split(',')
            .map(|entry| {
                entry
                    .trim()
                    .parse::<i64>()
                    .map_err(|error| error.to_string())
            })
            .collect::<Result<Vec<_>, _>>()?;
        if row.len() != 6 {
            return Err("parent descriptor needs six columns".into());
        }
        let auxiliary = row[0] == 0 && row[1] == 0;
        if auxiliary_seen && !auxiliary {
            return Err("physical slots must precede auxiliaries".into());
        }
        auxiliary_seen |= auxiliary;
        let mut momentum = Atom::Zero;
        let mut components = Vec::new();
        for (axis, coefficient) in row[2..].iter().enumerate() {
            if *coefficient != 0 {
                momentum += Atom::num(*coefficient) * function!(S.k, Atom::num(axis + 1));
                components.push(format!("({coefficient})*k{}", axis + 1));
            }
        }
        if momentum.is_zero() {
            return Err("zero descriptor momentum".into());
        }
        propagators.push(TextPropagator {
            id: format!("D{}", propagators.len() + 1),
            expression: format!("({})^2-1", components.join("+")),
            target_power: i64::from(!auxiliary),
            power_shift: None,
        });
        if !auxiliary {
            physical_momenta.push(momentum);
        }
    }
    if propagators.len() != 10 || physical_momenta.is_empty() {
        return Err("four-loop descriptor needs ten complete denominator slots".into());
    }
    let family = Compiler::new(Limits::default())
        .map_err(|error| error.to_string())?
        .compile_text(TextProject {
            name: None,
            parameters: None,
            loop_momenta: (1..=4).map(|axis| format!("k{axis}")).collect(),
            external_momenta: Vec::new(),
            dimension: "d".into(),
            propagators,
            external_gram: Vec::new(),
            numerator: None,
        })
        .map_err(|error| error.to_string())?
        .into_lowered(LoweringLimits::default())
        .map_err(|error| error.to_string())?
        .into_family();
    Ok(Descriptor {
        physical_momenta,
        family_fingerprint: family.fingerprint().to_owned(),
    })
}

fn load(index: usize) -> Result<ExperimentalRustRed, String> {
    let input = &INPUTS[index];
    let descriptor = &DESCRIPTORS.as_ref().map_err(Clone::clone)?[index];
    let mut bytes = Vec::new();
    GzDecoder::new(input.program)
        .take(MAX_PROGRAM_BYTES as u64 + 1)
        .read_to_end(&mut bytes)
        .map_err(|error| format!("decompress {}: {error}", input.label))?;
    if bytes.len() > MAX_PROGRAM_BYTES {
        return Err("shipped program exceeds byte limit".into());
    }
    let limits = CandidateBundleLimits {
        max_bundle_bytes: MAX_PROGRAM_BYTES,
        max_collection_entries: 8_000_000,
        max_total_coefficient_bytes: MAX_COEFFICIENT_BYTES,
        ..CandidateBundleLimits::default()
    };
    let (family, reducer) =
        load_generated_candidate_bundle::<10>(&bytes, limits, ReductionLimits::default())
            .map_err(|error| error.to_string())?;
    if family.fingerprint() != descriptor.family_fingerprint {
        return Err("shipped program differs from its physical/auxiliary descriptor".into());
    }
    let catalog =
        OfflineTerminalCatalog::decode_generated(input.catalog, family.fingerprint(), 10)?;
    catalog.require_complete()?;
    if catalog.terms().keys().ne(reducer.terminals().iter()) {
        return Err("offline catalog does not exactly cover declared program terminals".into());
    }
    let native = NativeCandidate::from_reducer_with_terminal_normalization(
        Arc::new(family),
        descriptor.physical_momenta.clone(),
        reducer,
        input.normalization,
    )?;
    ExperimentalRustRed::new(Arc::new(native), catalog.into_terms())
        .map_err(|error| error.to_string())
}

fn select(topology: &Topology) -> Result<usize, RustRedEvaluationError> {
    let integral = topology.get_integral();
    if !matches!(topology, Topology::FourLoop(_)) || integral.n_loops != 4 {
        return Err(invalid("not a matched four-loop family"));
    }
    let descriptors = DESCRIPTORS
        .as_ref()
        .map_err(|detail| invalid(detail.clone()))?;
    let index = descriptors
        .iter()
        .position(|descriptor| {
            integral
                .validate_parent_routing(&descriptor.physical_momenta)
                .is_ok()
        })
        .ok_or_else(|| invalid("no shipped program matches the retained parent witness"))?;
    let expression = integral
        .canonical_expression
        .as_ref()
        .ok_or_else(|| invalid("missing canonical expression"))?;
    let mut mass: Option<Atom> = None;
    for slot in 1..=integral.n_props {
        let Some(properties) = get_prop_with_id(expression.as_view(), slot) else {
            continue;
        };
        if properties
            .get(&vk_symbol!("pow_"))
            .and_then(|power| get_integer_from_atom(power.as_view()))
            .is_none()
        {
            return Err(invalid(
                "four-loop program requires integer propagator powers",
            ));
        }
        let current = properties
            .get(&vk_symbol!("mUVsq_"))
            .ok_or_else(|| invalid("missing mass squared"))?;
        if current.is_zero() || mass.as_ref().is_some_and(|value| value != current) {
            return Err(invalid(
                "four-loop program requires one nonzero common mass",
            ));
        }
        mass.get_or_insert_with(|| current.clone());
    }
    if mass.is_none() {
        return Err(invalid("matched family has no physical mass"));
    }
    Ok(index)
}

fn invalid(detail: impl Into<String>) -> RustRedEvaluationError {
    RustRedEvaluationError::InvalidMatchedFamily {
        detail: detail.into(),
    }
}

pub(super) fn supports(topology: &Topology) -> bool {
    select(topology).is_ok()
}

pub(super) fn evaluate(
    settings: &VakintSettings,
    numerator: AtomView,
    matched: &ReplacementRules,
    options: &RustRedEvaluationOptions,
) -> Result<Atom, VakintError> {
    let index = select(&matched.canonical_topology)?;
    let program = match index {
        0 => &*H,
        1 => &*FG,
        2 => &*BMW,
        3 => &*X,
        _ => unreachable!(),
    };
    let program = program
        .as_ref()
        .map_err(|detail| RustRedEvaluationError::ArtifactLoad {
            family: INPUTS[index].label,
            detail: detail.clone(),
        })?;
    program
        .evaluate_matched(settings, numerator, matched, options.substitute_masters)
        .map(|result| result.value)
}

#[cfg(test)]
mod tests {
    use super::*;
    use symbolica::atom::AtomCore;

    #[test]
    fn shipped_descriptors_have_distinct_family_bindings() {
        let descriptors = DESCRIPTORS.as_ref().unwrap();
        assert_eq!(
            descriptors
                .iter()
                .map(|d| d.physical_momenta.len())
                .collect::<Vec<_>>(),
            [9, 8, 8, 9]
        );
        for (i, left) in descriptors.iter().enumerate() {
            for right in &descriptors[i + 1..] {
                assert_ne!(left.family_fingerprint, right.family_fingerprint);
                assert_ne!(left.physical_momenta, right.physical_momenta);
            }
        }
    }

    #[test]
    fn shipped_programs_cold_load_with_exact_terminal_catalogs() {
        for index in 0..INPUTS.len() {
            load(index).unwrap_or_else(|error| panic!("{}: {error}", INPUTS[index].label));
        }
    }

    #[test]
    fn shipped_normalizations_keep_raw_keys_and_close_onto_positive_outputs() {
        use rustred::persistence::BinaryIoLimits;
        use rustred::reduction::terminal_normalization::{
            TerminalNormalizationLimits, TerminalNormalizationPlan,
        };

        // These are expectations for vendored input data, not engine limits or
        // a dispatch based on a topology name. Original rules/catalogs stay raw.
        for (input, (raw_count, output_count)) in
            INPUTS
                .iter()
                .zip([(386, 22), (145, 16), (179, 17), (445, 19)])
        {
            let mut bytes = Vec::new();
            GzDecoder::new(input.program)
                .read_to_end(&mut bytes)
                .unwrap();
            let (family, reducer) = load_generated_candidate_bundle::<10>(
                &bytes,
                CandidateBundleLimits {
                    max_bundle_bytes: MAX_PROGRAM_BYTES,
                    max_collection_entries: 8_000_000,
                    max_total_coefficient_bytes: MAX_COEFFICIENT_BYTES,
                    ..CandidateBundleLimits::default()
                },
                ReductionLimits::default(),
            )
            .unwrap();
            let plan = TerminalNormalizationPlan::decode_generated(
                input.normalization,
                &family,
                reducer.terminals(),
                reducer.ordering(),
                TerminalNormalizationLimits::default(),
                BinaryIoLimits::default(),
            )
            .unwrap();
            assert_eq!(plan.raw_terminals(), reducer.terminals());
            assert_eq!(plan.raw_terminals().len(), raw_count, "{}", input.label);
            assert_eq!(
                plan.canonical_terminals().len(),
                output_count,
                "{}",
                input.label
            );
            assert!(
                plan.canonical_terminals()
                    .iter()
                    .all(|key| key.powers().iter().all(|power| *power >= 0))
            );
            assert!(plan.canonical_terminals().is_subset(plan.raw_terminals()));
        }
    }

    #[test]
    fn public_four_loop_raw_master_mode_does_not_expand_the_finite_master_table() {
        use crate::utils::vakint_macros::vk_parse;
        use crate::{EvaluationMethod, EvaluationOrder, LoopNormalizationFactor};

        let vakint = Vakint::new().unwrap();
        let descriptors = DESCRIPTORS.as_ref().unwrap();
        let mut propagators = Atom::num(1);
        for (slot, line) in INPUTS[1]
            .csv
            .lines()
            .map(str::trim)
            .filter(|line| !line.is_empty() && !line.starts_with('#'))
            .enumerate()
        {
            let row = line
                .split(',')
                .map(|entry| entry.trim().parse::<i64>().unwrap())
                .collect::<Vec<_>>();
            if row[0] == 0 && row[1] == 0 {
                continue;
            }
            propagators *= vk_parse!(format!(
                "prop({},edge({},{}),{},muvsq,1)",
                slot + 1,
                row[0],
                row[1],
                descriptors[1].physical_momenta[slot].to_canonical_string()
            ))
            .unwrap();
        }
        let input = function!(S.topo, propagators);
        let settings = VakintSettings {
            form_exe_path: "/four-loop-raw-master-mode-must-not-use-form".into(),
            number_of_terms_in_epsilon_expansion: 7,
            integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD,
            evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
                RustRedEvaluationOptions {
                    substitute_masters: false,
                },
            )]),
            ..VakintSettings::default()
        };
        let result = vakint
            .evaluate_integral(&settings, input.as_view())
            .unwrap();
        let text = result.to_canonical_string();
        assert!(
            text.contains("PR"),
            "raw FMFT master basis expected: {text}"
        );
        assert!(
            !text.contains("Oep") && !text.contains("PRep"),
            "master table must stay unexpanded: {text}"
        );
    }
}
