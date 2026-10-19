//! Shipped, finite-tested four-loop candidate programs.
//!
//! These are not `ClosedArtifact`s: a requested point succeeds only when
//! RustRed's candidate applier finds guarded, descending rules ending in the
//! declared terminal catalog. Unknown leaves remain errors. Discovery and
//! FMFT terminal preparation are offline tasks, never runtime fallbacks.
//! Native Symbolica programs are trusted, build-time embedded data from the
//! pinned RustRed/Symbolica stack, not a public untrusted-file input boundary.

use std::sync::Arc;
pub(super) mod cache;

use rustred::input::{Compiler, Limits, LoweringLimits, TextProject, TextPropagator};
use rustred::persistence::BinaryIoLimits;
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
fn catalog_limits() -> BinaryIoLimits {
    BinaryIoLimits {
        max_program_bytes: 64 * 1024,
        max_collection_entries: 4096,
        max_state_bytes: 32 * 1024,
        max_atom_bytes: 32 * 1024,
        max_total_atom_bytes: 32 * 1024,
        ..BinaryIoLimits::default()
    }
}

struct Input {
    label: &'static str,
    csv: String,
    program: Vec<u8>,
    catalog: Vec<u8>,
    normalization: Vec<u8>,
}

struct Descriptor {
    physical_momenta: Vec<Atom>,
    family_fingerprint: String,
}

// Independent lazy cells avoid loading all large programs when a calculation
// only visits one parent. NativeCandidate owns one shared memoized applier.
struct LoadedProgram {
    native: Arc<NativeCandidate<10>>,
    evaluator: ExperimentalRustRed,
    cache_program: rustred::solver::CandidateUnitMassProgram<10>,
}

pub(super) struct FourLoopPrograms {
    #[cfg(test)]
    inputs: Vec<Input>,
    descriptors: Vec<Descriptor>,
    programs: Vec<LoadedProgram>,
}

impl FourLoopPrograms {
    pub(super) fn load(data: &super::data::LoopData) -> Result<Self, RustRedEvaluationError> {
        Vakint::initialize_vakint_symbols();
        let mut inputs = Vec::new();
        for (name, label) in [
            ("h", "four-loop-h"),
            ("fg", "four-loop-fg"),
            ("bmw", "four-loop-bmw"),
            ("x", "four-loop-x"),
        ] {
            inputs.push(Input {
                label,
                csv: String::from_utf8(data.read(&format!("ibps/{name}.csv"))?)
                    .map_err(super::data::failure)?,
                program: data.read(&format!("ibps/{name}.candidates.rrbin.gz"))?,
                catalog: data.read(&format!("ibps/{name}.rrcat.bin.gz"))?,
                normalization: data.read(&format!("ibps/{name}.rrnorm.bin"))?,
            });
        }
        let descriptors = inputs
            .iter()
            .map(|input| descriptor(&input.csv))
            .collect::<Result<Vec<_>, _>>()
            .map_err(super::data::failure)?;
        let programs = inputs
            .iter()
            .zip(&descriptors)
            .map(|(input, descriptor)| load(input, descriptor))
            .collect::<Result<Vec<_>, _>>()
            .map_err(super::data::failure)?;
        Ok(Self {
            #[cfg(test)]
            inputs,
            descriptors,
            programs,
        })
    }

    pub(super) fn native_program(&self, index: usize) -> Option<(Arc<NativeCandidate<10>>, usize)> {
        self.programs.get(index).map(|program| {
            (
                Arc::clone(&program.native),
                self.descriptors[index].physical_momenta.len(),
            )
        })
    }

    pub(super) fn physical_momenta(&self, index: usize) -> Option<&[Atom]> {
        self.descriptors
            .get(index)
            .map(|descriptor| descriptor.physical_momenta.as_slice())
    }
}

impl Vakint {
    pub(super) fn four_loop_data(
        &self,
    ) -> Result<Arc<super::data::RustRedData>, RustRedEvaluationError> {
        let data = self.rustred_data()?;
        if data.four_loop.is_none() {
            return Err(super::data::failure(
                "4L RustRed data is not loaded; call load_official_external_data(4, ...) or load_external_data(path, 4)",
            ));
        }
        Ok(data)
    }
}

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

fn load(input: &Input, descriptor: &Descriptor) -> Result<LoadedProgram, String> {
    let bytes = super::experimental::compressed::decode(&input.program, MAX_PROGRAM_BYTES)
        .map_err(|error| format!("decompress {}: {error}", input.label))?;
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
    let catalog = OfflineTerminalCatalog::decode_generated_compressed(
        &input.catalog,
        family.fingerprint(),
        family.denominator_count(),
        catalog_limits(),
    )?;
    let mut native = NativeCandidate::from_reducer_with_terminal_normalization(
        Arc::new(family),
        descriptor.physical_momenta.clone(),
        reducer,
        &input.normalization,
    )?;
    let cache_program =
        native.cache_program(cache::source_digest(&bytes, &input.normalization)?)?;
    catalog.require_complete()?;
    catalog.require_exact_keys(native.output_terminals())?;
    native.preserve_unresolved()?;
    let native = Arc::new(native);
    let evaluator = ExperimentalRustRed::new(native.clone(), catalog.into_terms())
        .map_err(|error| error.to_string())?;
    Ok(LoadedProgram {
        native,
        evaluator,
        cache_program,
    })
}

fn select(
    programs: &FourLoopPrograms,
    topology: &Topology,
) -> Result<usize, RustRedEvaluationError> {
    let integral = topology.get_integral();
    if !matches!(topology, Topology::FourLoop(_)) || integral.n_loops != 4 {
        return Err(invalid("not a matched four-loop family"));
    }
    let descriptors = &programs.descriptors;
    let index = descriptors
        .iter()
        .position(|descriptor| {
            integral
                .validate_parent_routing(&descriptor.physical_momenta)
                .is_ok()
        })
        .ok_or_else(|| invalid("no shipped program matches the retained parent witness"))?;
    validate_concrete(topology)?;
    Ok(index)
}

fn validate_concrete(topology: &Topology) -> Result<(), RustRedEvaluationError> {
    let integral = topology.get_integral();
    if !matches!(topology, Topology::FourLoop(_)) || integral.n_loops != 4 {
        return Err(invalid("not a matched four-loop family"));
    }
    let (parent, coordinates) = integral
        .parent_routing
        .as_ref()
        .ok_or_else(|| invalid("missing parent routing"))?;
    if coordinates.len() != 4 {
        return Err(invalid("invalid four-loop routing dimension"));
    }
    let momenta = (1..=integral.n_props)
        .map(|slot| {
            get_prop_with_id(parent.as_view(), slot)
                .and_then(|row| row.get(&vk_symbol!("q_")).cloned())
                .ok_or_else(|| invalid("missing defining-parent propagator"))
        })
        .collect::<Result<Vec<_>, _>>()?;
    integral.validate_parent_routing(&momenta)?;
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
    Ok(())
}

fn invalid(detail: impl Into<String>) -> RustRedEvaluationError {
    RustRedEvaluationError::InvalidMatchedFamily {
        detail: detail.into(),
    }
}

pub(super) fn supports(topology: &Topology) -> bool {
    validate_concrete(topology).is_ok()
}

pub(super) fn evaluate(
    vakint: &Vakint,
    settings: &VakintSettings,
    numerator: AtomView,
    matched: &ReplacementRules,
    _options: &RustRedEvaluationOptions,
) -> Result<Atom, VakintError> {
    let data = vakint.four_loop_data()?;
    let programs = data.four_loop.as_ref().unwrap();
    let index = select(programs, &matched.canonical_topology)?;
    let program = &programs.programs[index];
    log::warn!(
        "RustRed 4L uses finite-tested candidate rules; there is no guarantee of a complete reduction to the known masters outside the documented validation workload"
    );
    program.evaluator.evaluate_catalog(
        settings,
        numerator,
        matched,
        &programs.descriptors[index].family_fingerprint,
        &data.masters,
        |target| vakint.cached_candidate_reduction(index, target),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use symbolica::atom::AtomCore;

    #[test]
    fn shipped_descriptors_have_distinct_family_bindings() {
        let engine = super::super::cache_tests::external_engine();
        let data = engine.four_loop_data().unwrap();
        let descriptors = &data.four_loop.as_ref().unwrap().descriptors;
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
        let engine = super::super::cache_tests::external_engine();
        let first_data = engine.four_loop_data().unwrap();
        let second_data = engine.clone().four_loop_data().unwrap();
        for index in 0..4 {
            let first = &first_data.four_loop.as_ref().unwrap().programs[index];
            let second = &second_data.four_loop.as_ref().unwrap().programs[index];
            assert!(std::ptr::eq(first, second));
        }
    }

    #[test]
    fn shipped_normalizations_keep_raw_keys_and_close_onto_positive_outputs() {
        use rustred::persistence::BinaryIoLimits;
        use rustred::reduction::terminal_normalization::{
            TerminalNormalizationLimits, TerminalNormalizationPlan,
        };

        // These are expectations for vendored input data, not engine limits or
        // a dispatch based on a topology name. Original rules stay raw; catalog
        // values cover only the exact output convention installed by RustRed.
        let engine = super::super::cache_tests::external_engine();
        let data = engine.four_loop_data().unwrap();
        for (input, (raw_count, output_count)) in data
            .four_loop
            .as_ref()
            .unwrap()
            .inputs
            .iter()
            .zip([(386, 22), (145, 16), (179, 17), (445, 19)])
        {
            let bytes =
                super::super::experimental::compressed::decode(&input.program, MAX_PROGRAM_BYTES)
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
                &input.normalization,
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
            let catalog = OfflineTerminalCatalog::decode_generated_compressed(
                &input.catalog,
                family.fingerprint(),
                family.denominator_count(),
                catalog_limits(),
            )
            .unwrap();
            assert_eq!(catalog.terms().len(), output_count);
            assert!(catalog.terms().keys().eq(plan.canonical_terminals().iter()));
            assert!(catalog.terms().len() < raw_count);
        }
    }

    #[test]
    fn public_four_loop_raw_master_mode_does_not_expand_the_finite_master_table() {
        use crate::utils::vakint_macros::vk_parse;
        use crate::{EvaluationMethod, EvaluationOrder, LoopNormalizationFactor};

        let vakint = super::super::cache_tests::external_engine();
        let input = cache::test_input_mass(
            1,
            &[1, 1, 1, 1, 1, 1, 1, 1, 0, 0],
            &vk_parse!("muvsq").unwrap(),
        );
        let settings = VakintSettings {
            form_exe_path: "/four-loop-raw-master-mode-must-not-use-form".into(),
            number_of_terms_in_epsilon_expansion: 7,
            integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD,
            evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
                RustRedEvaluationOptions {
                    master_mode: crate::RustRedMasterMode::Retain,
                },
            )]),
            ..VakintSettings::default()
        };
        let result = vakint
            .evaluate_integral(&settings, input.as_view())
            .unwrap();
        let text = result.to_canonical_string();
        assert!(
            text.contains("RustRedMaster"),
            "raw family-qualified master basis expected: {text}"
        );
        assert!(
            !text.contains("Oep") && !text.contains("PRep"),
            "master table must stay unexpanded: {text}"
        );
    }
}
