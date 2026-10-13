//! Four-loop cache steering using authenticated shipped native programs.
//! Rule application, terminal normalization and cache policy stay in RustRed.
use super::*;
use rustred::solver::CandidateUnitMassProgram;
use sha2::{Digest, Sha256};
#[cfg(test)]
use symbolica::atom::AtomCore;

/// Pair the native program with the exact decoded bundle and installed
/// normalization bytes consumed by its loader. The framing excludes ambiguous
/// concatenations and does not involve numerical terminal-value tables.
pub(super) fn source_digest(program: &[u8], normalization: &[u8]) -> Result<[u8; 32], String> {
    let mut digest = Sha256::new();
    digest.update(b"vakint-native-candidate-and-normalization-v1\0");
    for bytes in [program, normalization] {
        let length = u64::try_from(bytes.len()).map_err(|error| error.to_string())?;
        digest.update(length.to_le_bytes());
        digest.update(bytes);
    }
    Ok(digest.finalize().into())
}

/// This is an explicit prepare/load boundary, not an ordinary statistics or
/// clear-cache helper: authenticating the exact bundle may load all programs.
pub(in crate::rustred_evaluation) fn programs()
-> Result<Vec<CandidateUnitMassProgram<10>>, RustRedEvaluationError> {
    (0..INPUTS.len())
        .map(|index| {
            loaded(index)
                .as_ref()
                .map(|loaded| loaded.cache_program.clone())
                .map_err(|detail| RustRedEvaluationError::ArtifactLoad {
                    family: INPUTS[index].label,
                    detail: detail.clone(),
                })
        })
        .collect()
}

pub(in crate::rustred_evaluation) fn parent_masks()
-> Result<Vec<(usize, Vec<bool>)>, RustRedEvaluationError> {
    DESCRIPTORS
        .as_ref()
        .map_err(|detail| invalid(detail.clone()))?
        .iter()
        .enumerate()
        .map(|(index, descriptor)| {
            // The existing descriptor validator authenticates this prefix;
            // auxiliary rows are explicitly excluded from positive support.
            Ok((
                index,
                (0..10)
                    .map(|slot| slot < descriptor.physical_momenta.len())
                    .collect(),
            ))
        })
        .collect()
}

#[cfg(test)]
pub(in crate::rustred_evaluation) fn test_input(
    namespace: usize,
    powers: &[i64],
    mass: i64,
) -> Atom {
    test_input_mass(namespace, powers, &Atom::num(mass))
}
#[cfg(test)]
pub(in crate::rustred_evaluation) fn test_input_mass(
    namespace: usize,
    powers: &[i64],
    mass: &Atom,
) -> Atom {
    use crate::utils::vakint_macros::vk_parse;
    assert_eq!(powers.len(), 10);
    Vakint::initialize_vakint_symbols();
    let descriptors = DESCRIPTORS.as_ref().unwrap();
    let mut propagators = Atom::num(1);
    for (slot, line) in INPUTS[namespace]
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
            assert_eq!(powers[slot], 0);
            continue;
        }
        propagators *= vk_parse!(format!(
            "prop({},edge({},{}),{},{},{})",
            slot + 1,
            row[0],
            row[1],
            descriptors[namespace].physical_momenta[slot].to_canonical_string(),
            mass.to_canonical_string(),
            powers[slot]
        ))
        .unwrap();
    }
    function!(S.topo, propagators)
}

#[cfg(test)]
pub(in crate::rustred_evaluation) fn test_denominator(
    namespace: usize,
    slot: usize,
    mass: i64,
) -> Atom {
    let q = &DESCRIPTORS.as_ref().unwrap()[namespace].physical_momenta[slot];
    function!(S.dot, q, q) - Atom::num(mass)
}

#[cfg(test)]
pub(in crate::rustred_evaluation) fn ordinary_rule_applications(namespace: usize) -> usize {
    loaded(namespace)
        .as_ref()
        .unwrap()
        .native
        .rule_applications()
        .unwrap()
}
#[cfg(test)]
pub(in crate::rustred_evaluation) fn ordinary_statistics(
    namespace: usize,
) -> rustred::solver::CandidateStatistics {
    loaded(namespace)
        .as_ref()
        .unwrap()
        .native
        .cache_test_statistics()
}

#[cfg(test)]
mod tests {
    use super::*;
    use rustred_app::{CachePreparationOptions, CachePreparationWorkload};

    #[test]
    fn actual_four_parent_primer_fits_default_but_refuses_explicit_smaller_cap() {
        let bindings = programs()
            .unwrap()
            .iter()
            .map(|p| p.binding().clone())
            .collect::<Vec<_>>();
        let masks = parent_masks().unwrap();
        assert_eq!(
            masks
                .iter()
                .map(|(_, mask)| mask.iter().filter(|&&p| p).count())
                .collect::<Vec<_>>(),
            [9, 8, 8, 9]
        );
        assert!(
            CachePreparationWorkload::parent_pinches(
                bindings.clone(),
                &masks,
                Default::default(),
                100_000,
            )
            .is_err()
        );
        let options = CachePreparationOptions::default();
        let workload = CachePreparationWorkload::parent_pinches(
            bindings,
            &masks,
            Default::default(),
            options.max_roots,
        )
        .unwrap();
        assert_eq!(workload.targets().len(), 164_352);
        #[cfg(target_pointer_width = "64")]
        assert_eq!(
            options.resolve(4, Default::default()).unwrap().ram_bytes,
            16_000_000_000
        );
        // Enumeration/admission shape only: no reduction or new rules.
    }

    #[test]
    fn program_and_normalization_identity_is_unambiguously_framed() {
        assert_ne!(
            source_digest(b"ab", b"c").unwrap(),
            source_digest(b"a", b"bc").unwrap()
        );
        assert_ne!(
            source_digest(b"program", b"norm1").unwrap(),
            source_digest(b"program", b"norm2").unwrap()
        );
        assert_eq!(
            source_digest(b"program", b"norm").unwrap(),
            source_digest(b"program", b"norm").unwrap()
        );
    }
}
