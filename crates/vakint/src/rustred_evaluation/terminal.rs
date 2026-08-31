//! Immutable terminal-value catalogs paired with authenticated RustRed artifacts.
//!
//! RustRed owns reduction closure and typed terminal keys. Vakint separately
//! owns the values used to materialize those terminals. Keeping that boundary
//! explicit lets a certified artifact close on a finite, nonminimal basis:
//! terminals may either project exactly to MATAD's symbolic basis or carry a
//! finite numerical Laurent record with an explicit unknown-tail sentinel.

use std::collections::BTreeMap;

use rustred::family::IntegralKey;
use rustred::foundry::artifact::{ArtifactSchemaVersion, ClosedArtifact};
use symbolica::atom::{Atom, AtomCore};

use crate::utils::vakint_macros::{vk_parse, vk_symbol};

use super::RustRedEvaluationError;
use super::artifact::ArtifactFamily;

/// Source record compiled once while the corresponding artifact is loaded.
pub(super) struct TerminalSource<'source> {
    powers: &'source [i64],
    value: TerminalValueSource<'source>,
}

/// Immutable consumer manifest for one exact artifact/catalog pair.
pub(super) struct TerminalManifest<'source> {
    expected_schema: ArtifactSchemaVersion,
    expected_algorithm_id: &'source str,
    expected_family_fingerprint: &'source str,
    sources: &'source [TerminalSource<'source>],
}

impl<'source> TerminalManifest<'source> {
    pub(super) const fn new(
        expected_schema: ArtifactSchemaVersion,
        expected_algorithm_id: &'source str,
        expected_family_fingerprint: &'source str,
        sources: &'source [TerminalSource<'source>],
    ) -> Self {
        Self {
            expected_schema,
            expected_algorithm_id,
            expected_family_fingerprint,
            sources,
        }
    }
}

enum TerminalValueSource<'source> {
    /// Exact unit-mass Minkowski value in Vakint's pre-normalization MATAD basis.
    ExactMatadBasis(&'source str),
    /// Consecutive coefficients beginning at `leading_power`.
    ///
    /// The compiler appends `ep^(last+1)*Oep(...)`. Vakint's existing master
    /// finalizer therefore fails closed if a reduction coefficient exposes an
    /// epsilon order that was not shipped, including through spurious poles.
    NumericalLaurent {
        leading_power: i64,
        coefficients: &'source [&'source str],
    },
}

impl<'source> TerminalSource<'source> {
    pub(super) const fn exact_matad_basis(powers: &'source [i64], value: &'source str) -> Self {
        Self {
            powers,
            value: TerminalValueSource::ExactMatadBasis(value),
        }
    }

    #[allow(dead_code)]
    pub(super) const fn numerical_laurent(
        powers: &'source [i64],
        leading_power: i64,
        coefficients: &'source [&'source str],
    ) -> Self {
        Self {
            powers,
            value: TerminalValueSource::NumericalLaurent {
                leading_power,
                coefficients,
            },
        }
    }
}

#[derive(Debug)]
enum TerminalValue {
    ExactMatadBasis(Atom),
    NumericalLaurent { symbolic: Atom, series: Atom },
}

impl TerminalValue {
    fn expression(&self, substitute_masters: bool) -> Atom {
        match self {
            Self::ExactMatadBasis(value) => value.clone(),
            Self::NumericalLaurent { symbolic, series } => {
                if substitute_masters {
                    series.clone()
                } else {
                    symbolic.clone()
                }
            }
        }
    }
}

/// Exact terminal-key cover paired with one authenticated artifact.
#[derive(Debug)]
pub(super) struct TerminalCatalog {
    values: BTreeMap<IntegralKey, TerminalValue>,
}

impl TerminalCatalog {
    pub(super) fn compile(
        artifact: &ClosedArtifact,
        manifest: &TerminalManifest<'_>,
    ) -> Result<Self, String> {
        if artifact.schema() != manifest.expected_schema {
            return Err(format!(
                "artifact schema {} does not match terminal manifest {}",
                artifact.schema().stable_id(),
                manifest.expected_schema.stable_id()
            ));
        }
        if artifact.algorithm_id() != manifest.expected_algorithm_id {
            return Err(format!(
                "artifact algorithm {:?} does not match terminal manifest {:?}",
                artifact.algorithm_id(),
                manifest.expected_algorithm_id
            ));
        }
        if artifact.family_fingerprint() != manifest.expected_family_fingerprint {
            return Err(format!(
                "artifact family fingerprint does not match terminal manifest (actual={:?}, expected={:?})",
                artifact.family_fingerprint(),
                manifest.expected_family_fingerprint
            ));
        }

        let mut values = BTreeMap::new();
        for source in manifest.sources {
            let key = IntegralKey::try_new(source.powers.iter().copied())
                .map_err(|error| format!("invalid terminal key {:?}: {error}", source.powers))?;
            if key.powers().len() != artifact.arity() {
                return Err(format!(
                    "terminal {:?} has arity {}, expected {}",
                    key.powers(),
                    key.powers().len(),
                    artifact.arity()
                ));
            }
            let value = compile_value(&key, &source.value)?;
            if values.insert(key.clone(), value).is_some() {
                return Err(format!("duplicate terminal value for {:?}", key.powers()));
            }
        }

        if values.len() != artifact.masters().len() || values.keys().ne(artifact.masters().iter()) {
            let missing = artifact
                .masters()
                .iter()
                .filter(|master| !values.contains_key(*master))
                .map(|master| master.powers().to_vec())
                .collect::<Vec<_>>();
            let foreign = values
                .keys()
                .filter(|master| !artifact.masters().contains(*master))
                .map(|master| master.powers().to_vec())
                .collect::<Vec<_>>();
            return Err(format!(
                "terminal catalog does not exactly cover the artifact masters (missing={missing:?}, foreign={foreign:?})"
            ));
        }

        Ok(Self { values })
    }

    pub(super) fn materialize(
        &self,
        family: ArtifactFamily,
        master: &IntegralKey,
        loop_count: usize,
        mass_squared: &Atom,
        substitute_masters: bool,
    ) -> Result<Atom, RustRedEvaluationError> {
        let value =
            self.values
                .get(master)
                .ok_or_else(|| RustRedEvaluationError::UnsupportedMaster {
                    family: family.name(),
                    powers: master.powers().to_vec(),
                })?;

        // The per-loop normalization later supplies (m^2)^(-L*ep). This
        // classical factor, combined with RustRed's homogeneous reduction
        // ratio (m^2)^(sum(master)-sum(target)), yields the target dimension.
        let twice_loop_count = i128::try_from(loop_count)
            .ok()
            .and_then(|loops| loops.checked_mul(2))
            .ok_or_else(|| RustRedEvaluationError::Reduction {
                detail: "loop count overflowed terminal mass homogeneity".to_owned(),
            })?;
        let master_power_sum = master
            .powers()
            .iter()
            .try_fold(0_i128, |sum, &power| sum.checked_add(i128::from(power)));
        let exponent = master_power_sum
            .and_then(|sum| twice_loop_count.checked_sub(sum))
            .ok_or_else(|| RustRedEvaluationError::Reduction {
                detail: "terminal powers overflowed common-mass homogeneity".to_owned(),
            })?;
        let exponent = i64::try_from(exponent)
            .map_err(|_| RustRedEvaluationError::MassExponentOverflow { exponent })?;

        Ok(mass_squared.clone().pow(Atom::num(exponent)) * value.expression(substitute_masters))
    }
}

fn compile_value(
    key: &IntegralKey,
    source: &TerminalValueSource<'_>,
) -> Result<TerminalValue, String> {
    match source {
        TerminalValueSource::ExactMatadBasis(value) => vk_parse!(value)
            .map(TerminalValue::ExactMatadBasis)
            .map_err(|error| format!("could not parse terminal {:?}: {error}", key.powers())),
        TerminalValueSource::NumericalLaurent {
            leading_power,
            coefficients,
        } => {
            if coefficients.is_empty() {
                return Err(format!(
                    "numerical terminal {:?} has no Laurent coefficients",
                    key.powers()
                ));
            }
            let epsilon = vk_parse!("ep").expect("Vakint's private epsilon symbol parses");
            let mut series = Atom::Zero;
            for (offset, coefficient) in coefficients.iter().enumerate() {
                let offset = i64::try_from(offset).map_err(|_| {
                    format!(
                        "too many Laurent coefficients for terminal {:?}",
                        key.powers()
                    )
                })?;
                let power = leading_power.checked_add(offset).ok_or_else(|| {
                    format!("Laurent power overflow for terminal {:?}", key.powers())
                })?;
                let coefficient = vk_parse!(coefficient).map_err(|error| {
                    format!(
                        "could not parse epsilon^{power} coefficient for terminal {:?}: {error}",
                        key.powers()
                    )
                })?;
                if coefficient.contains_symbol(vk_symbol!("ep"))
                    || coefficient.contains_symbol(vk_symbol!("Oep"))
                {
                    return Err(format!(
                        "epsilon^{power} coefficient for terminal {:?} must not contain ep or Oep",
                        key.powers()
                    ));
                }
                series += coefficient * epsilon.clone().pow(Atom::num(power));
            }
            let coefficient_count = i64::try_from(coefficients.len()).map_err(|_| {
                format!(
                    "too many Laurent coefficients for terminal {:?}",
                    key.powers()
                )
            })?;
            let unknown_power = leading_power
                .checked_add(coefficient_count)
                .ok_or_else(|| {
                    format!(
                        "Laurent tail power overflow for terminal {:?}",
                        key.powers()
                    )
                })?;
            let powers = key
                .powers()
                .iter()
                .map(i64::to_string)
                .collect::<Vec<_>>()
                .join(",");
            let symbolic =
                vk_parse!(format!("RustRedMaster({powers})").as_str()).map_err(|error| {
                    format!(
                        "could not construct symbolic terminal {:?}: {error}",
                        key.powers()
                    )
                })?;
            let unknown =
                vk_parse!(format!("Oep({unknown_power},RustRedMaster({powers}))").as_str())
                    .map_err(|error| {
                        format!(
                            "could not construct Laurent tail for terminal {:?}: {error}",
                            key.powers()
                        )
                    })?;
            series += epsilon.pow(Atom::num(unknown_power)) * unknown;
            Ok(TerminalValue::NumericalLaurent { symbolic, series })
        }
    }
}

#[cfg(test)]
mod tests {
    use rustred::foundry::artifact::derive_one_loop_unit_mass_tadpole;
    use symbolica::atom::AtomCore;

    use rustred::foundry::artifact::ArtifactSchemaVersion;

    use super::{TerminalCatalog, TerminalManifest, TerminalSource};
    use crate::utils::vakint_macros::vk_parse;

    #[test]
    fn catalog_requires_an_exact_cover_of_authenticated_terminals() {
        let artifact = derive_one_loop_unit_mass_tadpole().unwrap();
        let missing_manifest = TerminalManifest::new(
            ArtifactSchemaVersion::V3,
            artifact.algorithm_id(),
            artifact.family_fingerprint(),
            &[],
        );
        let missing = TerminalCatalog::compile(&artifact, &missing_manifest).unwrap_err();
        assert!(missing.contains("missing=[[1]]"), "{missing}");

        let foreign_sources = [
            TerminalSource::exact_matad_basis(&[1], "1"),
            TerminalSource::exact_matad_basis(&[2], "1"),
        ];
        let foreign_manifest = TerminalManifest::new(
            ArtifactSchemaVersion::V3,
            artifact.algorithm_id(),
            artifact.family_fingerprint(),
            &foreign_sources,
        );
        let foreign = TerminalCatalog::compile(&artifact, &foreign_manifest).unwrap_err();
        assert!(foreign.contains("foreign=[[2]]"), "{foreign}");
    }

    #[test]
    fn catalog_manifest_rejects_wrong_artifact_identity() {
        let artifact = derive_one_loop_unit_mass_tadpole().unwrap();
        let sources = [TerminalSource::exact_matad_basis(&[1], "1")];
        let wrong_algorithm = TerminalManifest::new(
            ArtifactSchemaVersion::V3,
            "rustred.generated.some-other-family.v1",
            artifact.family_fingerprint(),
            &sources,
        );
        let error = TerminalCatalog::compile(&artifact, &wrong_algorithm).unwrap_err();
        assert!(error.contains("algorithm"), "{error}");

        let wrong_family = TerminalManifest::new(
            ArtifactSchemaVersion::V3,
            artifact.algorithm_id(),
            "rustred-integral-family-v2;wrong",
            &sources,
        );
        let error = TerminalCatalog::compile(&artifact, &wrong_family).unwrap_err();
        assert!(error.contains("family fingerprint"), "{error}");
    }

    #[test]
    fn numerical_catalogs_expose_symbolic_raw_masters_and_guard_unknown_orders() {
        let artifact = derive_one_loop_unit_mass_tadpole().unwrap();
        let sources = [TerminalSource::numerical_laurent(
            &[1],
            -1,
            &["2", "3", "5"],
        )];
        let manifest = TerminalManifest::new(
            ArtifactSchemaVersion::V3,
            artifact.algorithm_id(),
            artifact.family_fingerprint(),
            &sources,
        );
        let catalog = TerminalCatalog::compile(&artifact, &manifest).unwrap();
        let master = artifact.masters().first().unwrap();
        let mass = vk_parse!("m2").unwrap();

        let raw = catalog
            .materialize(
                super::ArtifactFamily::UnitMassVacuumK1,
                master,
                1,
                &mass,
                false,
            )
            .unwrap();
        assert!(raw.to_canonical_string().contains("RustRedMaster(1)"));

        let numerical = catalog
            .materialize(
                super::ArtifactFamily::UnitMassVacuumK1,
                master,
                1,
                &mass,
                true,
            )
            .unwrap();
        let rendered = numerical.to_canonical_string();
        assert!(
            rendered.contains("Oep(2,") && rendered.contains("RustRedMaster(1)"),
            "{rendered}"
        );
        assert!(rendered.contains("m2"), "{rendered}");
    }

    #[test]
    fn numerical_catalog_coefficients_cannot_embed_epsilon_or_unknown_tails() {
        let artifact = derive_one_loop_unit_mass_tadpole().unwrap();
        for coefficient in ["1 + ep", "Oep(0,hidden)"] {
            let coefficients = [coefficient];
            let sources = [TerminalSource::numerical_laurent(&[1], 0, &coefficients)];
            let manifest = TerminalManifest::new(
                ArtifactSchemaVersion::V3,
                artifact.algorithm_id(),
                artifact.family_fingerprint(),
                &sources,
            );
            let error = TerminalCatalog::compile(&artifact, &manifest).unwrap_err();
            assert!(error.contains("must not contain ep or Oep"), "{error}");
        }
    }
}
