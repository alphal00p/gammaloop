use std::{
    collections::BTreeSet,
    fmt::{Display, Formatter},
};

use crate::utils::{
    GS,
    serde_utils::{
        IsDefault, is_default_form_path, is_default_pysecdec_relative_precision,
        is_default_python_path, is_default_vakint_evaluation_methods,
        is_default_vakint_normalization, is_false, is_minus_one_string, is_true, is_usize,
    },
};
use bincode_trait_derive::{Decode, Encode};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use vakint::{AlphaLoopOptions, EvaluationMethod, FMFTOptions, MATADOptions, PySecDecOptions};

#[derive(
    Debug,
    Clone,
    Copy,
    Serialize,
    Deserialize,
    Encode,
    Decode,
    PartialEq,
    Eq,
    Hash,
    JsonSchema,
    Default,
)]
#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
#[serde(deny_unknown_fields)]
pub enum ApproximationType {
    #[default]
    #[serde(rename = "MUV", alias = "muv")]
    MUV,
    #[serde(rename = "PolePart", alias = "pole_part")]
    PolePart,
    #[serde(rename = "OS", alias = "os")]
    OS,
    #[serde(rename = "IR", alias = "ir")]
    IR,
    #[serde(rename = "Unsubtracted", alias = "unsubtracted")]
    Unsubtracted,
    #[serde(rename = "VaccuumLimit", alias = "vaccuum_limit")]
    VaccuumLimit,
}

impl Display for ApproximationType {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            ApproximationType::MUV => write!(f, "MUV"),
            ApproximationType::PolePart => write!(f, "PolePart"),
            ApproximationType::OS => write!(f, "OS"),
            ApproximationType::IR => write!(f, "IR"),
            ApproximationType::Unsubtracted => write!(f, "Unsubtracted"),
            ApproximationType::VaccuumLimit => write!(f, "VaccuumLimit"),
        }
    }
}

#[derive(
    Debug,
    Clone,
    Copy,
    Serialize,
    Deserialize,
    Encode,
    Decode,
    PartialEq,
    Eq,
    Hash,
    JsonSchema,
    Default,
)]
#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
#[serde(deny_unknown_fields)]
pub enum UVOrchestrator {
    #[serde(rename = "legacy_dag_forest", alias = "LegacyDagForest")]
    LegacyDagForest,
    #[default]
    #[serde(rename = "hedge_poset", alias = "HedgePoset")]
    HedgePoset,
    #[serde(rename = "compare", alias = "Compare")]
    Compare,
}

impl Display for UVOrchestrator {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            UVOrchestrator::LegacyDagForest => write!(f, "legacy_dag_forest"),
            UVOrchestrator::HedgePoset => write!(f, "hedge_poset"),
            UVOrchestrator::Compare => write!(f, "compare"),
        }
    }
}

#[cfg_attr(
    feature = "python_api",
    pyo3::pyclass(from_py_object, get_all, set_all)
)]
#[derive(
    Debug,
    Clone,
    Serialize,
    Deserialize,
    Encode,
    Decode,
    Default,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Hash,
    JsonSchema,
)]
#[serde(default, deny_unknown_fields)]
pub struct CTIdentifier {
    #[serde(skip_serializing_if = "BTreeSet::is_empty")]
    pub external_pdg_set: BTreeSet<isize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub internal_pdg_set: Option<BTreeSet<isize>>,
}

impl CTIdentifier {
    pub fn new(
        external_pdg_set: BTreeSet<isize>,
        internal_pdg_set: Option<BTreeSet<isize>>,
    ) -> Self {
        Self {
            external_pdg_set,
            internal_pdg_set,
        }
    }

    pub fn matches(&self, candidate: &Self) -> bool {
        self.external_pdg_set == candidate.external_pdg_set
            && self
                .internal_pdg_set
                .as_ref()
                .is_none_or(|internal| candidate.internal_pdg_set.as_ref() == Some(internal))
    }
}

#[cfg_attr(
    feature = "python_api",
    pyo3::pyclass(from_py_object, get_all, set_all)
)]
#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, Eq, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CTRenormalizationRule {
    pub ct_identifier: CTIdentifier,
    pub prescription: ApproximationType,
}

impl CTRenormalizationRule {
    pub fn new(ct_identifier: CTIdentifier, prescription: ApproximationType) -> Self {
        Self {
            ct_identifier,
            prescription,
        }
    }
}

#[cfg_attr(
    feature = "python_api",
    pyo3::pyclass(from_py_object, get_all, set_all)
)]
#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct RenormalizationPrescriptionSettings {
    /// Approximation applied to logarithmically divergent ultraviolet counterterms.
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub log_divergent: ApproximationType,
    /// Approximation applied to power-divergent counterterms with massive external structure.
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub massive_power_divergent: ApproximationType,
    /// Approximation applied to power-divergent counterterms with massless external structure.
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub massless_power_divergent: ApproximationType,
    /// Counterterm-specific prescription rules evaluated before the divergence-class defaults.
    #[serde(
        default,
        skip_serializing_if = "IsDefault::is_default",
        with = "ct_renormalization_overrides_serde"
    )]
    #[schemars(with = "Vec<CTRenormalizationRule>")]
    pub overrides: Vec<CTRenormalizationRule>,
}

impl Default for RenormalizationPrescriptionSettings {
    fn default() -> Self {
        Self {
            log_divergent: ApproximationType::MUV,
            massive_power_divergent: ApproximationType::MUV,
            massless_power_divergent: ApproximationType::MUV,
            overrides: Vec::new(),
        }
    }
}

impl RenormalizationPrescriptionSettings {
    pub fn approximation_scheme_for(
        &self,
        ct_identifier: &CTIdentifier,
        dod: i32,
        has_massive_externals: bool,
    ) -> ApproximationType {
        if let Some(rule) = self
            .overrides
            .iter()
            .find(|rule| &rule.ct_identifier == ct_identifier)
        {
            return rule.prescription;
        }

        if let Some(rule) = self.overrides.iter().find(|rule| {
            rule.ct_identifier.internal_pdg_set.is_none()
                && rule.ct_identifier.matches(ct_identifier)
        }) {
            return rule.prescription;
        }

        if dod == 0 {
            self.log_divergent
        } else if has_massive_externals {
            self.massive_power_divergent
        } else {
            self.massless_power_divergent
        }
    }
}

#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct MATADSettings {
    /// Expand MATAD master-integral abbreviations before returning the result.
    #[serde(skip_serializing_if = "is_true")]
    pub expand_masters: bool,
    /// Substitute MATAD master integrals by their known analytic expressions.
    #[serde(skip_serializing_if = "is_true")]
    pub susbstitute_masters: bool,
    /// Substitute harmonic polylogarithms produced by MATAD where formulas are available.
    #[serde(skip_serializing_if = "is_true")]
    pub substitute_hpls: bool,
    /// Numerically substitute constants directly in the MATAD result.
    #[serde(skip_serializing_if = "is_true")]
    pub direct_numerical_substition: bool,
}

impl Default for MATADSettings {
    fn default() -> Self {
        Self {
            expand_masters: true,
            susbstitute_masters: true,
            substitute_hpls: true,
            direct_numerical_substition: true,
        }
    }
}

#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct AlphaLoopSettings {
    /// Substitute AlphaLoop master integrals by their registered expressions.
    #[serde(skip_serializing_if = "is_true")]
    pub susbstitute_masters: bool,
}

impl Default for AlphaLoopSettings {
    fn default() -> Self {
        Self {
            susbstitute_masters: true,
        }
    }
}

#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct FMFTSettings {
    /// Expand FMFT master-integral abbreviations before returning the result.
    #[serde(skip_serializing_if = "is_true")]
    pub expand_masters: bool,
    /// Substitute FMFT master integrals by their known analytic expressions.
    #[serde(skip_serializing_if = "is_true")]
    pub susbstitute_masters: bool,
}

impl Default for FMFTSettings {
    fn default() -> Self {
        Self {
            expand_masters: true,
            susbstitute_masters: true,
        }
    }
}
#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct PySecDecSettings {
    /// Suppress pySecDec progress output while Vakint evaluates an integral.
    #[serde(skip_serializing_if = "is_true")]
    pub quiet: bool,
    /// Requested relative precision for the pySecDec numerical integration.
    #[serde(skip_serializing_if = "is_default_pysecdec_relative_precision")]
    pub relative_precision: f64,
    /// Minimum pySecDec sample count before a precision-based stop is accepted.
    #[serde(skip_serializing_if = "is_usize::<10_000>")]
    pub min_n_evals: usize,
    /// Maximum pySecDec sample count before the numerical integration stops.
    #[serde(skip_serializing_if = "is_usize::<1_000_000_000_000>")]
    pub max_n_evals: usize,
    /// Existing pySecDec output directory to reuse instead of regenerating integration code.
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub reuse_existing_output: Option<String>,
}

impl Default for PySecDecSettings {
    fn default() -> Self {
        Self {
            quiet: true,
            relative_precision: 1.0e-7,
            min_n_evals: 10_000,
            max_n_evals: 1_000_000_000_000,
            reuse_existing_output: None,
        }
    }
}

#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct VakintSettings {
    /// FORM executable used by analytic Vakint evaluation backends.
    #[serde(skip_serializing_if = "is_default_form_path")]
    pub form_exe_path: String,
    /// Python interpreter used to invoke pySecDec.
    #[serde(skip_serializing_if = "is_default_python_path")]
    pub python_exe_path: String,
    /// Ordered Vakint evaluation backends tried for each matched vacuum integral.
    #[serde(skip_serializing_if = "is_default_vakint_evaluation_methods")]
    pub evaluation_methods: Vec<String>,
    /// MATAD-specific master-integral and substitution controls.
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub matad: MATADSettings,
    /// AlphaLoop-specific master-integral substitution controls.
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub alphaloop: AlphaLoopSettings,
    /// FMFT-specific master-integral expansion and substitution controls.
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub fmft: FMFTSettings,
    /// pySecDec precision, sample, and output-reuse controls.
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub pysecdec: PySecDecSettings,
    /// Decimal precision used when substituting backend results at run time.
    #[serde(skip_serializing_if = "is_usize::<16>")]
    pub run_time_decimal_precision: usize,
    /// Remove Vakint's temporary backend directory after a successful evaluation.
    #[serde(skip_serializing_if = "is_true")]
    pub clean_tmp_dir: bool,
    /// Parent directory for temporary Vakint backend files; the system default is used when absent.
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub temporary_directory: Option<String>,
    /// Named normalization convention applied to matched vacuum integrals.
    #[serde(skip_serializing_if = "is_default_vakint_normalization")]
    pub normalization: String,
    /// Extra symbolic normalization factor multiplied into every Vakint result.
    #[serde(skip_serializing_if = "is_minus_one_string")]
    pub additional_normalization: String,
}

impl VakintSettings {
    pub fn true_settings(&self) -> vakint::VakintSettings {
        vakint::VakintSettings {
            form_exe_path: self.form_exe_path.clone(),
            python_exe_path: self.python_exe_path.clone(),
            verify_numerator_identification: false,
            run_time_decimal_precision: self.run_time_decimal_precision as u32,
            allow_unknown_integrals: false,
            clean_tmp_dir: self.clean_tmp_dir,
            precision_for_input_float_rationalization:
                vakint::InputFloatRationalizationPrecision::FullPrecision,
            evaluation_order: vakint::EvaluationOrder(
                self.evaluation_methods
                    .iter()
                    .map(|a| match a.as_str() {
                        "alphaloop" => EvaluationMethod::AlphaLoop(AlphaLoopOptions {
                            susbstitute_masters: self.alphaloop.susbstitute_masters,
                        }),
                        "matad" => EvaluationMethod::MATAD(MATADOptions {
                            expand_masters: self.matad.expand_masters,
                            susbstitute_masters: self.matad.susbstitute_masters,
                            substitute_hpls: self.matad.substitute_hpls,
                            direct_numerical_substition: self.matad.direct_numerical_substition,
                        }),
                        "fmft" => EvaluationMethod::FMFT(FMFTOptions {
                            expand_masters: self.fmft.expand_masters,
                            susbstitute_masters: self.fmft.susbstitute_masters,
                        }),
                        "pysecdec" => EvaluationMethod::PySecDec(PySecDecOptions {
                            quiet: self.pysecdec.quiet,
                            relative_precision: self.pysecdec.relative_precision,
                            min_n_evals: self.pysecdec.min_n_evals as u64,
                            max_n_evals: self.pysecdec.max_n_evals as u64,
                            reuse_existing_output: self.pysecdec.reuse_existing_output.clone(),
                            ..Default::default()
                        }),
                        _ => panic!("Unknown vakint evaluation method: {}", a),
                    })
                    .collect(),
            ),
            use_dot_product_notation: false,
            temporary_directory: self.temporary_directory.clone(),
            epsilon_symbol: GS.dim_epsilon.get_name().into(),
            mu_r_sq_symbol: GS.mu_r_sq.get_name().into(),
            integral_normalization_factor: match self.normalization.as_str() {
                "MSbar" => vakint::LoopNormalizationFactor::MSbar,
                "FMFTandMATAD" => vakint::LoopNormalizationFactor::FMFTandMATAD,
                "pySecDec" => vakint::LoopNormalizationFactor::pySecDec,
                _ => vakint::LoopNormalizationFactor::Custom(self.normalization.clone()),
            },
            //Custom("1".to_string()),
            number_of_terms_in_epsilon_expansion: 5,
            // ..Default::default()
        }
    }
}

impl Default for VakintSettings {
    fn default() -> Self {
        Self {
            form_exe_path: "form".to_string(),
            python_exe_path: "python3".to_string(),
            // Supported evaluation methods: "alphaloop", "matad", "fmft" and "pysecdec"
            evaluation_methods: vec![
                "alphaloop".to_string(),
                "matad".to_string(),
                "fmft".to_string(),
            ],
            matad: MATADSettings::default(),
            alphaloop: AlphaLoopSettings::default(),
            fmft: FMFTSettings::default(),
            pysecdec: PySecDecSettings::default(),
            run_time_decimal_precision: 100,
            clean_tmp_dir: true,
            temporary_directory: None,
            normalization: "MSbar".to_string(),
            additional_normalization: "-1".to_string(),
        }
    }
}

#[derive(Debug, Clone, Default, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[cfg_attr(
    feature = "python_api",
    pyo3::pyclass(from_py_object, get_all, set_all)
)]
pub enum FinalIntegrandDimension {
    FourD,
    #[default]
    ThreeD,
}

impl Display for FinalIntegrandDimension {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FinalIntegrandDimension::FourD => write!(f, "4D"),
            FinalIntegrandDimension::ThreeD => write!(f, "3D"),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[cfg_attr(
    feature = "python_api",
    pyo3::pyclass(from_py_object, get_all, set_all)
)]
#[serde(default, deny_unknown_fields)]
pub struct UVgenerationSettings {
    /// Generate soft ultraviolet counterterms where required by the selected prescription.
    #[serde(skip_serializing_if = "is_true")]
    pub softct: bool,
    /// Generate analytically integrated ultraviolet counterterms.
    #[serde(skip_serializing_if = "is_true")]
    pub generate_integrated: bool,
    /// Subtract generated ultraviolet counterterms from the final integrand.
    #[serde(skip_serializing_if = "is_true")]
    pub subtract_uv: bool,
    /// Dimensional representation used for the final locally subtracted integrand.
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub final_integrand: FinalIntegrandDimension,
    /// Insert explicit marker functions around generated ultraviolet counterterms.
    #[serde(skip_serializing_if = "is_false")]
    pub add_marker: bool,
    /// Preserve ultraviolet marker functions after counterterm orchestration.
    #[serde(skip_serializing_if = "is_true")]
    pub keep_marker: bool,
    /// Rewrite Lorentz contractions into scalar inner products during ultraviolet processing.
    #[serde(skip_serializing_if = "is_true")]
    pub inner_products: bool,
    /// Algorithm used to enumerate and combine nested ultraviolet counterterms.
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub orchestrator: UVOrchestrator,
    /// Default and counterterm-specific renormalization approximation schemes.
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub renormalization_prescription: RenormalizationPrescriptionSettings,
    /// Vakint backend configuration for integrated ultraviolet counterterms.
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub vakint: VakintSettings,
}

impl Default for UVgenerationSettings {
    fn default() -> Self {
        UVgenerationSettings {
            softct: true,
            generate_integrated: true,
            subtract_uv: true,
            final_integrand: FinalIntegrandDimension::default(),
            inner_products: true,
            orchestrator: UVOrchestrator::default(),
            add_marker: false,
            keep_marker: true,
            renormalization_prescription: RenormalizationPrescriptionSettings::default(),
            vakint: VakintSettings::default(),
        }
    }
}

impl UVgenerationSettings {
    pub fn approximation_scheme_for(
        &self,
        ct_identifier: &CTIdentifier,
        dod: i32,
        has_massive_externals: bool,
    ) -> ApproximationType {
        self.renormalization_prescription.approximation_scheme_for(
            ct_identifier,
            dod,
            has_massive_externals,
        )
    }
}

mod ct_renormalization_overrides_serde {
    use super::CTRenormalizationRule;
    use serde::{Deserialize, Deserializer, Serialize, Serializer, de::Error};
    use std::collections::BTreeSet;

    pub fn serialize<S>(
        overrides: &Vec<CTRenormalizationRule>,
        serializer: S,
    ) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        overrides.serialize(serializer)
    }

    pub fn deserialize<'de, D>(deserializer: D) -> Result<Vec<CTRenormalizationRule>, D::Error>
    where
        D: Deserializer<'de>,
    {
        let rules = Vec::<CTRenormalizationRule>::deserialize(deserializer)?;
        let mut seen = BTreeSet::new();

        for rule in &rules {
            if !seen.insert(rule.ct_identifier.clone()) {
                return Err(D::Error::custom(format!(
                    "duplicate CTIdentifier in renormalization_prescription.overrides: {:?}",
                    rule.ct_identifier
                )));
            }
        }

        Ok(rules)
    }
}

#[cfg(test)]
mod tests {
    use super::{
        ApproximationType, CTIdentifier, CTRenormalizationRule,
        RenormalizationPrescriptionSettings, UVOrchestrator, UVgenerationSettings,
    };
    use std::collections::BTreeSet;

    fn pdg_set(values: impl IntoIterator<Item = isize>) -> BTreeSet<isize> {
        values.into_iter().collect()
    }

    #[test]
    fn renormalization_prescription_prefers_exact_internal_match_over_wildcard() {
        let exact_identifier = CTIdentifier::new(pdg_set([1]), Some(pdg_set([1, 22])));
        let wildcard_identifier = CTIdentifier::new(pdg_set([1]), None);
        let settings = UVgenerationSettings {
            renormalization_prescription: RenormalizationPrescriptionSettings {
                overrides: vec![
                    CTRenormalizationRule::new(
                        wildcard_identifier,
                        ApproximationType::Unsubtracted,
                    ),
                    CTRenormalizationRule::new(exact_identifier.clone(), ApproximationType::OS),
                ],
                ..Default::default()
            },
            ..Default::default()
        };

        assert_eq!(
            settings.approximation_scheme_for(&exact_identifier, 0, false),
            ApproximationType::OS
        );
    }

    #[test]
    fn renormalization_prescription_matches_wildcard_internal_rule() {
        let candidate_identifier = CTIdentifier::new(pdg_set([1]), Some(pdg_set([1, 22])));
        let settings = UVgenerationSettings {
            renormalization_prescription: RenormalizationPrescriptionSettings {
                overrides: vec![CTRenormalizationRule::new(
                    CTIdentifier::new(pdg_set([1]), None),
                    ApproximationType::OS,
                )],
                ..Default::default()
            },
            ..Default::default()
        };

        assert_eq!(
            settings.approximation_scheme_for(&candidate_identifier, 0, false),
            ApproximationType::OS
        );
    }

    #[test]
    fn renormalization_prescription_can_select_unsubtracted() {
        let candidate_identifier = CTIdentifier::new(pdg_set([1]), Some(pdg_set([1, 22])));
        let settings = UVgenerationSettings {
            renormalization_prescription: RenormalizationPrescriptionSettings {
                overrides: vec![CTRenormalizationRule::new(
                    CTIdentifier::new(pdg_set([1]), None),
                    ApproximationType::Unsubtracted,
                )],
                ..Default::default()
            },
            ..Default::default()
        };

        assert_eq!(
            settings.approximation_scheme_for(&candidate_identifier, 0, false),
            ApproximationType::Unsubtracted
        );
    }

    #[test]
    fn renormalization_prescription_defaults_to_msbar() {
        let settings = UVgenerationSettings::default();

        assert_eq!(
            settings.approximation_scheme_for(
                &CTIdentifier::new(pdg_set([1]), Some(pdg_set([1, 22]))),
                1,
                false
            ),
            ApproximationType::MUV
        );
    }

    #[test]
    fn renormalization_prescription_accepts_blanket_pole_part() {
        let settings: RenormalizationPrescriptionSettings = toml::from_str(
            r#"
log_divergent = "PolePart"
massive_power_divergent = "PolePart"
massless_power_divergent = "PolePart"
"#,
        )
        .unwrap();

        assert_eq!(settings.log_divergent, ApproximationType::PolePart);
        assert_eq!(
            settings.massive_power_divergent,
            ApproximationType::PolePart
        );
        assert_eq!(
            settings.massless_power_divergent,
            ApproximationType::PolePart
        );
    }

    #[test]
    fn renormalization_prescription_uses_log_divergent_bucket() {
        let settings = UVgenerationSettings {
            renormalization_prescription: RenormalizationPrescriptionSettings {
                log_divergent: ApproximationType::OS,
                massive_power_divergent: ApproximationType::IR,
                massless_power_divergent: ApproximationType::Unsubtracted,
                ..Default::default()
            },
            ..Default::default()
        };

        assert_eq!(
            settings.approximation_scheme_for(
                &CTIdentifier::new(pdg_set([1]), Some(pdg_set([1, 22]))),
                0,
                false
            ),
            ApproximationType::OS
        );
    }

    #[test]
    fn renormalization_prescription_uses_massive_power_divergent_bucket() {
        let settings = UVgenerationSettings {
            renormalization_prescription: RenormalizationPrescriptionSettings {
                massive_power_divergent: ApproximationType::OS,
                massless_power_divergent: ApproximationType::IR,
                ..Default::default()
            },
            ..Default::default()
        };

        assert_eq!(
            settings.approximation_scheme_for(
                &CTIdentifier::new(pdg_set([1]), Some(pdg_set([1, 22]))),
                1,
                true
            ),
            ApproximationType::OS
        );
    }

    #[test]
    fn renormalization_prescription_uses_massless_power_divergent_bucket() {
        let settings = UVgenerationSettings {
            renormalization_prescription: RenormalizationPrescriptionSettings {
                massive_power_divergent: ApproximationType::OS,
                massless_power_divergent: ApproximationType::IR,
                ..Default::default()
            },
            ..Default::default()
        };

        assert_eq!(
            settings.approximation_scheme_for(
                &CTIdentifier::new(pdg_set([1]), Some(pdg_set([1, 22]))),
                1,
                false
            ),
            ApproximationType::IR
        );
    }

    #[test]
    fn renormalization_prescription_rejects_duplicate_overrides() {
        let toml = r#"
[[overrides]]
prescription = "OS"

[overrides.ct_identifier]
external_pdg_set = [1]

[[overrides]]
prescription = "IR"

[overrides.ct_identifier]
external_pdg_set = [1]
"#;

        let err = toml::from_str::<RenormalizationPrescriptionSettings>(toml).unwrap_err();
        assert!(
            err.to_string()
                .contains("duplicate CTIdentifier in renormalization_prescription.overrides"),
            "{err}"
        );
    }

    #[test]
    fn orchestrator_defaults_to_hedge_poset_and_accepts_explicit_modes() {
        assert_eq!(UVOrchestrator::default(), UVOrchestrator::HedgePoset);
        assert_eq!(
            UVgenerationSettings::default().orchestrator,
            UVOrchestrator::HedgePoset
        );
        assert_eq!(
            toml::from_str::<UVgenerationSettings>("")
                .unwrap()
                .orchestrator,
            UVOrchestrator::HedgePoset
        );

        let legacy = UVgenerationSettings {
            orchestrator: UVOrchestrator::LegacyDagForest,
            ..Default::default()
        };
        assert_eq!(legacy.orchestrator, UVOrchestrator::LegacyDagForest);
        let legacy_toml = toml::to_string(&legacy).unwrap();
        assert!(legacy_toml.contains("orchestrator = \"legacy_dag_forest\""));
        assert_eq!(
            toml::from_str::<UVgenerationSettings>(&legacy_toml).unwrap(),
            legacy
        );

        let compare = UVgenerationSettings {
            orchestrator: UVOrchestrator::Compare,
            ..Default::default()
        };
        assert_eq!(compare.orchestrator, UVOrchestrator::Compare);
    }
}
