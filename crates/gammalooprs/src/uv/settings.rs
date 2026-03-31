use std::collections::{BTreeMap, BTreeSet};

use crate::utils::{
    GS,
    serde_utils::{
        IsDefault, is_default_form_path, is_default_pysecdec_relative_precision,
        is_default_python_path, is_default_vakint_evaluation_methods,
        is_default_vakint_normalization, is_false, is_one_string, is_true, is_usize,
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
pub enum RenormalizationScheme {
    #[default]
    #[serde(rename = "MSbar", alias = "msbar")]
    MSbar,
    #[serde(rename = "OS", alias = "os")]
    OS,
    #[serde(rename = "Unsubtracted", alias = "unsubtracted")]
    Unsubtracted,
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
    pub renormalization_scheme: RenormalizationScheme,
}

impl CTRenormalizationRule {
    pub fn new(ct_identifier: CTIdentifier, renormalization_scheme: RenormalizationScheme) -> Self {
        Self {
            ct_identifier,
            renormalization_scheme,
        }
    }
}

#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct MATADSettings {
    #[serde(skip_serializing_if = "is_true")]
    pub expand_masters: bool,
    #[serde(skip_serializing_if = "is_true")]
    pub susbstitute_masters: bool,
    #[serde(skip_serializing_if = "is_true")]
    pub substitute_hpls: bool,
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
    #[serde(skip_serializing_if = "is_true")]
    pub expand_masters: bool,
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
    #[serde(skip_serializing_if = "is_true")]
    pub quiet: bool,
    #[serde(skip_serializing_if = "is_default_pysecdec_relative_precision")]
    pub relative_precision: f64,
    #[serde(skip_serializing_if = "is_usize::<10_000>")]
    pub min_n_evals: usize,
    #[serde(skip_serializing_if = "is_usize::<1_000_000_000_000>")]
    pub max_n_evals: usize,
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
    #[serde(skip_serializing_if = "is_default_form_path")]
    pub form_exe_path: String,
    #[serde(skip_serializing_if = "is_default_python_path")]
    pub python_exe_path: String,
    #[serde(skip_serializing_if = "is_default_vakint_evaluation_methods")]
    pub evaluation_methods: Vec<String>,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub matad: MATADSettings,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub alphaloop: AlphaLoopSettings,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub fmft: FMFTSettings,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub pysecdec: PySecDecSettings,
    #[serde(skip_serializing_if = "is_usize::<16>")]
    pub run_time_decimal_precision: usize,
    #[serde(skip_serializing_if = "is_true")]
    pub clean_tmp_dir: bool,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub temporary_directory: Option<String>,
    #[serde(skip_serializing_if = "is_default_vakint_normalization")]
    pub normalization: String,
    #[serde(skip_serializing_if = "is_one_string")]
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
            additional_normalization: "1".to_string(),
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
    #[serde(skip_serializing_if = "is_true")]
    pub softct: bool,
    #[serde(skip_serializing_if = "is_true")]
    pub generate_integrated: bool,
    #[serde(skip_serializing_if = "is_true")]
    pub subtract_uv: bool,
    #[serde(skip_serializing_if = "is_false")]
    pub only_integrated: bool,
    #[serde(skip_serializing_if = "is_false")]
    pub pole_part: bool,
    #[serde(skip_serializing_if = "is_false")]
    pub add_sigma: bool,
    #[serde(skip_serializing_if = "is_true")]
    pub inner_products: bool,
    #[serde(skip_serializing_if = "is_true")]
    pub use_legacy: bool,
    #[serde(skip_serializing_if = "is_true")]
    pub cached_integrated: bool,
    #[serde(
        default,
        skip_serializing_if = "IsDefault::is_default",
        with = "ct_renormalization_map_serde"
    )]
    #[schemars(with = "Vec<CTRenormalizationRule>")]
    pub renormalization_schemes: BTreeMap<CTIdentifier, RenormalizationScheme>,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub vakint: VakintSettings,
}

impl Default for UVgenerationSettings {
    fn default() -> Self {
        UVgenerationSettings {
            softct: true,
            generate_integrated: true,
            subtract_uv: true,
            pole_part: false,
            only_integrated: false,
            inner_products: true,
            use_legacy: true,
            cached_integrated: true,
            add_sigma: false,
            renormalization_schemes: BTreeMap::default(),
            vakint: VakintSettings::default(),
        }
    }
}

impl UVgenerationSettings {
    pub fn renormalization_scheme_for(
        &self,
        ct_identifier: &CTIdentifier,
    ) -> RenormalizationScheme {
        if let Some(scheme) = self.renormalization_schemes.get(ct_identifier) {
            return *scheme;
        }

        self.renormalization_schemes
            .iter()
            .find_map(|(rule, scheme)| {
                (rule.internal_pdg_set.is_none() && rule.matches(ct_identifier)).then_some(*scheme)
            })
            .unwrap_or_default()
    }
}

mod ct_renormalization_map_serde {
    use super::{CTIdentifier, CTRenormalizationRule, RenormalizationScheme};
    use serde::{Deserialize, Deserializer, Serialize, Serializer, de::Error};
    use std::collections::BTreeMap;

    pub fn serialize<S>(
        map: &BTreeMap<CTIdentifier, RenormalizationScheme>,
        serializer: S,
    ) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        map.iter()
            .map(|(ct_identifier, renormalization_scheme)| {
                CTRenormalizationRule::new(ct_identifier.clone(), *renormalization_scheme)
            })
            .collect::<Vec<_>>()
            .serialize(serializer)
    }

    pub fn deserialize<'de, D>(
        deserializer: D,
    ) -> Result<BTreeMap<CTIdentifier, RenormalizationScheme>, D::Error>
    where
        D: Deserializer<'de>,
    {
        let rules = Vec::<CTRenormalizationRule>::deserialize(deserializer)?;
        let mut map = BTreeMap::new();

        for rule in rules {
            if map
                .insert(rule.ct_identifier.clone(), rule.renormalization_scheme)
                .is_some()
            {
                return Err(D::Error::custom(format!(
                    "duplicate CTIdentifier in renormalization_schemes: {:?}",
                    rule.ct_identifier
                )));
            }
        }

        Ok(map)
    }
}

#[cfg(test)]
mod tests {
    use super::{CTIdentifier, RenormalizationScheme, UVgenerationSettings};
    use std::collections::{BTreeMap, BTreeSet};

    fn pdg_set(values: impl IntoIterator<Item = isize>) -> BTreeSet<isize> {
        values.into_iter().collect()
    }

    #[test]
    fn renormalization_scheme_prefers_exact_internal_match_over_wildcard() {
        let exact_identifier = CTIdentifier::new(pdg_set([1]), Some(pdg_set([1, 22])));
        let wildcard_identifier = CTIdentifier::new(pdg_set([1]), None);
        let settings = UVgenerationSettings {
            renormalization_schemes: BTreeMap::from([
                (wildcard_identifier, RenormalizationScheme::Unsubtracted),
                (exact_identifier.clone(), RenormalizationScheme::OS),
            ]),
            ..Default::default()
        };

        assert_eq!(
            settings.renormalization_scheme_for(&exact_identifier),
            RenormalizationScheme::OS
        );
    }

    #[test]
    fn renormalization_scheme_matches_wildcard_internal_rule() {
        let candidate_identifier = CTIdentifier::new(pdg_set([1]), Some(pdg_set([1, 22])));
        let settings = UVgenerationSettings {
            renormalization_schemes: BTreeMap::from([(
                CTIdentifier::new(pdg_set([1]), None),
                RenormalizationScheme::OS,
            )]),
            ..Default::default()
        };

        assert_eq!(
            settings.renormalization_scheme_for(&candidate_identifier),
            RenormalizationScheme::OS
        );
    }

    #[test]
    fn renormalization_scheme_can_select_unsubtracted() {
        let candidate_identifier = CTIdentifier::new(pdg_set([1]), Some(pdg_set([1, 22])));
        let settings = UVgenerationSettings {
            renormalization_schemes: BTreeMap::from([(
                CTIdentifier::new(pdg_set([1]), None),
                RenormalizationScheme::Unsubtracted,
            )]),
            ..Default::default()
        };

        assert_eq!(
            settings.renormalization_scheme_for(&candidate_identifier),
            RenormalizationScheme::Unsubtracted
        );
    }

    #[test]
    fn renormalization_scheme_defaults_to_msbar() {
        let settings = UVgenerationSettings::default();

        assert_eq!(
            settings.renormalization_scheme_for(&CTIdentifier::new(
                pdg_set([1]),
                Some(pdg_set([1, 22]))
            )),
            RenormalizationScheme::MSbar
        );
    }
}
